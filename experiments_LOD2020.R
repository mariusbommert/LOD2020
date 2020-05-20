# In this file the experiments are described. 
# This is the second part for reproducing the results for LOD 2020.

# load all packages which are needed for the following code:
library(batchtools)
library(DiceKriging)
library(mvtnorm)
library(parallelMap)

# set a seed for reproducability:
set.seed(1)

# make a registry for the following jobs:
reg = makeExperimentRegistry(file.dir = "experiments_LOD2020")

# function for loading Kriging models: 
load_model = function(data, job, n) {
  # path to rds files from results generated with generate_models_LOD2020.R needed:
  model = readRDS(file = paste0("generate_models/results/", n, ".rds"))
  return(model)
}

# function for 2-dimensional integration for calculating E and V given a model:
integration2d = function(x1, x2, model, mean, sigma, lower1 = 0, upper1 = 100, 
                         lower2 = 20, upper2 = 120, subdivisions = 100L, 
                         rel.tol = 1e-2) {
  # integrand for expectation:
  integrand1 = function(x1, x2, c1, c2, model, mean, sigma) {
    DiceKriging::predict(model, newdata = data.frame(x1 = x1, x2 = x2, c1 = c1, 
                                                     c2 = c2, row.names = NULL), 
                         type = "UK")$mean * 
      mvtnorm::dmvnorm(cbind(c1, c2), mean = mean, sigma = sigma)
  }
  
  # integrand for variance:
  integrand2 = function(x1, x2, c1, c2, model, mean, sigma, e) {
    (DiceKriging::predict(model, newdata = data.frame(x1 = x1, x2 = x2, c1 = c1, 
                                                      c2 = c2, row.names = NULL), 
                          type = "UK")$mean - e)^2 * 
      mvtnorm::dmvnorm(cbind(c1, c2), mean = mean, sigma = sigma)
  }
  
  # integration for expectation:
  e = integrate(function(z) {
    sapply(z, function(b) {
      integrate(function(a) integrand1(x1, x2, a, b, model, mean, sigma), 
                lower1, upper1, subdivisions = subdivisions, rel.tol = rel.tol)$value
    })
  }, lower2, upper2, subdivisions = subdivisions, rel.tol = rel.tol)$value
  
  # integration for variance: 
  v = integrate(function(z) {
    sapply(z, function(b) {
      integrate(function(a) integrand2(x1, x2, a, b, model, mean, sigma, e), 
                lower1, upper1, subdivisions = subdivisions, rel.tol = rel.tol)$value
    })
  }, lower2, upper2, subdivisions = subdivisions, rel.tol = rel.tol)$value
  return(data.frame(e = e, v = v))
}

# function for parallel integration using integration2d:
integration2d_ev = function(data, job, instance, rel.tol = 1e-2, ncpus = 10) { 
  library(parallelMap)

  # grid for calculating E and V:
  x = expand.grid(x1 = seq(0, 100, length.out = 101), 
                  x2 = seq(20, 120, length.out = 101))
  # split grid for cpus:
  sx = split(x, rep(1:ncpus, each = ceiling(nrow(x) / ncpus)))
  
  # estimate mean and sigma:
  mean = colMeans(instance@X[, c("c1", "c2")])
  sigma = cov(instance@X[, c("c1", "c2")])

  # export for parallel calculation: 
  parallelExport("instance", "x", "sx", "mean", "sigma", "integration2d")

  # calculate expectation and variance using integration2d:
  res = parallelLapply(sx, function(y) {
    r = apply(y, 1, function(z) {
      tmp = integration2d(x1 = z[1], x2 = z[2], model = instance, mean = mean,
                          sigma = sigma, lower1 = 0, upper1 = 100, lower2 = 20, 
                          upper2 = 120, subdivisions = 100L, rel.tol = rel.tol)
      return(tmp)
    })
    do.call(rbind, r)
  })
  res = do.call(rbind, res)
  return(res)
}

# function for sampling for calculating E and V given a model: 
sampling2d = function (model, x1, x2, rn) {
  # calculate values for given x and C:
  value = DiceKriging::predict(object = model, 
                               newdata = data.frame(x1 = x1, x2 = x2, c1 = rn[, 1], 
                                                    c2 = rn[, 2], row.names = NULL), 
                               type = "UK")$mean
  
  # calculate expectation and variance:
  e = mean(value)
  v = var(value)
  return(data.frame(e = e, v = v))
}

# function for parallel sampling using sampling2d:
sampling2d_ev = function(data, job, instance, nsamples = 1e4, ncpus = 10) { 
  library(parallelMap)

  # grid for calculating E and V:
  x = expand.grid(x1 = seq(0, 100, length.out = 101), 
                  x2 = seq(20, 120, length.out = 101))
  # split grid for cpus:
  sx = split(x, rep(1:ncpus, each = ceiling(nrow(x) / ncpus)))
  
  # estimate mean and sigma and generate random numbers:
  mean = colMeans(instance@X[, c("c1", "c2")])
  sigma = cov(instance@X[, c("c1", "c2")])
  rn = mvtnorm::rmvnorm(nsamples, mean = mean, sigma = sigma)

  # export for parallel calculation: 
  parallelExport("instance", "x", "sx", "rn", "sampling2d")

  # calculate expectation and variance using sampling2d:
  res = parallelLapply(sx, function(y) {
    r = apply(y, 1, function(z) {
      tmp = sampling2d(model = instance, x1 = z[1], x2 = z[2], rn = rn)
      return(tmp)
    })
    do.call(rbind, r)
  })
  res = do.call(rbind, res)
  return(res)
}

# add problems and algorithms: 
addProblem(name = "load_model", fun = load_model, seed = 1)
addAlgorithm(name = "integration2d_ev", fun = integration2d_ev)
addAlgorithm(name = "sampling2d_ev", fun = sampling2d_ev)

# add experiments: 
pdes = expand.grid(n = 1:22500)
addExperiments(prob.designs = list(load_model = pdes),
               algo.designs = list(integration2d_ev = data.frame()), repls = 1)
addExperiments(prob.designs = list(load_model = pdes),
               algo.designs = list(sampling2d_ev = data.frame()), repls = 1)

# ids for experiments used in LOD 2020:
experimentids = c(seq(from = 22, to = 22500, by = 45), 
                  seq(from = 22, to = 22500, by = 45) + 22500)

# submit jobs for LOD 2020, maybe more resources are needed:
submitJobs(ids = experimentids, resources = list(walltime = 8 * 3600, memory = 6000,
                                                 ncpus = 10, pm.backend = "multicore"))
res_25_25_m1_3 = reduceResultsDataTable(experimentids, missing.val = list())
