# In this file the generation of the Kriging models is described.
# This is the first part for reproducing the results for LOD 2020.

# load all packages which are needed for the following code:
library(batchtools)
library(DiceKriging)
library(lhs)
library(mvtnorm)
library(ParamHelpers)

# set a seed for reproducability:
set.seed(1)

# make a registry for the following jobs:
reg = makeRegistry(file.dir = "generate_models")

# function for generating the x-values:
generate2dx = function(n, lower1, upper1, lower2, upper2, type) {
  if (is.null(type)) stop("type must be given.")
  lh = switch(type, 
              improved = lhs::improvedLHS, 
              random = lhs::randomLHS, 
              maximin = lhs::maximinLHS)
  design = ParamHelpers::generateDesign(n, par.set = ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam(id = "x1", lower = lower1, upper = upper1), 
    ParamHelpers::makeNumericParam(id = "x2", lower = lower2, upper = upper2)), fun = lh)
  return(design)
}

# objective function:
f_nv_2dx2dc = function(x1, x2, c1, c2, p = 5, q = 3, r = 4, s = 2) {
  # calculate sum of two newsvendor models and return it:
  return(p * pmin(x1, c1) - q * x1 + r * pmin(x2, c2) - s * x2)
}

# function for generating models:
generate_models = function(n, v1, v2, cv) {
  # generate data set with 2-dimensional x and 2-dimensional c:
  x2d = generate2dx(n = n, lower1 = 0, upper1 = 100, lower2 = 20, upper2 = 120,
                    type = "improved")
  c2d = mvtnorm::rmvnorm(n = n, mean = c(50, 70), 
                         sigma = matrix(data = c(v1, cv, cv, v2), nrow = 2))
  data = cbind(x2d, c2d)
  data = cbind(data, apply(X = data, MARGIN = 1, FUN = function(x) 
    f_nv_2dx2dc(x[1], x[2], x[3], x[4])))
  colnames(data) = c("x1", "x2", "c1", "c2", "value")
  
  # generate Kriging model:
  model = DiceKriging::km(formula = ~1, design = data[, 1:4], response = data[, 5],
                          covtype = "matern3_2", nugget.estim = TRUE, multistart = 100,
                          control = list(trace = 0, maxit = 1000,
                                         factr = 1e7, pop.size = 20))
  
  return(model)
}

# determine parameter combinations for models:
computecovv = function(v1, v2, corv) {
  return(corv * (sqrt(v1) * sqrt(v2)))
}
v1 = c(9, 25, 49)
v2 = c(9, 25, 49)
v1v2 = expand.grid(v1 = v1, v2 = v2)
corv = c(-2/3, -1/3, 0, 1/3, 2/3)
pars = data.frame(sapply(v1v2, rep, each = 5),
                  cv = as.numeric(do.call(rbind, lapply(corv, function(x)
                    apply(X = v1v2, MARGIN = 1, FUN = function(y) 
                      computecovv(y[1], y[2], x))))))
n = data.frame(n = rep(c(50, 100, 200, 500, 1000), each = 100))
pars = merge(x = pars, y = n)

# register jobs on batchsystem:
batchMap(fun = generate_models, args = pars)

# ids for models used in LOD:
modelids = seq(from = 22, to = 22500, by = 45)

# submit jobs needed for LOD:
# maybe it is necessary to use higher values for the resources on the used system:
submitJobs(ids = modelids, resources = list(walltime = 8 * 3600, memory = 7500))
