# In this file the evaluation of the experiments is described.
# This is the third part for reproducing the results for LOD 2020.

# load all packages which are needed for the following code:
library(ecr)
library(ggplot2)
library(data.table)

# set a seed for reproducability:
set.seed(1)

# function for estimating expectation and variance:
estimate_ev = function(f, x1, x2, rn, ...){
  fv = f(x1 = x1, x2 = x2, c1 = rn[, 1], c2 = rn[, 2], ...)
  e = mean(fv)
  v = var(fv)
  return(data.frame(e = e, v = v))
}

# objective function: 
f_nv_2dx2dc = function(x1, x2, c1, c2, p = 5, q = 3, r = 4, s = 2) {
  # calculate sum of two newsvendor models and return it:
  return(p * pmin(x1, c1) - q * x1 + r * pmin(x2, c2) - s * x2)
}

# calculating true values for expectation and variance using 1e7 random numbers:
rn = mvtnorm::rmvnorm(n = 1e7, mean = c(50, 70), 
                      sigma = matrix(c(25, -25/3, -25/3, 25), nrow = 2))
grid = expand.grid(x1 = seq(0, 100, length.out = 101), 
                   x2 = seq(20, 120, length.out = 101))
fv2 = apply(grid, 1, function(y) {
  estimate_ev(f = f_nv_2dx2dc, x1 = y[1], x2 = y[2], rn = rn)
})
fv_25_25_m1_3 = do.call(rbind, fv2)
true_25_25_m1_3 = cbind(-fv_25_25_m1_3[, 1], sqrt(fv_25_25_m1_3[, 2]))
colnames(true_25_25_m1_3) = c("E", "S")

# save(true_25_25_m1_3, file = "true_25_25_m1_3.RData")
# This object is available at GitHub.
# load(true_25_25_m1_3) can be used to load it. 

################################################################################

# determine which true values are non-dominated:
ind = which(true_25_25_m1_3[, 2] == 0)
ind = ind[-48]
ind_25_25_m1_3 = nondominated(t(true_25_25_m1_3))
ind2 = numeric(47)
unique(true_25_25_m1_3[ind, 1])
for(i in 1:47) {
  ind2[i] = which(true_25_25_m1_3[, 1] == unique(true_25_25_m1_3[ind, 1])[i])[1]
}
ind_25_25_m1_3[ind2] = TRUE

# combine all non-dominated true values:
t25_25_m1_3 = as.data.frame(true_25_25_m1_3)
t25_25_m1_3 = cbind(t25_25_m1_3, col = ifelse(ind_25_25_m1_3, "black", "red"))

# plot the truth:
pdf("true.pdf", width = 8, height = 4)
ggplot(t25_25_m1_3[!ind_25_25_m1_3, ], mapping = aes(x = E, y = S, col = col)) + geom_point() +
  geom_point(data = t25_25_m1_3[ind_25_25_m1_3, ], mapping = aes(x = E, y = S, col = col)) +
  theme_bw() + labs(title = "", x = "-E(f(x, C))", y = "S(f(x, C))") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.title = element_blank(), legend.position = "bottom",
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)) +
  scale_color_manual("", values = c("red", "black"),
                     labels = c("non dominated", "dominated")) +
  guides(colour = guide_legend(nrow = 1))
dev.off()

################################################################################

# use res_25_25_m1_3 which is the result from experiments_LOD2020:
# only some necessary columns are used for the evaluation: 
res_25_25_m1_3 = res_25_25_m1_3[, list(job.id, error, time.running, problem, 
                                       algorithm, n, rel.tol, ncpus, result)]
# save(res_25_25_m1_3, file = "res_25_25_m1_3.RData")
# This object is available at GitHub. 
# Because of the storage limit res_25_25_m1_3 is splitted in two parts:
# The following code can be used for generating res_25_25_m1_3:
# load(res_25_25_m1_3_integration)
# load(res_25_25_m1_3_sampling)
# res_25_25_m1_3 = rbind(res_25_25_m1_3_integration, res_25_25_m1_3_sampling)

# determine parameter combinations for models:
v1 = c(9, 25, 49)
v2 = c(9, 25, 49)
v1v2 = expand.grid(v1 = v1, v2 = v2)
corv = c(-2/3, -1/3, 0, 1/3, 2/3)
computecovv = function(v1, v2, corv) {
  return(corv * (sqrt(v1) * sqrt(v2)))
}
pars = data.frame(sapply(v1v2, rep, each = 5),
                  cv = as.numeric(do.call(rbind, lapply(corv, function(x)
                    apply(v1v2, 1, function(y) computecovv(y[1], y[2], x))))))
n = data.frame(n = rep(c(50, 100, 200, 500, 1000), each = 100))
pars = merge(x = pars, y = n)

# ids used for LOD 2020: 
ids = c(seq(22, 22500, 45), seq(22, 22500, 45) + 22500)

# combine res_25_25_m1_3 and corresponding parameters:
res_25_25_m1_3 = cbind(res_25_25_m1_3, rbind(pars[ids, ][1:500, ], pars[ids, ][1:500, ]))
# rename column 6: 
colnames(res_25_25_m1_3)[6] = "id"

# determine which jobs are done and which jobs are not done: 
notdone_25_25_m1_3 = which(!is.na(res_25_25_m1_3$error))
done_25_25_m1_3 = which(!(1:1000) %in% notdone_25_25_m1_3)

################################################################################

# function for calculating Euclidean distance: 
euclid = function (x, y) {
  return(dist(x = rbind(x, y), method = "euclidean"))
}

# function for calculating mean Euclidean distance:
calculate_distall = function (est, true) {
  est = as.data.table(est)
  est = est[, list(E, V)]
  est = apply(est, 2, as.numeric)
  true = apply(true, 2, as.numeric)
  distall = mean(vapply(X = 1:nrow(est), 
                        FUN = function(x) euclid(as.numeric(est[x, ]), 
                                                 as.numeric(true[x, ])), 
                        FUN.VALUE = numeric(1)))
  return(distall)
}

# function for calculating quality_measure mean Euclidean distance for given data:   
calculate_qual_2dx2dc = function(res, true, done, n = 1000) {
  qual = numeric(n)
  for(i in 1:n) {
    if(i %in% done) {
      # if job done calculate mean Euclidean distance:
      est = data.frame(E = -res$result[[i]]$e, V = sqrt(res$result[[i]]$v))
      qual[i] = calculate_distall(est, true)
    } else {
      # else use NA:
      qual[i] = NA
    }
  }
  return(qual)
}


# calculate quality measure mean Euclidean distance for given data: 
qual_25_25_m1_3 = calculate_qual_2dx2dc(res_25_25_m1_3, true_25_25_m1_3, done_25_25_m1_3)
# save(qual_25_25_m1_3, file = "qual_25_25_m1_3.RData")
# This object is available at GitHub.
# load(qual_25_25_m1_3) can be used to load it. 

# combine res_25_25_m1_3 and qual_25_25_m1_3:
qual_25_25_m1_3 = cbind(res_25_25_m1_3, distall = qual_25_25_m1_3)

# prepare quality for plotting: 
qual = qual_25_25_m1_3
qual$algorithm = rep(c("integration", "sampling"), each = 500)
qual$n = factor(paste0("n = ", qual$n), levels = paste0("n = ", unique(qual$n)))

# generate plot for mean Euclidean distance:
pdf("distall_2dx2dc.pdf", width = 9, height = 3)
ggplot(qual, aes(x = algorithm, y = distall)) + geom_boxplot() +
  facet_wrap(~ n, nrow = 1) +
  labs(x = "Algorithm", y = "Mean Euclidean distance")  +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 13))
dev.off()

# calculate difference in mean Euclidean distance:
diff_distall = data.frame(diff = qual$distall[1:500] - qual$distall[501:1000],
                          n = qual$n[1:500])

# generate plot for differences in mean Euclidean distance: 
pdf("diff_distall_2dx2dc.pdf", width = 9, height = 3)
ggplot(diff_distall, aes(x = n, y = diff)) + geom_boxplot() +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  labs(x = "Number of observations for building the metamodel", 
       y = "Difference in mean Eucl. distance") +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 13))
dev.off()

################################################################################

# function for calculating eaf: 
eaf_ev = function (data, percentiles) {
  eaf = eaf::eafs(points = data[, c("E", "S")], 
                  sets = data[, "group"], percentiles = percentiles)
  eaf = as.data.frame(eaf)
  colnames(eaf) = c("E", "S", "group")
  return(eaf)
}

# function for calculating different eafs: 
eaf_5_50_95 = function(res, algorithm, n) {
  ind = res$algorithm == algorithm & res$n == n
  tmp = res[ind, ]
  tmp = do.call(rbind, lapply(1:nrow(tmp), function(x) {
    if(!is.null(tmp$result[x][[1]])) {
      cbind(tmp$result[x][[1]], x)
    }
  }))
  tmp[, 1] = -tmp[, 1]
  tmp[, 2] = sqrt(tmp[, 2])
  colnames(tmp) = c("E", "S", "group")
  eaf = eaf_ev(data = tmp, percentiles = c(5, 50, 95))
  eaf$n = rep(n, nrow(eaf))
  return(eaf)
}

# function for plotting eafs: 
eaf_sequential_5_50_95_facet = function(res, true, n = c(50, 100, 200, 500, 1000)) {
  grid = expand.grid(algorithm = c("integration", "sampling"), n = n)
  
  tmp = vector("list", length(grid))
  for(i in seq_len(nrow(grid))) {
    d = eaf_5_50_95(res = res, algorithm = grid[i, 1], n = as.numeric(grid[i, 2]))
    tmp[[i]] = cbind(d, algorithm = grid[i, 1])
  }
  data = do.call(rbind, tmp)
  
  ind = ecr::nondominated(t(true))
  ind[ind2] = TRUE
  true = true[ind, ]
  true = as.data.frame(true)
  
  true = cbind(true[rep(1:nrow(true), length(n)), ],
               group = "true", n = rep(n, each = nrow(true)), algorithm = "true")
  data2 = rbind(data, true)
  data2 = as.data.table(data2)
  data2[, ga := paste0(algorithm, " ", group)]
  
  data2$ga[data2$ga == "true true"] = "true"
  data2$ga = factor(data2$ga, levels = c("integration 5", "sampling 5",
                                         "integration 50", "sampling 50",
                                         "integration 95", "sampling 95",
                                         "true"))
  data2$n = factor(paste0("n = ", data2$n), levels = paste0("n = ", unique(data2$n)))

  data = data2[data2$algorithm != "true", ]
  true = data2[data2$algorithm == "true", ]
  
  p = ggplot(data = data, aes(x = E, y = S, colour = ga, group = ga, linetype = ga)) +
    geom_step(size = 1) + geom_line(data = true, aes(x = E, y = S, colour = ga, 
                                                     group = ga), size = 1) +
    labs(title = "", x = "-E(f(x, C))", y = "S(f(x, C))") + facet_wrap(~ n, nrow = 2) +
    theme_bw() +
    theme(legend.key = element_rect(fill = "transparent", colour = NA),
          legend.key.width = unit(2, "line"),
          legend.title = element_blank(), legend.position = "bottom",
          legend.text = element_text(size = 16),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 15)) +
    scale_color_manual(values = c("red", "orange", "blue", "cyan", 
                                  "darkgreen", "chartreuse", "black"),
                       breaks = c("integration 5", "sampling 5",
                                  "integration 50", "sampling 50",
                                  "integration 95", "sampling 95",
                                  "true")) +
    scale_linetype_manual(values = c("dashed", "dotted", "dashed", "dotted", 
                                     "dashed", "dotted", "solid"),
                          breaks = c("integration 5", "sampling 5",
                                     "integration 50", "sampling 50",
                                     "integration 95", "sampling 95",
                                     "true"))
  print(p)
}

# preparing results for plot: 
res = res_25_25_m1_3
res$algorithm = rep(c("integration", "sampling"), each = 500)

# generating plot with eafs:
pdf("eaf_2dx2dc.pdf", width = 10, height = 6)
eaf_sequential_5_50_95_facet(res = res, true = true_25_25_m1_3)
dev.off()

################################################################################

# preparing results for plot: 
res$n = factor(paste0("n = ", res$n), levels = paste0("n = ", unique(res$n)))

# generating plot with runtimes: 
pdf("time_2dx2dc.pdf", width = 9, height = 3)
ggplot(res[done_25_25_m1_3, ], aes(x = algorithm, y = as.numeric(time.running))) + 
  geom_boxplot() +
  facet_wrap(~ n, nrow = 1) +
  labs(x = "Algorithm", y = "Runtime (in seconds)") + scale_y_log10() +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 13))
dev.off()
