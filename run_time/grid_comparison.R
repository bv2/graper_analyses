library(graper)
library(tidyverse)
library(SGL)
library(parallel)
library(ipflasso)
library(GRridge)
library(varbvs)
library(grpreg)
library(randomForest)
library(glmnet)
source("../utils_eval.R", chdir = TRUE)

# ## make grid
# n_base <- 100
# p_base <- 300
# g_base <- 6
# pi_low_base <- 0.2 
# rho_base <- 0
# tau_base <- 1
#  
# n <- seq(100,2000,200)
# p <- seq(100,5000,200)
# g <- seq(2,10,2)
# 
# grid <- expand.grid(n=n, p=p, g=g, pi_low=pi_low_base, rho = rho_base, tau =tau_base)
# grid_n <- dplyr::filter(grid,  pi_low==pi_low_base, p==p_base, rho==rho_base, tau==tau_base, g == g_base)
# grid_p <- dplyr::filter(grid,  pi_low==pi_low_base, n==n_base, rho==rho_base, tau==tau_base, g == g_base)
# grid_g <- dplyr::filter(grid,  pi_low==pi_low_base, p==p_base, n==n_base, rho==rho_base, tau==tau_base)
# grid <- rbind(grid_n,grid_p,grid_g)
# nrow(grid)
# 
# save(grid, file = "data/grid.RData")

load("data/grid.RData")
nrow(grid)

# get setting parameters
param <- grid[as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')),]
g <- param$g
true.gamma <- rep(1,g)
n <- param$n
pg <- rep(round(param$p/g),g)
true.pi <- rep(param$pi, g)

nm <- paste0("n",n,"_p", sum(pg),
             "_g",param$g)
print(nm)

# fit different models
detectCores()

res <- mclapply(1:50, function(it){
  set.seed(2468 + it*777)
  dat <- makeExampleDataWithUnequalGroups(n=n, pg=pg,
                                          gammas = true.gamma,
                                          pis = true.pi,
                                          rho = param$rho,
                                          tau=param$tau)
  Xtrain <- dat$X
  ytrain <- dat$y
  
  # fit on the full dat sets to evaluate estimation of parameters
  fits <- runMethods(Xtrain, ytrain, dat$annot,
                     include_graper_nonfacQ = TRUE, includeIPF = TRUE,
                     n_rep=1 , th=0.01, beta0 = dat$beta,
                     standardize = FALSE, includeGRridge = TRUE,
                     includeVarbvs = TRUE, includeRF = FALSE,
                     includeGroupLasso = FALSE, includeAdaLasso = FALSE,
                     includeSparseGroupLasso = TRUE)
  runtime <- sapply(fits$summaryList, function(l) l$runtime)

  # collect relevant performance meaasures
  data.frame(runtime= runtime, method = names(runtime), it =it)
}, mc.cores=10) %>% bind_rows()

# for debugging
print(res[[1]])

# save results
outdir <- as.character(Sys.Date())
if(!dir.exists(outdir)) dir.create(outdir)
save(res, file=file.path(outdir, paste0(nm, ".RData")))

print(sessionInfo())
