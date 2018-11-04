library(graper)
library(tidyverse)
library(SGL)
library(parallel)
library(ipflasso)
library(GRridge)
library(varbvs)
library(glmnet)
library(grpreg)
library(randomForest)
source("../utils_eval.R", chdir = TRUE)

## make grid
# n_base <- 100
# pg_base <- 50
# pi_low_base <- 0.2
# rho_base <- 0
# tau_base <- 1
# 
# n <- seq(20,500,20)
# pg <- seq(10,200,10)
# pi_low <- c(0.001, 0.01, seq(0.05,1,0.05))
# rho <- seq(0, 0.9, 0.1)
# tau <- c(0.01,0.1,1,10,100)
# 
# n_scenario <- length(n) + length(pg) + length(pi_low) +length(rho) +length(tau)
# 
# grid <- expand.grid(n=n, pg=pg, pi_low=pi_low, rho=rho, tau=tau)
# grid_rho <- dplyr::filter(grid, n==n_base, pg==pg_base, pi_low==pi_low_base, tau==tau_base)
# grid_n <- dplyr::filter(grid,  pi_low==pi_low_base, pg==pg_base, rho==rho_base, tau==tau_base)
# grid_pg <- dplyr::filter(grid, n==n_base, rho==rho_base, pi_low==pi_low_base, tau==tau_base)
# grid_pi <- dplyr::filter(grid, n==n_base, pg==pg_base, rho==rho_base, tau==tau_base)
# grid_tau <- dplyr::filter(grid, n==n_base, pg==pg_base, rho==rho_base , pi_low==pi_low_base)
# 
# grid <- rbind(grid_rho,grid_n,grid_pg,grid_pi, grid_tau)
# grid$pi_high <- pmin(1,grid$pi_low*1.5)
# save(grid, file = "data/grid.RData")
# 
load("data/grid.RData")
nrow(grid)

# get setting parameters
param <- grid[as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')),]
n_groups <- 6
true.gamma <- c(0.01, 0.01, 1, 1, 100, 100)
n <- param$n
pg <- rep(param$pg,n_groups)
true.pi <- rep(c(param$pi_low ,param$pi_high), n_groups/2)

nm <- paste0("n",n,"_p", sum(pg),
             "_pil",param$pi_low,
             "_rho",param$rho,
             "_tau", param$tau)
print(nm)

# fit different models
detectCores()

res <- mclapply(1:10, function(it){
  set.seed(9876 + it*777)
  dat <- makeExampleDataWithUnequalGroups(n=n+1000, pg=pg,
                                          gammas = true.gamma,
                                          pis = true.pi,
                                          rho = param$rho,
                                          tau=param$tau)
  Xtrain <- dat$X[1:n,]
  ytrain <- dat$y[1:n]
  # fit on the full dat sets to evaluate estimation of parameters
  fits <- runMethods(Xtrain, ytrain, dat$annot,
                     include_graper_nonfacQ = TRUE, includeIPF = TRUE,
                     n_rep = 1 , th = 0.01, beta0 = dat$beta,
                     standardize=FALSE, includeGRridge = TRUE,
                     includeSparseGroupLasso = TRUE,
                     includeGroupLasso = TRUE, includeRF = TRUE,
                     includeAdaLasso = TRUE, includeVarbvs = TRUE,
                     include_graper_SS_ungrouped = TRUE)
  runtime <- sapply(fits$summaryList, function(l) l$runtime)
  betas <- sapply(fits$summaryList[names(fits$summaryList)!="RandomForest"], function(l) l$beta)
  MSE_beta <- 1/length(beta)*colSums((betas - dat$beta)^2)

  # prediction performance evaluated using cross-validation
  Xtest <- dat$X[(n+1):(n+1000),]
  ytest <- dat$y[(n+1):(n+1000)]
  comp <- evaluateFits(fits, Xtest, ytest)
  RMSE <- sapply(comp$summaryList, function(l) l$RMSE)

  # collect relevant performance meaasures
  data.frame(runtime= t(runtime),
             RMSE = t(RMSE),
             gamma = t(fits$summaryList$graper_SS$pf),
             pi = t(fits$summaryList$graper_SS$sparsity),
             true.gamma = t(true.gamma),
             true.pi=t(true.pi),
             MSE_beta = t(MSE_beta))
}, mc.cores=10) %>% bind_rows()

# for debugging
print(res[[1]])

# save results
outdir <- paste0("all_",as.character(Sys.Date()))
if(!dir.exists(outdir)) dir.create(outdir)
save(res, file=file.path(outdir, paste0(nm, ".RData")))

print(sessionInfo())
