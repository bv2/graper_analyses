library(graper)
library(Matrix)
library(varbvs)
library(glmnet)
library(ipflasso)
library(SGL)
library(GRridge)
library(grpreg)
library(parallel)
library(randomForest)
source("../utils_eval.R", chdir = TRUE)

# load data as prepared by the scripts in the data directory
load("data/dataGTEx_t5_PCs.RData")
data <- dataGTEx

# compare regression methods in a cross-validation scheme
resultList <- compareMethodsCV(X=data$X, y=data$y, annot=data$annot,
                         seed=9876, family="gaussian", 
                         ncores=10, parallel=TRUE, plot_cv=FALSE,
                         include_graper_nonfacQ = TRUE, max_iter = 5000, nfolds = 10,
                         includeSparseGroupLasso = TRUE,
                         includeGroupLasso = TRUE, includeIPF = TRUE,
                         includeGRridge=TRUE, includeAdaLasso = TRUE,
                         includeVarbvs = TRUE, saveFits =TRUE,
                         th = 0.01, n_rep = 3, standardize = TRUE)

# for debugging:
print(resultList[[1]]) 

# save results
outdir <- as.character(Sys.Date())
if(!dir.exists(outdir)) dir.create(outdir)
save(resultList, file=file.path(outdir, "result_GTEx.RData"))

print(sessionInfo())