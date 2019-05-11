library(graper)
library(Matrix)
library(varbvs)
library(glmnet)
library(SGL)
library(ipflasso)
library(GRridge)
library(grpreg)
library(randomForest)
library(parallel)
source("../utils_eval.R", chdir = TRUE)

# load data from the CLL study as provided in MOFAdata
data("CLL_data", package="MOFAdata")

# use methylation data, gene expression data and drug responses as predictors
CLL_data <- CLL_data[1:3]
CLL_data <- lapply(CLL_data,t)
ngr <- sapply(CLL_data,ncol)
CLL_data <- Reduce(cbind, CLL_data)

#only include patient samples profiles in all three omics
CLL_data <- CLL_data[apply(CLL_data,1, function(p) !any(is.na(p))),]
dim(CLL_data)

# prepare design matrix and response
X <- CLL_data[,!grepl("D_002", colnames(CLL_data))]
y <- rowMeans(CLL_data[,grepl("D_002", colnames(CLL_data))])
annot <- rep(1:3, times = ngr-c(5,0,0)) # group annotations to drugs, meth and RNA


# compare regression methods in a cross-validation scheme
resultList <- compareMethodsCV(X, y, annot, family="gaussian", ncores=10,
                         parallel=TRUE, plot_cv=FALSE, max_iter=5000, nfolds = 10,
                         include_graper_nonfacQ = TRUE,
                         includeSparseGroupLasso = TRUE,
                         includeGroupLasso = TRUE, includeIPF = TRUE,
                         includeGRridge=TRUE, includeAdaLasso = TRUE,
                         includeVarbvs = TRUE, 
                         saveFits=TRUE, seed=9876, 
                         th=0.01, n_rep=3, standardize = TRUE)

# for debugging:
print(resultList[[1]]) 

# save results
outdir <- as.character(Sys.Date())
if(!dir.exists(outdir)) dir.create(outdir)
save(resultList, file=file.path(outdir, "result_CLL.RData"))

print(sessionInfo())