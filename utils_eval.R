
source("utils.R")
source("getFromSummary.R")
source("mygrridge.R") # not required in GRridge version 1.7.1 from GitHub, which can be used alternatively from install_github("markvdwiel/GRridge")

#' @title Run various regression methods
#' @name runMethods
#' @description Function to fit a linear or logistic regression model using several different methods.
#' @details This function fits a linear of logistic regression model to the data using various different methods.
#' graper is always included in its factorized form (both dense (graper_FF) and sparse (graper_SS)). In addition, a graper model is included with all
#' coefficients set to 0 whose posterior inclusion probability (s) is below 50% (graper_SScutoff). As a comparison, ridge regression,
#' Lasso and elastic net as well as a intercept-only model are fitted. Other methods can be included via the respective options, such as GRridge,
#'  non-factorized graper (graper), IPF-Lasso, sparse group Lasso, group Lasso, adaptive Lasso, varbvs, random forest and
#'  graper without group annotations (graper_SS_ungrouped) or without different slab precisions (graper_nogamma).
#' 
#'  The fitted methods can be evaluated on test data using the function \code{\link{evaluateFits}}.
#' @param Xtrain design matrix with samples in rows and features in columns (n x p)
#' @param ytrain response vector of length n
#' @param annot factor of length p indicating group membership of each feature
#' @param intercept  whether to include an intercept into the model
#' @param standardize  whether to standardize the predictors to unit variance.
#'  Note this does not affect GRridge and group Lasso  where standardization is default.
#' @param beta0 true coefficients in the linear model if known,
#'  NULL otherwise (default)
#' @param trueintercept true intercept in the linear model if known,
#'  NULL otherwise (default)
#' @param max_iter maximum number of iterations for graper methods (see also  \code{\link{graper}})
#' @param family likelihood model for the response,
#'  either "gaussian" for linear regression
#' or "binomial" for logistic regression
#' @param calcELB whether to calculate the evidence lower bound (ELB) for graper (see also  \code{\link{graper}})
#' @param freqELB frequency at which the evidence lower bound (ELB) is to be calculated for graper,
#'  i.e. each freqELB-th iteration (see also  \code{\link{graper}})
#' @param th convergence threshold for the evidence lower bound (ELB) in graper (see also  \code{\link{graper}})
#' @param n_rep number of repetitions with different random initializations to be fit in graper (see also  \code{\link{graper}})
#' @param verbose  whether to print out intermediate messages during fitting
#' @param includeGRridge  whether to fit GRridge
#' @param include_graper_nonfacQ  whether to fit graper method with multivariate variational distribution be fitted (can be slow for large data sets)
#' @param includeIPF  whether to fit IPF-Lasso be fitted
#' @param includeSparseGroupLasso  whether to fit sparse group Lasso
#' @param includeGroupLasso  whether to fit group Lasso
#' @param includeAdaLasso  whether to fit Lasso
#' @param includeRF whether to fit random forest
#' @param includeVarbvs  whether to fit varbvs
#' @param include_graper_SS_nogamma  whether to fit graper with same penalty factor but different sparsity levels  per group
#' @param include_graper_SS_ungrouped whether to fit graper without group annotations
#' @param verbose_progress  whether to print out details on the overall progress
#' @importFrom glmnet cv.glmnet
#' @importFrom varbvs varbvs
#' @importFrom randomForest randomForest
#' @importFrom grpreg cv.grpreg
#' @importFrom SGL cvSGL
#' @importFrom GRridge CreatePartition grridge
#' @importFrom ipflasso cvr2.ipflasso
#' @importFrom stats coef
#' @return a list containing
#' \describe{
#' \item{summaryList}{list with the fitted models for each method that was included in the comparison, such as
#' \itemize{
#' \item graper_SS (sparse graper model with factorized variational distribution)
#' \item graper_FF (dense graper model with factorized variational distribution)
#' \item graper (dense graper model with non-factorized variational distribution)
#' \item graper_SScutoff (as graper_SS where coefficients with posterior inclusion probabilities below 0.5 are set to zero)
#' \item graper_SS_nogamma (as graper_SS with a common slab precision across groups)
#' \item graper_SS_ungrouped (as graper_SS without group annotation)
#' \item Ridge (ridge regression)
#' \item Lasso (Lasso regression)
#' \item ElasticNet (elastic net regression)
#' \item varbvs (varbvs regression)
#' \item GroupLasso (group Lasso regression)
#' \item SparseGroupLasso (sparse group Lasso regression)
#' \item SparseGroupLasso (sparse group Lasso regression)
#' \item GRridge (GRridge regression)
#' \item NullModel (regression with the intercept term only)
#' \item TrueModel (regression with the true model coefficients (if known))
#' \item IPFLasso (IPF-Lasso regression)
#' \item adaptiveLasso (adaptive Lasso regression)
#' }
#' Each entry contains the
#' runtime, estimated penalty factors, coefficients, intercepts, sparsity level and the full output returned by the methods' call (out).
#' }
#' \item{groupnames}{groupnames from \code{annot}}
#' \item{varnames}{predictor names}
#' \item{family}{likelihoods model used}
#' \item{n}{number of samples}
#' \item{p}{number of features}
#' \item{G}{number of groups}
#' \item{annot}{annotation of features to groups as specified when calling \code{\link{runMethods}}}
#' }
#' @export
#' @examples
#' dat <- makeExampleData()
#' allFits <- runMethods(Xtrain=dat$X, ytrain=dat$y, annot=dat$annot)

runMethods <- function(Xtrain, ytrain, annot,
                       family = "gaussian", intercept = TRUE, standardize = TRUE,
                       beta0 = NULL, trueintercept = NULL,
                       max_iter = 5000, freqELB = 10, calcELB = TRUE, th = 0.01,
                       n_rep=1, verbose = FALSE, verbose_progress = TRUE,
                       includeGRridge = FALSE, include_graper_nonfacQ = FALSE,
                       includeSparseGroupLasso = FALSE, includeIPF = FALSE,
                       includeGroupLasso = FALSE, includeRF = FALSE,
                       includeAdaLasso = FALSE, includeVarbvs = FALSE,
                       include_graper_SS_nogamma = FALSE,
                       include_graper_SS_ungrouped = FALSE) {

  # sanity checks
  stopifnot(nrow(Xtrain) == length(ytrain))
  if (!is.null(beta0)) stopifnot(ncol(Xtrain) == length(beta0))

  # turn annot to factor
  annot <- as.factor(annot)

  # extract important parameters and names
  p <- ncol(Xtrain)
  n <- nrow(Xtrain)
  G <- length(unique(annot))
  groupnames <- as.character(annot)
  varnames <- colnames(Xtrain)
  if (length(varnames) == 0)
    varnames <- factor(paste("Feature", seq_len(ncol(Xtrain)), sep = ""))

  #### RUN DIFFERENT METHODS ####
  summaryList <- list()

  # graper: multivariate variational distribution, normal prior (dense)
  if (include_graper_nonfacQ) {
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting graper model...")
    graper <- graper(Xtrain, ytrain, annot = annot, factoriseQ = FALSE, spikeslab = FALSE,
                       max_iter = max_iter, intercept = intercept,
                       verbose = verbose, freqELB = freqELB, calcELB = calcELB,
                       family = family, th = th, standardize=standardize, n_rep=n_rep)
    timeNF <- difftime(Sys.time(), tmp, units = "secs")
    graper_summary <- list()
    graper_summary$runtime <- as.numeric(timeNF)
    graper_summary$pf <- as.numeric(graper$EW_gamma)
    graper_summary$beta <- graper$EW_beta
    graper_summary$intercept <- graper$intercept
    graper_summary$sparsity <- rep(1, G)  #dense - no sparsity per groups
    graper_summary$out <- graper
    rm(graper)
    summaryList$graper <- graper_summary
  }

  # graper_FF : fully factorized variational distribution, normal prior (dense)
  tmp <- Sys.time()
  if(verbose_progress) message(" ### Fitting graper_FF model...")
  graper_FF <- graper(Xtrain, ytrain, annot = annot, factoriseQ = TRUE, spikeslab = FALSE,
                        max_iter = max_iter, intercept = intercept,
                        verbose = verbose, freqELB = freqELB, calcELB = calcELB,
                        family = family, th = th, standardize=standardize, n_rep=n_rep)
  timeFF <- difftime(Sys.time(), tmp, units = "secs")
  graper_FF_summary <- list()
  graper_FF_summary$runtime <- as.numeric(timeFF)
  graper_FF_summary$pf <- as.numeric(graper_FF$EW_gamma)
  graper_FF_summary$beta <- graper_FF$EW_beta
  graper_FF_summary$intercept <- graper_FF$intercept
  graper_FF_summary$sparsity <- rep(1, G)  #dense - no sparsity per groups
  graper_FF_summary$out <- graper_FF
  rm(graper_FF)
  summaryList$graper_FF <- graper_FF_summary

  # graper_SS: fully factorized variational distribution, spike and slab prior (sparse)
  tmp <- Sys.time()
  if(verbose_progress) message(" ### Fitting graper_SS model...")
  graper_SS <- graper(Xtrain, ytrain, annot = annot, factoriseQ = TRUE, spikeslab = TRUE,
                        max_iter = max_iter, intercept = intercept,
                        verbose = verbose, freqELB = freqELB, calcELB = calcELB,
                        th = th, family = family, standardize=standardize, n_rep=n_rep)
  timeSS <- difftime(Sys.time(), tmp, units = "secs")
  graper_SS_summary <- list()
  graper_SS_summary$runtime <- as.numeric(timeSS)
  graper_SS_summary$pf <- as.numeric(graper_SS$EW_gamma)
  graper_SS_summary$beta <- graper_SS$EW_beta
  graper_SS_summary$intercept <- graper_SS$intercept
  graper_SS_summary$sparsity <- graper_SS$EW_pi
  graper_SS_summary$out <- graper_SS
  summaryList$graper_SS <- graper_SS_summary

  # set factors with a low posteriori inclusion probability to zero
  graper_SScutoff_summary <- list()
  graper_SScutoff_summary$runtime <- as.numeric(timeSS)
  graper_SScutoff_summary$pf <- as.numeric(graper_SS$EW_gamma)
  graper_SScutoff_summary$beta <- ifelse(graper_SS$EW_s < 0.5, 0, graper_SS$EW_beta)
  graper_SScutoff_summary$intercept <- graper_SS$intercept
  graper_SScutoff_summary$sparsity <- graper_SS$EW_pi
  graper_SScutoff_summary$out <- NULL
  summaryList$graper_SScutoff <- graper_SScutoff_summary
  rm(graper_SS)

  # graper_SS_nogamma: fully factorized variational distribution,
  # spike and slab prior (sparse) without different slab precisions
  if(include_graper_SS_nogamma){
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting graper_SS model without gamma...")
    graper_SS_nogamma <- graper(Xtrain, ytrain, annot = annot, factoriseQ = TRUE,
                                  spikeslab = TRUE, max_iter = max_iter, intercept = intercept,
                                  verbose = verbose, freqELB = freqELB, calcELB = calcELB, th = th,
                                  family = family,  nogamma = TRUE, standardize=standardize, n_rep=n_rep)
    timeSS_nogamma <- difftime(Sys.time(), tmp, units = "secs")
    graper_SS_nogamma_summary <- list()
    graper_SS_nogamma_summary$runtime <- as.numeric(timeSS_nogamma)
    graper_SS_nogamma_summary$pf <- rep(as.numeric(graper_SS_nogamma$EW_gamma), G)
    graper_SS_nogamma_summary$beta <- graper_SS_nogamma$EW_beta
    graper_SS_nogamma_summary$intercept <- graper_SS_nogamma$intercept
    graper_SS_nogamma_summary$sparsity <- graper_SS_nogamma$EW_pi
    graper_SS_nogamma_summary$out <- graper_SS_nogamma
    summaryList$graper_SS_nogamma <- graper_SS_nogamma_summary
    rm(graper_SS_nogamma)
  }

  # graper_SS_ungrouped: fully factorized variational distribution, spike and slab prior
  # (sparse) without group annotations
  if(include_graper_SS_ungrouped){
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting graper_SS model without group annotations...")
    graper_SS_ungrouped <- graper(Xtrain, ytrain, annot = rep(1,ncol(Xtrain)),
                                    factoriseQ = TRUE, spikeslab = TRUE,
                                    max_iter = max_iter, intercept = intercept,
                                    verbose = verbose, freqELB = freqELB,
                                    calcELB = calcELB, th = th,
                                    family = family,  nogamma = TRUE,
                                    standardize=standardize, n_rep=n_rep)
    timeSS_ungrouped <- difftime(Sys.time(), tmp, units = "secs")
    graper_SS_ungrouped_summary <- list()
    graper_SS_ungrouped_summary$runtime <- as.numeric(timeSS_ungrouped)
    graper_SS_ungrouped_summary$pf <- rep(as.numeric(graper_SS_ungrouped$EW_gamma), G)
    graper_SS_ungrouped_summary$beta <- graper_SS_ungrouped$EW_beta
    graper_SS_ungrouped_summary$intercept <- graper_SS_ungrouped$intercept
    graper_SS_ungrouped_summary$sparsity <- rep(graper_SS_ungrouped$EW_pi, G)
    graper_SS_ungrouped_summary$out <- graper_SS_ungrouped
    summaryList$graper_SS_ungrouped <- graper_SS_ungrouped_summary
    rm(graper_SS_ungrouped)
  }

  # ridge regression
  tmp <- Sys.time()
  if(verbose_progress) message(" ### Fitting ridge regression...")
  RidgeFit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = 0, intercept = intercept,
                                standardize = standardize, family = family)
  tmp <- difftime(Sys.time(), tmp, units = "secs")
  beta_ridge <- as.vector(stats::coef(RidgeFit, RidgeFit$lambda.min))[-1]
  Ridge_summary <- list()
  Ridge_summary$runtime <- as.numeric(tmp)
  Ridge_summary$pf <- rep(RidgeFit$lambda.min, G)
  Ridge_summary$beta <- beta_ridge
  if(intercept) {
    Ridge_summary$intercept <- as.vector(stats::coef(RidgeFit, RidgeFit$lambda.min))[1]
  }
  Ridge_summary$sparsity <- rep(1, G)  #dense - no sparsity per groups
  Ridge_summary$out <- RidgeFit
  rm(RidgeFit, beta_ridge)
  summaryList$Ridge <- Ridge_summary

  # Lasso
  tmp <- Sys.time()
  if(verbose_progress) message(" ### Fitting Lasso...")
  LassoFit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = 1, intercept = intercept,
                                standardize = standardize, family = family)
  tmp <- difftime(Sys.time(), tmp, units = "secs")
  beta_lasso <- as.vector(stats::coef(LassoFit, LassoFit$lambda.min))[-1]
  Lasso_summary <- list()
  Lasso_summary$runtime <- as.numeric(tmp)
  Lasso_summary$pf <- rep(LassoFit$lambda.min, G)
  Lasso_summary$beta <- beta_lasso
  if(intercept) {
    Lasso_summary$intercept <- as.vector(stats::coef(LassoFit, LassoFit$lambda.min))[1]
  }
  Lasso_summary$sparsity <- vapply(unique(annot),
                                   function(gr) sum(beta_lasso[annot == gr] != 0)/sum(annot == gr),
                                   numeric(1))
  Lasso_summary$out <- LassoFit
  rm(LassoFit, beta_lasso)
  summaryList$Lasso <- Lasso_summary

  # elastic net
  tmp <- Sys.time()
  if(verbose_progress) message(" ### Fitting elastic net...")
  ENFit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = 0.2, intercept = intercept,
                             standardize = standardize, family = family)
  tmp <- difftime(Sys.time(), tmp, units = "secs")
  beta_EN <- as.vector(stats::coef(ENFit, ENFit$lambda.min))[-1]
  ElasticNet_summary <- list()
  ElasticNet_summary$runtime <- as.numeric(tmp)
  ElasticNet_summary$pf <- rep(ENFit$lambda.min, G)
  ElasticNet_summary$beta <- beta_EN
  if(intercept) ElasticNet_summary$intercept <- as.vector(stats::coef(ENFit, ENFit$lambda.min))[1]
  ElasticNet_summary$sparsity <- vapply(unique(annot),
                                        function(gr) sum(beta_EN[annot == gr] != 0)/sum(annot == gr),
                                        numeric(1))
  ElasticNet_summary$out <- ENFit
  rm(ENFit, beta_EN)
  summaryList$ElasticNet <- ElasticNet_summary

  # varbvs
  if(includeVarbvs){
    if(!intercept | standardize) warning("varbvs always includes an intercept and does not standardize")
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting varbvs...")
    varbvsFit <- varbvs::varbvs(X=Xtrain, Z=NULL, y=ytrain, family = family)
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    #if (intercept)
    beta_varbvs <- apply(varbvsFit$alpha * varbvsFit$mu,1, function(b) sum(b * varbvsFit$w))
    varbvs_summary <- list()
    varbvs_summary$runtime <- as.numeric(tmp)
    varbvs_summary$pf <- rep(sum(varbvsFit$sa * varbvsFit$w),G)
    varbvs_summary$beta <- beta_varbvs
    varbvs_summary$intercept <- sum(varbvsFit$mu.cov * varbvsFit$w)
    varbvs_summary$sparsity <- vapply(unique(annot),
                                      function(gr) mean(varbvsFit$pip[annot == gr]),
                                      numeric(1))
    varbvs_summary$out <- varbvsFit
    rm(varbvsFit, beta_varbvs)
    summaryList$varbvs <- varbvs_summary
  }

  # Random Forest
  if (includeRF) {
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting random forest...")
    if(family=="gaussian") {
      rf.out <- randomForest::randomForest(x = Xtrain, y = ytrain)
    } else if(family=="binomial") {
      rf.out <- randomForest::randomForest(x = Xtrain, y = as.factor(ytrain))
    }
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    RandomForest_summary <- list()
    RandomForest_summary$runtime <- as.numeric(tmp)
    RandomForest_summary$pf <- NULL
    RandomForest_summary$beta <- NULL
    RandomForest_summary$intercept <- NULL
    RandomForest_summary$sparsity <- NULL
    RandomForest_summary$out <- rf.out
    rm(rf.out)
    summaryList$RandomForest <- RandomForest_summary
  }

  # group lasso
  if (includeGroupLasso) {
    if (!standardize | !intercept)
      warning("Group Lasso implementation always includes an intercept and standardizes the predictors.")

    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting group Lasso...")
    GrpLassoFit <- try(grpreg::cv.grpreg(Xtrain, ytrain, group = as.factor(annot),
                                         penalty = "grLasso", family = family))
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    if (class(GrpLassoFit) == "try-error") {
      warning("Group Lasso encountered errors, not included in the comparison!")
    } else {
      beta_GrpLasso <- as.vector(stats::coef(GrpLassoFit, GrpLassoFit$lambda.min))[-1]
      GroupLasso_summary <- list()
      GroupLasso_summary$runtime <- as.numeric(tmp)
      GroupLasso_summary$pf <- rep(GrpLassoFit$lambda.min, G)
      GroupLasso_summary$beta <- beta_GrpLasso
      GroupLasso_summary$intercept <- as.vector(stats::coef(GrpLassoFit, GrpLassoFit$lambda.min))[1]
      GroupLasso_summary$sparsity <- vapply(unique(annot),
                                            function(gr) {
                                              sum(beta_GrpLasso[annot == gr] != 0)/sum(annot == gr)
                                              }, numeric(1))
      GroupLasso_summary$out <- GrpLassoFit
      rm(GrpLassoFit, beta_GrpLasso)
      summaryList$GroupLasso <- GroupLasso_summary
    }
  }

  # sparse group lasso
  if (includeSparseGroupLasso) {
    if(!intercept & standardize){
      warning("Sparse group Lasso includes an intercept when standardize = TRUE")
    }
    if(intercept & !standardize & family == "gaussian") {
      # cvSGL uses the standardize option for both including an intercept and scaling
      # here we rmove the intercept beforehand if standardize is FALSE
        meansX <- apply(Xtrain, 2, mean)
        XtrainSGL <- t(t(Xtrain) - meansX)
        meany <- mean(ytrain)
        ytrainSGL <- ytrain - meany
      } else {
        ytrainSGL <- ytrain
        XtrainSGL <- Xtrain
      }
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting sparse group Lasso...")
    SpGrpLassoFit <- try(SGL::cvSGL(list(x=XtrainSGL, y=ytrainSGL),
                                    index = as.factor(annot), standardize = standardize,
                                    type = ifelse(family=="gaussian", "linear", "logit")))
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    if (class(SpGrpLassoFit) == "try-error") {
      warning("Sparse Group Lasso encountered errors, not included in the comparison!")
    } else {
      beta_SpGrpLasso <- SpGrpLassoFit$fit$beta[, which.min(SpGrpLassoFit$lldiff)]
      SpGroupLasso_summary <- list()
      SpGroupLasso_summary$runtime <- as.numeric(tmp)
      SpGroupLasso_summary$pf <- rep(SpGrpLassoFit$lambdas[which.min(SpGrpLassoFit$lldiff)], G)
      SpGroupLasso_summary$beta <- beta_SpGrpLasso
      # cvSGL implementation does not return coefficients on origingal scale when standardizing but for X/v: b = b_SGL/v
      if(standardize) SpGroupLasso_summary$beta <- 1/apply(XtrainSGL,2,function(x) sqrt(sum((x-mean(x))^2))) * SpGroupLasso_summary$beta
      if(intercept & !standardize){
        if(family == "gaussian"){
        intercept4SGL <- meany - sum(meansX*beta_SpGrpLasso)
        } else intercept4SGL <- 0  # not implemented for unstd. logistic regression in SGL
      }
      if(intercept) SpGroupLasso_summary$intercept <- ifelse(standardize, SpGrpLassoFit$fit$intercept,intercept4SGL)
      # adapt intercept to X offset (SGL returns for centered (X-m)/v) b0 = b_SGL0 - m/v * b_SGL = b_SGL0 - m*b
      if(standardize & intercept) {
        SpGroupLasso_summary$intercept <-  SpGroupLasso_summary$intercept - sum(apply(XtrainSGL,2,mean) * SpGroupLasso_summary$beta)
      }
        
      SpGroupLasso_summary$sparsity <- vapply(unique(annot), function(gr){
        sum(beta_SpGrpLasso[annot == gr] != 0)/sum(annot == gr)
      }, numeric(1))
      SpGroupLasso_summary$out <- SpGrpLassoFit
      rm(SpGrpLassoFit, beta_SpGrpLasso)
      summaryList$SparseGroupLasso <- SpGroupLasso_summary
    }
  }

  # GRridge
  if (includeGRridge) {
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting GRridge...")
    partition <- GRridge::CreatePartition(as.factor(annot))
    if (intercept) {
      #notes itself which type of response
      if(packageVersion("GRridge") < "1.7.1") {
        MessagesGR <- utils::capture.output(GRfit <- try(mygrridge(highdimdata=t(Xtrain), response=as.numeric(ytrain),
                                                                        partitions=list(partition), unpenal = ~1,
                                                                        standardizeX = standardize)))
      } else {
        MessagesGR <- utils::capture.output(GRfit <- try(GRridge::grridge(highdimdata=t(Xtrain), response=as.numeric(ytrain),
                                                                   partitions=list(partition), unpenal = ~1,
                                                                   standardizeX = standardize)))
      }
    } else {
      #notes itself which type of response
      if(packageVersion("GRridge") < "1.7.1") {
      MessagesGR <- utils::capture.output(GRfit <- try(mygrridge(highdimdata=t(Xtrain), response=as.numeric(ytrain),
                                                                        partitions=list(partition), unpenal = ~0,
                                                                        standardizeX = standardize)))
      } else {
        MessagesGR <- utils::capture.output(GRfit <- try(GRridge::grridge(highdimdata=t(Xtrain), response=as.numeric(ytrain),
                                                                   partitions=list(partition), unpenal = ~0,
                                                                   standardizeX = standardize)))
      }
    }
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    if (verbose)
      MessagesGR
    if (class(GRfit) == "try-error") {
      warning("GRridge encountered errors, not included in the comparison!")
    } else {
      GRridge_summary <- list()
      GRridge_summary$runtime <- as.numeric(tmp)
      GRridge_summary$pf <- as.numeric(GRfit$lambdamults[[1]])
      GRridge_summary$beta <- GRfit$betas  
      if(standardize) GRridge_summary$beta <- GRridge_summary$beta * (1/apply(Xtrain,2,sd)) # return on original scale
      if(intercept) {
        if(family=="gaussian") GRridge_summary$intercept <- coef(GRfit$predobj$GroupRegul)[1]
        else GRridge_summary$intercept <- GRfit$predobj$GroupRegul@unpenalized
      }
      # adapt intercept to X offset (grridge returns for centered (X-m)/v) b0 = b_gr - m/v * b_gr = b0_gr - m*b
      if(standardize & intercept) {
        GRridge_summary$intercept <-  GRridge_summary$intercept - sum(apply(Xtrain,2,mean) * GRridge_summary$beta)
      }
      GRridge_summary$sparsity <- rep(1, G)
      GRridge_summary$out <- GRfit
      rm(GRfit)
      summaryList$GRridge <- GRridge_summary
    }
  }

  # zero model
  tmp <- Sys.time()
  if(verbose_progress) message(" ### Fitting zero model...")
  if (intercept) intercept_zeromodel <- mean(ytrain) else intercept_zeromodel <- 0
  tmp <- difftime(Sys.time(), tmp, units = "secs")
  NullModel_summary <- list()
  NullModel_summary$runtime <- as.numeric(tmp)
  NullModel_summary$pf <- rep(Inf, G)
  NullModel_summary$beta <- rep(0, p)
  NullModel_summary$intercept <- intercept_zeromodel
  NullModel_summary$sparsity <- rep(0, G)
  NullModel_summary$out <- NULL
  summaryList$NullModel <- NullModel_summary

  # True Model
  if (!is.null(beta0)) {
    if(verbose_progress) message(" ### Including true model...")
    TrueModel_summary <- list()
    TrueModel_summary$runtime <- 0
    TrueModel_summary$pf <- 1/vapply(unique(annot),
                                     function(gr) mean((beta0[annot == gr])^2),
                                     numeric(1))
    TrueModel_summary$beta <- beta0
    TrueModel_summary$intercept <- trueintercept
    TrueModel_summary$sparsity <- vapply(unique(annot),
                                         function(gr) sum(beta0[annot == gr] != 0)/sum(annot == gr),
                                         numeric(1))
    TrueModel_summary$out <- NULL
    summaryList$TrueModel <- TrueModel_summary
  }

  # IPF -Lasso
  if (includeIPF) {
    if(!intercept){
      warning("IPF-Lasso always includes an intercept.")
    }
    # penalty factors to consider for cross-validation
    # (unclear how to choose, take a very rough grid here to make it applicable to larger number of groups)
    lambda_1d <- c(0.1, 0.5, 1, 2, 10)
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting IPF-Lasso...")
    pfgrid <- expand.grid(rep(list(lambda_1d), G))
    pflist <- lapply(seq_len(nrow(pfgrid)), function(i) pfgrid[i, ])
    type.measure <- ifelse(family == "gaussian", "mse", "class")
    ipf.out <- try(ipflasso::cvr2.ipflasso(Xtrain, ytrain, alpha = 1, standardize = standardize,
                                           family = family, type.measure = type.measure,
                                           blocks = lapply(unique(annot), function(gr) which(annot == gr)),
                                           pflist = pflist, nfolds = 10, ncv = 3))
    #using same cv parameter as standard glmnet leads to errors, needs to be ncv>1
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    if (class(ipf.out) == "try-error") {
      warning("IPF-Lasso encountered errors, not included in the comparison!")
    } else {
      IPFLasso_summary <- list()
      IPFLasso_summary$runtime <- as.numeric(tmp)
      IPFLasso_summary$pf <- pflist[[ipf.out$ind.bestpf]]
      IPFLasso_summary$beta <- ipf.out$coeff[-1,ipf.out$ind.bestlambda]
      IPFLasso_summary$intercept <- ipf.out$coeff[1,ipf.out$ind.bestlambda]
      IPFLasso_summary$sparsity <-  vapply(unique(annot), function(gr) {
        sum(IPFLasso_summary$beta[annot == gr] != 0)/sum(annot == gr)
        }, numeric(1))
      IPFLasso_summary$out <- ipf.out
      rm(ipf.out)
      summaryList$IPFLasso <- IPFLasso_summary
    }

  }

  # Adaptive Lasso
  if (includeAdaLasso) {
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting adaptive Lasso...")
    ## Ridge Regression to create the Adaptive Weights Vector
    RidgeFit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = 0, intercept = intercept,
                                  standardize = standardize, family = family)
    wRidge <- pmin(1/abs((stats::coef(RidgeFit, s = RidgeFit$lambda.min))), 1e+300)
    wRidge <- wRidge[-1]
    adaLassoFit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = 1, intercept = intercept,
                                     standardize = standardize, family = family,
                                     penalty.factor = wRidge)
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    beta_adalasso <- as.vector(stats::coef(adaLassoFit, adaLassoFit$lambda.min))[-1]
    adaLasso_summary <- list()
    adaLasso_summary$runtime <- as.numeric(tmp)
    adaLasso_summary$pf <- vapply(unique(annot), function(gr){
      mean(adaLassoFit$lambda.min * wRidge[annot == gr])
    }, numeric(1))
    adaLasso_summary$beta <- beta_adalasso
    if(intercept) {
      adaLasso_summary$intercept <- as.vector(stats::coef(adaLassoFit, adaLassoFit$lambda.min))[1]
    }
    adaLasso_summary$sparsity <- vapply(unique(annot), function(gr) {
      sum(beta_adalasso[annot == gr] != 0)/sum(annot == gr)
      }, numeric(1))
    adaLasso_summary$out <- adaLassoFit
    rm(adaLassoFit, beta_adalasso)
    summaryList$adaptiveLasso <- adaLasso_summary
  }

  return(list(summaryList = summaryList,
              groupnames = groupnames,
              varnames = varnames,
              family = family,
              n = n, p = p,
              G = G, annot = annot))
}


#' @title Evaluate fits from various regression methods
#' @name evaluateFits
#' @description Function to evaluate results from \code{\link{runMethods}} on test data.
#' @param allFits List as produced by \code{\link{runMethods}}
#' @param Xtest Design matrix of size n' x p
#' (same feature structure as used in \code{\link{runMethods}})
#' @param ytest Response vector of size n'
#' @return a list as produced by \code{\link{runMethods}} with
#'  additional predicition performance slots in the summaryList, i.e. the list contains
#' \describe{
#' \item{summaryList}{list with the fitted models for each method that was included in the comparison, such as
#' \itemize{
#' \item graper_SS (sparse graper model with factorized variational distribution)
#' \item graper_FF (dense graper model with factorized variational distribution)
#' \item graper (dense graper model with non-factorized variational distribution)
#' \item graper_SScutoff (as graper_SS where coefficients with posterior inclusion probabilities below 0.5 are set to zero)
#' \item graper_SS_nogamma (as graper_SS with a common slab precision across groups)
#' \item graper_SS_ungrouped (as graper_SS without group annotation)
#' \item Ridge (ridge regression)
#' \item Lasso (Lasso regression)
#' \item ElasticNet (elastic net regression)
#' \item varbvs (varbvs regression)
#' \item GroupLasso (group Lasso regression)
#' \item SparseGroupLasso (sparse group Lasso regression)
#' \item SparseGroupLasso (sparse group Lasso regression)
#' \item GRridge (GRridge regression)
#' \item NullModel (regression with the intercept term only)
#' \item TrueModel (regression with the true model coefficients (if known))
#' \item IPFLasso (IPF-Lasso regression)
#' \item adaptiveLasso (adaptive Lasso regression)
#' }
#'Each entry contains the
#' runtime, estimated penalty factors, coefficients, intercepts, sparsity level and the full output returned by the methods' call (out) as in
#' the output produced by \code{\link{runMethods}}.
#'
#' In addition, it contains performance statistics such as root mean squared error (RMSE) and
#' feautre selection properties (if \code{beta0} was specified)}
#' \item{groupnames}{groupnames from \code{annot}}
#' \item{varnames}{predictor names}
#' \item{family}{likelihoods model used}
#' \item{n}{number of samples}
#' \item{p}{number of features}
#' \item{G}{number of groups}
#' \item{annot}{annotation of features to groups as specified when calling \code{\link{runMethods}}}
#' }
#' @export
#' @details This function can be used to evaluate the fits of various method as produced by \code{\link{runMethods}} on a test dataset.
#' If the true coefficients of the model are known they can be specified via \code{beta0} and \code{trueintercept}.
#' Then, additionally the error on the estimates is evaluates as well as the feature selection performance.
#' Note that for graper the selected features are determined by the posterior inclusion probabilities with
#' a feature being called active for s>0.5 (The method graper_cutoff sets inactive with s<=0.5 coefficients to exactly zero).
#' @importFrom stats predict
#' @examples
#' dat <- makeExampleData()
#' ntrain <- dat$n/2
#' fit <- runMethods(dat$X[seq_len(ntrain),],
#'  dat$y[seq_len(ntrain)], dat$annot,
#'   beta0 = dat$beta)
#' evalFit <- evaluateFits(fit, dat$X[seq_len(ntrain) + dat$n/2,],
#'  dat$y[seq_len(ntrain)+ dat$n/2])

evaluateFits <- function(allFits, Xtest, ytest) {

  stopifnot(nrow(Xtest) == length(ytest))
  stopifnot(ncol(Xtest) == length(allFits$summaryList[[1]]$beta))

  family <- allFits$family
  summaryList <- allFits$summaryList
  ytest <- as.vector(ytest)
  ntest <- length(ytest)
  beta0 <- summaryList$TrueModel$beta
  intercept0 <- summaryList$TrueModel$intercept

  ###### Prediction Performance

  # For gaussian family calculate RMSE as measure of prediciton performance
  if (family == "gaussian") {
    summaryList <- lapply(summaryList, function(summary) {
      beta <- summary$beta
      intercept <- summary$intercept
      if (!is.null(beta)) {
        # for cases with linear coefficients
        RMSE <- evaluateModel(beta, intercept = intercept, Xtest, ytest,
                              beta0 = beta0, family = "gaussian")$RMSE_test
        summary$RMSE <- RMSE
      }
      summary
    })
    if ("RandomForest" %in% names(summaryList))
      summaryList$RandomForest$RMSE <- sqrt(1/length(ytest) * sum((predict(summaryList$RandomForest$out, Xtest) - ytest)^2))
    if ("varbvs" %in% names(summaryList))
      summaryList$varbvs$RMSE <- sqrt(1/length(ytest) * sum((predict(summaryList$varbvs$out, Xtest) - ytest)^2))
  }

  # For binomial family calculate AUC and Brier Score
  if (family == "binomial") {
    summaryList <- lapply(summaryList, function(summary) {
      beta <- summary$beta
      intercept <- summary$intercept
      if (!is.null(beta)) {
        # for cases with linear coeeficients
        eval.out <- evaluateModel(beta, intercept = intercept, Xtest,
                                  ytest, beta0 = beta0, family = "binomial")
        summary$AUC <- eval.out$AUC
        summary$BS <- eval.out$BrierScore
      }
      # TODO add varbvs and RandomForest
      summary
    })
  }

  ###### Feature Recovery
  if (!is.null(beta0)) {
    summaryList <- lapply(summaryList, function(summary) {
      beta <- summary$beta
      intercept <- summary$intercept
      if (!is.null(beta)) {
        # for cases without linear coeeficients e.g. Random Forest
        eval.out <- evaluateModel(beta, intercept = intercept, Xtest,
                                  ytest, beta0 = beta0, family = family)
        summary$FPR <- eval.out$FPR
        summary$FNR <- eval.out$FNR
        summary$precision <- eval.out$precision
        summary$recall <- eval.out$recall
        summary$F1score <- eval.out$F1score
      }
      summary
    })
  }


  ###### Error on estimate

  # l1-Error in estimation of coeffcients
  if (!is.null(beta0)) {
    summaryList <- lapply(summaryList, function(summary) {
      beta <- summary$beta
      if (!is.null(beta))
        summary$l1error_beta <- sum(abs(beta - beta0))
      summary
    })
  }

  # l1-Error in estimation of intercept
  if (!is.null(intercept0)) {
    summaryList <- lapply(summaryList, function(summary) {
      intercept <- summary$intercept
      if (!is.null(intercept))
        summary$l1error_intercept <- sum(abs(intercept - intercept0))
      summary
    })
  }

  allFits$summaryList <- summaryList

  return(allFits)
}

#' @title Compare various regression method via cross-validation
#' @name compareMethodsCV
#' @description  Function to fit a regression  model using several different methods
#' and evaluate them in a cross-validated fashion in terms of prediction and estimation performance.
#' @param X design matrix with samples in rows and features in columns (n x p)
#' @param y response vector of length n
#' @param annot factor of length p indicating group membership of each feature
#' @param family likelihood model for the response,
#'  either "gaussian" for linear regression
#' or "binomial" for logistic regression
#' @param nfolds number of folds for evaluation
#' @param ncores number of cores to use
#' @param plot_cv whether to produce summary plots from evaluation
#' @param seed optional seed for the choice of folds
#' @param parallel whether to run cross-validation in parallel
#' @param saveFits whether to save the fit of each fold
#' @param ... other parameters that can be passed to \code{\link{runMethods}}
#' @return List with one entry per fold containing a list with
#' \describe{
#' \item{FPR}{False positive rate (requires \code{beta0} to be specified)}
#' \item{FNR}{False negative rate (requires \code{beta0} to be specified)}
#' \item{RMSE}{Root mean squared error on the left-out fold}
#' \item{pf_mat}{matrix with learnt penalty factors per group (rows) and method (column)}
#' \item{beta_mat}{matrix with estimated coefficients per feature (rows) and method (column)}
#' \item{intercepts}{estimated intercepts per method}
#' \item{sparsity_mat}{matrix with learnt sparsity levels per group (rows) and method (column) (0=sparse, 1=dense)}
#' \item{annot}{annotation of features to groups as specified when calling \code{\link{compareMethodsCV}}}
#' \item{runtime}{vector of runtimes for the different methods}
#' \item{l1error_intercept}{absolute error on the estimated intercept  (requires \code{trueintercept} to be specified)}
#' \item{l1error_beta}{absolute error on the estimated coefficients (requires \code{beta0} to be specified)}
#' }
#' @import ggplot2
#' @import parallel
#' @details This function can be used to test various method for regression on a dataset in a cross-validated fashion.
#' It fits the methods on all except one fold and evaluates their prediction performance on the remaining fold.
#' See the documentation of \code{\link{runMethods}} on the methods that can be included in the comparison.

#' If the true coefficients of the model are known they can be specified via \code{beta0} and \code{trueintercept}.
#' Then, additionally the error on the estimates is evaluates as well as the feature selection performance.
#' Note that for graper the selected features are determined by the posterior inclusion probabilities with
#' a feature being called active for s>0.5 (The method graper_cutoff sets inactive with s<=0.5 coefficients to exactly zero).
#' @export
#' @examples
#' dat <- makeExampleData()
#' cv.out <- compareMethodsCV(dat$X, dat$y, dat$annot, nfolds=3)

compareMethodsCV <- function(X, y, annot, family="gaussian",
                       ncores=1, nfolds=10, plot_cv=FALSE,
                       seed=NULL, parallel=FALSE, saveFits=FALSE, ...){

  # split observations into folds
  if(!is.null(seed)) set.seed(seed)
  foldid <- sample(rep(seq(nfolds), length=nrow(X)))

  # function for each fold in cross-validation: train method on all
  # excpet one fold, evaluate on the hold-out fold
  runPerFold <- function(foldidx){
    # split in train and test data
    use4test <- foldid==foldidx
    ytrain <- y[ !use4test]
    ytest <- y[use4test]
    Xtrain <- X[ !use4test,]
    Xtest <- X[use4test,]

    #fit models
    AllFits <- runMethods(Xtrain = as.matrix(Xtrain),
                          ytrain =as.vector(ytrain),
                          annot = annot, family = family,
                          ...)

    # save the fit
    if(saveFits) save(AllFits, file=paste0("AllFits",foldidx,".RData"))

    # extract relevant parameter from the model
    pf_mat <- getPenaltyFactors(AllFits)
    sparsity_mat <- getSparsityLevel(AllFits)
    beta_mat <- getCoefficients(AllFits)
    intercepts <- getIntercept(AllFits)

    # evaluate prediciton performance
    AllFits <- evaluateFits(AllFits, Xtest=as.matrix(Xtest), ytest=ytest)
    runtime <- getRunTime(AllFits)
    FNR <- getFNR(AllFits)
    FPR <- getFPR(AllFits)
    l1error_intercept <- getl1error_intercept(AllFits)
    l1error_beta <- getl1error_beta(AllFits)

    if (family=="gaussian"){
      RMSE <- getRMSE(AllFits)
      l <- list(FPR=FPR, FNR=FNR, RMSE=RMSE,
                pf_mat=pf_mat, beta_mat=beta_mat,
                intercepts=intercepts, sparsity_mat=sparsity_mat,
                annot=AllFits$annot, runtime=runtime,
                l1error_intercept=l1error_intercept,
                l1error_beta=l1error_beta)
    } else if(family=="binomial"){
      BS <- getBS(AllFits)
      AUC <- getAUC(AllFits)
      l <- list(FPR=FPR, FNR=FNR, BS=BS, AUC=AUC,
                pf_mat=pf_mat, beta_mat=beta_mat,
                intercepts=intercepts, sparsity_mat=sparsity_mat,
                annot=AllFits$annot, runtime=runtime,
                l1error_intercept=l1error_intercept,
                l1error_beta=l1error_beta)
    }
    else stop("Family not implemented")
    return(l)
  }


  if(parallel){
    resultList <- parallel::mclapply(seq_len(nfolds), runPerFold,
                                     mc.cores = ncores)
  } else resultList <- lapply(seq_len(nfolds), runPerFold)

  # plot results
  if(plot_cv) plotMethodComparison(resultList, family = family)

  return(resultList)
}



