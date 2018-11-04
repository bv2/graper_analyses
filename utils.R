library(cowplot)
library(reshape2)
library(dplyr)

#' @title Plot performance measures for various methods
#' @name plotMethodComparison
#' @description Function to plot performance measures of the fits obtained from various regression methods that were evaluated using to \code{\link{compareMethodsCV}}.
#' @param resultList list as created by \code{\link{compareMethodsCV}}
#' @param family gaussian or binomial (same as used in \code{\link{compareMethodsCV}})
#' @param methods2plot which methods to include into the fit
#' @importFrom dplyr filter bind_rows
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @export
#' @return a ggplot object
#' @examples
#' dat <- makeExampleData()
#' cv.out <- compareMethodsCV(dat$X, dat$y, dat$annot, nfolds=3)
#' plotMethodComparison(cv.out)

plotMethodComparison <- function(resultList, family = "gaussian", methods2plot="all") {
  # avoid notes on global varibale binding in check
  runtime <- method <- group <- penalty_factor <- NULL
  sparsity_level <- value <- runtime <-  NULL
  
  # get results in dataframe format
  if(family=="gaussian"){
    if(!all(is.na(resultList[[1]]$FPR))) {
      eval_summary <- reshape2::melt(lapply(resultList,function(l) {
        rbind(FPR = l$FPR, FNR = l$FNR,F1score = l$F1score,
              RMSE = l$RMSE, l1error_beta = l$l1error_beta)
      }), varnames = c("measure", "method"), level = "run")
    } else {
      eval_summary <- reshape2::melt(lapply(resultList, function(l){
        rbind(RMSE = l$RMSE)
      }), varnames = c("measure", "method"), level = "run")
    }
  } else {
    if(!all(is.na(resultList[[1]]$FPR))) {
      eval_summary <- reshape2::melt(lapply(resultList, function(l) {
        rbind(FPR = l$FPR, FNR = l$FNR,  F1score = l$F1score,
              BS = l$BS, AUC = l$AUC, l1error_beta = l$l1error_beta)
      }), varnames = c("measure", "method"), level = "run")
    } else{
      eval_summary <- reshape2::melt(lapply(resultList, function(l) {
        rbind(BS = l$BS, AUC = l$AUC)
      }), varnames = c("measure", "method"), level = "run")
    }
  }
  
  pf_summary <- lapply(seq_along(resultList), function(i){
    cbind(reshape2::melt(resultList[[i]]$pf_mat,
                         varnames = c("group", "method"),
                         value.name = "penalty_factor"), Lrun = i)
  })
  pf_summary <- dplyr::bind_rows(pf_summary)
  
  sparsity_summary <- lapply(seq_along(resultList), function(i) {
    cbind(reshape2::melt(resultList[[i]]$sparsity_mat,
                         varnames = c("group","method"),
                         value.name = "sparsity_level"),
          Lrun = i)
  })
  sparsity_summary <- dplyr::bind_rows(sparsity_summary)
  
  runtime_summary <- reshape2::melt(lapply(resultList, function(l){
    t(l$runtime)
  }),  varnames = c("const", "method"),
  value.name = "runtime", level = "run")[, 2:4]
  
  if(!any(methods2plot=="all")) {
    eval_summary <- dplyr::filter(eval_summary, method %in% methods2plot)
    pf_summary <- dplyr::filter(pf_summary, method %in% methods2plot)
    sparsity_summary <- dplyr::filter(sparsity_summary, method %in% methods2plot)
    runtime_summary <- dplyr::filter(runtime_summary, method %in% methods2plot)
  }
  
  gg_pf <- ggplot(pf_summary, aes(x = as.factor(group), y = penalty_factor,
                                  fill = as.factor(group), group = as.factor(group))) +
    geom_boxplot() + facet_wrap(~method, scales = "free_y") +
    ggtitle("Penalty Factors per group") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  gg_sparse <- ggplot(sparsity_summary, aes(x = as.factor(group), y = sparsity_level,
                                            fill = as.factor(group), group = as.factor(group))) +
    geom_boxplot() + facet_wrap(~method, scales = "free_y") +
    ggtitle("Sparsity Level per group (1=dense)") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  gg_perf <- ggplot(dplyr::filter(eval_summary, method != "TrueModel"),
                    aes(x = method, y = value, fill = method)) +
    geom_boxplot() + ggtitle("Method comparison") +
    facet_wrap(~measure, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    ggtitle("Performance measures")
  
  gg_runtime <- ggplot(runtime_summary,
                       aes(x = method, y = runtime, group = method, fill = method)) +
    geom_boxplot() + theme(axis.text.x = element_text(angle = 60,
                                                      hjust = 1)) + ggtitle("Runtime") + ylab("secs")
  
  print(cowplot::plot_grid(gg_perf,gg_runtime, gg_pf, gg_sparse))
  
}


#' @title Evaluate a fitted regression model on test data
#' @name evaluateModel
#' @description Function to evaluate a fitted regression model on test data in terms of their prediction and feature selection properties.
#' @param beta_est estimated model coefficients
#' @param intercept estimated intercept
#' @param X_test predictor matrix of test data of size
#'  number of test samples (n_test) times features (p)
#' @param y_test response vector of length number of test samples (n_test)
#' @param beta0 true model coefficients (if known)
#' @param family likelihood model for the response, either
#'  "gaussian" for linear regression or "binomial" for logistic regression
#' @importFrom GRridge auc roc
#' @export
#' @return A list containing various performance measures such as RMSE, FNR, FPR etc.
#' @examples
#' dat <- makeExampleData()
#' # take a null model with all coefficeints set to zero for evaluation
#' beta <- rep(0, dat$p)
#' eval.out <- evaluateModel(beta, intercept = 0, X_test = dat$X,
#'  y_test = dat$y, beta0 = dat$beta, family = "gaussian")

evaluateModel <- function(beta_est, intercept, X_test,
                          y_test, beta0 = NULL, family) {
  #sanity check
  if(length(beta_est) != ncol(X_test)) stop("Number of columns in X_test has to agree with length of beta_est.")

  if (is.null(intercept))
        intercept <- 0
    intercept <- as.numeric(intercept)
    p <- ncol(X_test)
    n_test <- nrow(X_test)
    stopifnot(family %in% c("gaussian", "binomial"))
    RMSE_test <- FNR <- sensitivity <- specificity <- FPR <- BrierScore <- NULL
    ROC <- AUC <- predprob <- pred_gauss <- precision <- recall <- F1score <- NULL

    # if 'true features' known evaluate sensitivity and specificity
    if (!is.null(beta0)) {
        stopifnot(length(beta0) == length(beta_est))
        # Feature selection performance
        FPR <- sum(beta0 == 0 & beta_est != 0)/sum(beta0 == 0)
        specificity <- 1 - FPR
        FNR <- sum(beta0 != 0 & beta_est == 0)/sum(beta0 != 0)
        sensitivity <- 1 - FNR
        precision <- sum(beta0 != 0 & beta_est != 0)/sum(beta_est != 0)
        recall <- sum(beta0 != 0 & beta_est != 0)/sum(beta0 != 0)
        F1score <- 2*precision*recall/(precision + recall)
    } else NULL

    # Prediction performance on test set
    if (family == "gaussian") {
        # RMSE
        pred_gauss <- intercept + X_test %*% beta_est
        RMSE_test <- sqrt(sum((y_test - pred_gauss)^2)/length(y_test))
    } else if (family == "binomial") {
        # Brier score
      # to avoid NaN for too large numbers
        predexp <- pmin(intercept + X_test %*% beta_est, 500)
        predprob <- exp(predexp)/(1 + exp(predexp))
        BrierScore <- sum((y_test - predprob)^2)/length(y_test)

        # ROC (use function from GRridge)
        cutoffs <- rev(seq(0, 1, by = 0.01))
        ROC <- GRridge::roc(predprob, y_test, cutoffs)

        # AUC (use function from GRridge)
        AUC <- GRridge::auc(ROC)
    }


    return(list(beta_est = beta_est, RMSE_test = RMSE_test, BrierScore = BrierScore,
                ROC = ROC, AUC = AUC, specificity = specificity,
                sensitivity = sensitivity, FNR = FNR, FPR = FPR, precision=precision,
                recall=recall, F1score=F1score, beta0 = beta0, y_test = y_test,
                predprob = predprob, pred_gauss = pred_gauss, intercept = intercept))
}
