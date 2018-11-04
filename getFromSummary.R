#' @title get penalty factors
#' @name getPenaltyFactors
#' @description Function to obtain the sparsity levels (1=dense, 0=sparse) estimated by various methods from the fitted regression models.
#' @param AllFits object as produced by \code{\link{runMethods}}
#' @export
#' @return A matrix of penalty factors per group (rows) and method (columns).
#' @examples
#' dat <- makeExampleData()
#' fit <- runMethods(dat$X, dat$y, dat$annot)
#' getPenaltyFactors(fit)

getPenaltyFactors <- function(AllFits) {
    G <- AllFits$G
    pfmat <- vapply(AllFits$summaryList, function(l) {
        pfs <- l$pf
        if (!is.null(pfs))
            as.numeric(pfs) else rep(NA, G)
    }, numeric(G))
    rownames(pfmat) <- unique(AllFits$groupnames)
    return(pfmat)
}

#' @title get sparsity levels
#' @name getSparsityLevel
#' @description Function to obtain the sparsity levels (1=dense, 0=sparse) estimated by various methods from the fitted regression models.
#' @param AllFits object as produced by \code{\link{runMethods}}
#' @export
#' @return A matrix of sparsity levels per group (rows) and method (columns).
#' @examples
#' dat <- makeExampleData()
#' fit <- runMethods(dat$X, dat$y, dat$annot)
#' getSparsityLevel(fit)

getSparsityLevel <- function(AllFits) {
    G <- AllFits$G
    sparsity_mat <- vapply(AllFits$summaryList, function(l) {
        sparsity <- l$sparsity
        if (!is.null(sparsity))
            as.numeric(sparsity) else rep(NA, G)
    }, numeric(G))
    rownames(sparsity_mat) <- unique(AllFits$groupnames)
    return(sparsity_mat)
}

#' @title get coefficients
#' @name getCoefficients
#' @description Function to obtain the coefficients estimated by various methods from the fitted regression models.
#' @param AllFits object as produced by \code{\link{runMethods}}
#' @export
#' @return A matrix of coefficients per feature (rows) and method (columns).
#' @examples
#' dat <- makeExampleData()
#' fit <- runMethods(dat$X, dat$y, dat$annot)
#' getCoefficients(fit)

getCoefficients <- function(AllFits) {
    p <- AllFits$p
    coefmat <- vapply(AllFits$summaryList, function(l) {
        coef <- l$beta
        if (!is.null(coef))
            as.numeric(coef) else rep(NA, p)
    }, numeric(p))
    rownames(coefmat) <- AllFits$varnames
    return(coefmat)
}

#' @title get intercept
#' @name getIntercept
#' @description Function to obtain the intercepts estimated by various methods from the fitted regression models.
#' @param AllFits object as produced by \code{\link{runMethods}}
#' @export
#' @return A vector of intercept values for each method.
#' @examples
#' dat <- makeExampleData()
#' fit <- runMethods(dat$X, dat$y, dat$annot)
#' getIntercept(fit)

getIntercept <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        intercept <- l$intercept
        if (is.null(intercept))
            intercept <- NA
        intercept
    }, numeric(1))
}

#' @title get run times
#' @name getRunTime
#' @description Function to obtain the runtimes of various methods from the fitted regression models.
#' @param AllFits object as produced by \code{\link{runMethods}}
#' @export
#' @return A vector of run times for each method.
#' @examples
#' dat <- makeExampleData()
#' fit <- runMethods(dat$X, dat$y, dat$annot)
#' getRunTime(fit)

getRunTime <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        runtime <- l$runtime
        if (is.null(runtime))
            runtime <- NA
        runtime
    }, numeric(1))
}

#' @title get RMSE
#' @name getRMSE
#' @description Function to obtain the root mean squaed error (RMSE) of various methods in a regression model.
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export
#' @return A vector of the root mean squared error on the test data for each method.
#' @examples
#' dat <- makeExampleData()
#' ntrain <- dat$n/2
#' fit <- runMethods(dat$X[seq_len(ntrain),],
#'  dat$y[seq_len(ntrain)], dat$annot)
#' evalFit <- evaluateFits(fit, dat$X[seq_len(ntrain) + dat$n/2,],
#'  dat$y[seq_len(ntrain)+ dat$n/2])
#' getRMSE(evalFit)

getRMSE <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        RMSE <- l$RMSE
        if (is.null(RMSE))
            RMSE <- NA
        RMSE
    }, numeric(1))
}

#' @title get Brier Score
#' @name getBS
#' @description Function to obtain the Brier Score of various methods in a logistic regression model.
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export
#' @return A vector of the Brier score on the test data for each method.
#' @examples
#' dat <- makeExampleData(response = "bernoulli")
#' ntrain <- dat$n/2
#' fit <- runMethods(dat$X[seq_len(ntrain),], dat$y[seq_len(ntrain)], dat$annot, family = "binomial")
#' evalFit <- evaluateFits(fit, dat$X[seq_len(ntrain) + dat$n/2,],
#'  dat$y[seq_len(ntrain)+ dat$n/2])
#' getBS(evalFit)

getBS <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        BS <- l$BS
        if (is.null(BS))
            BS <- NA
        BS
    }, numeric(1))
}



#' @title get AUC
#' @name getAUC
#' @description  Function to obtain the AUC of various methods in a logistic regression model.
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export
#' @return A vector of the AUC values on the test data for each method.
#' @examples
#' dat <- makeExampleData(response = "bernoulli")
#' ntrain <- dat$n/2
#' fit <- runMethods(dat$X[seq_len(ntrain),], dat$y[seq_len(ntrain)], dat$annot, family = "binomial")
#' evalFit <- evaluateFits(fit, dat$X[seq_len(ntrain) + dat$n/2,],
#' dat$y[seq_len(ntrain)+ dat$n/2])
#' getAUC(evalFit)

getAUC <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        AUC <- l$AUC
        if (is.null(AUC))
            AUC <- NA
        AUC
    }, numeric(1))
}


#' @title get MSE
#' @name getMSE
#' @description get MSE of various methods
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export
#' @return A vector of the mean squared error on the test data for each method.
#' @examples
#' dat <- makeExampleData()
#' ntrain <- dat$n/2
#' fit <- runMethods(dat$X[seq_len(ntrain),],
#'  dat$y[seq_len(ntrain)], dat$annot)
#' evalFit <- evaluateFits(fit, dat$X[seq_len(ntrain) + dat$n/2,],
#'  dat$y[seq_len(ntrain)+ dat$n/2])
#' getMSE(evalFit)

getMSE <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        MSE <- l$RMSE^2
        if (is.null(MSE))
            MSE <- NA
        MSE
    }, numeric(1))
}

#' @export

#' @title get FNR
#' @name getFNR
#' @description get false negative rates of various methods
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export
#' @return A vector of the FNR on the selected coefficients for each method.
#' @examples
#' dat <- makeExampleData()
#' ntrain <- dat$n/2
#' fit <- runMethods(dat$X[seq_len(ntrain),],
#'  dat$y[seq_len(ntrain)], dat$annot, beta0 = dat$beta)
#' evalFit <- evaluateFits(fit, dat$X[seq_len(ntrain) + dat$n/2,],
#'  dat$y[seq_len(ntrain)+ dat$n/2])
#' getFNR(evalFit)

getFNR <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        FNR <- l$FNR
        if (is.null(FNR))
            FNR <- NA
        FNR
    }, numeric(1))
}


#' @title get FPR
#' @name getFPR
#' @description get false positive rates of various methods
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export
#' @return A vector of the FPR on the selected coefficients for each method.
#' @examples
#' dat <- makeExampleData()
#' ntrain <- dat$n/2
#' fit <- runMethods(dat$X[seq_len(ntrain),],
#'  dat$y[seq_len(ntrain)], dat$annot, beta0 = dat$beta)
#' evalFit <- evaluateFits(fit, dat$X[seq_len(ntrain) + dat$n/2,],
#'  dat$y[seq_len(ntrain)+ dat$n/2])
#' getFPR(evalFit)

getFPR <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        fpr <- l$FPR
        if (is.null(fpr))
            fpr <- NA
        fpr
    }, numeric(1))
}

#' @title get l1 error on beta
#' @name getl1error_beta
#' @description get absolute error on beta
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export
#' @return A vector of the absolute error on the estimated coefficients for each method.
#' @examples
#' dat <- makeExampleData()
#' ntrain <- dat$n/2
#' fit <- runMethods(dat$X[seq_len(ntrain),],
#'  dat$y[seq_len(ntrain)], dat$annot, beta0 = dat$beta)
#' evalFit <- evaluateFits(fit, dat$X[seq_len(ntrain) + dat$n/2,],
#'  dat$y[seq_len(ntrain)+ dat$n/2])
#' getl1error_beta(evalFit)

getl1error_beta <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        l1error_beta <- l$l1error_beta
        if (is.null(l1error_beta))
            l1error_beta <- NA
        l1error_beta
    }, numeric(1))
}

#' @title get l1 error on intercept
#' @name getl1error_intercept
#' @description get absolute error on intercept
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export
#' @return A vector of the absolute error on the estimated intercept for each method.
#' @examples
#' dat <- makeExampleData(intercept = 1)
#' ntrain <- dat$n/2
#' fit <- runMethods(dat$X[seq_len(ntrain),],
#'  dat$y[seq_len(ntrain)], dat$annot,
#'   beta0 = dat$beta, trueintercept = 1)
#' evalFit <- evaluateFits(fit, dat$X[seq_len(ntrain) + dat$n/2,],
#'  dat$y[seq_len(ntrain)+ dat$n/2])
#' getl1error_intercept(evalFit)

getl1error_intercept <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
      l1error_intercept <- l$l1error_intercept
        if (is.null(l1error_intercept))
          l1error_intercept <- NA
        l1error_intercept
    }, numeric(1))
}
