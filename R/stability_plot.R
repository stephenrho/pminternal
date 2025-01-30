#' Plot prediction stability across bootstrap replicates
#'
#' @description
#' A prediction (in)stability plot shows estimated risk probabilities from
#' models developed on resampled data evaluated on the original development data
#' as a function of the 'apparent' prediction (prediction from original/development
#' model evaluated on original data). A stable model should produce points that
#' exhibit minimal dispersion. See Riley and Collins (2023).
#'
#' @param x an object produced by \code{\link{validate}} with method = "boot_\*" (or \code{\link{boot_optimism}} with method="boot")
#' @param bounds width of the 'stability interval' (percentiles of the bootstrap model predictions). NULL = do not add bounds to plot.
#' @param smooth_bounds if TRUE, use \code{loess} to smooth the bounds (default = FALSE)
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param pch plotting character (default = 16)
#' @param cex controls point size (default = 0.4)
#' @param col color of points (default = grDevices::grey(.5, .5))
#' @param lty line type for bounds (default = 2)
#' @param span controls the degree of smoothing (see \code{loess}; default = 0.75)
#' @param subset vector of observations to include (row indices). If dataset is large plotting N points for B bootstrap resamples is demanding. This can be used to select a random subset of observations.
#' @param plot if FALSE just returns stability matrix
#'
#' @return plots prediction (in)stability. The stability bounds are not smoothed.
#' Invisibly returns stability matrix (where column 1 are original predictions)
#' that can be used for creating plots with other packages/software.
#'
#' @references Riley RD, Collins GS. (2023). Stability of clinical prediction models developed using statistical or machine learning methods. Biom J. doi:10.1002/bimj.202200302. Epub ahead of print.
#'
#' @export
#'
#' @examples
#' set.seed(456)
#' # simulate data with two predictors that interact
#' dat <- pmcalibration::sim_dat(N = 2000, a1 = -2, a3 = -.3)
#' mean(dat$y)
#' dat$LP <- NULL # remove linear predictor
#'
#' # fit a (misspecified) logistic regression model
#' m1 <- glm(y ~ ., data=dat, family="binomial")
#'
#' # internal validation of m1 via bootstrap optimism with 10 resamples
#' # B = 10 for example but should be >= 200 in practice
#' m1_iv <- validate(m1, method="boot_optimism", B=10)
#'
#' prediction_stability(m1_iv)
#'
prediction_stability <- function(x, bounds=.95, smooth_bounds=FALSE,
                                 xlab, ylab, pch, cex, col, lty, span, subset, plot=TRUE){

  stabil <- get_stability(x)
  stabil <- stabil$stability

  b <- ncol(stabil)-1
  # settings
  if (missing(pch)) pch <- 16
  if (missing(cex)) cex <- .4
  if (missing(lty)) lty <- 2
  if (missing(span)) span <- .75
  if (missing(col)) col <- grey(.5, .3)
  if (missing(xlab)) xlab <- "Estimated risk from development model"
  if (missing(ylab)) ylab <- sprintf("Estimated risk from bootstrap models (n = %i)", b)
  if (!missing(subset)) stabil <- stabil[subset, ]

  # make plot
  stabil <- stabil[order(stabil[,1]), ] # [, 1] = original predictions
  if (plot){
    matplot(stabil[, 1], stabil[, -1], type = "p", pch=pch, cex=cex,
            col=col, xlab = xlab, ylab = ylab)
  }

  if (!is.null(bounds)){
    if (length(bounds) > 1 || !is.numeric(bounds) || bounds <= 0 || bounds >= 1) stop("bounds should be a single number (0,1)")

    probs <- c((1 - bounds)/2, 1 - (1 - bounds)/2)
    lims <- t(apply(stabil[, -1], 1, quantile, probs=probs))
    if (smooth_bounds){
      lims <- apply(lims, 2, function(x) predict(loess(x ~ stabil[, 1], span = span)))
    }
    if (plot){
      matplot(stabil[, 1], lims, type='l', lty=lty, add = T, col="black")
    }
  }

  colnames(stabil)[-1] <- paste0("b", seq(b))
  rownames(stabil) <- NULL
  return(invisible(stabil))
}


#' Plot calibration stability across bootstrap replicates
#'
#' @description
#' A calibration (in)stability plot shows calibration curves for bootstrap
#' models evaluated on original outcome. A stable model should produce
#' boot calibration curves that differ minimally from the 'apparent' curve.
#' See Riley and Collins (2023).
#'
#' @param x an object produced by \code{\link{validate}} with method = "boot_\*" (or \code{\link{boot_optimism}} with method="boot")
#' @param calib_args settings for calibration curve (see \code{pmcalibration::pmcalibration}). If unspecified settings
#' are given by \code{\link{cal_defaults}} with 'eval' set to 100 (evaluate each curve at 100 points between min and max prediction).
#' @param xlim x limits (default = c(0,1))
#' @param ylim y limits (default = c(0,1))
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param col color of lines for bootstrap models (default = grDevices::grey(.5, .3))
#' @param subset vector of observations to include (row indices). If dataset is large fitting B curves is demanding. This can be used to select a random subset of observations.
#' @param plot if FALSE just returns curves (see value)
#'
#' @return plots calibration (in)stability.
#' Invisibly returns a list containing data for each curve (p=x-axis, pc=y-axis).
#' The first element of this list is the apparent curve (original model on original outcome).
#'
#' @references Riley RD, Collins GS. (2023). Stability of clinical prediction models developed using statistical or machine learning methods. Biom J. doi:10.1002/bimj.202200302. Epub ahead of print.
#'
#' @export
#'
#' @examples
#' set.seed(456)
#' # simulate data with two predictors that interact
#' dat <- pmcalibration::sim_dat(N = 2000, a1 = -2, a3 = -.3)
#' mean(dat$y)
#' dat$LP <- NULL # remove linear predictor
#'
#' # fit a (misspecified) logistic regression model
#' m1 <- glm(y ~ ., data=dat, family="binomial")
#'
#' # internal validation of m1 via bootstrap optimism with 10 resamples
#' # B = 10 for example but should be >= 200 in practice
#' m1_iv <- validate(m1, method="boot_optimism", B=10)
#'
#' calibration_stability(m1_iv)
#'
calibration_stability <- function(x, calib_args,
                                  xlim, ylim, xlab, ylab, col, subset, plot=TRUE){

  stabil <- get_stability(x)
  y <- stabil$y # original outcome
  stabil <- stabil$stability

  if (!missing(subset)){
    stabil <- stabil[subset, ]
    y <- y[subset]
  }

  if (missing(calib_args)){
    calib_args <- cal_defaults()
    calib_args[["eval"]] <- 100
  } else{
    if (!"smooth" %in% names(calib_args)){
      calib_args[['smooth']] <- 'gam'
    }
    if (!"transf" %in% names(calib_args)){
      calib_args[['transf']] <- 'logit'
    }
    if ("ci" %in% names(calib_args)){
      message("ci in calib_args is ignored")
    }
    calib_args["ci"] <- "none"
    if ("eval" %in% names(calib_args)){
      if (length(calib_args[["eval"]]) > 1){
        message("replacing vector eval argument with eval = 100. ",
                "To evalulate calibration curves at more points ",
                "please set a higher value for eval.")
        calib_args["eval"] <- 100
      } else{
        calib_args[["eval"]] <- 100
      }
    }
  }

  # make calibration curves
  curves <- apply(stabil, 2, function(p){
    calib_args[["p"]] <- p
    calib_args[["y"]] <- y
    cc <- do.call(pmcalibration::pmcalibration, calib_args)
    data.frame(p = cc$plot$p, pc = cc$plot$p_c_plot)
  })

  # settings
  if (missing(xlim)) xlim <- c(0,1)
  if (missing(ylim)) ylim <- c(0,1)
  if (missing(xlab)) xlab <- "Estimated risk"
  if (missing(ylab)) ylab <- "Observed risk"
  if (missing(col)) col <- grey(.5, .3)

  # max(unlist(curves))
  if (plot){
    plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
    l <- lapply(rev(seq(length(curves))), function(cc){
      with(curves[[cc]], lines(p, pc, lty=if(cc>1) 1 else 2,
                               col=if(cc>1) col else "black", lwd=if(cc>1) 1 else 2))
    })
    legend(x=xlim[1], y=ylim[2], legend = c("Original", sprintf("Bootstrap (n = %i)", ncol(stabil)-1)),
           lty=2:1, col=c("black", col),
           lwd=2:1)
  }

  return(invisible(curves))
}

#' Mean absolute predictor error (MAPE) stability plot
#'
#' @description
#' A MAPE (in)stability plot shows mean absolute predictor error (average absolute difference between original
#' estimated risk and risk from B bootstrap models) as a function of apparent
#' estimated risk (prediction from original/development model). See Riley and Collins (2023).
#'
#' @param x an object produced by \code{\link{validate}} with method = "boot_\*" (or \code{\link{boot_optimism}} with method="boot")
#' @param xlim x limits (default = range of estimated risks)
#' @param ylim y limits (default = c(0, maximum mape))
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param pch plotting character (default = 16)
#' @param cex controls point size (default = 1)
#' @param col color of points (default = grDevices::grey(.5, .5))
#' @param subset vector of observations to include (row indices). This can be used to select a random subset of observations.
#' @param plot if FALSE just returns MAPE values (see value)
#'
#' @return plots calibration (in)stability.
#' Invisibly returns a list containing individual and average MAPE.
#'
#' @references Riley RD, Collins GS. (2023). Stability of clinical prediction models developed using statistical or machine learning methods. Biom J. doi:10.1002/bimj.202200302. Epub ahead of print.
#'
#' @export
#'
#' @examples
#' set.seed(456)
#' # simulate data with two predictors that interact
#' dat <- pmcalibration::sim_dat(N = 2000, a1 = -2, a3 = -.3)
#' mean(dat$y)
#' dat$LP <- NULL # remove linear predictor
#'
#' # fit a (misspecified) logistic regression model
#' m1 <- glm(y ~ ., data=dat, family="binomial")
#'
#' # internal validation of m1 via bootstrap optimism with 10 resamples
#' # B = 10 for example but should be >= 200 in practice
#' m1_iv <- validate(m1, method="boot_optimism", B=10)
#'
#' mape_stability(m1_iv)
#'
mape_stability <- function(x, xlim, ylim, xlab, ylab, pch, cex,
                           col, subset, plot=TRUE){
  stabil <- get_stability(x)
  stabil <- stabil$stability

  if (!missing(subset)){
    stabil <- stabil[subset, ]
  }

  p_orig <- stabil[, 1]
  p_boot <- stabil[, -1]

  pe <- apply(p_boot, 2, function(pb) abs(p_orig - pb))
  individual_mape <- apply(pe, 1, mean)
  average_mape <- mean(individual_mape)

  # settings
  if (missing(pch)) pch <- 16
  if (missing(cex)) cex <- 1
  if (missing(col)) col <- grey(.5, .5)
  if (missing(xlim)) xlim <- range(p_orig)
  if (missing(ylim)) ylim <- c(0, max(individual_mape))
  if (missing(xlab)) xlab <- "Estimated risk from development model"
  if (missing(ylab)) ylab <- "MAPE"

  if (plot){
    matplot(p_orig, individual_mape, type = "p", pch = pch,
            col = col, cex = cex, xlim=xlim, ylim=ylim,
            xlab=xlab, ylab=ylab)
  }

  out <- list('individual_mape' = individual_mape, 'average_mape' = average_mape)

  return(invisible(out))
}

#' Classification instability plot
#'
#' @description
#' Classification instability plot shows the relationship between original model estimated risk
#' and the classification instability index (CII). The CII is the proportion of bootstrap replicates
#' where the predicted class (0 if p <= threshold; 1 if p > threshold) is different to that
#' obtained from the original model. Those with risk predictions around the threshold will exhibit
#' elevated CII but an unstable model will exhibit high CII across a range of risk predictions.
#' See Riley and Collins (2023).
#'
#' @param x an object produced by \code{\link{validate}} with method = "boot_\*" (or \code{\link{boot_optimism}} with method="boot")
#' @param threshold estimated risks above the threshold get a predicted 'class' of 1, otherwise 0.
#' @param xlim x limits (default = range of estimated risks)
#' @param ylim y limits (default = c(0, maximum CII))
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param pch plotting character (default = 16)
#' @param cex controls point size (default = 1)
#' @param col color of points (default = grDevices::grey(.5, .5))
#' @param subset vector of observations to include (row indices). This can be used to select a random subset of observations.
#' @param plot if FALSE just returns CII values (see value)
#'
#' @return plots classification (in)stability.
#' Invisibly returns estimates of CII for each observation.
#'
#' @references Riley RD, Collins GS. (2023). Stability of clinical prediction models developed using statistical or machine learning methods. Biom J. doi:10.1002/bimj.202200302. Epub ahead of print.
#'
#' @export
#'
#' @examples
#' set.seed(456)
#' # simulate data with two predictors that interact
#' dat <- pmcalibration::sim_dat(N = 2000, a1 = -2, a3 = -.3)
#' mean(dat$y)
#' dat$LP <- NULL # remove linear predictor
#'
#' # fit a (misspecified) logistic regression model
#' m1 <- glm(y ~ ., data=dat, family="binomial")
#'
#' # internal validation of m1 via bootstrap optimism with 10 resamples
#' # B = 10 for example but should be >= 200 in practice
#' m1_iv <- validate(m1, method="boot_optimism", B=10)
#'
#' classification_stability(m1_iv, threshold=.2)
#'
classification_stability <- function(x, threshold, xlim, ylim,
                                     xlab, ylab, pch, cex, col, subset, plot=TRUE){
  stabil <- get_stability(x)
  stabil <- stabil$stability

  if (!missing(subset)){
    stabil <- stabil[subset, ]
  }

  p_orig <- stabil[, 1]
  c_orig <- p_orig > threshold
  c_boot <- stabil[, -1] > threshold

  # classification instability
  ci <- apply(c_boot, 2, function(pcla) c_orig == pcla)
  cii <- 1-apply(ci, 1, mean)

  # settings
  if (missing(pch)) pch <- 16
  if (missing(cex)) cex <- 1
  if (missing(col)) col <- grey(.5, .5)
  if (missing(xlim)) xlim <- range(p_orig)
  if (missing(ylim)) ylim <- c(0, max(cii))
  if (missing(xlab)) xlab <- "Estimated risk from development model"
  if (missing(ylab)) ylab <- "Classification Instability Index"

  if (plot){
    plot(p_orig, cii, type = "p", pch=pch, col=col, cex=cex,
         xlim=xlim, ylim=ylim, xlab=xlab,
         ylab=ylab)
    abline(v = threshold, lty=2)
  }

  return(invisible(cii))
}

#' Plot decision curve stability across bootstrap replicates
#'
#' @description
#' A decision curve (in)stability plot shows decision curves for bootstrap
#' models evaluated on original outcome. A stable model should produce
#' curves that differ minimally from the 'apparent' curve.
#' See Riley and Collins (2023).
#'
#' @param x an object produced by \code{\link{validate}} with method = "boot_\*" (or \code{\link{boot_optimism}} with method="boot")
#' @param thresholds points at which to evaluate the decision curves (see \code{dcurves::dca})
#' @param xlim x limits (default = range of thresholds)
#' @param ylim y limits (default = range of net benefit)
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param col color of points (default = grDevices::grey(.5, .5))
#' @param subset vector of observations to include (row indices). This can be used to select a random subset of observations.
#' @param plot if FALSE just returns curves (see value)
#'
#' @return plots decision curve (in)stability.
#' Invisibly returns a list containing data for each curve. These are returned from \code{dcurves::dca}.
#' The first element of this list is the apparent curve (original model on original outcome).
#'
#' @references Riley RD, Collins GS. (2023). Stability of clinical prediction models developed using statistical or machine learning methods. Biom J. doi:10.1002/bimj.202200302. Epub ahead of print.
#'
#' @export
#'
#' @examples
#' set.seed(456)
#' # simulate data with two predictors that interact
#' dat <- pmcalibration::sim_dat(N = 2000, a1 = -2, a3 = -.3)
#' mean(dat$y)
#' dat$LP <- NULL # remove linear predictor
#'
#' # fit a (misspecified) logistic regression model
#' m1 <- glm(y ~ ., data=dat, family="binomial")
#'
#' # internal validation of m1 via bootstrap optimism with 10 resamples
#' # B = 10 for example but should be >= 200 in practice
#' m1_iv <- validate(m1, method="boot_optimism", B=10)
#'
#' dcurve_stability(m1_iv)
#'
dcurve_stability <- function(x, thresholds = seq(0, .99, by=0.01),
                             xlim, ylim, xlab, ylab, col, subset, plot=TRUE){

  stabil <- get_stability(x)
  y <- stabil$y # original outcome
  stabil <- stabil$stability

  if (!missing(subset)){
    stabil <- stabil[subset, ]
    y <- y[subset]
  }

  # make calibration curves
  curves <- apply(stabil, 2, function(p){
    dc <- suppressMessages(dcurves::dca(y ~ p, data=data.frame(y=y, p=p),
                                        thresholds = thresholds))
    curve <- dc$dca
    #subset(curve, label == "p")
    curve[curve$label == "p",]
  })

  # settings
  if (missing(col)) col <- grey(.5, .5)
  if (missing(xlim)) xlim <- range(thresholds)
  if (missing(ylim)) ylim <- range(unlist(lapply(curves, function(x) x$net_benefit)))
  if (missing(xlab)) xlab <- "Estimated risk"
  if (missing(ylab)) ylab <- "Net Benefit"

  # max(unlist(curves))
  if (plot){
    plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
    l <- lapply(rev(seq(length(curves))), function(dc){
      with(curves[[dc]], lines(threshold, net_benefit, lty=if(dc>1) 1 else 2,
                               col=if(dc>1) col else "black", lwd=if(dc>1) 1 else 2))
    })
    legend(x=xlim[2], y=ylim[2], legend = c("Original", sprintf("Bootstrap (n = %i)", ncol(stabil)-1)),
           lty=2:1, col=c("black", col),
           lwd=2:1, xjust = 1)
    abline(h=0, lty=2)
  }

  return(invisible(curves))
}
