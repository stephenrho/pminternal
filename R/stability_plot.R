#' Plot prediction stability across bootstrap replicates
#'
#' @description
#' A prediction (in)stability plot shows predicted risk probabilities from
#' models devloped on resampled data evaluated on the original development data
#' as a function of the 'apparent' prediction (prediction from original/development
#' model evalulated on original data). A stable model should produce points that
#' exhibit minimal dispersion. See Riley and Collins (2023).
#'
#' @param x an object produced by \code{\link{validate}} with method = "boot_\*" (or \code{\link{boot_optimism}} with method="boot")
#' @param bounds width of the 'stability interval' (percentiles of the bootstrap model predictions). NULL = do not add bounds to plot.
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
#' library(pminternal)
#' set.seed(456)
#' # simulate data with two predictors that interact
#' dat <- pmcalibration::sim_dat(N = 2000, a1 = -2, a3 = -.3)
#' mean(dat$y)
#' dat$LP <- NULL # remove linear predictor
#'
#' # fit a (misspecified) logistic regression model
#' m1 <- glm(y ~ ., data=dat, family="binomial")
#'
#' # internal validation of m1 via bootstrap optimism with 20 resamples
#' # B = 20 for example but should be >= 200 in practice
#' m1_iv <- validate(m1, method="boot_optimism", B=20)
#'
#' prediction_stability(m1_iv)
#'
prediction_stability <- function(x, bounds=.95){

  stabil <- get_stability(x)
  stabil <- stabil$stability

  b <- ncol(stabil)-1
  # make plot
  stabil = stabil[order(stabil[,1]), ] # [, 1] = original predictions
  matplot(stabil[, 1], stabil[, -1], type = "p", pch=16, cex=.4,
          col=grey(.5, .3), xlab = "Estimated risk from development model",
          ylab = sprintf("Estimated risk from bootstrap models (n = %i)", b))

  if (!is.null(bounds)){
    if (length(bounds) > 1 || !is.numeric(bounds) || bounds <= 0 || bounds >= 1) stop("bounds should be a single number (0,1)")

    probs <- c((1 - bounds)/2, 1 - (1 - bounds)/2)
    lims <- t(apply(stabil[, -1], 1, quantile, probs=probs))
    matplot(stabil[, 1], lims, type='l', lty=2, add = T, col="black")
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
#' are given by \code{\link{cal_defaults}} with 'eval' set to 100 (evalulate each curve at 100 points between min and max prediction).
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
#' library(pminternal)
#' set.seed(456)
#' # simulate data with two predictors that interact
#' dat <- pmcalibration::sim_dat(N = 2000, a1 = -2, a3 = -.3)
#' mean(dat$y)
#' dat$LP <- NULL # remove linear predictor
#'
#' # fit a (misspecified) logistic regression model
#' m1 <- glm(y ~ ., data=dat, family="binomial")
#'
#' # internal validation of m1 via bootstrap optimism with 20 resamples
#' # B = 20 for example but should be >= 200 in practice
#' m1_iv <- validate(m1, method="boot_optimism", B=20)
#'
#' calibration_stability(m1_iv)
#'
calibration_stability <- function(x, calib_args){

  stabil <- get_stability(x)
  y <- stabil$y # original outcome
  stabil <- stabil$stability

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

  # max(unlist(curves))
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="Predicted risk", ylab="Estimated risk")
  l <- lapply(rev(seq(length(curves))), function(cc){
    with(curves[[cc]], lines(p, pc, lty=if(cc>1) 1 else 2,
         col=if(cc>1) grey(.5, .3) else "black", lwd=if(cc>1) 1 else 2))
  })
  legend(x=0, y=1, legend = c("Original", "Bootstrap"), lty=2:1, col=c("black", grey(.5, 1)),
         lwd=2:1)

  return(invisible(curves))
}

#' Mean absolute predictor error (MAPE) stability across bootstrap replicates
#'
#' @description
#' A MAPE (in)stability plot shows mean absolute predictor error (average absolute difference between original
#' predicted risk and predicted risk from B bootstrap models) as a function of apparent
#' predicted risk (prediction from original/development model). See Riley and Collins (2023).
#'
#' @param x an object produced by \code{\link{validate}} with method = "boot_\*" (or \code{\link{boot_optimism}} with method="boot")
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
#' library(pminternal)
#' set.seed(456)
#' # simulate data with two predictors that interact
#' dat <- pmcalibration::sim_dat(N = 2000, a1 = -2, a3 = -.3)
#' mean(dat$y)
#' dat$LP <- NULL # remove linear predictor
#'
#' # fit a (misspecified) logistic regression model
#' m1 <- glm(y ~ ., data=dat, family="binomial")
#'
#' # internal validation of m1 via bootstrap optimism with 20 resamples
#' # B = 20 for example but should be >= 200 in practice
#' m1_iv <- validate(m1, method="boot_optimism", B=20)
#'
#' mape_stability(m1_iv)
#'
mape_stability <- function(x){
  stabil <- get_stability(x)
  stabil <- stabil$stability

  p_orig <- stabil[, 1]
  p_boot <- stabil[, -1]

  pe <- apply(p_boot, 2, function(pb) abs(p_orig - pb))
  individual_mape <- apply(pe, 1, mean)
  average_mape <- mean(individual_mape)

  xlim <- range(p_orig)
  ylim <- c(0, max(individual_mape))
  matplot(p_orig, individual_mape, type = "p", pch = 16, col=grey(.5, .5),
          xlim=xlim, ylim=ylim, xlab="Predicted risk from development model",
          ylab="MAPE")

  out <- list('individual_mape' = individual_mape, 'average_mape' = average_mape)

  return(invisible(out))
}
