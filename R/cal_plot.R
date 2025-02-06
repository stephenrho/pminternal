#' Plot apparent and bias-corrected calibration curves
#'
#' @param x an object returned from \code{\link{validate}}. Original call should
#' have specified 'eval' argument. See \code{\link{score_binary}}.
#' @param xlim x limits (default = c(0, max of either curve))
#' @param ylim y limits (default = c(0, max of either curve))
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param app_col color of the apparent calibration curve (default = 'black')
#' @param bc_col color of the bias-corrected calibration curve (default = 'black')
#' @param app_lty line type of the apparent calibration curve (default = 1)
#' @param bc_lty line type of the bias-corrected calibration curve (default = 2)
#' @param plotci plot confidence intervals ('yes') or not ('no'). If 'yes' x should have confidence intervals
#' added by \code{\link{confint.internal_validate}}. 'if' (default) plots CIs if they are
#' available.
#'
#' @return plots apparent and bias-corrected curves. Silently returns a data.frame
#' that can be used to produce a more 'publication ready' plot. Columns are as
#' follows: predicted = values for the x-axis, apparent = value of apparent curve,
#' bias_corrected = value of bias-corrected curve. Confidence intervals are included
#' if available.
#' @export
#' @examples
#' library(pminternal)
#' set.seed(456)
#' # simulate data with two predictors that interact
#' dat <- pmcalibration::sim_dat(N = 2000, a1 = -2, a3 = -.3)
#' mean(dat$y)
#' dat$LP <- NULL # remove linear predictor
#'
#' # fit a (misspecified) logistic regression model
#' m1 <- glm(y ~ x1 + x2, data=dat, family="binomial")
#'
#' # to get a plot of bias-corrected calibration we need
#' # to specify 'eval' argument via 'calib_args'
#' # this argument specifies at what points to evalulate the
#' # calibration curve for plotting. The example below uses
#' # 100 equally spaced points between the min and max
#' # original prediction.
#'
#' p <- predict(m1, type="response")
#' p100 <- seq(min(p), max(p), length.out=100)
#'
#' m1_iv <- validate(m1, method="cv_optimism", B=10,
#'                   calib_args = list(eval=p100))
#' # calib_ags can be used to set other calibration curve
#' # settings: see pmcalibration::pmcalibration
#'
#' cal_plot(m1_iv)
#'
cal_plot <- function(x, xlim, ylim, xlab, ylab,
                     app_col, bc_col, app_lty, bc_lty, plotci=c("if", "yes", "no")){

  plotci <- match.arg(plotci)

  if (isFALSE(is(x, "internal_validate"))){
    stop("x should be an object returned from call to pminternal::validate")
  }

  if ("calib_args" %in% names(x$dots) && "eval" %in% names(x$dots$calib_args)){
    p <- x$dots$calib_args$eval
  } else{
    stop("original call to pminternal::validate should specify calib_args and specifically ",
         "eval (i.e., calib_args = list(eval = p); where p is a vector of points at which ",
         "to evalulate the calibration curve). See ?pminternal::validate")
  }

  i <- grep(pattern = "cal_plot", x = names(x$apparent), value = TRUE)
  if (length(i) == 0){
    stop("no scores starting 'cal_plot' found. In call to validate please sspecify calib_args and specifically ",
         "eval (i.e., calib_args = list(eval = p); where p is a vector of points at which ",
         "to evalulate the calibration curve). See ?pminternal::validate")
  }
  app <- x$apparent[i]
  bcor <- x$corrected[i]

  ci <- FALSE
  if ("confint" %in% names(x) & plotci %in% c("if", "yes")){
    if (all(i %in% rownames(x$confint$ci$Apparent))){
      ci <- TRUE
      ci_app <- x$confint$ci$Apparent[i, ]
      ci_bcor <- x$confint$ci$Corrected[i, ]
    }
  }

  if (isFALSE(ci) & plotci == "yes") message("Confidence intervals requested but not available in x. Please run confint to add CIs")

  shade <- function(x, ymin, ymax, color){
    graphics::polygon(x=c(x, rev(x)), y = c(ymin, rev(ymax)),
                      col = grDevices::adjustcolor(color, alpha.f = .2), border = NA)
  }

  lim <- c(0, max(c(bcor, app))*1.1)

  # settings
  if (missing(app_col)) app_col <- "black"
  if (missing(bc_col)) bc_col <- "black"
  if (missing(app_lty)) app_lty <- 1
  if (missing(bc_lty)) bc_lty <- 2
  if (missing(xlab)) xlab <- "Estimated risk"
  if (missing(ylab)) ylab <- "Observed risk"
  if (missing(xlim)) xlim <- lim
  if (missing(ylim)) ylim <- lim

  plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
  abline(a = 0, b = 1, col="lightgrey", lwd=2, lty=3)
  if (ci){
    shade(x = p, ymin = ci_app[, 1], ymax = ci_app[, 2],
          color = app_col)

    shade(x = p, ymin = ci_bcor[, 1], ymax = ci_bcor[, 2],
          color = bc_col)
  }
  lines(p, app, lwd=2, lty=app_lty, col=app_col)
  lines(p, bcor, lty=bc_lty, lwd=2, col=bc_col)
  legend(x = xlim[1], y = ylim[2], legend = c("Apparent", "Bias Corrected"),
         col = c(app_col, bc_col), lty = c(app_lty,bc_lty), lwd=2)

  if (ci){
    colnames(ci_app) <- c("app_lower", "app_upper")
    colnames(ci_bcor) <- c("corr_lower", "corr_upper")

    out <- data.frame(predicted = p, apparent = app, ci_app,
                      bias_corrected = bcor, ci_bcor,
                      row.names = NULL)

  } else{
    out <- data.frame(predicted = p, apparent = app,
                      bias_corrected = bcor,
                      row.names = NULL)
  }
  return(invisible(out))
}
