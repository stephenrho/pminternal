#' Score predictions for binary events
#'
#' @description
#' Calculate scores summarizing discrimination/calibration of predictions
#' against observed binary events. If score_fun is not defined when calling
#' \code{\link{validate}} this function is used.
#'
#' @param y vector containing a binary outcome
#' @param p vector of predictions
#' @param ... additional arguments. This function only supports calib_args as
#' an optional argument. calib_args should contain arguments for pmcalibration::pmcalibration.
#' If a calibration plot (apparent vs bias corrected calibration curves via \code{\link{cal_plot}})
#' is desired the argument 'eval' should be provided. This should be the points at which to evaluate
#' the calibration curve on each boot resample or crossvalidation fold. A good option would be
#' calib_args = list(eval = seq(min(p), max(p), length.out=100)); where p are predictions from the
#' original model evaluated on the original data. Dots can be used to supply additional arguments to
#' user-defined functions.
#'
#' @return a named vector of scores
#' @export
#'
#' @examples
#' p <- runif(100)
#' y <- rbinom(length(p), 1, p)
#' score_binary(y = y, p = p)
score_binary <- function(y, p, ...){

  # TODO: add logistic intercept/slope

  dots <- list(...)

  # pmcalib_defaults <- list('smooth' = 'gam', 'k'=10,
  #                          'ci' = 'none', 'transf' = 'logit',
  #                          'eval' = 0)

  if ("calib_args" %in% names(dots)){
    calib_args <- dots[["calib_args"]]
    if (isFALSE(is.list(calib_args))){
      stop("calib_args should be a list")
    }
    if (!"smooth" %in% names(calib_args)){
      calib_args[['smooth']] <- 'gam'
      calib_args[['k']] <- 10
    }
    if (!"transf" %in% names(calib_args)){
      calib_args[['transf']] <- 'logit'
    }
    if ("ci" %in% names(calib_args)){
      message("ci in calib_args is ignored")
    }
    calib_args["ci"] <- "none"
    if ("eval" %in% names(calib_args)){
      if (length(calib_args[["eval"]]) == 1){
        warning("if specifying 'eval' in calib_args we strongly suggest specifying ",
                "a vector of probabilities at which to evaluate the calibration curve. ",
                "e.g., eval = seq(min(p), max(p), length.out=100), where p are predictions from ",
                "original model evaluated on the original data. Only specifying only the number ",
                "of points to eval (e.g., eval = 100) will likely cause issues/errors.")
      }
    }
  } else{
    calib_args <- cal_defaults() # pmcalib_defaults
  }

  calib_args[["y"]] <- y
  calib_args[["p"]] <- p

  cal <- do.call(pmcalibration::pmcalibration, calib_args)

  brier <- mean((y - p)^2)
  names(brier) <- "brier"
  auc <- as.numeric(suppressMessages(pROC::auc(response = y, predictor = p)))
  names(auc) <- "C"

  if (!is.null(cal$plot$p_c_plot)){
    cal_plot <- cal$plot$p_c_plot
    # names(cal_plot) <- paste0( "cal_plot_", calib_args[["eval"]] )
    names(cal_plot) <- paste0( "cal_plot_", cal$plot$p )
  } else{
    cal_plot <- c()
  }

  scores <- c(brier, auc, cal$metrics, cal_plot)
  #names(scores) <- c("brier", "C")

  return(scores)
}