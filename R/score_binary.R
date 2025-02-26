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
#' @details
#' The following measures are returned in a named vector.
#'
#' \describe{
#' \item{C}{the c-statistic (aka area under the ROC curve). Probability that randomly selected
#' observation with y = 1 with have higher p compared to randomly selected y = 0.}
#' \item{Brier}{mean squared error - mean((y - p)^2)}
#' \item{Intercept}{Intercept from a logistic calibration model: glm(y ~ 1 + qlogis(p), family="binomial")}
#' \item{Slope}{Slope from a logistic calibration model: glm(y ~ 1 + qlogis(p), family="binomial")}
#' \item{Eavg}{average absolute difference between p and calibration curve
#' (aka integrated calibration index or ICI).}
#' \item{E50}{median absolute difference between p and calibration curve}
#' \item{E90}{90th percentile absolute difference between p and calibration curve}
#' \item{Emax}{maximum absolute difference between p and calibration curve}
#' \item{ECI}{average squared difference between p and calibration curve. Estimated
#' calibration index (Van Hoorde et al. 2015)}
#' \item{cal_plot}{if eval is specified (via calib_args), values for
#' plotting apparent and bias-corrected calibration curves are returned (see \code{\link{cal_plot}}).
#' By default these are omitted from the summary printed (see \code{\link{summary.internal_validate}}).}
#' }
#'
#' Logistic calibration and other calibration metrics from non-linear calibration curves
#' assessing 'moderate-calibration' (Eavg, E50, E90, Emax, ECI; see references) are calculated
#' via the \code{pmcalibration} package. The default settings can be modified by passing
#' calib_args to \code{\link{validate}} call. calib_args should be a named list corresponding to
#' arguments to \code{pmcalibration::pmcalibration}.
#'
#' @return a named vector of scores (see Details)
#'
#' @references Austin PC, Steyerberg EW. (2019) The Integrated Calibration Index (ICI) and related metrics for quantifying the calibration of logistic regression models. \emph{Statistics in Medicine}. 38, pp. 1–15. https://doi.org/10.1002/sim.8281
#' @references Van Hoorde, K., Van Huffel, S., Timmerman, D., Bourne, T., Van Calster, B. (2015). A spline-based tool to assess and visualize the calibration of multiclass risk predictions. \emph{Journal of Biomedical Informatics}, 54, pp. 283-93
#' @references Van Calster, B., Nieboer, D., Vergouwe, Y., De Cock, B., Pencina M., Steyerberg E.W. (2016). A calibration hierarchy for risk models was defined: from utopia to empirical data. \emph{Journal of Clinical Epidemiology}, 74, pp. 167-176
#'
#' @export
#'
#' @examples
#' p <- runif(100)
#' y <- rbinom(length(p), 1, p)
#' score_binary(y = y, p = p)
score_binary <- function(y, p, ...){

  # remove missing
  missy <- is.na(y); missp <- is.na(p)
  miss <- missy | missp
  if (any(miss)){
    message("score_binary: ", sum(miss), " cases removed for missing values.\n",
            sum(missy), " missing y; ", sum(missp), " missing p")
    y <- y[!miss]; p <- p[!miss]
  }

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
      if (length(calib_args[["eval"]]) == 1 && !calib_args[["eval"]] == 0){
        # warning("If specifying 'eval' in calib_args we strongly suggest specifying ",
        #         "a vector of probabilities at which to evaluate the calibration curve. ",
        #         "e.g., eval = seq(min(p), max(p), length.out=100), where p are predictions from ",
        #         "original model evaluated on the original data. Only specifying only the number ",
        #         "of points to eval (e.g., eval = 100) will likely cause issues/errors.")
        NULL # need to figure out somewhere else to put this warning as prints too many times
      }
    } else{
      calib_args[["eval"]] <- 0
    }
  } else{
    calib_args <- cal_defaults() # pmcalib_defaults
  }
  calib_args[["plot"]] <- FALSE

  # avoid error if resample returns constant y or p
  if (length(unique(y)) == 1 | length(unique(p)) == 1){
    warning("Returning NA scores as y and/or p is constant")

    logcalcoef <- c("Intercept" = NA_real_, "Slope" = NA_real_)
    brier <- c("Brier" = NA_real_)
    auc <- c("C" = NA_real_)
    # hacky(?) way to get current metrics from pmc
    # cal <- pmcalibration::pmcalibration(y=rbinom(100, 1, prob = .5),
    #                              p=runif(100), smooth = "none")
    cal <- pmcalibration::pmcalibration(y=c(0,0,0,1,1,1),
                                        p=c(.2,.2,.2,.6,.6,.6),
                                        smooth = "none", ci = "none", plot = FALSE)

    cal$metrics[1:(length(cal$metrics))] <- NA_real_

    if (length(calib_args[["eval"]]) > 1 || calib_args[["eval"]] != 0){
      if (length(calib_args[["eval"]]) > 1){
        cal_plot <- rep(NA_real_, times = length(calib_args[["eval"]]))
        names(cal_plot) <- paste0( "cal_plot_", calib_args[["eval"]] )
      } else{
        cal_plot <- rep(NA_real_, times = calib_args[["eval"]]) # not having names might cause issues?
      }
    } else{
      cal_plot <- c()
    }

  } else {
    calib_args[["y"]] <- y
    calib_args[["p"]] <- p

    cal <- do.call(pmcalibration::pmcalibration, calib_args)

    logcal <- pmcalibration::logistic_cal(y = y, p = p)
    CIL <- logcal$calibration_intercept$coefficients
    # return CIL?
    logcalcoef <- logcal$calibration_slope$coefficients
    names(logcalcoef) <- c("Intercept", "Slope")

    brier <- mean((y - p)^2)
    names(brier) <- "Brier"
    auc <- as.numeric(suppressMessages(pROC::auc(response = y, predictor = p)))
    names(auc) <- "C"

    if (!is.null(cal$plot$p_c_plot)){
      cal_plot <- cal$plot$p_c_plot
      # names(cal_plot) <- paste0( "cal_plot_", calib_args[["eval"]] )
      names(cal_plot) <- paste0( "cal_plot_", cal$plot$p )
    } else{
      cal_plot <- c()
    }
  }

  scores <- c(auc, brier, logcalcoef, cal$metrics, cal_plot)
  #names(scores) <- c("brier", "C")

  return(scores)
}
