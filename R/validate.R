#' Get bias-corrected performance measures via bootstrapping or cross-validation
#'
#' @description
#' Performs internal validation of a prediction model development procedure via bootstrapping
#' or cross-validation. Many model types are supported via the \code{insight} and \code{marginaleffects}
#' packages or users can supply user-defined functions that implement the model development
#' procedure and retrieve predictions. Bias-corrected scores and estimates of optimism (where applicable)
#' are provided. See \code{\link{confint.internal_validate}} for calculation of confidence intervals.
#'
#' @param fit a model object. If fit is given the \code{insight} package is
#' used to extract data, outcome, and original model call. Therefore, it is important
#' that fit be supported by \code{insight} and implements the entire model development
#' process (see Harrell 2015). A fit given after selection of variables by some method will
#' not give accurate bias-correction. Model predictions are obtained via
#' \code{marginaleffects::get_predict} with type = "response" so fit should be
#' compatible with this function. If fit is provided the arguments data, outcome,
#' model_fun, and pred_fun are all ignored.
#' @param method bias-correction method. Valid options are "boot_optimism", "boot_simple",
#' ".632", "cv_optimism", "cv_average", or "none" (return apparent performance). See details.
#' @param data a data.frame containing data used to fit development model
#' @param outcome character denoting the column name of the outcome in data
#' @param model_fun for models that cannot be supplied via fit this should be a function
#' that takes one named argument: 'data' (function should include ... among arguments).
#' This function should implement the entire model development
#' procedure (hyperparameter tuning, variable selection, imputation etc) and return an object
#' that can be used by pred_fun. Additional arguments can be supplied by ...
#' @param pred_fun for models that cannot be supplied via fit this should be a function
#' that takes two named arguments: 'model' and 'data' (function should include ... among arguments).
#' 'model' is an object returned by model_fun.
#' The function should return a vector of predicted risk probabilities of the same length as the number
#' of rows in data. Additional arguments can be supplied by ...
#' @param score_fun function used to produce performance measures from predicted risks
#' and observed binary outcome. Should take two named arguments: 'y' and 'p' (function should include ... among arguments).
#' This function should return a named vector of scores. If unspecified \code{\link{score_binary}}
#' is used and this should be good for most purposes.
#' @param B number of bootstrap replicates or crossvalidation folds. If unspecified B is set to
#' 200 for method = "boot_\*"/".632", or is set to 10 for method = "cv_\*".
#' @param ... additional arguments for user-defined functions. Arguments for
#' producing calibration curves can be set via 'calib_args' which should be
#' a named list (see \code{\link{cal_plot}} and \code{\link{score_binary}}).
#' For method = "boot_optimism", "boot_simple", or ".632" users can specify a
#' \code{cores} argument (e.g., \code{cores = 4}) to run bootstrap samples in parallel.
#'
#' @details
#' Internal validation can provide bias-corrected estimates of performance (e.g., C-statistic/AUC)
#' for a model development procedure (i.e., expected performance if the same procedure were applied
#' to another sample of the same size from the same population; see references). There are several approaches to producing
#' bias-corrected estimates (see below). It is important that the fit or model_fun provided implement
#' the entire model development procedure, including any hyperparameter tuning and/or variable selection.
#'
#' Note that \code{\link{validate}} does very little to check for missing values in predictors/features. If \code{fit} is
#' supplied \code{insight::get_data} will extract the data used to fit the model and usually
#' this will result in complete cases being used. User-defined model and predict functions can
#' be specified to handle missing values among predictor variables. Currently any user supplied data will
#' have rows with missing outcome values removed.
#'
#' \bold{method}
#' Different options for the method argument are described below:
#' \describe{
#' \item{boot_optimism}{ (default) estimates optimism for each score and subtracts from apparent score (score calculated
#' with the original/development model evaluated on the original sample). A new model is fit using the same procedure
#' using each bootstrap resample. Scores are calculated when applying the boot model to the boot sample (\eqn{S_{boot}})
#' and the original sample (\eqn{S_{orig}}) and the difference gives an estimate of optimism for a given resample (\eqn{S_{boot} - S_{orig}}).
#' The average optimism across the B resamples is subtracted from the apparent score to produce the bias corrected score.}
#' \item{boot_simple}{implements the simple bootstrap. B bootstrap models are fit and evaluated on the original data.
#' The average score across the B replicates is the bias-corrected score.}
#' \item{.632}{implements Harrell's adaption of Efron's .632 estimator for binary outcomes
#' (see rms::predab.resample and rms::validate). In this case the estimate of optimism is
#' \eqn{0.632 \times (S_{app} - mean(S_{omit} \times w))} where \eqn{S_{app}} is the apparent performance
#' score and \eqn{S_{omit}} is the score estimated using the bootstrap model evaluated on the out-of-sample
#' observations and \eqn{w} weights for the proportion of observations omitted (see Harrell 2015, p. 115).}
#' \item{cv_optimism}{estimate optimism via B-fold crossvalidation. Optimism is the average of the difference
#' in performance measure between predictions made on the training vs test (held out fold) data. This is the approach
#' implemented in \code{rms::validate} with method="crossvalidation".}
#' \item{cv_average}{bias corrected scores are the average of scores calculated by assessing the model developed on each
#' fold evaluated on the test/held out data. This approach is described and compared to "boot_optimism" and ".632" in
#' Steyerberg et al. (2001).}}
#'
#' \bold{Calibration curves}
#' To make calibration curves and calculate the associated estimates (ICI, ECI, etc - see \code{\link{score_binary}})
#' \code{validate} uses the default arguments in \code{\link{cal_defaults}}. These arguments are passed to the \code{pmcalibration} package
#' (see \code{?pmcalibration::pmcalibration} for options).
#'
#' If a calibration plot (apparent vs bias corrected calibration curves via \code{\link{cal_plot}})
#' is desired, the argument 'eval' should be provided. This should be the points at which to evaluate
#' the calibration curve on each boot resample or crossvalidation fold. A good option would be
#' \code{calib_args = list(eval = seq(min(p), max(p), length.out=100))}; where p are predictions from the
#' original model evaluated on the original data.
#'
#' \bold{Number of resamples/folds is less than requested}
#' If the \code{model_fun} produces an error or if \code{score_binary} is supplied with constant predictions
#' or outcomes (e.g. all(y == 0)) the returned scores will all be NA. These will be omitted from the calculation
#' of optimism or other bias-corrected estimates (cv_average, boot_simple) and the number of successful resamples/folds
#' will be < B. \code{validate} collects errors and will produce a warning summarizing them. The number of successful
#' samples is given in the 'n' column in the printed summary of an 'internal_validate' object.
#'
#' It is important to understand what is causing the loss of resamples/folds. Some potential sources (which will need to be added to) are that
#' for rare events the resamples/folds may be resulting in samples that have zero outcomes. For 'cv_*' this will especially
#' be the case if B (n folds) is set high. There may be problems with factor/binary predictor variables with rare levels, which could be dealt with
#' by specifying a \code{model_fun} that omits variables for the model formula if only one level is present. The issue may be related to the construction
#' of calibration curves and may be addressed by more carefully selecting settings (see section above).
#'
#' @return an object of class internal_validate containing apparent and bias-corrected
#' estimates of performance scores. If method = "boot_*" it also contains results pertaining
#' to stability of predictions across bootstrapped models (see Riley and Collins, 2023).
#'
#' @references Steyerberg, E. W., Harrell Jr, F. E., Borsboom, G. J., Eijkemans, M. J. C., Vergouwe, Y., & Habbema, J. D. F. (2001). Internal validation of predictive models: efficiency of some procedures for logistic regression analysis. Journal of clinical epidemiology, 54(8), 774-781.
#' @references Harrell Jr F. E. (2015). Regression Modeling Strategies: with applications to linear models, logistic and ordinal regression, and survival analysis. New York: Springer Science, LLC.
#' @references Efron (1983). “Estimating the error rate of a prediction rule: improvement on cross-validation”. Journal of the American Statistical Association, 78(382):316-331
#' @references Van Calster, B., Steyerberg, E. W., Wynants, L., and van Smeden, M. (2023). There is no such thing as a validated prediction model. BMC medicine, 21(1), 70.
#' @references Riley, R. D., & Collins, G. S. (2023). Stability of clinical prediction models developed using statistical or machine learning methods. Biometrical Journal, 65(8), 2200302. doi:10.1002/bimj.202200302
#'
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
#' m1 <- glm(y ~ ., data=dat, family="binomial")
#'
#' # internal validation of m1 via bootstrap optimism with 10 resamples
#' # B = 10 for example but should be >= 200 in practice
#' m1_iv <- validate(m1, method="boot_optimism", B=10)
#' m1_iv
#'
validate <- function(fit,
                     method=c("boot_optimism", "boot_simple",
                              ".632", "cv_optimism", "cv_average", "none"),
                     data,
                     outcome,
                     model_fun,
                     pred_fun,
                     score_fun,
                     B,
                     ...){

  call <- match.call()
  method <- match.arg(method)

  if (missing(B)){
    if (method %in% c("boot_optimism", "boot_simple", ".632")){
      B <- 200
    } else{
      B <- 10
    }
  }

  if (method == "none") B <- NA

  if (missing(fit) & any(missing(data), missing(outcome),
                         missing(model_fun), missing(pred_fun))){
    stop("if fit not provided, data, outcome, model_fun, and pred_fun must be specified")
  }

  if (!missing(fit)){
    if (any(#!missing(data), !missing(outcome),
      !missing(model_fun), !missing(pred_fun))){
      warning("If fit is specified, model_fun and pred_fun are ignored")
    }
    # check class
    if (isFALSE(class(fit)[1] %in% insight::supported_models())){
      stop("fit provided is not a class supported - see insight::supported_models()")
    }
    if (missing(data)){
      data <- insight::get_data(fit)
    }
    if (missing(outcome)){
      outcome <- insight::find_response(fit)
    }

    mcall <- insight::get_call(fit)

    model_fun <- function(data, ...){
      # see marginaleffects:::bootstrap_boot
      mcall[["data"]] <- data
      out <- eval(mcall)
      eval(out)
    }

    pred_fun <- function(model, data, ...){
      if (is(model, "lrm")){
        type <- "fitted"
      } else{
        type <- "response"
      }

      marginaleffects::get_predict(model = model,
                                   newdata = data,
                                   type = type)$estimate
    }
  }

  if (missing(score_fun)){
    score_fun <- score_binary
  }

  missy <- is.na(data[[outcome]])
  if ( any(missy) ){
    message("Removing ", sum(missy), " rows with missing outcome")
    data <- data[which(!missy), ]
  }

  if (method != "none"){
    if ( method %in% c("boot_optimism", "boot_simple", ".632") ){
      if (B < 200) message("It is recommended that B >= 200 for bootstrap validation")
      m <- if (method == ".632") ".632" else "boot"
      res <- boot_optimism(data = data,
                           outcome = outcome,
                           model_fun = model_fun,
                           pred_fun = pred_fun,
                           score_fun = score_fun,
                           method = m,
                           B = B, ...)
    } else if ( method %in% c("cv_optimism", "cv_average") ){
      res <- crossval(data = data,
                      outcome = outcome,
                      model_fun = model_fun,
                      pred_fun = pred_fun,
                      score_fun = score_fun,
                      k = B, ...)
    }

    apparent <- res$apparent
    if (method %in% c("boot_optimism", ".632", "cv_optimism")){
      optimism <- res$optimism
      corrected <- res$corrected
      B_actual <- res$n_optimism
    } else{
      optimism <- NA
      if (method == "boot_simple"){
        corrected <- res$simple
        B_actual <- res$n_simple
      } else if (method == "cv_average"){
        corrected <- res$cv_average
        B_actual <- res$n_cv_average
      }
    }

    if (any(B_actual < B)) {
      message(strwrap("Note there were some unsuccessful resamples (see n column in
                    summary). It is worth trying to understand why these unsuccessful
                    runs are happening. validate will print the warnings and messages
                    encountered.",
                      prefix = " ", initial = ""),
              "\n\nThe greater the proportion of samples lost the more important it is to figure out the source and try to address. See ?validate details for potential sources.")
    }
  } else{
    fit <- model_fun(data=data, ...)
    p_app <- pred_fun(model = fit, data = data, ...)
    apparent <- score_fun(y = data[[outcome]], p = p_app, ...)
    optimism <- B_actual <- corrected <-  rep(NA_real_, times = length(apparent))

    res <- NULL
  }

  out <- list(apparent = apparent,
              optimism = optimism,
              corrected = corrected,
              scores = res,
              method = method,
              call = call,
              B = B,
              B_actual = B_actual,
              fit = if (missing(fit)) NULL else fit,
              data = data,
              outcome = outcome,
              model_fun = model_fun,
              pred_fun = pred_fun,
              score_fun = score_fun,
              dots = list(...))

  class(out) <- "internal_validate"
  return(out)
}

#' Summarize a internal_validate object
#'
#' @param object created by call to \code{\link{validate}}
#' @param ignore_scores a string used to identify scores to omit from summary.
#' \code{\link{score_binary}} produces scores with prefix 'cal_plot' when a calibration plot
#' is desired (see \code{\link{cal_plot}}) and these are ignored by default.
#' @param ... ignored
#'
#' @return a data.frame with 4 columns (apparent score, optimism, bias-corrected score, number of successful resamples/folds)
#' and one row per score. Not all methods produce an optimism estimate so this row may be all NA. If confidence intervals
#' have been added for all measures via \code{\link{confint.internal_validate}}, two additional columns containing lower and upper bounds for
#' bias-corrected performance.
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
#' # internal validation of m1 via bootstrap optimism with 10 resamples
#' # B = 10 for example but should be >= 200 in practice
#' m1_iv <- validate(m1, method="boot_optimism", B=10)
#' summary(m1_iv)
#'
summary.internal_validate <- function(object, ignore_scores="^cal_plot", ...){
  i <- !grepl(ignore_scores, names(object$apparent))

  if ("confint" %in% names(object)){
    #if ("confint" %in% names(object) & all( names(object$corrected[i]) %in% rownames(object$confint$ci$Corrected) )){

    out <- data.frame(cbind(object$apparent[i], object$optimism[i],
                            object$corrected[i], object$B_actual[i]))
    out <- merge(x=out, y=data.frame(object$confint$ci$Corrected),
                 by="row.names", all.x=TRUE)
    row.names(out) <- out$Row.names
    out$Row.names <- NULL
    colnames(out) <- c("apparent", "optimism", "corrected", "n", "corrected_lower", "corrected_upper")

  } else{
    out <- data.frame(cbind(object$apparent[i], object$optimism[i],
                            object$corrected[i], object$B_actual[i]))
    colnames(out) <- c("apparent", "optimism", "corrected", "n")
  }

  class(out) <- c("internal_validatesummary", "data.frame")

  return(out)
}

#' Print summary of internal_validate object
#'
#' @param x a \code{internal_validatesummary} object
#' @param digits number of digits to print
#' @param ... ignored
#'
#' @returns invisible(x) - prints a summary
#' @rdname print.internal_validatesummary
#' @export
print.internal_validatesummary <- function(x, digits = 2, ...){
  #print.data.frame(x, digits = digits, scientific=F)
  print(format(round(x, 6), digits = digits, scientific=FALSE))
  return(invisible(x))
}

#' print a internal_validate object
#'
#' @param x a \code{internal_validate} object
#' @param digits number of digits to print
#' @param ... optional arguments passed to print
#'
#' @returns prints a summary
#'
#' @rdname print.internal_validate
#' @export
print.internal_validate <- function(x, digits = 2, ...) {
  print(summary(x), digits = digits, ...)
}
