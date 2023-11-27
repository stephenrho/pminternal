#' Get bias-corrected performance measures via bootstrapping or cross-validation
#'
#' @description
#' Performs internal validation of a prediction model development procedure via bootstrapping
#' or cross-validation. Many model types are supported via the \code{insight} and \code{marginaleffects}
#' packages or users can supply user-defined functions that implement the model development
#' procedure and retrieve predictions. Bias-corrected scores and estimates of optimism (where applicable)
#' are provided.
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
#' ".632", "cv_optimism", or "cv_average". See details.
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
#' Note that \code{validate} does very little to check for missing values. If \code{fit} is
#' supplied \code{insight::get_data} will extract the data used to fit the model and usually
#' this will result in complete cases being used. User-defined model and predict functions can
#' be specified to handle missing values among predictor variables. Currently any user supplied data will
#' have rows with missing outcome values removed.
#'
#' \bold{method}
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
#' Steyerberg et al. (2001).}
#' }
#'
#' @return an object of class internal_validate containing apparent and bias-corrected
#' estimates of performance scores. If method = "boot_*" it also contains results pertaining
#' to stability of predictions across bootstrapped models (see Riley and Collins, 2023).
#'
#' @references Steyerberg, E. W., Harrell Jr, F. E., Borsboom, G. J., Eijkemans, M. J. C., Vergouwe, Y., & Habbema, J. D. F. (2001). Internal validation of predictive models: efficiency of some procedures for logistic regression analysis. Journal of clinical epidemiology, 54(8), 774-781.
#' @references Harrell Jr F. E. (2015). Regression Modeling Strategies: with applications to linear models, logistic and ordinal regression, and survival analysis. New York: Springer Science, LLC.
#' @references Efron (1983). “Estimating the error rate of a prediction rule: improvement on cross-validation”. Journal of the American Statistical Association, 78(382):316-331
#' @references Van Calster, B., Steyerberg, E. W., Wynants, L., and van Smeden, M. (2023). There is no such thing as a validated prediction model. BMC medicine, 21(1), 70.
#' @references Riley RD, Collins GS. (2023). Stability of clinical prediction models developed using statistical or machine learning methods. Biom J. doi:10.1002/bimj.202200302. Epub ahead of print.
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
                              ".632", "cv_optimism", "cv_average"),
                     data,
                     outcome,
                     model_fun,
                     pred_fun,
                     score_fun,
                     B,
                     ...){
  # TODO:
  # pbapply
  # add plot arg to validate to automatically include stuff to make cal_plot
  # make stability plots more flexible?

  call <- match.call()
  method <- match.arg(method)

  if (missing(B)){
    if (method %in% c("boot_optimism", "boot_simple", ".632")){
      B <- 200
    } else{
      B <- 10
    }
  }

  if (missing(fit) & any(missing(data), missing(outcome),
                         missing(model_fun), missing(pred_fun))){
    stop("if fit not provided, data, outcome, model_fun, and pred_fun must be specified")
  }

  if (!missing(fit)){
    if (any(!missing(data), !missing(outcome),
            !missing(model_fun), !missing(pred_fun))){
      warning("if fit is specified, data, outcome, model_fun, and pred_fun are ignored")
    }
    # check class
    if (isFALSE(class(fit)[1] %in% insight::supported_models())){
      stop("fit provided is not a class supported - see insight::supported_models()")
    }
    data <- insight::get_data(fit)
    outcome <- insight::find_response(fit)
    mcall <- insight::get_call(fit)

    model_fun <- function(data, ...){
      # see marginaleffects:::bootstrap_boot
      mcall[["data"]] <- data
      out <- eval(mcall)
      eval(out)
    }

    pred_fun <- function(model, data, ...){
      #insight::get_predicted(x = model, data = data, ci=NULL)
      # metd <- marginaleffects:::type_dictionary_build()
      # metd <- metd[!duplicated(metd$class), ]
      # if (class(model)[1] %in% metd$class){
      #   type <- metd$type[which(class(model)[1] == metd$class)]
      # } else{
      #   type <- "response"
      # }
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

  # test
  fit <- model_fun(data=data, ...)
  p_app <- pred_fun(model = fit, data = data, ...)
  score_app <- score_fun(y = data[[outcome]], p = p_app, ...)

  if ( method %in% c("boot_optimism", "boot_simple", ".632") ){
    if (B < 200) message("It is recommended that B >= 200 for bootstrap validation")
    m <- if (method == ".632") ".632" else "boot"
    res <- boot_optimism(data = data, outcome = outcome,
                         model_fun = model_fun, pred_fun = pred_fun, score_fun = score_fun,
                         method = m, B = B, ...)
  } else if ( method %in% c("cv_optimism", "cv_average") ){
    res <- crossval(data = data, outcome = outcome,
                    model_fun = model_fun, pred_fun = pred_fun,
                    score_fun = score_fun, k = B, ...)
  }

  apparent <- res$apparent
  if (method %in% c("boot_optimism", ".632", "cv_optimism")){
    optimism <- res$optimism
    corrected <- res$corrected
  } else{
    optimism <- NA
    if (method == "boot_simple"){
      corrected <- res$simple
    } else if (method == "cv_average"){
      corrected <- res$cv_average
    }
  }

  out <- list(apparent = apparent,
              optimism = optimism,
              corrected = corrected,
              scores = res,
              method = method,
              call = call,
              B = B,
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
#' @return a data.frame with 3 rows (apparent score, optimism, bias-corrected score)
#' and one column per score. Not all methods produce an optimism estimate so this row may be all NA.
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

  out <- rbind(object$apparent[i], object$optimism[i], object$corrected[i])
  rownames(out) <- c("Apparent", "Optimism", "Corrected")

  out <- data.frame(out)
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
  print(format(round(x, 6), digits = digits, scientific=F))
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
