#' Calculate bias-corrected scores via cross-validation
#'
#' @description
#' Estimate bias-corrected scores via cross-validation. CV is used to calculate optimism
#' which is then subtracted from apparent scores and to calculate average performance in the
#' out of sample (held out) data.
#' This function is called by \code{\link{validate}}.
#'
#' @param data the data used in developing the model. Should contain all variables considered (i.e., even those excluded by variable selection in the development sample)
#' @param outcome character denoting the column name of the outcome in \code{data}.
#' @param model_fun a function that takes at least one argument, \code{data}. This function should implement the entire model development procedure (i.e., hyperparameter tuning, variable selection, imputation). Additional arguments can be provided via \code{...}. This function should return an object that works with \code{pred_fun}.
#' @param pred_fun function that takes at least two arguments, \code{model} and \code{data}. This function should return a numeric vector of predicted probabilities of the outcome with the same length as the number of rows in \code{data} so it is important to take into account how missing data is treated (e.g., \code{predict.glm} omits predictions for rows with missing values).
#' @param score_fun a function to calculate the metrics of interest. If this is not specified \code{\link{score_binary}} is used.
#' @param k number of folds. Typically scores need >> 2 observations to be calculated so folds should be chosen with this in mind.
#' @param ... additional arguments for \code{model_fun}, \code{pred_fun}, and/or \code{score_fun}.
#'
#' @return a list of class \code{internal_cv} containing:
#' \itemize{
#' \item{\code{apparent} - scores calculated on the original data using the original model.}
#' \item{\code{optimism} - estimates of optimism for each score (average difference in score for training data vs test data on each fold) which can be subtracted from 'apparent' performance calculated using the original model on the original data.}
#' \item{\code{cv_optimism_corrected} - 'bias corrected' scores (apparent - optimism). This is what is produced by \code{rms::validate}, \code{rms::predab.resample}.}
#' \item{\code{cv_average} - average of scores calculated on the test (held out) data. This is the metric described in Steyerberg et al. (2001).}
#' \item{\code{indices} - indices used to define test set on each fold.}
#' }
#'
#' @references Steyerberg, E. W., Harrell Jr, F. E., Borsboom, G. J., Eijkemans, M. J. C., Vergouwe, Y., & Habbema, J. D. F. (2001). Internal validation of predictive models: efficiency of some procedures for logistic regression analysis. Journal of clinical epidemiology, 54(8), 774-781.
#'
#' @export
#' @examples
#' library(pminternal)
#' set.seed(456)
#' # simulate data with two predictors that interact
#' dat <- pmcalibration::sim_dat(N = 1000, a1 = -2, a3 = -.3)
#' mean(dat$y)
#' dat$LP <- NULL # remove linear predictor
#'
#' # fit a (misspecified) logistic regression model
#' #m1 <- glm(y ~ x1 + x2, data=dat, family="binomial")
#'
#' model_fun <- function(data, ...){
#'   glm(y ~ x1 + x2, data=data, family="binomial")
#' }
#'
#' pred_fun <- function(model, data, ...){
#'   predict(model, newdata=data, type="response")
#' }
#'
#' # CV Corrected = Apparent - CV Optimism
#' # CV Average = average score in held out fold
#' crossval(data=dat, outcome="y", model_fun=model_fun, pred_fun=pred_fun, k=10)
#'
crossval <- function(data, outcome,
                     model_fun, pred_fun, score_fun,
                     k=10, ...){

  dots <- list(...)

  if (missing(score_fun)){
    score_fun <- score_binary
  }

  # apparent
  fit <- model_fun(data=data, ...)
  p_app <- pred_fun(model = fit, data = data, ...)
  y <- data[[outcome]]
  score_app <- score_fun(y = y, p = p_app, ...)

  n <- nrow(data)

  # k=25
  # n = 1234
  nperfold <- n/k
  if (nperfold < 2) stop("Number of holdouts per fold is less than 2. Decrease the number of folds.")

  indices <- lapply(seq(k), function(i){
    start <- 1 + round((i - 1)*nperfold)
    stop <- min(n, round(i*nperfold))
    #stop <- min(n, round(start + nperfold - 1)) # results in some overlapping indices
    start:stop
  })

  # (unlist(indices) |> length()) == n
  # (unlist(indices) |> max()) == n
  # (unlist(indices) |> unique() |> length()) == n
  # (sapply(indices, length) |> mean()) == nperfold

  S <- lapply(indices, function(i){
    data_train <- data[-i, ]
    data_test <- data[i, ]

    if (all(data_train[[outcome]] == 0, na.rm = TRUE) |
        all(data_train[[outcome]] == 1, na.rm = TRUE) |
        all(data_test[[outcome]] == 0, na.rm = TRUE) |
        all(data_test[[outcome]] == 1, na.rm = TRUE)
        ){
      stop("crossval error: one fold had all 0 or 1 in outcomes. ",
           "Try fewer folds.")
    }

    fit <- model_fun(data = data_train, ...)

    p_test <- pred_fun(model = fit, data = data_test, ...)
    score_test <- score_fun(y = data_test[[outcome]], p = p_test, ...)

    p_train <- pred_fun(model = fit, data = data_train, ...)
    score_train <- score_fun(y = data_train[[outcome]], p = p_train, ...)

    list("score_test" = score_test,
         "optimism" = score_train - score_test)
  })

  cv_average <- do.call(rbind, lapply(S, function(x) x$score_test))
  cv_average <- apply(cv_average, 2, mean)

  opt <- do.call(rbind, lapply(S, function(x) x$optimism))
  opt <- apply(opt, 2, mean)

  bcorr <- score_app - opt

  out <- list("apparent" = score_app,
              "optimism" = opt,
              "corrected" = bcorr,
              "cv_average" = cv_average,
              "indices" = indices)

  class(out) <- "internal_cv"

  return(out)
}

#' Print a internal_cv object
#'
#' @param x an object created with \code{crossval}
#' @param digits number of digits to print (default = 2)
#' @param ... additional arguments to print
#'
#' @return invisibly returns \code{x} and prints estimates to console
#' @export
print.internal_cv <- function(x, digits=2, ...){
  out <- rbind(x$apparent, x$optimism, x$corrected, x$cv_average)
  rownames(out) <- c("Apparent", "CV Optimism", "CV Corrected", "CV Average")
  print(out, digits = digits, ...)
  return(invisible(x))
}
