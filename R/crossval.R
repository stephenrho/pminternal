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
#' @param k number of folds. Typically scores need greater than 2 observations to be calculated so folds should be chosen with this in mind.
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

  # - deal with model fit fails
  mf <- function(data, ...){
    tryCatch(expr = model_fun(data=data, ...),
             error = function(e) {
               message("model fit failed on resample: ", e)
               NULL
             }
    )
  }

  ef <- function(model, data, outcome, ...){
    # eval fit and return NA scores if model fit failed
    if (is.null(model)){
      # hacky way to get names...
      #score <- score_fun(y = sample(x = 0:1, size = 1000, replace = TRUE), p = runif(1000, .01, .99), ...)
      score <- suppressWarnings(score_fun(y = rep(0, 4), p = rep(0, 4), ...))

      score[1:length(score)] <- NA_real_
      p <- rep(NA_real_, times = length(data[[outcome]]))
    } else{
      p <- pred_fun(model = fit, data = data, ...)
      y <- data[[outcome]]
      score <- score_fun(y = y, p = p, ...)
    }
    return(list(score=score, p=p))
  }

  # apparent
  fit <- mf(data=data, ...) #model_fun(data=data, ...)
  if (is.null(fit)) stop("Model fit failed on assessment of apparent performance")

  # p_app <- pred_fun(model = fit, data = data, ...)
  y <- data[[outcome]]
  score_app <- ef(model = fit, data = data, outcome = outcome, ...) #score_fun(y = y, p = p_app, ...)
  p_app <- score_app$p
  score_app <- score_app$score

  # fit <- model_fun(data=data, ...)
  # p_app <- pred_fun(model = fit, data = data, ...)
  # y <- data[[outcome]]
  # score_app <- score_fun(y = y, p = p_app, ...)

  n <- nrow(data)

  # set up folds
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

  fold <- function(i){
    data_train <- data[-i, ]
    data_test <- data[i, ]

    fit <- mf(data = data_train, ...)
    score_test <- ef(model = fit, data = data_test, outcome = outcome, ...)$score
    score_train <- ef(model = fit, data = data_train, outcome = outcome, ...)$score

    list("score_test" = score_test,
         "optimism" = score_train - score_test)
  }

  qfold <- purrr::quietly(fold)

  S <- lapply(indices, FUN = qfold)

  # extrct warnings/messages
  warns <- lapply(S, function(x) x$warnings)
  mess <- lapply(S, function(x) x$messages)

  if (any(sapply(warns, length) > 0)){
    w <- unique(unlist(warns))
    nw <- sapply(w, function(x){
      sum(sapply(warns, function(xx) any(xx == x)))
    })
    warning(paste0("The following warnings occurred during the call to validate\n\n",
                   paste0(seq(w), ": ", w, " (", nw, " occurrences)", collapse = "\n") ))
  }

  if (any(sapply(mess, length) > 0)){
    w <- unique(unlist(mess))
    nw <- sapply(w, function(x){
      sum(sapply(mess, function(xx) any(xx == x)))
    })
    message(paste0("The following messages occurred during the call to validate\n\n",
                   paste0(seq(w), ": ", w, " (", nw, " occurrences)", collapse = "\n") ))
  }

  # output
  cv_average <- do.call(rbind, lapply(S, function(x) x$result$score_test))
  n_cv_average <- apply(cv_average, 2, function(x) sum(!is.na(x)))
  cv_average <- apply(cv_average, 2, mean, na.rm=TRUE)

  opt <- do.call(rbind, lapply(S, function(x) x$result$optimism))
  nopt <- apply(opt, 2, function(x) sum(!is.na(x)))
  opt <- apply(opt, 2, mean, na.rm=TRUE)

  bcorr <- score_app - opt

  out <- list("apparent" = score_app,
              "optimism" = opt,
              "n_optimism" = nopt,
              "corrected" = bcorr,
              "cv_average" = cv_average,
              "n_cv_average" = n_cv_average,
              "indices" = indices,
              "warnings" = warns,
              "messages" = mess)

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
  out <- rbind(x$apparent, x$optimism, x$n_optimism, x$corrected,
               x$cv_average, x$n_cv_average)
  rownames(out) <- c("Apparent", "CV Optimism", "B Optimism", "Optimism Corrected",
                     "CV Average", "B Average")
  print(out, digits = digits, ...)
  return(invisible(x))
}
