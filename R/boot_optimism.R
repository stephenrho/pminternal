#' Calculate optimism and bias-corrected scores via bootstrap resampling
#'
#' @description
#' Estimate bias-corrected scores via calculation of bootstrap optimism (standard or .632).
#' Can also produce estimates for assessing the stability of prediction model predictions.
#' This function is called by \code{\link{validate}}.
#'
#' @param data the data used in developing the model. Should contain all variables considered (i.e., even those excluded by variable selection in the development sample)
#' @param outcome character denoting the column name of the outcome in \code{data}.
#' @param model_fun a function that takes at least one argument, \code{data}. This function should implement the entire model development procedure (i.e., hyperparameter tuning, variable selection, imputation). Additional arguments can be provided via \code{...}. This function should return an object that works with \code{pred_fun}.
#' @param pred_fun function that takes at least two arguments, \code{model} and \code{data}. This function should return a numeric vector of predicted probabilities of the outcome with the same length as the number of rows in \code{data} so it is important to take into account how missing data is treated (e.g., \code{predict.glm} omits predictions for rows with missing values). see \code{vignette("missing-data", package = "pminternal")}.
#' @param score_fun a function to calculate the metrics of interest. If this is not specified \code{\link{score_binary}} is used.
#' @param method "boot" or ".632". The former estimates bootstrap optimism for each score and subtracts
#' from apparent scores (simple bootstrap estimates are also produced as a by-product).
#' The latter estimates ".632" optimism as described in Harrell (2015). See \code{\link{validate}} details.
#' @param B number of bootstrap resamples to run
#' @param ... additional arguments for \code{model_fun}, \code{pred_fun}, and/or \code{score_fun}.
#'
#' @return a list of class \code{internal_boot} containing:
#' \itemize{
#' \item{\code{apparent} - scores calculated on the original data using the original model.}
#' \item{\code{optimism} - estimates of optimism for each score (average difference in score for bootstrap models evaluated on bootstrap vs original sample) which can be subtracted from 'apparent' performance calculated using the original model on the original data.}
#' \item{\code{corrected} - 'bias corrected' scores (apparent - optimism)}
#' \item{\code{simple} - if method = "boot", estimates of scores derived from the 'simple bootstrap'. This is the average of each score calculated from the bootstrap models evaluated on the original outcome data. NULL if method = ".632"}
#' \item{\code{stability} - if method = "boot", a N,(B+1) matrix where N is the number of observations in \code{data} and \code{B} is the number of bootstrap samples. The first column contains the original predictions and each of subsequent B columns contain the predicted probabilities of the outcome from each bootstrap model evaluated on the original data. There may be fewer than B+1 columns if errors occur during resamples (when model_fun throws an error all scores are NA). NULL if method = ".632"}
#' }
#'
#' @references Steyerberg, E. W., Harrell Jr, F. E., Borsboom, G. J., Eijkemans, M. J. C., Vergouwe, Y., & Habbema, J. D. F. (2001). Internal validation of predictive models: efficiency of some procedures for logistic regression analysis. Journal of clinical epidemiology, 54(8), 774-781.
#' @references Harrell Jr F. E. (2015). Regression Modeling Strategies: with applications to linear models, logistic and ordinal regression, and survival analysis. New York: Springer Science, LLC.
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
#' model_fun <- function(data, ...){
#'   glm(y ~ x1 + x2, data=data, family="binomial")
#' }
#'
#' pred_fun <- function(model, data, ...){
#'   predict(model, newdata=data, type="response")
#' }
#'
#' boot_optimism(data=dat, outcome="y", model_fun=model_fun, pred_fun=pred_fun,
#'               method="boot", B=20) # B set to 20 for example but should be >= 200
#'
boot_optimism <- function(data, outcome,
                          model_fun, pred_fun, score_fun,
                          method=c("boot", ".632"), B=200, ...){

  dots <- list(...)

  method <- match.arg(method)

  if (missing(score_fun)){
    score_fun <- score_binary
  }

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
      # score <- score_fun(y = sample(x = 0:1, size = 1000, replace = TRUE), p = runif(1000, .01, .99), ...)
      score <- suppressWarnings(score_fun(y = rep(0, 4), p = rep(0, 4), ...))

      score[1:length(score)] <- NA_real_
      p <- rep(NA_real_, times = length(data[[outcome]]))
    } else{
      p <- pred_fun(model = model, data = data, ...)
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

  # see rms::predab.resample
  n <- nrow(data)
  indices <- matrix(integer(1), nrow=n, ncol = B)
  W <- matrix(TRUE, nrow=n, ncol = B)
  for (i in seq(B)){
    indices[, i] <- s <- sample(n, replace = TRUE)
    W[s, i] <- FALSE # used for method = .632
  }

  if (method == ".632"){
    nomit <- apply(W, 1, sum)
    if (any(nomit == 0)){
      stop("not every observation omitted at least once ",
           "in bootstrap samples.\nRe--run with larger B")
    }
    wt <- 1 - (1/B - apply(W/nomit, 2, sum)/n)
    #wt <- 1 + (apply(W/nomit, 2, sum)/n - 1/B)
  } else wt <- NULL

  # get cores
  if ("cores" %in% names(dots)){
    cores <- dots[["cores"]]
  } else{
    cores <- 1
  }

  boot <- function(i){
    # resample
    #data_i <- data[sample(x = nrow(data), replace = T), ]
    data_i <- data[indices[, i], ]
    # fit model on bootstrap resample
    model_i <- mf(data = data_i, ...) # model_fun

    if (method == "boot"){
      # # eval bootstrap model on boot and original data
      # p_orig <- pred_fun(model = model_i, data = data, ...)
      # p_boot <- pred_fun(model = model_i, data = data_i, ...)
      # # calculate scores for boot model eval'd on orig and boot data
      # score_orig <- score_fun(y = data[[outcome]], p = p_orig, ...)
      # score_boot <- score_fun(y = data_i[[outcome]], p = p_boot, ...)

      score_orig <- ef(model = model_i, data = data, outcome = outcome, ...)
      p_orig <- score_orig$p
      score_orig <- score_orig$score
      score_boot <- ef(model = model_i, data = data_i, outcome = outcome, ...)$score
      optimism <- score_boot - score_orig
    } else {
      p_orig <- score_orig <-  NULL
      # evaluate model on the left out indices
      data_omit_i <- data[-indices[, i], ]
      # p_omit <- pred_fun(model = model_i, data = data_omit_i, ...)
      # score_omit <- score_fun(y = data_omit_i[[outcome]], p = p_omit, ...)
      score_omit <- ef(model = model_i, data = data_omit_i, outcome = outcome, ...)$score
      optimism <- .632*(score_app - score_omit*wt[i])
    }

    list("optimism" = optimism,
         "p_orig" = p_orig,
         #"p_boot" = p_boot ,
         "score_orig" = score_orig#,
         #"score_boot" = score_boot
    )
  }

  qboot <- purrr::quietly(boot)

  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    # make all loaded packages available
    pkgs <- .packages()

    varlist <- c("B", "data", "indices", "wt", "method",
                 "model_fun", "pred_fun", "score_fun",
                 "mf", "ef", "boot", "qboot", "pkgs")

    if ("export" %in% names(dots)) varlist <- c(varlist, dots[["export"]])

    parallel::clusterExport(cl, varlist = varlist,
                            envir = environment())
    # load packages
    parallel::clusterEvalQ(cl, expr = sapply(pkgs, library, character.only = TRUE))
  } else {
    cl <- NULL # sequential eval
  }
  S <- pbapply::pblapply(seq(B), FUN = qboot, cl = cl)
  if (!is.null(cl)) {
    parallel::stopCluster(cl)
  }
  #closeAllConnections()

  # extrct warnings/messages
  # extract_warningsmessages(S, "warnings")
  # extract_warningsmessages(S, "messages")
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

  # make output
  opt <- do.call(rbind, lapply(S, function(x) x$result$optimism))
  nopt <- apply(opt, 2, function(x) sum(!is.na(x)))
  opt <- apply(opt, 2, mean, na.rm=TRUE)

  bcorr <- score_app - opt

  if (method == "boot"){
    simple_boot <- do.call(rbind, lapply(S, function(x) x$result$score_orig))
    nsimple_boot <- apply(simple_boot, 2, function(x) sum(!is.na(x)))
    simple_boot <- apply(simple_boot, 2, mean, na.rm=TRUE)
    stability <- do.call(cbind, lapply(S, function(x) x$result$p_orig))
    stability <- cbind(p_app = p_app, stability)
  } else{
    simple_boot <- stability <- nsimple_boot <- NULL
  }

  out <- list("apparent" = score_app,
              "optimism" = opt,
              "n_optimism" = nopt,
              "corrected" = bcorr,
              "simple" = simple_boot,
              "n_simple" = nsimple_boot,
              "stability" = stability,
              "y" = y,
              "method" = method,
              "B"=B,
              "warnings" = warns,
              "messages" = mess)

  class(out) <- "internal_boot"

  return(out)
}

#' Print a internal_boot object
#'
#' @param x an object created with \code{boot_optimism}
#' @param digits number of digits to print (default = 2)
#' @param ... additional arguments to print
#'
#' @return invisibly returns \code{x} and prints estimates to console
#' @export
print.internal_boot <- function(x, digits=2, ...){
  out <- rbind(x$apparent, x$optimism, x$n_optimism, x$corrected,
               x$simple, x$n_simple)
  rownames(out) <- c("Apparent", "Optimism", "B Optimism", "Optimism Corrected",
                     "Simple Corrected", "B Simple")
  print(out, digits = digits, ...)
  return(invisible(x))
}
