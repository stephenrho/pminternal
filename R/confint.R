
#' Confidence intervals for bias-corrected performance measures
#'
#' @description
#' Implements the methods discussed in Noma et al. (2021), plus some others that have not been tested.
#' Specifically, Noma et al. discuss bootstrap optimism correction ("boot_optimism" and ".632") and the percentile
#' bootstrap (\code{ci_type = "perc"}). Their paper contains some simulation results on coverage properties of these
#' CIs. If you used \code{\link{validate}} to do something other than bootstrap optimism correction or if you request
#' normal approximation CIs please note that these approaches have (to my knowledge) not been thoroughly tested.
#' \code{ci_type = "norm"} is included as it might be able to reduce the number of runs needed for "twostage" CIs.
#' See details for the difference between "shifted" and "twostage". "norm" CIs are likely to perform poorly for some
#' performance measures, such as calibration Intercept and Slope, which for regular glms are always 0 and 1, respectively,
#' on assessment of apparent performance. As "shifted" CIs are based on apparent performance they will be meaningless for these measures.
#' Use the untested methods with caution!
#'
#' @param object created by call to \code{\link{validate}}
#' @param parm a specification of which performance measures are
#' to be given confidence intervals, either a vector of numbers
#' or a vector of names. If missing, all scores are considered.
#' @param level the confidence level required
#' @param method "shifted" or "twostage" (see details)
#' @param ci_type percentile ("perc") or normal approximation ("norm") bootstrap CIs
#' @param R number of replicates
#' @param add return the object with an additional slot containing CIs (default) or
#' just return the CIs
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' The two methods are as follows (see Noma et al. (2021) for more details):
#' \describe{
#' \item{shifted}{ (default) This approach is based on shifting bootstrap CIs for apparent performance
#' by optimism. This makes it the faster option as only the calculation of apparent performance is needed for
#' each replicate. If the CI for apparent performance is [lower, upper], the resulting CI for bias-corrected performance
#' is [lower - optimism, upper - optimism]. Note this method is only available when using an optimism based approach
#' (and "cv_optimism" was untested in Noma et al).}
#' \item{twostage}{This approach creates a bootstrap resample of the data and runs the entire
#' validation procedure on the resample (with the same number of 'inner' replicates, determined by B in the original
#' validate call). The CI is then constructed using the corrected estimates from the R 'outer' replicates.
#' As this involves R*B replicates, this could take a long time. Note \code{\link{validate}}
#' takes a \code{cores} argument that can allow the inner samples to run in parallel. }
#' }
#'
#' @returns A list with two elements, each a matrix with columns giving lower and upper confidence limits for each measure. One for apparent and one for bias-corrected measures.
#' Columns will be labelled as (1-level)/2 and 1 - (1-level)/2 in \% (by default 2.5\% and 97.5\%).
#'
#' @references Noma, H., Shinozaki, T., Iba, K., Teramukai, S., & Furukawa, T. A. (2021). Confidence intervals of prediction accuracy measures for multivariable prediction models based on the bootstrap‐based optimism correction methods. Statistics in Medicine, 40(26), 5691-5701.
#'
#' @rdname confint.internal_validatesummary
#' @export
confint.internal_validate <- function(object, parm, level = 0.95,
                                      method = c("shifted", "twostage"),
                                      ci_type = c("perc", "norm"),
                                      R = 1000, add = TRUE, ...){

  method <- match.arg(method)
  ci_type <- match.arg(ci_type)

  a <- (1 - level)/2
  a <- c(a, 1 - a)

  if (object$method == "none") stop("No bias-corrected scores found. Re-run validate with method != 'none'")
  if (method == "shifted" & !object$method %in% c("boot_optimism", ".632", "cv_optimism")) stop("'shifted' CIs are only for optimism based validation methods")

  app <- object$apparent
  pnames <- names(app)
  if (missing(parm))
    parm <- pnames
  else {
    if (length(parm) == 0){
      parm <- pnames
    }
    if (is.numeric(parm)) {
      parm <- pnames[parm]
    }
  }

  # make function to bootstrap
  vcall <- object$call
  if (method == "shifted") vcall[["method"]] <- "none"

  vfun <- function(i, data){
    vcall[["data"]] <- data[i, ]
    out <- eval(vcall)
    out <- eval(out)
    return(list(app = out$apparent[parm],
                corr = out$corrected[parm],
                b = out$B_actual[parm]))
  }

  # resample
  data <- object$data
  n <- nrow(data)
  indices <- lapply(seq(R), function(x) sample(x = n, size = n, replace = TRUE))

  S <- lapply(indices, FUN = vfun, data = data)

  # extract things
  app <- do.call(rbind, lapply(S, function(x) x$app))
  corr <- do.call(rbind, lapply(S, function(x) x$corr))

  if (ci_type == "norm"){
    q <- stats::qnorm(a)
    if (method == "shifted"){
      ses <- apply(app, MARGIN = 2, stats::sd)
      app_ci <- object$apparent[parm] + ses %o% q
      corr_ci <- object$corrected[parm] + ses %o% q
      # for method = "shifted" this is the same as
      # finding ci for apparent then subtracting optimism from bounds
    } else{
      app_ses <- apply(app, MARGIN = 2, stats::sd)
      corr_ses <- apply(corr, MARGIN = 2, stats::sd)
      app_ci <- object$apparent[parm] + app_ses %o% q
      corr_ci <- object$corrected[parm] + corr_ses %o% q
    }
  } else if (ci_type == "perc"){
    app_ci <- t(apply(app, MARGIN = 2, quantile, probs = a))
    if (method == "shifted"){
      corr_ci <- app_ci - object$optimism[parm]
    } else{
      corr_ci <- t(apply(corr, MARGIN = 2, quantile, probs = a))
    }
  }

  colnames(app_ci) <- sprintf("%0.1f%%", a * 100)
  colnames(corr_ci) <- sprintf("%0.1f%%", a * 100)

  ci <- list(Apparent = app_ci, Corrected = corr_ci)

  if (add){
    object$confint <- list(ci = ci, level = level, method = method, ci_type = ci_type)
    return(object)
  } else{
    # add average b?
    return(ci)
  }
}
