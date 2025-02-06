#' Get stability from internal_validate or internal_boot object
#'
#' @param x a internal_validate (see \code{\link{validate}}) or internal_boot
#' (see \code{\link{boot_optimism}}) object. The former should have been
#' created with method = "boot_optimism" or "boot_simple".
#'
#' @return a list containing the stability matrix (an n x B+1 matrix.
#' Each column contains predictions for p observations. First column contains
#' predictions from original model. Other columns are from bootstrap models) and
#' the original binary outcome (y). Unsuccessful resamples are omitted.
#' @export
#'
#' @keywords internal
get_stability <- function(x){
  if (is(x, "internal_validate")){
    if (!"stability" %in% names(x$scores) || is.null(x$scores$stability) ){
      stop("stability table not found in internal_validate object. ",
           "Please rerun and set method = 'boot_simple' or 'boot_optimism'.")
    }
    stabil <- x$scores$stability
    y <- x$scores$y
  } else if (is(x, "internal_boot")){
    if (!"stability" %in% names(x)){
      stop("stability table not found in internal_boot object. ",
           "Please rerun and set method = 'boot'.")
    }
    stabil <- x$stability
    y <- x$y
  } else{
    stop("x should be an object returned by pminternal::validate ",
         "(with method = 'boot_simple' or 'boot_optimism') or by ",
         "pminternal::boot_optimism (with method = 'boot')")
  }

  # omit NAs
  stabil <- stabil[, apply(stabil, 2, function(x) !any(is.na(x)))]

  return(list(stability = stabil, y = y))
}

#' Get default settings for calibration curves
#'
#' @param x ignored
#'
#' @return a list containing default arguments to supply to \code{pmcalibration::pmcalibration}
#' @export
#'
#' @keywords internal
cal_defaults <- function(x=NULL){
  return(list('smooth' = 'gam', 'k'=10,
              'ci' = 'none', 'transf' = 'logit',
              'eval' = 0))
}

