
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pminternal: Internal Validation of Clinical Prediction Models

<!-- badges: start -->
<!-- badges: end -->

This is a work in progress.

The goal is to offer a package that can produce bias-corrected
performance measures for binary outcomes for a range of model
development approaches available in R (similar to `rms::validate`). Also
contains functions for assessing prediction stability as described here
<https://doi.org/10.1002/bimj.202200302>.

To install development version:

``` r
# install.packages("devtools")
devtools::install_github("stephenrho/pminternal", build_vignettes = TRUE)
```

Please send feedback to <steverho89@gmail.com> or open an issue.

## Example

In the example below we use bootstrapping to correct performance
measures for a `glm` via calculation of ‘optimism’ (see
`vignette("pminternal")` and `vignette("validate-examples")` for more
examples):

``` r
library(pminternal)

# make some data
set.seed(2345)
n <- 800
p <- 10

X <- matrix(rnorm(n*p), nrow = n, ncol = p)
LP <- -1 + apply(X[, 1:5], 1, sum) # first 5 variables predict outcome
y <- rbinom(n, 1, plogis(LP))

dat <- data.frame(y, X)

# fit a model
mod <- glm(y ~ ., data = dat, family = "binomial")

# calculate bootstrap optimism corrected performance measures
(val <- validate(fit = mod, method = "boot_optimism", B = 100))
#> It is recommended that B > 200 for bootstrap validation
#>                C   Brier Intercept Slope    Eavg     E50     E90    Emax
#> Apparent  0.8567  0.1423   4.4e-12 1.000  0.0045  0.0039  0.0081  0.0109
#> Optimism  0.0093 -0.0054   1.7e-02 0.053 -0.0048 -0.0050 -0.0107 -0.0057
#> Corrected 0.8474  0.1477  -1.7e-02 0.947  0.0093  0.0089  0.0187  0.0165
#>               ECI
#> Apparent   0.0027
#> Optimism  -0.0038
#> Corrected  0.0065
```

The other available methods for calculating bias corrected performance
are the simple bootstrap (`boot_simple`), 0.632 bootstrap optimism
(`.632`), optimism via cross-validation (`cv_optimism`), and regular
cross-validation (`cv_average`). Please see `?pminternal::validate` and
the references therein. Bias corrected calibration curves can also be
produced (see `cal_plot`).

For models that cannot be supported via `fit`, users are able to specify
their own model (`model_fun`) and prediction (`pred_fun`) functions as
shown below. Note that when specifying user-defined model and prediction
functions the data and outcome must also be provided. It is crucial that
`model_fun` implements the entire model development procedure (variable
selection, hyperparameter tuning, etc). For more examples, see
`vignette("pminternal")` and `vignette("validate-examples")`.

``` r
# fit a glm with lasso penalty
library(glmnet)
#> Loading required package: Matrix
#> Loaded glmnet 4.1-7

lasso_fun <- function(data, ...){
  y <- data$y
  x <- as.matrix(data[, which(colnames(data) != "y")])
  
  cv <- cv.glmnet(x=x, y=y, alpha=1, nfolds = 10, family="binomial")
  lambda <- cv$lambda.min
  
  glmnet(x=x, y=y, alpha = 1, lambda = lambda, family="binomial")
}

lasso_predict <- function(model, data, ...){
  y <- data$y
  x <- as.matrix(data[, which(colnames(data) != "y")])
  
  predict(model, newx = x, type = "response")[,1]
}

(val <- validate(data = dat, outcome = "y", 
                 model_fun = lasso_fun, pred_fun = lasso_predict, 
                 method = "boot_optimism", B = 100))
#> It is recommended that B > 200 for bootstrap validation
#>                C   Brier Intercept Slope   Eavg    E50    E90  Emax   ECI
#> Apparent  0.8558  0.1427     0.073  1.14 0.0184 0.0178 0.0366 0.040 0.044
#> Optimism  0.0062 -0.0037     0.015  0.04 0.0025 0.0033 0.0039 0.014 0.016
#> Corrected 0.8496  0.1463     0.057  1.10 0.0159 0.0144 0.0326 0.026 0.028
```

The output of `validate` (with `method = "boot_*"`) can be used to
produce plots for assessing the stability of model predictions (across
models developed on bootstrap resamples).

A prediction (in)stability plot shows predictions from the `B` (in this
case 100) bootstrap models applied to the development data.

``` r
prediction_stability(val)
```

<img src="man/figures/README-stability-1.png" width="90%" />

A MAPE plot shows the mean absolute prediction error, which is the
difference between the predicted risk from the development model and
each of the `B` bootstrap models.

``` r
mape_stability(val)
```

<img src="man/figures/README-mape-1.png" width="90%" />

A calibration (in)stability plot depict the original calibration curve
along with `B` calibration curves from the bootstrap models applied to
the original data (`y`).

``` r
calibration_stability(val)
```

<img src="man/figures/README-cal-1.png" width="90%" />

The classification instability index (CII) is the proportion of
individuals that change predicted class (present/absent, 1/0) when
predicted risk is compared to some threshold. For example, a patient
predicted to be in class 1 would receive a CII of 0.3 if 30% of the
bootstrap models led to a predicted class of 0.

``` r
classification_stability(val, threshold = .4)
```

<img src="man/figures/README-CII-1.png" width="90%" />

Decision curves implied by the original and bootstrap models can also be
plotted.

``` r
dcurve_stability(val)
```

<img src="man/figures/README-dcurve-1.png" width="90%" />
