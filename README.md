
# pSplineLocationScale

*pSplineLocationScale* is an R package which deals with the Frequentist
Estimation of Non-Linear Covariate Effects using p-Splines

The *pSplineLocationScale* package provides functions for the implementation of
penelized b-splines. It was written for the “Advanced Statistical
Programming” course at the University of Göttingen. Use the
package as an introduction to location-scale regression, or as a basis
for the implementation of additional inference algorithms and model
extensions.

## Installation

You can use the development version of *pSplineLocationScale*from Gitlab with your email/password:

``` r
install.packages("getPass")
devtools::install_git(
"[https://gitlab.gwdg.de/fabian.lukassen/psplineslocationscale.git](https://gitlab.gwdg.de/fabian.lukassen/psplineslocationscale.git)",
credentials = git2r::cred_user_pass("gitlab_mail", getPass::getPass())
)
```

## Example

The *pSplineLocationScale* package comes with the functions for generating test datasets, here is how you can use it:

``` r
library("pSplineLocationScale")
n <- 100
lower_boundary <- -5
upper_boundary <- 5
noise_sd <- 0.4
degree <- 3
noise_sd_min <- 0.1
noise_sd_max <- 0.6


################# fan data ###################################
data <- generate_fan_shaped_data(n)
test_run <- lmls_bspline(
    as.numeric(data$x), as.numeric(data$y)
)

summary_lmls_bspline(test_run)
plot_lmls_bspline(test_run, data)

################ poly data ################################
data <- generate_polynomial_data(
  n, degree, noise_sd, lower_boundary, upper_boundary
)

model <- lmls_bspline(
  data$x, data$y
)
summary_lmls_bspline(model)
plot_lmls_bspline(model, data)

################# poly sd oscilating ################################
data <- generate_polynomial_data_with_changing_sd_oscilating(
  n, degree, noise_sd_min, noise_sd_max, lower_boundary, upper_boundary
)

model <- lmls_bspline(
  data$x, data$y
)
summary_lmls_bspline(model)
plot_lmls_bspline(model, data)

################ poly changing sd linearly #############################
noise_sd_min <- 0.1
noise_sd_max <- 1

data <- generate_polynomial_data_with_changing_sd_linearly(
  n, degree, noise_sd_min, noise_sd_max, lower_boundary, upper_boundary
)

model <- lmls_bspline(
  data$x, data$y
)
summary_lmls_bspline(model)
plot_lmls_bspline(model, data)
```

![Example](./man/figures/example.png)
