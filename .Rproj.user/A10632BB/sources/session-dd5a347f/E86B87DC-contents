#' Generate polynomial data
#'
#' Generates a dataset with polynomial relationship between x and y, with additive Gaussian noise.
#'
#' @param n The number of data points to generate.
#' @param degree The degree of the polynomial.
#' @param noise_sd The standard deviation of the Gaussian noise.
#' @param lower_boundary The lower boundary of the x values.
#' @param upper_boundary The upper boundary of the x values.
#'
#' @return A data frame with columns 'x' and 'y'.
#'
#' @examples
#' generate_polynomial_data(100, 3, 0.4, -5, 5)
#' @export
#' @importFrom stats poly rnorm runif
generate_polynomial_data <- function(n, degree, noise_sd, lower_boundary, upper_boundary) {
  x <- seq(lower_boundary, upper_boundary, length.out = n)
  predictors_vector <- runif(degree, min = 1, max = 6)
  y <- poly(x, degree) %*% predictors_vector + rnorm(n, mean = 0, sd = noise_sd)

  data <- data.frame(x, y)
  return(data)
}

#' Generate polynomial data with oscillating noise standard deviation
#'
#' Generates a dataset with polynomial relationship between x and y, where the noise standard deviation
#' oscillates based on the sine function.
#'
#' @param n The number of data points to generate.
#' @param degree The degree of the polynomial.
#' @param noise_sd_min The minimum value of the noise standard deviation.
#' @param noise_sd_max The maximum value of the noise standard deviation.
#' @param lower_boundary The lower boundary of the x values.
#' @param upper_boundary The upper boundary of the x values.
#'
#' @return A data frame with columns 'x' and 'y'.
#'
#' @examples
#' generate_polynomial_data_with_changing_sd_oscilating(100, 3, 0.1, 0.6, -5, 5)
#' @export
#' @importFrom stats poly rnorm runif sd
generate_polynomial_data_with_changing_sd_oscilating <- function(n, degree, noise_sd_min, noise_sd_max, lower_boundary, upper_boundary) {
  x <- seq(lower_boundary, upper_boundary, length.out = n)
  noise_sd <- (sin(x) + 1) / 2 * (noise_sd_max - noise_sd_min) + noise_sd_min
  predictors_vector <- runif(degree, min = 1, max = 6)
  y <- poly(x, degree) %*% predictors_vector + rnorm(n, mean = 0, sd = noise_sd)

  data <- data.frame(x, y)
  return(data)
}

#' Generate polynomial data with linearly increasing noise standard deviation
#'
#' Generates a dataset with polynomial relationship between x and y, where the noise standard deviation
#' linearly increases from a start value to an end value.
#'
#' @param n The number of data points to generate.
#' @param degree The degree of the polynomial.
#' @param noise_sd_start The starting value of the noise standard deviation.
#' @param noise_sd_end The ending value of the noise standard deviation.
#' @param lower_boundary The lower boundary of the x values.
#' @param upper_boundary The upper boundary of the x values.
#'
#' @return A data frame with columns 'x' and 'y'.
#'
#' @examples
#' generate_polynomial_data_with_changing_sd_linearly(100, 3, 0.1, 1, -5, 5)
#' @export
#' @importFrom stats poly rnorm runif sd
generate_polynomial_data_with_changing_sd_linearly <- function(n, degree, noise_sd_start, noise_sd_end, lower_boundary, upper_boundary) {
  x <- seq(lower_boundary, upper_boundary, length.out = n)
  noise_sd <- seq(noise_sd_start, noise_sd_end, length.out = n)
  predictors_vector <- runif(degree, min = 1, max = 6)
  y <- poly(x, degree) %*% predictors_vector + rnorm(n, mean = 0, sd = noise_sd)

  data <- data.frame(x, y)
  return(data)
}

#' Generate fan-shaped data
#'
#' Generates a dataset with a fan-shaped relationship between x and y, with increasing variability.
#'
#' @param n The number of data points to generate.
#'
#' @return A data frame with columns 'x' and 'y'.
#'
#' @examples
#' generate_fan_shaped_data(100)
#' @export
#' @importFrom stats rnorm sd
generate_fan_shaped_data <- function(n) {
  x <- 1:n
  y <- x + rnorm(n, mean = 0, sd = sqrt(x)*2)

  data <- data.frame(x, y)
  return(data)
}
