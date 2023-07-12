############################ RS-Algorithm #####################################

######################## Auxilary functions ###################################

#' Initialize Starting Values
#'
#' This function initializes the starting values for the RS algorithm.
#' It uses a simple linear model to estimate the initial values of mu_hat and sigma_hat.
#' The function returns a list containing the initial values.
#'
#' @param vector_x The input vector of x values.
#' @param response_vector_y The response vector of y values.
#'
#' @return A list containing the initial values of mu_hat, sigma_hat, beta_hat, gamma_hat, and the design matrix Z.
#'
init <- function(vector_x, response_vector_y) {
  mu_hat = mean(response_vector_y)
  sigma_hat = sd(response_vector_y)

  mu_hat_vector = rep(mu_hat, length(response_vector_y))
  sigma_hat_vector = rep(sigma_hat, length(response_vector_y))

  p_splines <- PLS_b_Splines(
    knot_count = 40, lower_bound = min(vector_x),
    upper_bound = max(vector_x), poly_degree = 3, vector_x,
    response_vector_y, hyperparameter_lambda = 3
  )

  y_hat <- p_splines$Z %*% p_splines$gamma_hat
  residuals <- response_vector_y - y_hat

  p_splines_for_sigma <- PLS_b_Splines(
    knot_count = 40, lower_bound = min(vector_x),
    upper_bound = max(vector_x), poly_degree = 3, vector_x,
    residuals, hyperparameter_lambda = 500
  )

  return(
    list(
      mu_hat = mu_hat_vector,
      sigma_hat = sigma_hat_vector,
      beta_hat = p_splines$gamma_hat,
      gamma_hat = p_splines_for_sigma$gamma_hat,
      design_matrix_Z = p_splines$Z
    )
  )
}


#' Global Deviance
#'
#' This function calculates the twice negative log-likelihood of the model.
#' It is used to check for convergence within the RS algorithm.
#'
#' @param response_vector_y The response vector of y values.
#' @param mu_hat The current estimate of mu_hat.
#' @param sigma_hat The current estimate of sigma_hat.
#'
#' @return The twice negative log-likelihood.
#' @importFrom stats dnorm optim
global_deviance <- function(response_vector_y, mu_hat, sigma_hat) {
  log_likelihood <- sum(
    dnorm(response_vector_y, mean = mu_hat, sd = sigma_hat, log = TRUE)
  )

  return(-2 * log_likelihood)
}

############################# Finding Lambda #################################

#' Compute S Matrix
#'
#' This function computes the S matrix used in the RS algorithm.
#' It calculates S as Z * inv(t(Z) * W * Z + lambda * K) * t(Z) * W.
#'
#' @param Z The B-spline design matrix.
#' @param W The weight matrix.
#' @param lambda The current value of lambda.
#' @param K The penalty matrix.
#'
#' @return The S matrix.
compute_S <- function(Z, W, lambda, K) {
  Z %*% solve(t(Z) %*% W %*% Z + lambda * K) %*% t(Z) %*% W
}

#' Compute GCV
#'
#' This function computes the Generalized Cross Validation (GCV) criterion used in the RS algorithm.
#' It calculates GCV as n * ||epsilon||^2 / (n - tr(S))^2, where epsilon is the residuals and S is the S matrix.
#'
#' @param n The number of observations.
#' @param inner_euclidean The inner Euclidean vector.
#' @param trace_S The trace of the S matrix.
#'
#' @return The GCV value.
compute_GCV <- function(n, inner_euclidean, trace_S) {
  n * norm(inner_euclidean, type = "F")^2 / (n - trace_S)^2
}

#' Compute GAIC
#'
#' This function computes the Generalized Akaike Information Criterion (GAIC) used in the RS algorithm.
#' It calculates GAIC as ||epsilon||^2 + nu * tr(S), where epsilon is the residuals, S is the S matrix, and nu is a parameter.
#'
#' @param inner_euclidean The inner Euclidean vector.
#' @param trace_S The trace of the S matrix.
#' @param nu The parameter for penalization.
#'
#' @return The GAIC value.
compute_GAIC <- function(inner_euclidean, trace_S, nu = 2) {
  norm(inner_euclidean, type = "F")^2 + nu * trace_S
}

#' Optimise GCV
#'
#' This function optimizes the GCV criterion to find the optimal value of lambda.
#' It uses the Brent method for optimization.
#'
#' @param lambda_mu_flag A flag indicating whether lambda is for mu or sigma estimation.
#' @param Z The B-spline design matrix.
#' @param W The weight matrix.
#' @param w The weights.
#' @param lambda The current value of lambda.
#' @param K The penalty matrix.
#' @param epsilon The epsilon vector.
#' @param beta_hat The current estimate of beta_hat or gamma_hat.
#' @param n The number of observations.
#'
#' @return The optimal value of lambda for GCV.
#' @importFrom stats dnorm optim
optimise_GCV <- function(lambda_mu_flag, Z, W, w, lambda, K, epsilon, beta_hat, n) {
  lower_bound <- 0
  upper_bound <- if (lambda_mu_flag) 50 else 5000
  optimise_func <- function(lambda) {
    S = compute_S(Z, W, lambda, K)
    trace_S = sum(diag(S))
    inner_euclidean = sqrt(w) * (epsilon - Z %*% beta_hat)
    compute_GCV(n, inner_euclidean, trace_S)
  }
  optim_result <- optim(lambda, optimise_func, method = "Brent", lower = lower_bound, upper = upper_bound)
  return(optim_result$par)
}

#' Optimise GAIC
#'
#' This function optimizes the GAIC criterion to find the optimal value of lambda.
#' It uses the Brent method for optimization.
#'
#' @param lambda_mu_flag A flag indicating whether lambda is for mu or sigma estimation.
#' @param Z The B-spline design matrix.
#' @param W The weight matrix.
#' @param lambda The current value of lambda.
#' @param K The penalty matrix.
#' @param epsilon The epsilon vector.
#' @param beta_hat The current estimate of beta_hat or gamma_hat.
#' @param nu The parameter for penalization.
#'
#' @return The optimal value of lambda for GAIC.
#' @importFrom stats dnorm optim
optimise_GAIC <- function(lambda_mu_flag, Z, W, lambda, K, epsilon, beta_hat, nu = 2) {
  lower_bound <- 0
  upper_bound <- if (lambda_mu_flag) 50 else 5000
  optimise_func <- function(lambda) {
    S = compute_S(Z, W, lambda, K)
    trace_S = sum(diag(S))
    inner_euclidean = sqrt(W) * (epsilon - Z %*% beta_hat)
    compute_GAIC(inner_euclidean, trace_S, nu)
  }
  optim_result <- optim(lambda, optimise_func, method = "Brent", lower = lower_bound, upper = upper_bound)
  return(optim_result$par)
}

########################  # End finding  lambda ################################

################## Main Part of the RS Algo ##################################
#' Inner function for mu estimation
#'
#' This function performs the inner loop for estimating the mean parameter (mu_hat) in the Local Maximum Likelihood (LML) regression model.
#' It updates the beta_hat coefficients using Penalized Weighted Least Squares (PWLS) and iteratively updates mu_hat until convergence is reached.
#'
#' @param response_vector_y A vector of response values.
#' @param mu_hat The initial estimate of the mean parameter.
#' @param sigma_hat The estimate of the variance parameter.
#' @param design_matrix_Z The design matrix of B-spline basis functions.
#' @param beta_hat The current estimate of the B-spline coefficients.
#' @param lambda_mu The tuning parameter (lambda) for updating the mean parameter.
#' @param K_mu The penalty matrix for the mean parameter.
#'
#' @return A list containing the updated mean parameter (mu_hat), updated B-spline coefficients (beta_hat),
#' and the optimal lambda_mu value.
#'
inner_mu <- function(
    response_vector_y, mu_hat, sigma_hat, design_matrix_Z, beta_hat, lambda_mu, K_mu
) {
  prev_deviance <- global_deviance(response_vector_y, mu_hat, sigma_hat)

  for(i in 1:10) {
    u_mu <- (response_vector_y - mu_hat) / sigma_hat^2
    w_mu <- u_mu ^ 2
    z_mu <- mu_hat + (1 / u_mu)

    # Update beta_hat using PWLS
    W_mu <- diag(c(w_mu))
    beta_hat <- solve(
      t(design_matrix_Z) %*% W_mu %*% design_matrix_Z + lambda_mu * K_mu
    ) %*%

      t(design_matrix_Z) %*% W_mu %*% z_mu

    mu_hat <- design_matrix_Z %*% beta_hat

    current_deviance <- global_deviance(response_vector_y, mu_hat, sigma_hat)

    if (round(current_deviance, 2) == round(prev_deviance, 2)) {
      print("Convergence reached in inner mu.")
      break
    }

    prev_deviance <- current_deviance

  }

  opt_lambda_mu <- optimise_GCV(
    TRUE, design_matrix_Z, W_mu, w_mu, lambda_mu, K_mu, z_mu, beta_hat,
    length(response_vector_y)
  )

  return(list(mu_hat = mu_hat, beta_hat = beta_hat, lambda_mu = opt_lambda_mu))
}


#' Inner function for sigma estimation
#'
#' This function performs the inner loop for estimating the variance parameter (sigma_hat) in the Local Maximum Likelihood (LML) regression model.
#' It updates the gamma_hat coefficients using Penalized Weighted Least Squares (PWLS) and iteratively updates sigma_hat until convergence is reached.
#'
#' @param mu_hat The current estimate of the mean parameter.
#' @param sigma_hat The initial estimate of the variance parameter.
#' @param design_matrix_Z The design matrix of B-spline basis functions.
#' @param response_vector_y A vector of response values.
#' @param gamma_hat The current estimate of the gamma coefficients for variance.
#' @param lambda_sigma The tuning parameter (lambda) for updating the variance parameter.
#' @param K_sigma The penalty matrix for the variance parameter.
#'
#' @return A list containing the updated variance parameter (sigma_hat), updated gamma coefficients (gamma_hat),
#' and the optimal lambda_sigma value.
#'
inner_sigma <- function(
    mu_hat, sigma_hat, design_matrix_Z, response_vector_y, gamma_hat,
    lambda_sigma, K_sigma
) {

  prev_deviance <- global_deviance(response_vector_y, mu_hat, sigma_hat)

  for(i in 1:10) {

    u_sigma <- (-1/sigma_hat) + (response_vector_y - mu_hat)/sigma_hat^3
    w_sigma <- u_sigma ^ 2
    z_sigma <- log(sigma_hat) + 1 / u_sigma

    # Update gamma_hat using PWLS
    W_sigma <- diag(c(w_sigma))
    gamma_hat <- solve(t(design_matrix_Z) %*% W_sigma %*%
                         design_matrix_Z + lambda_sigma * K_sigma) %*%
      t(design_matrix_Z) %*% W_sigma %*% z_sigma

    sigma_hat <- exp(design_matrix_Z %*% gamma_hat)
    current_deviance <- global_deviance(response_vector_y, mu_hat, sigma_hat)

    if (round(current_deviance, 2) == round(prev_deviance, 2)) {
      print("Convergence reached in inner sigma.")
      break
    }

    prev_deviance <- current_deviance
  }

  opt_lambda_sigma <- optimise_GCV(
    FALSE, design_matrix_Z, W_sigma, w_sigma, lambda_sigma, K_sigma, z_sigma,
    gamma_hat, length(response_vector_y)
  )

  return(
    list(
      sigma_hat = sigma_hat,
      gamma_hat = gamma_hat,
      lambda_sigma = opt_lambda_sigma
    )
  )
}


#' Successively updates the mean and variance parameters
#'
#' This function performs the outer loop of the Local Maximum Likelihood (LML) regression model,
#' successively updating the mean (mu_hat) and variance (sigma_hat) parameters until convergence is reached.
#'
#' @param mu_hat The initial estimate of the mean parameter.
#' @param sigma_hat The initial estimate of the variance parameter.
#' @param design_matrix_Z The design matrix of B-spline basis functions.
#' @param vector_x A vector of predictor values.
#' @param response_vector_y A vector of response values.
#' @param beta_hat The initial estimate of the B-spline coefficients.
#' @param gamma_hat The initial estimate of the gamma coefficients for variance.
#' @param lambda_mu The tuning parameter (lambda) for updating the mean parameter.
#' @param lambda_sigma The tuning parameter (lambda) for updating the variance parameter.
#' @param max_iterations The maximum number of iterations for the algorithm.
#' @param K_mu The penalty matrix for the mean parameter.
#' @param K_sigma The penalty matrix for the variance parameter.
#' @param knot_count The number of knots in the B-spline basis.
#' @param poly_degree The degree of the B-spline polynomials.
#' @param opt_method The optimization method to use for tuning parameter selection ('GCV' or 'GAIC').
#' @param lambda_acc The tuning parameter for the optimization method.
#'
#' @return A list containing the final estimates of the mean (mu_hat) and variance (sigma_hat),
#' along with the estimates of the B-spline coefficients (beta_hat) and gamma coefficients for variance (gamma_hat).
#'
outer <- function(mu_hat, sigma_hat, design_matrix_Z, vector_x, response_vector_y,
                  beta_hat, gamma_hat, lambda_mu, lambda_sigma, max_iterations,
                  K_mu, K_sigma, knot_count, poly_degree, opt_method, lambda_acc) {

  prev_deviance <- global_deviance(response_vector_y, mu_hat, sigma_hat)

  for(i in 1:max_iterations) {

    mu_result <- inner_mu(
      response_vector_y, mu_hat, sigma_hat, design_matrix_Z, beta_hat, lambda_mu, K_mu
    )
    mu_hat <- mu_result$mu_hat
    beta_hat <- mu_result$beta_hat
    lambda_mu <- mu_result$lambda_mu

    sigma_result <- inner_sigma(
      mu_hat, sigma_hat, design_matrix_Z, response_vector_y, gamma_hat, lambda_sigma, K_sigma
    )
    sigma_hat <- sigma_result$sigma_hat
    gamma_hat <- sigma_result$gamma_hat
    lambda_sigma <- sigma_result$lambda_sigma

    current_deviance <- global_deviance(response_vector_y, mu_hat, sigma_hat)

    if (round(current_deviance, 2) == round(prev_deviance, 2)) {
      print("Convergence reached in outer.")
      break
    }

    prev_deviance <- current_deviance
  }

  return(
    list(
      mu_hat = mu_hat,
      sigma_hat = sigma_hat,
      beta_hat = beta_hat,
      gamma_hat = gamma_hat
    )
  )
}
