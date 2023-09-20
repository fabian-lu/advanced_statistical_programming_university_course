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
#' @param knot_count The number of knots to use in the B-spline basis.
#' @param poly_degree The degree of the B-spline polynomials.
#' @param hyperparameter_lambda_mu The mean of the hyperparameter lamda.
#' @param hyperparameter_lambda_sigma The variance of the hyperparameter lambda.
#' @param order_difference_matrix_r The order of the difference matrix used for penalty.
#'
#' @return A list containing the initial values of mu_hat, sigma_hat, beta_hat, gamma_hat, and the design matrix Z.
#' @importFrom stats fitted.values lm sigma
#'
init <- function(vector_x, response_vector_y, knot_count, poly_degree,
                 hyperparameter_lambda_mu, hyperparameter_lambda_sigma,
                 order_difference_matrix_r) {
  initial_model <- lm(response_vector_y ~ vector_x)
  y_hat <- fitted.values(initial_model)
  mu_hat_vector <- c(rep(mean(y_hat), times = length(response_vector_y)))
  sigma_hat_vector <- c(rep(sigma(initial_model), times = length(response_vector_y)))

  p_splines <- PLS_b_Splines(
    knot_count, lower_bound = min(vector_x),
    upper_bound = max(vector_x), poly_degree, vector_x,
    response_vector_y, hyperparameter_lambda_mu
  )

  penalty_matrix_K <- create_penalty_matrix_K(
    order_difference_matrix_r, poly_degree, knot_count
  )

  return(
    list(
      mu_hat = mu_hat_vector,
      sigma_hat = sigma_hat_vector,
      beta_hat = p_splines$gamma_hat,
      gamma_hat = p_splines$gamma_hat,
      design_matrix_Z = p_splines$Z,
      penalty_matrix_K = penalty_matrix_K
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
  log_likelihood <- - length(mu_hat)/2*log(2*pi) - sum(log(sigma_hat)) -
    1/2*sum(((response_vector_y-mu_hat)/sigma_hat)^2)

  return(-2 * log_likelihood)
}

############################# Finding Lambda #################################

#' Optimise GCV
#'
#' This function optimizes the GCV criterion to find the optimal value of lambda.
#' It uses the Brent method for optimization.
#'
#' @param lambda_mu_flag A flag indicating whether lambda is for mu or sigma estimation.
#' @param design_matrix_Z The B-spline design matrix.
#' @param W The weight matrix.
#' @param w A vector containing weights.
#' @param z A vector of values.
#' @param K The penalty matrix.
#' @param lambda_acc The tuning parameter for the optimization method.
#' @param n The number of observations.
#'
#' @return The optimal value of lambda for GCV.
compute_lambda_gcv <- function(lambda_mu_flag, design_matrix_Z, W, w, z, K, lambda_acc, n) {

  if(lambda_mu_flag) {
    lambda_gcv <- seq(0, 50, length.out = lambda_acc)
  } else {
    lambda_gcv <- seq(1, 5000, length.out = lambda_acc)
  }

  score_gcv <- rep(0, lambda_acc)

  for(j in 1:lambda_acc) {

    if(lambda_mu_flag) {
      parameter_hat <- solve(t(design_matrix_Z) %*% W %*% design_matrix_Z + lambda_gcv[j] * K, tol = 1e-30) %*% t(design_matrix_Z) %*% W %*% z
      first_element <- sqrt(w) * (z - design_matrix_Z %*% parameter_hat)
    } else {
      parameter_hat <- solve(t(design_matrix_Z) %*% W %*% design_matrix_Z + lambda_gcv[j] * K, tol = 1e-30) %*% t(design_matrix_Z) %*% W %*% z
      first_element <- sqrt(w) * (z - exp(design_matrix_Z %*% parameter_hat))
    }

    S <- design_matrix_Z %*% solve(t(design_matrix_Z) %*% W %*% design_matrix_Z + lambda_gcv[j] * K, tol = 1e-30) %*% t(design_matrix_Z) %*% W
    score_gcv[j] <- (n * (t(first_element) %*% first_element)) / (n - sum(diag(S)))^2

  }

  lambda_optimal <- lambda_gcv[which.min(score_gcv)]

  return(lambda_optimal)
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
#' @param lambda_acc The tuning parameter for the optimization method.
#'
#' @return A list containing the updated mean parameter (mu_hat), updated B-spline coefficients (beta_hat),
#' and the optimal lambda_mu value.
#'
inner_mu <- function(
    response_vector_y, mu_hat, sigma_hat, design_matrix_Z, beta_hat, lambda_mu,
    K_mu, lambda_acc
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
    , tol = 1e-30) %*%

      t(design_matrix_Z) %*% W_mu %*% z_mu

    mu_hat <- design_matrix_Z %*% beta_hat

    current_deviance <- global_deviance(response_vector_y, mu_hat, sigma_hat)

    if (round(current_deviance, 2) == round(prev_deviance, 2)) {
      print("Convergence reached in inner mu.")
      break
    }

    prev_deviance <- current_deviance

  }
  opt_lambda_mu <- compute_lambda_gcv(
    TRUE, design_matrix_Z, W_mu, w_mu, z_mu, K_mu, lambda_acc, length(response_vector_y)
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
#' @param lambda_acc The tuning parameter for the optimization method.
#'
#' @return A list containing the updated variance parameter (sigma_hat), updated gamma coefficients (gamma_hat),
#' and the optimal lambda_sigma value.
#'
inner_sigma <- function(
    mu_hat, sigma_hat, design_matrix_Z, response_vector_y, gamma_hat,
    lambda_sigma, K_sigma, lambda_acc
) {

  prev_deviance <- global_deviance(response_vector_y, mu_hat, sigma_hat)

  for(i in 1:10) {

    u_sigma <- ((response_vector_y-mu_hat)^2-sigma_hat^2)/sigma_hat^2
    w_sigma <- u_sigma^2
    z_sigma <- log(sigma_hat) + 1/u_sigma

    # Update gamma_hat using PWLS
    W_sigma <- diag(c(w_sigma))
    gamma_hat <- solve(t(design_matrix_Z) %*% W_sigma %*%
                         design_matrix_Z + lambda_sigma * K_sigma, tol = 1e-30) %*%
      t(design_matrix_Z) %*% W_sigma %*% z_sigma

    sigma_hat <- exp(design_matrix_Z %*% gamma_hat)
    current_deviance <- global_deviance(response_vector_y, mu_hat, sigma_hat)

    if (round(current_deviance, 2) == round(prev_deviance, 2)) {
      print("Convergence reached in inner sigma.")
      break
    }

    prev_deviance <- current_deviance
  }

  opt_lambda_sigma <- compute_lambda_gcv(
    FALSE, design_matrix_Z, W_sigma, w_sigma, z_sigma, K_sigma, lambda_acc, length(response_vector_y)
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
  iterations <- 0
  for(i in 1:max_iterations) {
    iterations <- iterations + 1
    mu_result <- inner_mu(
      response_vector_y, mu_hat, sigma_hat, design_matrix_Z, beta_hat, lambda_mu, K_mu,
      lambda_acc
    )
    mu_hat <- mu_result$mu_hat
    beta_hat <- mu_result$beta_hat
    lambda_mu <- mu_result$lambda_mu

    sigma_result <- inner_sigma(
      mu_hat, sigma_hat, design_matrix_Z, response_vector_y, gamma_hat, lambda_sigma,
      K_sigma, lambda_acc
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
      gamma_hat = gamma_hat,
      GD_conv = current_deviance,
      iterations = iterations

    )
  )
}
