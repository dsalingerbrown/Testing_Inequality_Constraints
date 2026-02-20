# Install required packages
# install.packages(c("quadprog", "MASS"))

library(quadprog)
library(MASS)

#' Simulate Chi-Bar-Square Weights
#'
#' @param omega Matrix. The covariance matrix of the constraints.
#' @param draws Integer. Number of Monte Carlo simulations (default 500).
#' @return A vector containing the estimated weights.
simulate_weights <- function(omega, draws = 500) {
  
  P <- ncol(omega)
  omega_inv <- solve(omega) 
  Dmat <- omega_inv
  
  Amat <- diag(P)
  bvec <- rep(0, P)
  
  counts <- rep(0, P + 1)
  names(counts) <- paste0(0:P, "_positives")
  
  set.seed(123) 
  Z_draws <- mvrnorm(n = draws, mu = rep(0, P), Sigma = omega)
  
  for (i in 1:draws) {
    Z <- Z_draws[i, ]
    dvec <- omega_inv %*% Z
    
    sol <- solve.QP(Dmat = Dmat, dvec = as.vector(dvec), Amat = Amat, bvec = bvec)
    theta <- sol$solution
    
    num_positive <- sum(theta > 1e-6)
    
    counts[num_positive + 1] <- counts[num_positive + 1] + 1
  }
  
  weights <- counts / draws
  return(weights)
}

#' Calculate the Exact P-Value using Chi-Square Mixture
#'
#' @param test_stat Numeric. The observed test statistic (e.g., Wald, LR, or KT).
#' @param weights Numeric vector. The simulated weights from `simulate_weights`.
#' @return Numeric. The exact p-value.
calculate_p_value <- function(test_stat, weights) {
  
  # The number of dimensions (P) is the length of the weights vector minus 1
  # (Because the vector includes 0 positives up to P positives)
  P <- length(weights) - 1
  
  # Initialize the p-value
  p_value <- 0
  
  # Loop through every possible number of positive elements (from 0 to P)
  for (positives in 0:P) {
    
    # Calculate the Degrees of Freedom. 
    # DF = Number of Binding Constraints = (Total Constraints - Positive Elements)
    df <- P - positives
    
    # We only calculate the probability if DF > 0.
    # A Chi-square with 0 DF is just a point mass at zero, so Pr(ChiSq >= c) = 0 for any c > 0.
    if (df > 0) {
      
      # Calculate the probability from the standard Chi-Square distribution
      # lower.tail = FALSE gives us the probability of getting a value GREATER than our test stat
      chi_prob <- pchisq(test_stat, df = df, lower.tail = FALSE)
      
      # Grab the corresponding weight (Adding +1 because R indices start at 1, not 0)
      weight <- weights[positives + 1]
      
      # Multiply the probability by the weight and add it to our running total
      p_value <- p_value + (weight * chi_prob)
    }
  }
  
  return(p_value)
}


# ==========================================
# --- FULL TEST RUN ---
# ==========================================

# 1. Provide the Covariance Matrix (Omega)
dummy_omega <- matrix(c(1.0, 0.5, 0.2,
                        0.5, 1.0, 0.3,
                        0.2, 0.3, 1.0), 
                      nrow = 3)

# 2. Simulate the Weights
cat("Simulating weights...\n")
estimated_weights <- simulate_weights(dummy_omega, draws = 500)

# 3. Provide an Observed Test Statistic 
# (In a real scenario, you calculate this using your actual data. 
# Here we use 5.581, which is the Likelihood Ratio statistic from the 1989 paper example)
dummy_test_stat <- 5.581

# 4. Calculate Final P-Value
final_p_value <- calculate_p_value(test_stat = dummy_test_stat, weights = estimated_weights)

cat("\n--- Final Results ---\n")
cat("Observed Test Statistic:", dummy_test_stat, "\n")
cat("Calculated P-Value:     ", round(final_p_value, 4), "\n")
