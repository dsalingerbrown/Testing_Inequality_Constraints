# ============================================================
# Wolak (1989) Inequality Constraints Test
# ============================================================

library(quadprog)
library(MASS)

# ------------------------------------------------------------
# Core function
#
# Arguments:
#   y      : (n x 1) response vector
#   X      : (n x K) regressor matrix
#   R      : (P x K) constraint matrix
#   r      : (P x 1) constraint vector  (H: R*beta >= r)
#   sigma2 : known error variance (scalar); Sigma = sigma2 * I_n
#   N_sim  : number of Monte Carlo draws for weight estimation
#
# Returns a list with: b_hat, b_tilde, IU, weights, p_value
# ------------------------------------------------------------
wolak_test <- function(y, X, R, r, sigma2 = 1, N_sim = 500) {

  P <- nrow(R)

  # --- Step 1: Unconstrained GLS (= OLS when Sigma = sigma2 * I) ---
  XtX     <- t(X) %*% X
  XtX_inv <- solve(XtX)
  b_hat   <- XtX_inv %*% t(X) %*% y

  # Covariance matrix of R*b_hat: A = sigma2 * R (X'X)^{-1} R'
  A <- sigma2 * R %*% XtX_inv %*% t(R)

  # --- Step 2: Inequality-constrained estimation (ICLS) via QP ---
  # min  (y - Xb)'(y - Xb) / sigma2   s.t.  R b >= r
  # solve.QP: min  -d'b + 0.5 b'D b   s.t.  Amat' b >= bvec
  Dmat <- XtX / sigma2
  dvec <- as.vector(t(X) %*% y) / sigma2
  qp   <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(R), bvec = r)
  b_tilde <- qp$solution

  # --- Step 3: Three equivalent test statistics ---

  # Wald (W): distance in constraint space
  diff_Rb <- as.vector(R %*% b_hat - R %*% b_tilde)
  W       <- as.numeric(t(diff_Rb) %*% solve(A) %*% diff_Rb)

  # Likelihood Ratio (LR): difference in weighted residual SS
  resid_hat   <- as.vector(y - X %*% b_hat)
  resid_tilde <- as.vector(y - X %*% b_tilde)
  LR <- (sum(resid_tilde^2) - sum(resid_hat^2)) / sigma2

  # Kuhn-Tucker (KT): via Lagrange multipliers
  # From KT conditions: lambda_bar = 2 * A^{-1} * (R*b_hat - R*b_tilde)
  lambda_bar <- as.vector(2 * solve(A) %*% diff_Rb)
  KT         <- as.numeric(t(lambda_bar) %*% A %*% lambda_bar / 4)

  # When Sigma is known, LR = W = KT (up to floating point)
  IU <- W

  # --- Step 4a: Monte Carlo weights ---
  # Draw z_i ~ N(0, A), project onto {mu >= 0} under A^{-1} metric:
  #   min_mu  (z_i - mu)' A^{-1} (z_i - mu)   s.t.  mu >= 0
  # w(P, k, A) = fraction of draws whose projection has exactly k > 0
  A_inv  <- solve(A)
  Z      <- mvrnorm(N_sim, mu = rep(0, P), Sigma = A)

  n_pos <- integer(N_sim)
  for (i in seq_len(N_sim)) {
    z_i    <- Z[i, ]
    qp_mc  <- solve.QP(Dmat = A_inv, dvec = as.vector(A_inv %*% z_i),
                       Amat = diag(P), bvec = rep(0, P))
    n_pos[i] <- sum(qp_mc$solution > 1e-8)
  }

  # weights[k+1] = w(P, k, A) for k = 0, 1, ..., P
  weights <- sapply(0:P, function(k) mean(n_pos == k))

  # --- Step 4b: P-value ---
  # Pr[IU >= c] = sum_{k=0}^{P}  w(P, P-k, A) * Pr(chi2_k >= c)
  p_value <- sum(sapply(0:P, function(k) {
    weights[P - k + 1] * pchisq(IU, df = k, lower.tail = FALSE)
  }))

  # --- Report ---
  cat(sprintf("  Unconstrained (b_hat)  : %s\n",
              paste(sprintf("%.4f", b_hat), collapse = "  ")))
  cat(sprintf("  Constrained   (b_tilde): %s\n",
              paste(sprintf("%.4f", b_tilde), collapse = "  ")))
  cat(sprintf("  LR statistic           : %.6f\n", LR))
  cat(sprintf("  Wald (W) statistic     : %.6f\n", W))
  cat(sprintf("  KT statistic           : %.6f\n", KT))
  cat(sprintf("  Weights w(P,k,A)       : %s\n",
              paste(sprintf("k=%d: %.4f", 0:P, weights), collapse = "  ")))
  cat(sprintf("  P-value (using W)      : %.4f\n", p_value))
  cat(sprintf("  Decision at 5%%        : %s\n\n",
              ifelse(p_value < 0.05,
                     "Reject H — constraints violated.",
                     "Fail to reject H — constraints consistent with data.")))

  invisible(list(b_hat = b_hat, b_tilde = b_tilde,
                 LR = LR, W = W, KT = KT,
                 weights = weights, p_value = p_value))
}

# ============================================================
# Shared setup
# ============================================================
set.seed(42)

n      <- 100
sigma2 <- 1
X      <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
R      <- diag(2)       # H: beta1 >= 0 and beta2 >= 0
r      <- c(0, 0)

# ============================================================
# Case 1: Null is correct  (beta = (0.5, 0.3)')
# ============================================================
cat("========================================\n")
cat(" Case 1: True beta = (0.5, 0.3)  [H is true]\n")
cat("========================================\n")

beta1 <- c(0.5, 0.3)
y1    <- X %*% beta1 + rnorm(n, 0, sqrt(sigma2))
res1  <- wolak_test(y1, X, R, r, sigma2)

# ============================================================
# Case 2: Null is violated  (beta = (-0.5, -0.3)')
# ============================================================
cat("========================================\n")
cat(" Case 2: True beta = (-0.5, -0.3)  [H is false]\n")
cat("========================================\n")

beta2 <- c(-0.5, -0.3)
y2    <- X %*% beta2 + rnorm(n, 0, sqrt(sigma2))
res2  <- wolak_test(y2, X, R, r, sigma2)

# ============================================================
# Example 2: Monotone ordering constraints (K=3, P=2)
# H: beta1 >= beta2 and beta2 >= beta3
# R = [[1,-1,0],[0,1,-1]],  r = (0,0)'
# ============================================================
X3    <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
R_ord <- matrix(c(1, -1,  0,
                  0,  1, -1), nrow = 2, byrow = TRUE)
r_ord <- c(0, 0)

# Case 3: Null is correct — effects are decreasing (0.8 > 0.5 > 0.2)
cat("========================================\n")
cat(" Case 3: True beta = (0.8, 0.5, 0.2)  [H is true]\n")
cat("         H: beta1 >= beta2 >= beta3\n")
cat("========================================\n")

beta3 <- c(0.8, 0.5, 0.2)
y3    <- X3 %*% beta3 + rnorm(n, 0, sqrt(sigma2))
res3  <- wolak_test(y3, X3, R_ord, r_ord, sigma2)

# Case 4: Null is violated — effects are increasing (0.2 < 0.5 < 0.8)
cat("========================================\n")
cat(" Case 4: True beta = (0.2, 0.5, 0.8)  [H is false]\n")
cat("         H: beta1 >= beta2 >= beta3\n")
cat("========================================\n")

beta4 <- c(0.2, 0.5, 0.8)
y4    <- X3 %*% beta4 + rnorm(n, 0, sqrt(sigma2))
res4  <- wolak_test(y4, X3, R_ord, r_ord, sigma2)

# ============================================================
# Extension: Unknown sigma^2 case
#
# wolak_test_unknown_sigma()
#
# Same structure as wolak_test() but sigma^2 is estimated from
# the data rather than supplied by the user.
#
# Key differences from wolak_test():
#   - W  uses sigma2_hat  (unconstrained residuals) -> largest
#   - KT uses sigma2_tilde (constrained residuals)  -> smallest
#   - LR = n * log(RSS_C / RSS_U)  (log-likelihood ratio)
#   - Finite-sample ordering: W >= LR >= KT
#   - Monte Carlo weights use A_hat (consistent estimate)
#   - A separate p-value is reported for each statistic
#
# Arguments:
#   y     : (n x 1) response vector
#   X     : (n x K) regressor matrix
#   R     : (P x K) constraint matrix
#   r     : (P x 1) constraint vector  (H: R*beta >= r)
#   N_sim : number of Monte Carlo draws for weight estimation
# ============================================================
wolak_test_unknown_sigma <- function(y, X, R, r, N_sim = 500) {

  n <- nrow(X)
  K <- ncol(X)
  P <- nrow(R)

  # --- Step 1: Unconstrained OLS + sigma2_hat ---
  XtX        <- t(X) %*% X
  XtX_inv    <- solve(XtX)
  b_hat      <- XtX_inv %*% t(X) %*% y
  resid_hat  <- as.vector(y - X %*% b_hat)
  RSS_U      <- sum(resid_hat^2)
  sigma2_hat <- RSS_U / (n - K)               # unbiased; used in W and A_hat
  A_hat      <- sigma2_hat * R %*% XtX_inv %*% t(R)

  # --- Step 2: Constrained OLS + sigma2_tilde ---
  # QP objective is scale-free in sigma2, so b_tilde is identical to
  # the known-sigma2 case; only the variance estimate differs
  qp           <- solve.QP(Dmat = XtX, dvec = as.vector(t(X) %*% y),
                            Amat = t(R), bvec = r)
  b_tilde      <- qp$solution
  resid_tilde  <- as.vector(y - X %*% b_tilde)
  RSS_C        <- sum(resid_tilde^2)
  sigma2_tilde <- RSS_C / (n - K)             # unbiased; used in KT and A_tilde
  A_tilde      <- sigma2_tilde * R %*% XtX_inv %*% t(R)

  # --- Step 3: Three statistics (no longer equal in finite samples) ---
  diff_Rb <- as.vector(R %*% b_hat - R %*% b_tilde)

  # Wald: uses A_hat (unconstrained sigma2) -> largest
  W <- as.numeric(t(diff_Rb) %*% solve(A_hat) %*% diff_Rb)

  # LR: log-likelihood ratio using MLE residual variances (divide by n)
  LR <- n * log(RSS_C / RSS_U)

  # KT: uses A_tilde (constrained sigma2) -> smallest
  lambda_tilde <- as.vector(2 * solve(A_tilde) %*% diff_Rb)
  KT           <- as.numeric(t(lambda_tilde) %*% A_tilde %*% lambda_tilde / 4)

  # --- Step 4a: Monte Carlo weights using A_hat ---
  A_hat_inv <- solve(A_hat)
  Z         <- mvrnorm(N_sim, mu = rep(0, P), Sigma = A_hat)

  n_pos <- integer(N_sim)
  for (i in seq_len(N_sim)) {
    z_i      <- Z[i, ]
    qp_mc    <- solve.QP(Dmat = A_hat_inv, dvec = as.vector(A_hat_inv %*% z_i),
                         Amat = diag(P), bvec = rep(0, P))
    n_pos[i] <- sum(qp_mc$solution > 1e-8)
  }
  weights <- sapply(0:P, function(k) mean(n_pos == k))

  # --- Step 4b: Separate p-value for each statistic ---
  pval <- function(stat) {
    sum(sapply(0:P, function(k) {
      weights[P - k + 1] * pchisq(stat, df = k, lower.tail = FALSE)
    }))
  }
  pval_W  <- pval(W)
  pval_LR <- pval(LR)
  pval_KT <- pval(KT)

  # --- Report ---
  cat(sprintf("  sigma2_hat   (unconstrained): %.4f\n", sigma2_hat))
  cat(sprintf("  sigma2_tilde (constrained)  : %.4f\n", sigma2_tilde))
  cat(sprintf("  Unconstrained (b_hat)       : %s\n",
              paste(sprintf("%.4f", b_hat), collapse = "  ")))
  cat(sprintf("  Constrained   (b_tilde)     : %s\n",
              paste(sprintf("%.4f", b_tilde), collapse = "  ")))
  cat(sprintf("  LR  : %.6f  (p = %.4f)\n", LR,  pval_LR))
  cat(sprintf("  W   : %.6f  (p = %.4f)\n", W,   pval_W))
  cat(sprintf("  KT  : %.6f  (p = %.4f)\n", KT,  pval_KT))
  cat(sprintf("  Ordering W >= LR >= KT: %s\n",
              ifelse(W >= LR - 1e-10 & LR >= KT - 1e-10, "TRUE", "FALSE")))
  cat(sprintf("  Weights w(P,k,A_hat): %s\n",
              paste(sprintf("k=%d: %.4f", 0:P, weights), collapse = "  ")))
  cat(sprintf("  Decision at 5%% — W : %s\n",
              ifelse(pval_W  < 0.05, "Reject H.", "Fail to reject H.")))
  cat(sprintf("  Decision at 5%% — LR: %s\n",
              ifelse(pval_LR < 0.05, "Reject H.", "Fail to reject H.")))
  cat(sprintf("  Decision at 5%% — KT: %s\n\n",
              ifelse(pval_KT < 0.05, "Reject H.", "Fail to reject H.")))

  invisible(list(b_hat = b_hat, b_tilde = b_tilde,
                 sigma2_hat = sigma2_hat, sigma2_tilde = sigma2_tilde,
                 LR = LR, W = W, KT = KT,
                 weights = weights,
                 pval_W = pval_W, pval_LR = pval_LR, pval_KT = pval_KT))
}

# ============================================================
# Redo all four cases with unknown sigma^2
# (same y1, y2, y3, y4 generated above — no new data)
# ============================================================

cat("\n\n")
cat("############################################\n")
cat("  EXTENSION: UNKNOWN SIGMA^2\n")
cat("############################################\n\n")

# --- Example 1: Non-negativity constraints (K=2, P=2) ---

cat("========================================\n")
cat(" Case 1 (unknown sigma): beta = (0.5, 0.3)  [H is true]\n")
cat("========================================\n")
res1_u <- wolak_test_unknown_sigma(y1, X, R, r)

cat("========================================\n")
cat(" Case 2 (unknown sigma): beta = (-0.5, -0.3)  [H is false]\n")
cat("========================================\n")
res2_u <- wolak_test_unknown_sigma(y2, X, R, r)

# --- Example 2: Monotone ordering constraints (K=3, P=2) ---

cat("========================================\n")
cat(" Case 3 (unknown sigma): beta = (0.8, 0.5, 0.2)  [H is true]\n")
cat("         H: beta1 >= beta2 >= beta3\n")
cat("========================================\n")
res3_u <- wolak_test_unknown_sigma(y3, X3, R_ord, r_ord)

cat("========================================\n")
cat(" Case 4 (unknown sigma): beta = (0.2, 0.5, 0.8)  [H is false]\n")
cat("         H: beta1 >= beta2 >= beta3\n")
cat("========================================\n")
res4_u <- wolak_test_unknown_sigma(y4, X3, R_ord, r_ord)
