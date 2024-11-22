# STEP A: Data simulation by using Euler-Maruyama method
# ================================================
# We simulate the trajectory of the Jacobi process:
# dX_t = -a(X_t - b) dt + sqrt(2 * X_t * (1 - X_t)) dB_t,
# where a = alpha + beta + 2 and b = (beta + 1) / (alpha + beta + 2).


# 1. Parameter values

alpha <-   # Enter a value
beta  <-  # Enter a value
a <- alpha + beta + 2
b <- (beta + 1) / (alpha + beta + 2)

# 2. Setup

set.seed(123)  # For reproducibility
T <-        # Enter time horizon
dt <- 0.01     # Time step of the process discretization
M <- T/dt      # Total number of time steps simulated
X <- numeric(M)  # Creates a vector of dimension N to keep the X_t's generated


# [For a stationary Jacobi process, set the initial condition as ~ Beta(a*b, a*(1-b))]:
 X[1] <- rbeta(1, a*b, a*(1-b)) 

# 3.  Euler-Maruyama scheme

for (t in 2:M) {
  dW <- rnorm(1, mean = 0, sd = sqrt(dt))  # generation of the noises
  X[t] <- X[t-1] + (-a * (X[t-1] - b) * dt) + sqrt(2 * X[t-1] * (1 - X[t-1])) * dW
  
# To keel the X_t's values in (0, 1)
  
  if (X[t] >= 1) X[t] <- 0.99
  if (X[t] <= 0) X[t] <- 0.01
}

# 4. Visualization of the simulated trajectory

plot(seq(1, M) * dt, X, type = "l", col = "green4", 
     xlab = "t", ylab = "X(t)", main = "Jacobi process",
     ylim = c(0, 1))  
abline(h = 1, col = "black", lwd = 2, lty = 2)  # Line at 1 to show the superior limit    

# Step B: A ``pre-estimation'' of  alpha and beta
# ================================================
# Fitting the Beta distribution parameters with the elements of the simulated 
# trajectory by using the method of moments (MM), which, for large times, can 
# give good approximations by the ergodicity of the process


mu <- mean(X)

var <- var(X)

s1 <- ((1 - mu) / var - 1 / mu) * mu ^ 2 

s2 <- s1 * (1 / mu - 1)

beta_prest  <- s1-1
alpha_prest <- s2-1

alpha_prest
beta_prest 


# Step C: Quasi-likelihood implementation

# ================================================

# \( p_N(x_t | x_{t-1}; \theta) = p(x_t; \theta) \left( 1 + \sum_{n=1}^N \exp(\lambda_n(\theta)) P_n(x_t; \theta) P_n(x_{t-1}; \theta) \right) \)

# We define \(\lambda_n\) y \(P_n\)
lambda_n <- function(n, alpha_o, beta_o) {
  -((alpha_o + beta_o + 2) * n + n * (n - 1))
}

P_n <- function(x, n, alpha_o, beta_o) {
  scalar <- sqrt(
    gamma(beta_o + n + 1) * (2 * n + alpha_o + beta_o + 1) * gamma(alpha_o + 1) * gamma(beta_o + 1) /
      (factorial(n) * gamma(alpha_o + beta_o + n + 1) * gamma(alpha_o + beta_o + 2) * gamma(alpha_o + n + 1))
  )
  sum_val <- 0
  for (m in 0:n) {
    binomial_coeff <- choose(n, m)
    gamma_term <- gamma(alpha_o + beta_o + n + m + 1) / gamma(beta_o + m + 1)
    sum_val <- sum_val + (-1)^m * binomial_coeff * gamma_term * x^m
  }
  return(scalar * sum_val)
}

p_N <- function(x_t, x_t_1, alpha_o, beta_o, N) {
  p_xt <- dbeta(x_t, beta_o + 1, alpha_o +1)
  sum_val <- 0
  for (n in 1:N) {
    lambda <- lambda_n(n, alpha_o, beta_o)
    Pn_xt <- P_n(x_t, n, alpha_o, beta_o)
    Pn_xt_1 <- P_n(x_t_1, n, alpha_o, beta_o)
    sum_val <- sum_val + exp(lambda) * Pn_xt * Pn_xt_1
  }
  return(p_xt * (1 + sum_val))
}


cat("p_N for N = 1, 2, 3, 4, 5\n")
for (N in 1:5) {
  p_val <- p_N(X[2], X[1], alpha_prest, beta_prest, N)
  cat("N =", N, ", p_N =", p_val, "\n")
}
cat("\n")


cat("Individual terms of the sum\n")
for (n in 1:5) {
  lambda <- lambda_n(n, alpha_prest, beta_prest)
  Pn_xt <- P_n(X[2], n, alpha_prest, beta_prest)
  Pn_xt_1 <- P_n(X[1], n, alpha_prest, beta_prest)
  cat("n =", n, ", lambda =", lambda, ", Pn_xt =", Pn_xt, ", Pn_xt_1 =", Pn_xt_1, "\n")
}
cat("\n")

#  lambda_n and Jacobi polinomials (P_n) for n = 1, ...,5
for (n in 1:5) {
  lambda <- lambda_n(n, alpha_prest, beta_prest)
  Pn_xt <- P_n(X[2], n, alpha_prest, beta_prest)
  Pn_xt_1 <- P_n(X[1], n, alpha_prest, beta_prest)
  cat("n =", n, ", lambda =", lambda, ", Pn_xt =", Pn_xt, ", Pn_xt_1 =", Pn_xt_1, "\n")
}

# Quasi-likelihood assessment

log_likelihood <- function(X, alpha_o, beta_o, N) {
  log_lik <- 0  # Initialization
  
 
  for (t in 2:length(X)) {
    
    p_xt_xt_1 <- max(p_N(X[t], X[t-1], alpha_o, beta_o, N), 1e-10)  # To avoid log(0)
    
    
    log_lik <- log_lik + log(p_xt_xt_1)
    
    
    if (t <= 5) {
      cat("t =", t, ", X[t] =", X[t], ", X[t-1] =", X[t-1], 
          ", p_N =", p_xt_xt_1, ", log(p_N) =", log(p_xt_xt_1), "\n")
    }
  }
  return(log_lik)
}


cat("Log-likelihood calculation for N = 5\n")
log_lik_val <- log_likelihood(X, alpha_prest, beta_prest, 5)
cat("Log-likelihood=", log_lik_val, "\n\n")




optim_likelihood <- function(params, X, N) {
  alpha_o <- params[1]
  beta_o <- params[2]
  log_lik_value <- tryCatch(log_likelihood(X, alpha_o, beta_o, N), error = function(e) -1e10)
  if (!is.finite(log_lik_value)) return(1e10)  # To avoid no-finite values
  return(-log_lik_value)  
}


result <- optim(par = c(alpha_prest, beta_prest),
                fn = optim_likelihood,
                X = X,
                N = 5,
                method = "L-BFGS-B",
                lower = c(1, 1),
                upper = c(50, 50))  # Constraints


cat("Optimization completed.\n")
cat("Alpha optimal =", result$par[1], ", Beta Optimal =", result$par[2], "\n")
cat("Optimal log-likelihood =", -result$value, "\n")
cat("Convergence =", result$convergence, "(0 means success)\n\n")





