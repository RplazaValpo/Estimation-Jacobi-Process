# STEP A: Data simulation by using Euler-Maruyama method
# ================================================
# We simulate the trajectory of the Jacobi process:
# dX_t = -a(X_t - b) dt + sqrt(2 * X_t * (1 - X_t)) dB_t,
# where a = alpha + beta + 2 and b = (beta + 1) / (alpha + beta + 2).


# 1. Parameter values

alpha <- 5  # Enter a value
beta  <- 4 # Enter a value
a <- alpha + beta + 2
b <- (beta + 1) / (alpha + beta + 2)

# 2. Setup

set.seed(123)  # For reproducibility
T <- 100       # Time horizon
dt <- 0.01     # Time step of the process discretization
M <- T/dt      # Total number of time steps simulated
X <- numeric(M)  # Creates a vector of dimension N to keep the X_t's generated

# Initial condition ~ Beta(1, 1)=Unif(0,1)

X[1] <- rbeta(1, 1, 1)  

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

# PASO B: A ``pre-estimation'' of  alpha and beta
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

