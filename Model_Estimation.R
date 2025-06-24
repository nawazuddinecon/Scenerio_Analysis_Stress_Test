# ---------------------------------------------------------
# Model Estimation
# ---------------------------------------------------------

rm(list = ls())

# -------------------------------
# Define parameters
# -------------------------------
library(vars)

lags <- 5  # Based on BIC criterion (see below)
mu1 <- 1     # Tightness on AR(1) own lag
mu5 <- 1     # SOC prior tightness τ
mu6 <- 1     # Weight on dummy observations

c1 <- 0.2      # Scale of prior covariance matrix (s scale)
c2 <- 0.2    # Not used currently (it is there in the E-views setup)
c3 <- 20     # Degrees of freedom for IW prior

# -------------------------------
# Prepare data
# -------------------------------
library(readxl)
library(dplyr)

data <- read_excel("data_for_model_stationary.xlsx", sheet = "Sheet1")
data <- na.omit(data)
data_matrix <- as.matrix(data)
ts_data_stat <- ts(data_matrix, start = c(2007, 5), frequency = 12)

n_vars <- ncol(ts_data_stat)
T_obs <- nrow(ts_data_stat)
# -------------------------------
# Rename variables
# -------------------------------
colnames(ts_data_stat) <- c(
  "Markup on Short-term loans",
  "Loans to NFCs",
  "Loans to Other Financial Intermediaries",
  "Bank Time Deposit",
  "NPA",
  "Real GDP",
  "Yield Spread",
  "WPI",
  "Short Term Interest Rate",
  "5-Year Bond Yield",
  "10-Year Bond Yield",
  "Lending Rate",
  "Capital Adequacy Ratio",
  "Interest Rate on G-Secs",
  "Stock Market Index",
  "Loans to HHs",
  "Stock of Bonds",
  "Bank Demand Deposits",
  "Stock of Indian Debt Bonds Held By Banks"
)

# -------------------------------
# Step 1: Sum-of-coefficients (SOC) dummy
# -------------------------------
soc_strength <- mu5 * mu6
Y_soc <- matrix(0, nrow = n_vars, ncol = n_vars)
X_soc <- matrix(0, nrow = n_vars, ncol = n_vars * lags + 1)
for (i in 1:n_vars) {
  for (j in 1:lags) {
    X_soc[i, (i - 1) * lags + j + 1] <- soc_strength
  }
}

# -------------------------------
# Step 2: Dummy-Initial-Observation
# -------------------------------
Y_minn <- matrix(0, nrow = n_vars * lags + 1, ncol = n_vars)
X_minn <- matrix(0, nrow = n_vars * lags + 1, ncol = n_vars * lags + 1)
for (i in 1:n_vars) {
  for (j in 1:lags) {
    row_idx <- (i - 1) * lags + j
    col_idx <- (i - 1) * lags + j + 1
    tightness <- mu1 / j
    X_minn[row_idx, col_idx] <- mu6 * tightness
    if (j == 1) Y_minn[row_idx, i] <- mu6 * tightness
  }
}
X_minn[n_vars * lags + 1, 1] <- mu6 * mu1 * 100

# -------------------------------
# Step 3: Real data matrices
# -------------------------------
X <- embed(ts_data_stat, lags + 1)[, -c(1:n_vars)]
Y_response <- embed(ts_data_stat, lags + 1)[, 1:n_vars]
X <- cbind(1, X)

# -------------------------------
# Step 4: Stack data
# -------------------------------
X_stack <- rbind(X_soc, X_minn, X)
Y_stack <- rbind(Y_soc, Y_minn, Y_response)

# -------------------------------
# Step 5: Estimate VAR with conjugate priors
# -------------------------------
B_hat <- solve(t(X_stack) %*% X_stack) %*% t(X_stack) %*% Y_stack
E_hat <- Y_stack - X_stack %*% B_hat
Sigma_hat <- (t(E_hat) %*% E_hat) / (nrow(Y_stack) - ncol(X_stack))
Sigma_hat <- Sigma_hat / c1

# -------------------------------
# Step 6: VAR Coefficients as Phi matrices
# -------------------------------
p <- lags
k <- ncol(X)
B_mat <- matrix(B_hat[-1, ], nrow = n_vars * p, byrow = FALSE)
Phi <- array(0, dim = c(n_vars, n_vars, p))
for (i in 1:p) {
  Phi[, , i] <- t(B_mat[((i - 1) * n_vars + 1):(i * n_vars), ])
}

# -------------------------------
# Step 7: Companion form & GIRFs
# -------------------------------
A_comp <- matrix(0, nrow = n_vars * p, ncol = n_vars * p)
for (i in 1:p) {
  A_comp[1:n_vars, ((i - 1) * n_vars + 1):(i * n_vars)] <- Phi[, , i]
}
if (p > 1) {
  A_comp[(n_vars + 1):(n_vars * p), 1:(n_vars * (p - 1))] <- diag(n_vars * (p - 1))
}

horizon <- 36

# -------------------------------
# Step 8: Posterior Sampling using NIW prior
# -------------------------------
library(MASS)
df_post <- nrow(Y_stack) - ncol(X_stack)
post_df <- c3 + df_post
S_0 <- Sigma_hat
rwish <- function(v, S) {
  return(rWishart(1, df = v, Sigma = solve(S))[, , 1])
}

n_draws <- 1000
shock_var <- which(colnames(ts_data_stat) == "Stock Market Index")
response_var <- which(colnames(ts_data_stat) == "NPA")
girf_matrix <- array(0, dim = c(n_draws, horizon + 1))

for (d in 1:n_draws) {
  Sigma_draw <- solve(rwish(post_df, solve(S_0)))
  shock_vec <- Sigma_draw[, shock_var] / sqrt(Sigma_draw[shock_var, shock_var])
  
  state <- matrix(0, nrow = n_vars * lags, ncol = 1)
  state[1:n_vars, 1] <- shock_vec
  response <- numeric(horizon + 1)
  response[1] <- shock_vec[response_var]
  
  for (t in 1:horizon) {
    state <- A_comp %*% state
    response[t + 1] <- state[response_var, 1]
  }
  
  girf_matrix[d, ] <- response
}
