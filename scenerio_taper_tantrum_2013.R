rm(list = ls())

# -------------------------------
# User-defined parameters
# -------------------------------
library(vars)

lags <- 5  # Based on BIC criterion (see below)
mu1 <- 1     # Tightness on AR(1) own lag
mu5 <- 1     # SOC prior tightness τ
mu6 <- 1     # Weight on dummy observations

c1 <- 0.2      # Scale of prior covariance matrix (s scale)
c2 <- 0.2    # Not used currently
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

pre_taper_data <- window(ts_data_stat, end = c(2013, 4))  # Ends at Apr 2013
tsp(pre_taper_data)
print(pre_taper_data)

n_vars <- ncol(pre_taper_data)
T_obs <- nrow(pre_taper_data)
# -------------------------------
# Rename variables
# -------------------------------
colnames(pre_taper_data) <- c(
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
# Step 2: Single-unit-root dummy-observation
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
X <- embed(pre_taper_data, lags + 1)[, -c(1:n_vars)]
Y_response <- embed(pre_taper_data, lags + 1)[, 1:n_vars]
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
#____________________*Plot one IRF: log_loanhh and log_GDP_detrended*_________________

#**skip this segment of final_g

# ----------------------------------------------------------
# Step 1: Define 2013 taper tantrum shock vector
# ----------------------------------------------------------
shock_2013 <- c(
  "Real GDP" = -0.02467,
  "WPI" = 2.65361,
  "Short Term Interest Rate" = 2.6162,
  "10-Year Bond Yield" = 0.644,
  "Stock Market Index" = -0.07312
)

shock_vector <- rep(0, n_vars)
for (name in names(shock_2013)) {
  idx <- which(colnames(pre_taper_data) == name)
  if (length(idx) == 1) shock_vector[idx] <- shock_2013[name]
}

# ----------------------------------------------------------
# Step 2: Conditional Forecasting Parameters
# ----------------------------------------------------------
forecast_horizon <- 36  # or 36 for 3 years, adjust as needed
n_draws <- 1000
cond_forecast_array <- array(NA, dim = c(forecast_horizon + 1, n_vars, n_draws))

# ----------------------------------------------------------
# Step 3: Simulate forecast paths under the taper shock
# ----------------------------------------------------------
for (d in 1:n_draws) {
  Sigma_draw <- solve(rwish(df_post, solve(Sigma_hat)))
  chol_Sigma <- t(chol(Sigma_draw))
  shock_draw <- chol_Sigma %*% shock_vector
  
  state <- matrix(0, nrow = n_vars * lags, ncol = 1)
  Y_path <- matrix(0, nrow = forecast_horizon + 1, ncol = n_vars)
  Y_path[1, ] <- shock_draw
  
  for (h in 1:forecast_horizon) {
    state[1:n_vars, 1] <- Y_path[h, ]
    state <- A_comp %*% state
    Y_path[h + 1, ] <- state[1:n_vars, 1]
  }
  
  cond_forecast_array[, , d] <- Y_path
}

# ----------------------------------------------------------
# Step 4: Compute forecast mean and credible intervals
# ----------------------------------------------------------
forecast_mean <- apply(cond_forecast_array, c(1, 2), mean)
forecast_lower <- apply(cond_forecast_array, c(1, 2), quantile, probs = 0.05)
forecast_upper <- apply(cond_forecast_array, c(1, 2), quantile, probs = 0.95)

# ----------------------------------------------------------
# Step 5: Plot the conditional forecasts for banking variables
# ----------------------------------------------------------
bank_vars <- c("Real GDP", "WPI","Short Term Interest Rate",
               "5-Year Bond Yield", "10-Year Bond Yield", "Yield Spread", "Stock Market Index",
               "Loans to NFCs", "Loans to HHs", "NPA", "Stock of Bonds", "Lending Rate", "Capital Adequacy Ratio" )

par(mfrow = c(3, 2))
for (vname in bank_vars) {
  var_idx <- which(colnames(pre_taper_data) == vname)
  plot(0:forecast_horizon, forecast_mean[, var_idx], type = "l", col = "black",
       ylim = range(forecast_lower[, var_idx], forecast_upper[, var_idx]),
       main = paste("Forecast:", vname),
       xlab = "Horizon (Months)", ylab = "Response")
  lines(0:forecast_horizon, forecast_lower[, var_idx], col = "blue", lty = 2)
  lines(0:forecast_horizon, forecast_upper[, var_idx], col = "blue", lty = 2)
  abline(h = 0, lty = 2, col = "gray")
}
par(mfrow = c(1, 1))


# ----------------------------------------------------------
# Step 6: Compute realisation probabilities (vs zero)
# ----------------------------------------------------------
realisation_probs <- matrix(NA, nrow = forecast_horizon + 1, ncol = length(bank_vars))
colnames(realisation_probs) <- bank_vars
rownames(realisation_probs) <- 0:forecast_horizon

for (i in seq_along(bank_vars)) {
  var_idx <- which(colnames(pre_taper_data) == bank_vars[i])
  for (h in 0:forecast_horizon) {
    realisation_probs[h + 1, i] <- mean(cond_forecast_array[h + 1, var_idx, ] < 0) * 100
  }
}

realisation_probs

# ----------------------------------------------------------
# Step 7: Combine mean and probability for tabular display
# ----------------------------------------------------------
forecast_summary <- list()

for (i in seq_along(bank_vars)) {
  var_idx <- which(colnames(pre_taper_data) == bank_vars[i])
  mean_path <- round(forecast_mean[, var_idx], 4)
  prob_path <- round(realisation_probs[, i], 1)
  forecast_summary[[bank_vars[i]]] <- paste0(mean_path, " (", prob_path, "%)")
}

# Convert to data.frame for nicer display
summary_df <- as.data.frame(forecast_summary)
rownames(summary_df) <- paste0("Month_", 0:forecast_horizon)
summary_df


############Year_wise_results:
# Define the month ranges
horizons <- list(
  `1yr` = 1:12,
  `2yr` = 13:24,
  `3yr` = 25:36
)

# Initialize results
summary_years <- data.frame(Horizon = names(horizons))

for (vname in bank_vars) {
  means <- sapply(horizons, function(h) mean(forecast_mean[h, which(colnames(pre_taper_data) == vname)]))
  probs <- sapply(horizons, function(h) mean(realisation_probs[h, which(bank_vars == vname)]))
  
  summary_years[[paste0(vname, "_summary")]] <- paste0(
    round(means, 4), " (", round(probs, 1), "%)"
  )
}
summary_years
str(summary_years)
write.csv(summary_years, "shock_realization_probabilities_2013.csv", row.names = FALSE)
