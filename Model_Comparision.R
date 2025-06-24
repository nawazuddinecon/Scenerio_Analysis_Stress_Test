# ---------------------------------------------------------
# Forecast Comparison for Different τ value
# ---------------------------------------------------------

# ---- PARAMETERS ----
library(forecast)
library(Metrics)

lags <- 5
tau_vals <- c(0.1, 0.5, 10, 20, 30)
horizon <- 12
n_vars <- ncol(ts_data_stat)
results <- data.frame()

mu1 <- 1
mu6 <- 1
c1 <- 0.2

for (tau in tau_vals) {
  mu5 <- tau
  soc_strength <- mu5 * mu6
  
  T_total <- nrow(ts_data_stat)
  test_start <- T_total - horizon + 1
  actuals <- ts_data_stat[test_start:T_total, ]
  
  preds_bvar <- preds_rw <- preds_ar1 <- matrix(NA, nrow = horizon, ncol = n_vars)
  
  for (h in 1:horizon) {
    train <- ts_data_stat[1:(test_start + h - 2), , drop = FALSE]
    
    embed_data <- embed(train, lags + 1)
    Y <- embed_data[, 1:n_vars]
    X <- embed_data[, -(1:n_vars)]
    X <- cbind(1, X)
    
    # --- SOC Dummies ---
    Y_soc <- matrix(0, nrow = n_vars, ncol = n_vars)
    X_soc <- matrix(0, nrow = n_vars, ncol = ncol(X))
    for (i in 1:n_vars) {
      for (j in 1:lags) {
        X_soc[i, (i - 1) * lags + j + 1] <- soc_strength
      }
    }
    
    # ---  Dummy-Initial-Observation ---
    Y_minn <- matrix(0, nrow = n_vars * lags + 1, ncol = n_vars)
    X_minn <- matrix(0, nrow = n_vars * lags + 1, ncol = ncol(X))
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
    
    # --- Stack all ---
    X_stack <- rbind(X_soc, X_minn, X)
    Y_stack <- rbind(Y_soc, Y_minn, Y)
    
    # --- Estimate with scaled residuals (NIW style) ---
    B_hat <- solve(t(X_stack) %*% X_stack) %*% t(X_stack) %*% Y_stack
    E_hat <- Y_stack - X_stack %*% B_hat
    Sigma_hat <- (t(E_hat) %*% E_hat) / (nrow(Y_stack) - ncol(X_stack))
    Sigma_hat <- Sigma_hat / c1  # scaling like inverse-Wishart shrinkage
    
    # --- Prediction using point estimates ---
    last_lags <- c(1, as.vector(t(train[(nrow(train) - lags + 1):nrow(train), ])))
    preds_bvar[h, ] <- last_lags %*% B_hat
    
    # --- AR(1) forecasts ---
    preds_ar1[h, ] <- apply(train, 2, function(y) {
      model <- Arima(y, order = c(1, 0, 0))
      forecast(model, h = 1)$mean[1]
    })
    
    # --- RW forecast ---
    preds_rw[h, ] <- train[nrow(train), ]
  }
  
  # --- Metrics ---
  for (v in 1:n_vars) {
    varname <- colnames(ts_data_stat)[v]
    
    for (model in c("BVAR", "AR1", "RW")) {
      pred <- switch(model,
                     BVAR = preds_bvar[, v],
                     AR1 = preds_ar1[, v],
                     RW = preds_rw[, v])
      
      mae_val <- mae(actuals[, v], pred)
      rmse_val <- rmse(actuals[, v], pred)
      direction_acc <- mean(sign(diff(actuals[, v])) == sign(diff(pred))) * 100
      
      # --- Stability for BVAR only ---
      stable_pct <- NA
      if (model == "BVAR") {
        Phi_flat <- B_hat[-1, ]  # remove intercept
        A <- matrix(0, nrow = n_vars * lags, ncol = n_vars * lags)
        A[1:n_vars, ] <- t(matrix(Phi_flat, ncol = n_vars))
        if (lags > 1) A[(n_vars + 1):(n_vars * lags), 1:(n_vars * (lags - 1))] <- diag(n_vars * (lags - 1))
        eigs <- eigen(A)$values
        stable_pct <- mean(Mod(eigs) < 1) * 100
      }
      
      results <- rbind(results, data.frame(
        Tau = tau,
        Variable = varname,
        Model = model,
        MAE = round(mae_val, 4),
        RMSE = round(rmse_val, 4),
        DirectionAcc = round(direction_acc, 2),
        Stability = ifelse(is.na(stable_pct), "", round(stable_pct, 2))
      ))
    }
  }
}

# --- Output ---
print(results)
write.csv(results, "bvar_model_comparison_results.csv", row.names = FALSE)

########Summary_table
library(dplyr)

# Ensure Stability is numeric (handle empty strings if any)
results$Stability <- as.numeric(as.character(results$Stability))

# Summarise: average by Tau and Model
model_summary <- results %>%
  group_by(Tau, Model) %>%
  summarise(
    Avg_MAE = round(mean(MAE, na.rm = TRUE), 4),
    Avg_RMSE = round(mean(RMSE, na.rm = TRUE), 4),
    Avg_DirectionAcc = round(mean(DirectionAcc, na.rm = TRUE), 2),
    Avg_Stability = round(mean(Stability, na.rm = TRUE), 2),
    .groups = "drop"
  )

# Print summary
print(model_summary)

# Optional: save to CSV
write.csv(model_summary, "bvar_model_comparision_by_tau.csv", row.names = FALSE)





