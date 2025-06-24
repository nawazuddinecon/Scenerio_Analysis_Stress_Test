# ---------------------------------------------------------
# Direction Accuracey for Different τ=10
# ---------------------------------------------------------

library(Metrics)

# --- Parameters ---
lags <- 5
tau <- 10
mu5 <- tau
mu1 <- 1
mu6 <- 1
c1 <- 0.2
soc_strength <- mu5 * mu6
horizons <- c(3, 12, 24, 36)

n_vars <- ncol(ts_data_stat)
T_total <- nrow(ts_data_stat)

direction_results <- data.frame()

for (H in horizons) {
  test_start <- T_total - H + 1
  actuals <- ts_data_stat[test_start:T_total, ]
  preds_bvar <- matrix(NA, nrow = H, ncol = n_vars)
  
  for (h in 1:H) {
    train <- ts_data_stat[1:(test_start + h - 2), ]
    
    # Design matrix
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
    
    # --- Minnesota Dummies ---
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
    
    # --- Stack all data ---
    X_stack <- rbind(X_soc, X_minn, X)
    Y_stack <- rbind(Y_soc, Y_minn, Y)
    
    # --- OLS estimation ---
    B_hat <- solve(t(X_stack) %*% X_stack) %*% t(X_stack) %*% Y_stack
    
    # --- Forecast next period ---
    x_pred <- c()
    for (j in 1:lags) {
      x_pred <- c(x_pred, train[nrow(train) - j + 1, ])
    }
    x_pred <- c(1, x_pred)
    preds_bvar[h, ] <- x_pred %*% B_hat
  }
  
  # --- Direction Accuracy ---
  for (v in 1:n_vars) {
    act_diff <- diff(actuals[, v])
    pred_diff <- diff(preds_bvar[, v])
    valid_idx <- !is.na(act_diff) & !is.na(pred_diff)
    dir_acc <- mean(sign(act_diff[valid_idx]) == sign(pred_diff[valid_idx])) * 100
    
    direction_results <- rbind(direction_results, data.frame(
      Variable = colnames(ts_data_stat)[v],
      Horizon = H,
      DirectionAccuracy = round(dir_acc, 2)
    ))
  }
}

# --- Output ---
print(direction_results)
write.csv(direction_results, "bvar_model_direction_results.csv", row.names = FALSE)