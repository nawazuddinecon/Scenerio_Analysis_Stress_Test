# GIRF Confidence Bands and Plot
# -------------------------------


##Graph 1: Shock in Stock Market Index: Create Graph GIRF_Shock_StockMarket.png
# Define shock variable (same for all responses)
shock_var <- which(colnames(ts_data_stat) == "Stock Market Index")

# Define multiple response variables
response_vars <- which(colnames(ts_data_stat) %in% c(
  "Stock Market Index", "Yield Spread", "Loans to HHs",
  "Loans to NFCs", "NPA", "Lending Rate", "Markup on Short-term loans", "Interest Rate on G-Secs"
))

# Loop over response variables and store GIRFs
girf_list <- list()

for (rv in response_vars) {
  girf_matrix <- array(0, dim = c(n_draws, horizon + 1))
  
  for (d in 1:n_draws) {
    Sigma_draw <- solve(rwish(df_post, solve(Sigma_hat)))
    shock_vec <- Sigma_draw[, shock_var] / sqrt(Sigma_draw[shock_var, shock_var])
    
    state <- matrix(0, nrow = n_vars * lags, ncol = 1)
    state[1:n_vars, 1] <- shock_vec
    response <- numeric(horizon + 1)
    response[1] <- shock_vec[rv]
    
    for (t in 1:horizon) {
      state <- A_comp %*% state
      response[t + 1] <- state[rv, 1]
    }
    
    girf_matrix[d, ] <- response
  }
  
  # Store mean + CI in list
  girf_list[[colnames(ts_data_stat)[rv]]] <- list(
    mean = apply(girf_matrix, 2, median),
    lower = apply(girf_matrix, 2, quantile, probs = 0.05),
    upper = apply(girf_matrix, 2, quantile, probs = 0.95)
  )
}

# Set layout: 4 rows x 2 columns with outer title space
png("GIRF_Shock_StockMarket.png", width = 664, height = 664)

# 1. Reduce vertical spacing between graphs using smaller outer margin (oma)
# 2. Use tighter plot margins (mar)
# 3. Ensure x-axis goes up to 36
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 4, 0))  

for (varname in names(girf_list)) {
  girf_mean <- girf_list[[varname]]$mean
  girf_lower <- girf_list[[varname]]$lower
  girf_upper <- girf_list[[varname]]$upper
  
  plot(0:horizon, girf_mean, type = "l", lwd = 2, col = "black",
       ylim = range(girf_lower, girf_upper),
       ylab = "", 
       xlab = "", 
       main = paste(varname),
       xaxt = "n", yaxt = "n")
  
  axis(1)  
  axis(2, las = 1)                       # Horizontal y-axis labels
  
  lines(0:horizon, girf_lower, col = "blue", lty = 2)
  lines(0:horizon, girf_upper, col = "blue", lty = 2)
  abline(h = 0, lty = 2, col = "gray")
  
  # Bold x-axis label "Horizon"
  mtext("Horizon", side = 1, line = 2.5, font = 2)
}

# Grand title
mtext("Shock in Stock Market Index", outer = TRUE, cex = 1, font = 2)

dev.off()
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
##Graph 2: Shock in 5-Year Bond Yield: Create Graph GIRF_Shock_5YearYield.png
# Define shock variable (same for all responses)
shock_var <- which(colnames(ts_data_stat) == "5-Year Bond Yield")

# Define multiple response variables
response_vars <- which(colnames(ts_data_stat) %in% c(
  "5-Year Bond Yield", "Yield Spread", "Loans to HHs",
  "Loans to NFCs","Loans to Other Financial Intermediaries", "NPA", "Lending Rate", "Interest Rate on G-Secs"
))

# Loop over response variables and store GIRFs
girf_list <- list()

for (rv in response_vars) {
  girf_matrix <- array(0, dim = c(n_draws, horizon + 1))
  
  for (d in 1:n_draws) {
    Sigma_draw <- solve(rwish(df_post, solve(Sigma_hat)))
    shock_vec <- Sigma_draw[, shock_var] / sqrt(Sigma_draw[shock_var, shock_var])
    
    state <- matrix(0, nrow = n_vars * lags, ncol = 1)
    state[1:n_vars, 1] <- shock_vec
    response <- numeric(horizon + 1)
    response[1] <- shock_vec[rv]
    
    for (t in 1:horizon) {
      state <- A_comp %*% state
      response[t + 1] <- state[rv, 1]
    }
    
    girf_matrix[d, ] <- response
  }
  
  # Store mean + CI in list
  girf_list[[colnames(ts_data_stat)[rv]]] <- list(
    mean = apply(girf_matrix, 2, median),
    lower = apply(girf_matrix, 2, quantile, probs = 0.05),
    upper = apply(girf_matrix, 2, quantile, probs = 0.95)
  )
}

# Set layout: 4 rows x 2 columns with outer title space
png("GIRF_Shock_5YearYield.png", width = 664, height = 664)

# 1. Reduce vertical spacing between graphs using smaller outer margin (oma)
# 2. Use tighter plot margins (mar)
# 3. Ensure x-axis goes up to 36
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 4, 0))  

for (varname in names(girf_list)) {
  girf_mean <- girf_list[[varname]]$mean
  girf_lower <- girf_list[[varname]]$lower
  girf_upper <- girf_list[[varname]]$upper
  
  plot(0:horizon, girf_mean, type = "l", lwd = 2, col = "black",
       ylim = range(girf_lower, girf_upper),
       ylab = "", 
       xlab = "", 
       main = paste(varname),
       xaxt = "n", yaxt = "n")
  
  axis(1, at = seq(0, horizon, by = 3))  # Show x-axis ticks every 3 months up to 36
  axis(2, las = 1)                       # Horizontal y-axis labels
  
  lines(0:horizon, girf_lower, col = "blue", lty = 2)
  lines(0:horizon, girf_upper, col = "blue", lty = 2)
  abline(h = 0, lty = 2, col = "gray")
  
  # Bold x-axis label "Horizon"
  mtext("Horizon", side = 1, line = 2.5, font = 2)
}

# Grand title
mtext("Shock in 5-Year Bond Yield", outer = TRUE, cex = 1, font = 2)

dev.off()

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
##Graph 3:: Shock in Real GDP: Create Graph GIRF_Shock_RealGDP.png
# Define shock variable (same for all responses)
shock_var <- which(colnames(ts_data_stat) == "Real GDP")

# Define multiple response variables
response_vars <- which(colnames(ts_data_stat) %in% c(
  "Real GDP", "Yield Spread", "Loans to HHs",
  "Loans to NFCs","Loans to Other Financial Intermediaries", "NPA", "Lending Rate", "Interest Rate on G-Secs"
))

# Loop over response variables and store GIRFs
girf_list <- list()

for (rv in response_vars) {
  girf_matrix <- array(0, dim = c(n_draws, horizon + 1))
  
  for (d in 1:n_draws) {
    Sigma_draw <- solve(rwish(df_post, solve(Sigma_hat)))
    shock_vec <- Sigma_draw[, shock_var] / sqrt(Sigma_draw[shock_var, shock_var])
    
    state <- matrix(0, nrow = n_vars * lags, ncol = 1)
    state[1:n_vars, 1] <- shock_vec
    response <- numeric(horizon + 1)
    response[1] <- shock_vec[rv]
    
    for (t in 1:horizon) {
      state <- A_comp %*% state
      response[t + 1] <- state[rv, 1]
    }
    
    girf_matrix[d, ] <- response
  }
  
  # Store mean + CI in list
  girf_list[[colnames(ts_data_stat)[rv]]] <- list(
    mean = apply(girf_matrix, 2, median),
    lower = apply(girf_matrix, 2, quantile, probs = 0.05),
    upper = apply(girf_matrix, 2, quantile, probs = 0.95)
  )
}

# Set layout: 4 rows x 2 columns with outer title space
png("GIRF_Shock_RealGDP.png", width = 664, height = 664)

# 1. Reduce vertical spacing between graphs using smaller outer margin (oma)
# 2. Use tighter plot margins (mar)
# 3. Ensure x-axis goes up to 36
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 4, 0))  

for (varname in names(girf_list)) {
  girf_mean <- girf_list[[varname]]$mean
  girf_lower <- girf_list[[varname]]$lower
  girf_upper <- girf_list[[varname]]$upper
  
  plot(0:horizon, girf_mean, type = "l", lwd = 2, col = "black",
       ylim = range(girf_lower, girf_upper),
       ylab = "", 
       xlab = "", 
       main = paste(varname),
       xaxt = "n", yaxt = "n")
  
  axis(1, at = seq(0, horizon, by = 3))  # Show x-axis ticks every 3 months up to 36
  axis(2, las = 1)                       # Horizontal y-axis labels
  
  lines(0:horizon, girf_lower, col = "blue", lty = 2)
  lines(0:horizon, girf_upper, col = "blue", lty = 2)
  abline(h = 0, lty = 2, col = "gray")
  
  # Bold x-axis label "Horizon"
  mtext("Horizon", side = 1, line = 2.5, font = 2)
}

# Grand title
mtext("Shock in Real GDP", outer = TRUE, cex = 1, font = 2)

dev.off()

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
##Graph 4: Shock in Yield Spread: Create Graph GIRF_Shock_Yield Spread.png
# Define shock variable (same for all responses)
shock_var <- which(colnames(ts_data_stat) == "Yield Spread")

# Define multiple response variables
response_vars <- which(colnames(ts_data_stat) %in% c(
  "Yield Spread", "Real GDP", "Loans to HHs",
  "Loans to NFCs","Loans to Other Financial Intermediaries", "NPA", "Lending Rate", "Interest Rate on G-Secs"
))

# Loop over response variables and store GIRFs
girf_list <- list()

for (rv in response_vars) {
  girf_matrix <- array(0, dim = c(n_draws, horizon + 1))
  
  for (d in 1:n_draws) {
    Sigma_draw <- solve(rwish(df_post, solve(Sigma_hat)))
    shock_vec <- Sigma_draw[, shock_var] / sqrt(Sigma_draw[shock_var, shock_var])
    
    state <- matrix(0, nrow = n_vars * lags, ncol = 1)
    state[1:n_vars, 1] <- shock_vec
    response <- numeric(horizon + 1)
    response[1] <- shock_vec[rv]
    
    for (t in 1:horizon) {
      state <- A_comp %*% state
      response[t + 1] <- state[rv, 1]
    }
    
    girf_matrix[d, ] <- response
  }
  
  # Store mean + CI in list
  girf_list[[colnames(ts_data_stat)[rv]]] <- list(
    mean = apply(girf_matrix, 2, median),
    lower = apply(girf_matrix, 2, quantile, probs = 0.05),
    upper = apply(girf_matrix, 2, quantile, probs = 0.95)
  )
}

# Set layout: 4 rows x 2 columns with outer title space
png("GIRF_Shock_Yield Spread.png", width = 664, height = 664)

# 1. Reduce vertical spacing between graphs using smaller outer margin (oma)
# 2. Use tighter plot margins (mar)
# 3. Ensure x-axis goes up to 36
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 4, 0))  

for (varname in names(girf_list)) {
  girf_mean <- girf_list[[varname]]$mean
  girf_lower <- girf_list[[varname]]$lower
  girf_upper <- girf_list[[varname]]$upper
  
  plot(0:horizon, girf_mean, type = "l", lwd = 2, col = "black",
       ylim = range(girf_lower, girf_upper),
       ylab = "", 
       xlab = "", 
       main = paste(varname),
       xaxt = "n", yaxt = "n")
  
  axis(1, at = seq(0, horizon, by = 3))  # Show x-axis ticks every 3 months up to 36
  axis(2, las = 1)                       # Horizontal y-axis labels
  
  lines(0:horizon, girf_lower, col = "blue", lty = 2)
  lines(0:horizon, girf_upper, col = "blue", lty = 2)
  abline(h = 0, lty = 2, col = "gray")
  
  # Bold x-axis label "Horizon"
  mtext("Horizon", side = 1, line = 2.5, font = 2)
}

# Grand title
mtext("Shock in Yield Spread", outer = TRUE, cex = 1, font = 2)

dev.off()
