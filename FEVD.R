# ---------------------------------------------------------
# FEVD for 12, 24 and 36 Months
# ---------------------------------------------------------
gfevd_compute <- function(Phi, Sigma, H) {
  n <- dim(Sigma)[1]
  p <- dim(Phi)[3]
  
  # Companion matrix
  A_comp <- matrix(0, n * p, n * p)
  for (j in 1:p) {
    A_comp[1:n, ((j - 1) * n + 1):(j * n)] <- Phi[, , j]
  }
  if (p > 1) {
    A_comp[(n + 1):(n * p), 1:(n * (p - 1))] <- diag(n * (p - 1))
  }
  
  # Impulse response matrices
  Theta <- array(0, dim = c(n, n, H + 1))  # [response, shock, horizon]
  Theta[, , 1] <- diag(n)
  A_power <- diag(n * p)
  
  for (h in 1:H) {
    A_power <- A_power %*% A_comp
    Theta[, , h + 1] <- A_power[1:n, 1:n]
  }
  
  # GFEVD matrix
  gfevd <- matrix(0, nrow = n, ncol = n)
  
  for (j in 1:n) {  # loop over shock
    for (i in 1:n) {  # loop over response
      num <- 0
      denom <- 0
      for (h in 0:H) {
        th <- Theta[, , h + 1]
        num <- num + (Sigma[j, ] %*% th[i, ])^2 / Sigma[j, j]
        denom <- denom + th[i, ] %*% Sigma %*% matrix(th[i, ], ncol = 1)
      }
      gfevd[i, j] <- num / denom
    }
  }
  
  # Normalize rows to sum to 1
  row_sums <- rowSums(gfevd)
  gfevd <- gfevd / row_sums
  
  return(gfevd)
}


gfevd_12 <- gfevd_compute(Phi, Sigma_hat, 12)
gfevd_24 <- gfevd_compute(Phi, Sigma_hat, 24)
gfevd_36 <- gfevd_compute(Phi, Sigma_hat, 36)

# Pretty tables
round(gfevd_12 * 100, 2)
round(gfevd_24 * 100, 2)
round(gfevd_36 * 100, 2)

var_names <- colnames(ts_data_stat)

colnames(gfevd_12) <- paste("Shock:", var_names)
rownames(gfevd_12) <- var_names

colnames(gfevd_24) <- paste("Shock:", var_names)
rownames(gfevd_24) <- var_names

colnames(gfevd_36) <- paste("Shock:", var_names)
rownames(gfevd_36) <- var_names


write.csv(round(gfevd_12 * 100, 2), "gfevd_12.csv")
write.csv(round(gfevd_24 * 100, 2), "gfevd_24.csv")
write.csv(round(gfevd_36 * 100, 2), "gfevd_36.csv")

kable(round(gfevd_12 * 100, 2), caption = "GFEVD at Horizon = 12")