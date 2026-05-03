# exponential_smoothing.R
library(forecast)

train <- read.csv("train.csv")
test <- read.csv("test.csv")

calc_accuracy <- function(actual, forecast) {
  data.frame(
    RMSE = round(sqrt(mean((actual - forecast)^2, na.rm = TRUE)), 2),
    MAE = round(mean(abs(actual - forecast), na.rm = TRUE), 2)
  )
}

ses_roll <- function(train, test, alpha, level) {
  n <- length(test)
  fc <- numeric(n)
  lev <- level
  fc[1] <- lev
  for (i in 2:n) {
    lev <- alpha * test[i - 1] + (1 - alpha) * lev
    fc[i] <- lev
  }
  return(fc)
}

holt_roll <- function(train, test, alpha, beta, level, slope) {
  n <- length(test)
  fc <- numeric(n)
  lev <- level
  slp <- slope
  fc[1] <- lev + slp
  for (i in 2:n) {
    lev_new <- alpha * test[i - 1] + (1 - alpha) * (lev + slp)
    slp <- beta * (lev_new - lev) + (1 - beta) * slp
    lev <- lev_new
    fc[i] <- lev + slp
  }
  return(fc)
}

hw_roll <- function(train, test, alpha, beta, gamma, level, slope, season, period = 7) {
  n <- length(test)
  fc <- numeric(n)
  lev <- level
  slp <- slope
  sea <- as.numeric(season)
  for (i in 1:n) {
    s_idx <- ((i - 1) %% period) + 1
    fc[i] <- lev + slp + sea[s_idx]
    y <- test[i]
    lev_new <- alpha * (y - sea[s_idx]) + (1 - alpha) * (lev + slp)
    slp <- beta * (lev_new - lev) + (1 - beta) * slp
    sea[s_idx] <- gamma * (y - lev_new) + (1 - gamma) * sea[s_idx]
    lev <- lev_new
  }
  return(fc)
}

run_es <- function(train_vec, test_vec, name) {
  ts_train <- ts(train_vec, frequency = 7)
  
  fit_ses <- HoltWinters(ts_train, beta = FALSE, gamma = FALSE)
  fit_holt <- HoltWinters(ts_train, gamma = FALSE)
  fit_hw <- HoltWinters(ts_train, seasonal = "additive")
  
  fc_ses <- ses_roll(train_vec, test_vec, fit_ses$alpha, fit_ses$coefficients["a"])
  fc_holt <- holt_roll(train_vec, test_vec, fit_holt$alpha, fit_holt$beta, fit_holt$coefficients["a"], fit_holt$coefficients["b"])
  fc_hw <- hw_roll(train_vec, test_vec, fit_hw$alpha, fit_hw$beta, fit_hw$gamma, fit_hw$coefficients["a"], fit_hw$coefficients["b"], fit_hw$coefficients[3:9], 7)
  
  acc_ses <- calc_accuracy(test_vec, fc_ses)
  acc_holt <- calc_accuracy(test_vec, fc_holt)
  acc_hw <- calc_accuracy(test_vec, fc_hw)
  
  params <- data.frame(
    Series = name,
    Model = c("SES", "Holt", "HW (m=7)"),
    Alpha = c(round(fit_ses$alpha, 4), round(fit_holt$alpha, 4), round(fit_hw$alpha, 4)),
    Beta = c(NA, round(fit_holt$beta, 4), round(fit_hw$beta, 4)),
    Gamma = c(NA, NA, round(fit_hw$gamma, 4)),
    RMSE = c(acc_ses$RMSE, acc_holt$RMSE, acc_hw$RMSE),
    MAE = c(acc_ses$MAE, acc_holt$MAE, acc_hw$MAE)
  )
  
  print(params)
  return(list(params = params, fc_ses = fc_ses, fc_holt = fc_holt, fc_hw = fc_hw, fit_hw = fit_hw, fit_ses = fit_ses))
}

res_total <- run_es(train$total, test$total, "Total")
res_A <- run_es(train$influenza_a, test$influenza_a, "Influenza A")
res_B <- run_es(train$influenza_b, test$influenza_b, "Influenza B")

all_results <- rbind(res_total$params, res_A$params, res_B$params)
write.csv(all_results, "es_all_results.csv", row.names = FALSE)

test_dates <- as.Date(test$date)

plot_forecast <- function(dates, actual, forecast, model_name, series_name, rmse, mae, filename) {
  png(filename, width = 10, height = 6, units = "in", res = 150)
  plot(dates, actual, type = "l", col = "black", lwd = 1.5, xlab = "Date", ylab = "Daily Cases", main = paste(series_name, "-", model_name), sub = paste("RMSE =", rmse, "| MAE =", mae))
  lines(dates, forecast, col = "red", lwd = 1.5)
  legend("topleft", legend = c("Actual", "Forecast"), col = c("black", "red"), lty = 1, cex = 0.8)
  dev.off()
}

plot_forecast(test_dates, test$total, res_total$fc_hw, "Holt-Winters (m=7)", "Total Influenza", 85.09, 39.06, "plot_total_hw.png")
plot_forecast(test_dates, test$influenza_a, res_A$fc_hw, "Holt-Winters (m=7)", "Influenza A", 76.18, 32.92, "plot_A_hw.png")
plot_forecast(test_dates, test$influenza_b, res_B$fc_ses, "SES", "Influenza B", 38.67, 10.03, "plot_B_ses.png")

png("plot_all_best.png", width = 12, height = 8, units = "in", res = 150)
par(mfrow = c(3, 1), mar = c(4, 4, 3, 2))
plot(test_dates, test$total, type = "l", col = "black", lwd = 1.5, xlab = "", ylab = "Cases", main = "Total Influenza - Holt-Winters (m=7)")
lines(test_dates, res_total$fc_hw, col = "red", lwd = 1.5)
legend("topleft", legend = c("Actual", "HW Forecast"), col = c("black", "red"), lty = 1, cex = 0.7)
plot(test_dates, test$influenza_a, type = "l", col = "black", lwd = 1.5, xlab = "", ylab = "Cases", main = "Influenza A - Holt-Winters (m=7)")
lines(test_dates, res_A$fc_hw, col = "red", lwd = 1.5)
legend("topleft", legend = c("Actual", "HW Forecast"), col = c("black", "red"), lty = 1, cex = 0.7)
plot(test_dates, test$influenza_b, type = "l", col = "black", lwd = 1.5, xlab = "Date", ylab = "Cases", main = "Influenza B - SES")
lines(test_dates, res_B$fc_ses, col = "red", lwd = 1.5)
legend("topleft", legend = c("Actual", "SES Forecast"), col = c("black", "red"), lty = 1, cex = 0.7)
par(mfrow = c(1, 1))
dev.off()

cat("\n✓ Results saved to es_all_results.csv\n")
cat("✓ Plots saved: plot_total_hw.png, plot_A_hw.png, plot_B_ses.png, plot_all_best.png\n")