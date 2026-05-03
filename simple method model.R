# baseline_models.R
library(forecast)

train <- read.csv("train.csv")
test <- read.csv("test.csv")

calc_accuracy <- function(actual, forecast) {
  data.frame(
    RMSE = round(sqrt(mean((actual - forecast)^2, na.rm = TRUE)), 2),
    MAE = round(mean(abs(actual - forecast), na.rm = TRUE), 2)
  )
}

persist_fc <- function(train, test) {
  n_test <- length(test)
  fc <- numeric(n_test)
  fc[1] <- tail(train, 1)
  for (i in 2:n_test) fc[i] <- test[i - 1]
  return(fc)
}

mean_fc <- function(train, test) {
  return(rep(mean(train), length(test)))
}

snaive7_fc <- function(train, test) {
  n_train <- length(train)
  n_test <- length(test)
  fc <- numeric(n_test)
  for (i in 1:n_test) {
    idx <- n_train + i - 7
    if (idx <= n_train) fc[i] <- train[idx]
    else fc[i] <- test[idx - n_train]
  }
  return(fc)
}

snaive365_fc <- function(train, test) {
  n_train <- length(train)
  n_test <- length(test)
  fc <- numeric(n_test)
  for (i in 1:n_test) {
    idx <- n_train + i - 365
    if (idx <= n_train && idx >= 1) fc[i] <- train[idx]
    else if (idx > n_train) fc[i] <- test[idx - n_train]
    else fc[i] <- train[1]
  }
  return(fc)
}

ar1_fc <- function(train, test) {
  n_train <- length(train)
  ar1_model <- arima(train, order = c(1, 0, 0))
  coef_ar1 <- ar1_model$coef[1]
  coef_intercept <- ar1_model$coef[2]
  
  n_test <- length(test)
  fc <- numeric(n_test)
  fc[1] <- coef_intercept + coef_ar1 * train[n_train]
  for (i in 2:n_test) fc[i] <- coef_intercept + coef_ar1 * test[i - 1]
  
  return(list(fc = fc, ar1 = coef_ar1))
}

run_baselines <- function(train_vec, test_vec, name) {
  cat("\n=== ", name, " ===\n")
  
  fc_persist <- persist_fc(train_vec, test_vec)
  fc_mean <- mean_fc(train_vec, test_vec)
  fc_snaive7 <- snaive7_fc(train_vec, test_vec)
  fc_snaive365 <- snaive365_fc(train_vec, test_vec)
  ar1_result <- ar1_fc(train_vec, test_vec)
  
  results <- data.frame(
    Series = name,
    Model = c("Persistence", "Mean", "Seasonal Naive (7)", "Seasonal Naive (365)", "AR(1)"),
    RMSE = c(
      calc_accuracy(test_vec, fc_persist)$RMSE,
      calc_accuracy(test_vec, fc_mean)$RMSE,
      calc_accuracy(test_vec, fc_snaive7)$RMSE,
      calc_accuracy(test_vec, fc_snaive365)$RMSE,
      calc_accuracy(test_vec, ar1_result$fc)$RMSE
    ),
    MAE = c(
      calc_accuracy(test_vec, fc_persist)$MAE,
      calc_accuracy(test_vec, fc_mean)$MAE,
      calc_accuracy(test_vec, fc_snaive7)$MAE,
      calc_accuracy(test_vec, fc_snaive365)$MAE,
      calc_accuracy(test_vec, ar1_result$fc)$MAE
    )
  )
  print(results)
  
  cat("\nAR(1) coefficient for", name, ":", round(ar1_result$ar1, 4), "\n")
  
  return(list(results = results, ar1 = ar1_result$ar1))
}

res_total <- run_baselines(train$total, test$total, "Total")
res_A <- run_baselines(train$influenza_a, test$influenza_a, "Influenza A")
res_B <- run_baselines(train$influenza_b, test$influenza_b, "Influenza B")

ar1_table <- data.frame(
  Series = c("Total Influenza", "Influenza A", "Influenza B"),
  AR1_Coefficient = c(res_total$ar1, res_A$ar1, res_B$ar1)
)

cat("\n=== AR(1) Coefficient Table ===\n")
print(ar1_table)
write.csv(ar1_table, "ar1_coefficients.csv", row.names = FALSE)

all_results <- rbind(res_total$results, res_A$results, res_B$results)
write.csv(all_results, "baseline_results_all.csv", row.names = FALSE)

cat("\n✓ Files saved:\n")
cat("  - baseline_results_all.csv\n")
cat("  - ar1_coefficients.csv\n")