# arimax.R
library(forecast)
library(dplyr)
library(zoo)

raw <- read.csv("C:/Users/User/Downloads/Data set 1 (1).csv", stringsAsFactors = FALSE)
names(raw) <- trimws(names(raw))
names(raw) <- tolower(gsub(" ", "_", names(raw)))

data_ccf <- raw %>%
  mutate(date = as.Date(date)) %>%
  arrange(date) %>%
  transmute(
    date = date,
    total = as.numeric(total),
    temperature = as.numeric(temperature),
    relative_humidity = as.numeric(relative.humidity),
    wind_speed = as.numeric(wind.speed),
    sunshine_duration = as.numeric(x.sunshine.duration)
  )

met_vars <- c("temperature", "relative_humidity", "wind_speed", "sunshine_duration")
for (v in met_vars) data_ccf[[v]] <- na.locf(data_ccf[[v]], na.rm = FALSE)

n_ccf <- floor(nrow(data_ccf) * 0.8)
train_ccf <- data_ccf[1:n_ccf, ]

png("ccf_temperature.png", width = 10, height = 6, units = "in", res = 150)
ccf(train_ccf$temperature, train_ccf$total, lag.max = 30, main = "Temperature vs Total Influenza", xlab = "Lag (days)", ylab = "CCF")
abline(v = c(-7, -1, 1, 7), col = "red", lty = 2)
dev.off()

png("ccf_humidity.png", width = 10, height = 6, units = "in", res = 150)
ccf(train_ccf$relative_humidity, train_ccf$total, lag.max = 30, main = "Relative Humidity vs Total Influenza", xlab = "Lag (days)", ylab = "CCF")
abline(v = c(-7, -1, 1, 7), col = "red", lty = 2)
dev.off()

png("ccf_wind.png", width = 10, height = 6, units = "in", res = 150)
ccf(train_ccf$wind_speed, train_ccf$total, lag.max = 30, main = "Wind Speed vs Total Influenza", xlab = "Lag (days)", ylab = "CCF")
dev.off()

png("ccf_sunshine.png", width = 10, height = 6, units = "in", res = 150)
ccf(train_ccf$sunshine_duration, train_ccf$total, lag.max = 30, main = "Sunshine Duration vs Total Influenza", xlab = "Lag (days)", ylab = "CCF")
dev.off()

armax_fc <- function(tsdata, ntrain, order, seasonal, method, traincoef, include.mean, xreg) {
  obj <- arima(tsdata, order = order, seasonal = seasonal, init = traincoef, fixed = traincoef,
               method = method, include.mean = include.mean, xreg = xreg, optim.control = list(maxit = 0))
  fc <- tsdata - obj$residuals
  ntotal <- length(tsdata)
  holdout_fc <- fc[(ntrain + 1):ntotal]
  holdout <- tsdata[(ntrain + 1):ntotal]
  list(rmse = sqrt(mean((holdout - holdout_fc)^2)), mae = mean(abs(holdout - holdout_fc)), fc = holdout_fc)
}

raw <- read.csv("C:/Users/User/Downloads/Data set 1 (1).csv", stringsAsFactors = FALSE)
names(raw) <- trimws(names(raw))
names(raw) <- tolower(gsub(" ", "_", names(raw)))

data <- raw %>%
  mutate(date = as.Date(date)) %>%
  arrange(date) %>%
  transmute(
    date = date,
    total = as.numeric(total),
    flu_a = as.numeric(influneza.a),
    flu_b = as.numeric(influneza.b),
    temperature = as.numeric(temperature),
    relative_humidity = as.numeric(relative.humidity),
    wind_speed = as.numeric(wind.speed),
    sunshine_duration = as.numeric(x.sunshine.duration),
    vapour_pressure = as.numeric(vapour.pressure)
  )

met_vars <- c("temperature", "relative_humidity", "wind_speed", "sunshine_duration", "vapour_pressure")
for (v in met_vars) data[[v]] <- na.locf(data[[v]], na.rm = FALSE)

lag_names <- c()
for (v in met_vars) {
  for (k in c(1, 7)) {
    col <- paste0(v, "_lag", k)
    data[[col]] <- lag(data[[v]], k)
    lag_names <- c(lag_names, col)
  }
}
data <- na.omit(data)

ntotal <- nrow(data)
ntrain <- floor(ntotal * 0.8)
train <- data[1:ntrain, ]
test <- data[(ntrain + 1):ntotal, ]

xreg_set_A <- c("temperature_lag1", "temperature_lag7", "relative_humidity_lag1", "relative_humidity_lag7")
xreg_set_B <- lag_names

series_configs <- list(
  Total = list(y_col = "total", order = c(3, 1, 2), seasonal = list(order = c(2, 0, 2), period = 7)),
  FluA = list(y_col = "flu_a", order = c(5, 1, 1), seasonal = list(order = c(2, 0, 1), period = 7)),
  FluB = list(y_col = "flu_b", order = c(2, 0, 3), seasonal = list(order = c(2, 1, 0), period = 7))
)

results <- data.frame()

for (series_name in names(series_configs)) {
  cfg <- series_configs[[series_name]]
  y_train <- train[[cfg$y_col]]
  y_all <- data[[cfg$y_col]]
  
  fit_pure <- arima(y_train, order = cfg$order, seasonal = cfg$seasonal, method = "CSS")
  obj_pure <- arima(y_all, order = cfg$order, seasonal = cfg$seasonal,
                    init = fit_pure$coef, fixed = fit_pure$coef, method = "CSS",
                    optim.control = list(maxit = 0))
  fc_pure <- y_all - obj_pure$residuals
  holdout_pure <- fc_pure[(ntrain + 1):ntotal]
  holdout_y <- y_all[(ntrain + 1):ntotal]
  rmse_pure <- sqrt(mean((holdout_y - holdout_pure)^2))
  mae_pure <- mean(abs(holdout_y - holdout_pure))
  
  results <- rbind(results, data.frame(Series = series_name, XregSet = "Pure SARIMA",
                                       RMSE = round(rmse_pure, 2), MAE = round(mae_pure, 2)))
  
  for (set_label in c("temp+rh", "all_5vars")) {
    cols <- if (set_label == "temp+rh") xreg_set_A else xreg_set_B
    xreg_train <- as.matrix(train[, cols])
    xreg_all <- as.matrix(data[, cols])
    
    fit_arimax <- arima(y_train, order = cfg$order, seasonal = cfg$seasonal, xreg = xreg_train, method = "CSS")
    res_arimax <- armax_fc(y_all, ntrain, cfg$order, cfg$seasonal, "CSS", fit_arimax$coef, FALSE, xreg_all)
    
    results <- rbind(results, data.frame(Series = series_name, XregSet = set_label,
                                         RMSE = round(res_arimax$rmse, 2), MAE = round(res_arimax$mae, 2)))
  }
}

print(results, row.names = FALSE)
write.csv(results, "arimax_results.csv", row.names = FALSE)
save(data, train, test, ntrain, results, file = "arimax_workspace.RData")

cat("\n✓ Saved: arimax_results.csv, arimax_workspace.RData\n")