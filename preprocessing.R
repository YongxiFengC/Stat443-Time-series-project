# preprocessing.R
library(readr)
library(dplyr)
library(lubridate)
library(forecast)


data <- read_csv("C:/Users/User/Downloads/Data set 1 (1).csv") %>%
  mutate(date = as.Date(Date)) %>%
  arrange(date) %>%
  select(date, Total, `Influneza A`, `Influneza B`)

colnames(data) <- c("date", "total", "influenza_a", "influenza_b")
data <- na.omit(data)


split_idx <- floor(nrow(data) * 0.8)
train <- data[1:split_idx, ]
test <- data[(split_idx + 1):nrow(data), ]


write.csv(train, "train.csv", row.names = FALSE)
write.csv(test, "test.csv", row.names = FALSE)
saveRDS(list(train = train, test = test), file = "data_split.RData")

stats <- data.frame(
  Series = c("Total", "Influenza A", "Influenza B"),
  Mean = round(c(mean(train$total), mean(train$influenza_a), mean(train$influenza_b)), 1),
  SD = round(c(sd(train$total), sd(train$influenza_a), sd(train$influenza_b)), 1),
  Min = c(min(train$total), min(train$influenza_a), min(train$influenza_b)),
  Max = c(max(train$total), max(train$influenza_a), max(train$influenza_b))
)
print(stats)
write.csv(stats, "summary_stats.csv", row.names = FALSE)


par(mfrow = c(3, 1), mar = c(4, 4, 2, 2))

plot(train$date, train$total, type = "l", col = "black",
     xlab = "", ylab = "Daily Cases", main = "Total")
plot(train$date, train$influenza_a, type = "l", col = "red",
     xlab = "", ylab = "Daily Cases", main = "Influenza A")
plot(train$date, train$influenza_b, type = "l", col = "blue",
     xlab = "Year", ylab = "Daily Cases", main = "Influenza B")

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))  # reset


par(mfrow = c(1, 2))
Acf(train$total, lag.max = 30, main = "Total: ACF (lags 1-30 days)")
Pacf(train$total, lag.max = 30, main = "Total: PACF (lags 1-30 days)")
par(mfrow = c(1, 1))

par(mfrow = c(1, 2))
Acf(train$influenza_a, lag.max = 30, main = "Influenza A: ACF (lags 1-30 days)")
Acf(train$influenza_b, lag.max = 30, main = "Influenza B: ACF (lags 1-30 days)")
par(mfrow = c(1, 1))

par(mfrow = c(3, 1), mar = c(4, 4, 3, 2)) 

acf_total <- Acf(train$total, lag.max = 12*365, plot = FALSE)
plot(acf_total$lag / 365, acf_total$acf, type = "h",
     xlab = "Lag (Years)", ylab = "ACF",
     main = "Total Influenza",
     lwd = 1.5, col = "steelblue")
abline(h = 0)
abline(h = c(-1.96/sqrt(length(train$total)), 1.96/sqrt(length(train$total))),
       lty = 2, col = "red")

acf_A <- Acf(train$influenza_a, lag.max = 12*365, plot = FALSE)
plot(acf_A$lag / 365, acf_A$acf, type = "h",
     xlab = "Lag (Years)", ylab = "ACF",
     main = "Influenza A",
     lwd = 1.5, col = "steelblue")
abline(h = 0)
abline(h = c(-1.96/sqrt(length(train$influenza_a)), 1.96/sqrt(length(train$influenza_a))),
       lty = 2, col = "red")

acf_B <- Acf(train$influenza_b, lag.max = 12*365, plot = FALSE)
plot(acf_B$lag / 365, acf_B$acf, type = "h",
     xlab = "Lag (Years)", ylab = "ACF",
     main = "Influenza B",
     lwd = 1.5, col = "steelblue")
abline(h = 0)
abline(h = c(-1.96/sqrt(length(train$influenza_b)), 1.96/sqrt(length(train$influenza_b))),
       lty = 2, col = "red")

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))  # reset margins

ts_365 <- ts(train$total, frequency = 365)
stl_result <- stl(ts_365, s.window = "periodic", robust = TRUE)
plot(stl_result, main = "STL Decomposition: Total Influenza (Annual Seasonality)")

cat("\nDone. Saved: train.csv, test.csv, data_split.RData, summary_stats.csv\n")
cat("Files saved to:", getwd(), "\n")

# DIFFERENCING ANALYSIS

library(tseries)

cat("\n=== ADF Tests ===\n")
for (s in c("total", "influenza_a", "influenza_b")) {
  adf_orig <- adf.test(train[[s]])
  adf_diff <- adf.test(diff(train[[s]]))
  cat(sprintf("\n%s:\n", s))
  cat(sprintf("  Original:    ADF = %.4f,  p = %.4f\n",
              adf_orig$statistic, adf_orig$p.value))
  cat(sprintf("  Differenced: ADF = %.4f,  p = %.4f\n",
              adf_diff$statistic, adf_diff$p.value))
}

par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

for (s in c("total", "influenza_a", "influenza_b")) {
  label <- switch(s,
                  total       = "Total",
                  influenza_a = "Influenza A",
                  influenza_b = "Influenza B")
  d <- diff(train[[s]])
  Acf(d,  lag.max = 300, main = paste(label, "— Differenced ACF"))
  Pacf(d, lag.max = 300, main = paste(label, "— Differenced PACF"))
}

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

par(mfrow = c(3, 1), mar = c(4, 4, 3, 2))

for (s in c("total", "influenza_a", "influenza_b")) {
  label <- switch(s,
                  total       = "Total Influenza",
                  influenza_a = "Influenza A",
                  influenza_b = "Influenza B")
  d <- diff(train[[s]])
  acf_d <- Acf(d, lag.max = 10* 365, plot = FALSE)
  plot(acf_d$lag / 365, acf_d$acf, type = "h",
       xlab = "Lag (Years)", ylab = "ACF",
       main = paste(label, "— Differenced, ACF in Years"),
       lwd = 1.2, col = "steelblue")
  abline(h = 0)
  abline(h = c(-1.96, 1.96) / sqrt(length(d)),
         lty = 2, col = "red")
}

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

cat("\nDifferencing analysis done.\n")