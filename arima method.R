# arima.R
library(forecast)

train <- read.csv("train.csv")
test  <- read.csv("test.csv")

calc_acc <- function(actual, forecast) {
  data.frame(
    RMSE = round(sqrt(mean((actual - forecast)^2, na.rm = TRUE)), 2),
    MAE  = round(mean(abs(actual - forecast),     na.rm = TRUE), 2)
  )
}

arma_fc <- function(y_all, ntrain, order, seasonal = list(order=c(0,0,0), period=7),
                    method = "CSS", traincoef) {
  obj <- arima(y_all,
               order    = order,
               seasonal = seasonal,
               init     = traincoef,
               fixed    = traincoef,
               method   = method,
               optim.control = list(maxit = 0))
  fc     <- y_all - obj$residuals
  ntotal <- length(y_all)
  list(
    fc   = fc[(ntrain + 1):ntotal],
    rmse = sqrt(mean((y_all[(ntrain+1):ntotal] - fc[(ntrain+1):ntotal])^2)),
    mae  = mean(abs(y_all[(ntrain+1):ntotal]   - fc[(ntrain+1):ntotal]))
  )
}

run_arima <- function(train_vec, test_vec, name) {
  cat(sprintf("\n%s\n%s\n", strrep("=",50), name))
  
  y_all  <- c(train_vec, test_vec)
  ntrain <- length(train_vec)
  ts_train <- ts(train_vec, frequency = 7)
  
  if (name == "Influenza A") {
    fit_auto <- auto.arima(ts_train,
                           seasonal = TRUE,
                           lambda = "auto",
                           stepwise = FALSE,        # 关闭逐步搜索，做全搜索
                           approximation = FALSE,    # 关闭近似，精确计算
                           max.p = 5, max.q = 5,     # 允许更大的 p,q
                           max.P = 2, max.Q = 2)     # 允许更大的 P,Q
  } else {
    fit_auto <- auto.arima(ts_train, 
                           seasonal = TRUE, 
                           lambda = "auto")
  }
  cat("\nauto.arima selected:\n"); print(fit_auto)
  
  p <- fit_auto$arma[1]; d <- fit_auto$arma[6]; q <- fit_auto$arma[2]
  P <- fit_auto$arma[3]; D <- fit_auto$arma[7]; Q <- fit_auto$arma[4]
  m <- fit_auto$arma[5]
  auto_order    <- c(p, d, q)
  auto_seasonal <- list(order = c(P, D, Q), period = m)
  
  candidates <- list(
    list(label = "ARIMA(1,1,0)",
         order = c(1,1,0), seas = list(order=c(0,0,0), period=7)),
    list(label = "ARIMA(1,0,1)",
         order = c(1,0,1), seas = list(order=c(0,0,0), period=7)),
    list(label = "SARIMA(1,0,0)(1,0,0)[7]",
         order = c(1,0,0), seas = list(order=c(1,0,0), period=7)),
    list(label = sprintf("auto.arima: ARIMA(%d,%d,%d)(%d,%d,%d)[%d]",
                         p,d,q,P,D,Q,m),
         order = auto_order, seas = auto_seasonal)
  )
  
  results <- data.frame()
  
  for (cand in candidates) {
    # Fit on training set
    fit_tr <- tryCatch(
      arima(train_vec, order=cand$order, seasonal=cand$seas, method="CSS"),
      error = function(e) NULL
    )
    if (is.null(fit_tr)) {
      cat("  FAILED:", cand$label, "\n"); next
    }
    
    res <- arma_fc(y_all, ntrain,
                   order    = cand$order,
                   seasonal = cand$seas,
                   method   = "CSS",
                   traincoef = fit_tr$coef)
    
    cat(sprintf("  %-45s  RMSE=%7.2f  MAE=%7.2f\n",
                cand$label, res$rmse, res$mae))
    
    results <- rbind(results, data.frame(
      Series = name,
      Model  = cand$label,
      RMSE   = round(res$rmse, 2),
      MAE    = round(res$mae,  2)
    ))
  }
  
  cat("\n=== Residual Diagnostics (auto.arima) ===\n")
  checkresiduals(fit_auto)
  lb <- Box.test(residuals(fit_auto), lag=30, type="Ljung-Box",
                 fitdf=length(coef(fit_auto)))
  cat(sprintf("Ljung-Box (lag 30): Q* = %.2f  p = %.4f\n",
              lb$statistic, lb$p.value))
  
  print(results)
  return(results)
}

res_total <- run_arima(train$total,       test$total,       "Total")
res_A     <- run_arima(train$influenza_a, test$influenza_a, "Influenza A")
res_B     <- run_arima(train$influenza_b, test$influenza_b, "Influenza B")

all_results <- rbind(res_total, res_A, res_B)
write.csv(all_results, "arima_results.csv", row.names = FALSE)
save(res_total, res_A, res_B, file = "arima_workspace.RData")
cat("\n✓ Saved: arima_results.csv, arima_workspace.RData\n")