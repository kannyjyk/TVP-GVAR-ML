rm(list = ls())
library(vars)
library(rpart)
library(ranger)
library(glmnet)
library(Rcpp)
library(reticulate)
library(lubridate)

source('~/GitHub/GVAR/code/PaperCode/LoadData.R')
source('~/GitHub/GVAR/code/PaperCode/Functions.R')

file_date <- "20220720"

# ------------------ Preprocess Data --------------------

h_pred <- 6
time_all <- unique(dat_all$Time)
len_time <- length(time_all)
time_ind_train <- time_all[1:(len_time - h_pred)]
time_ind_test <- time_all[(len_time - h_pred + 1):len_time]
dat_train <- dat_all[dat_all$Time %in% time_ind_train, ]
dat_test <- dat_all[dat_all$Time %in% time_ind_test, ]
dat_cb_train <- dat_cb[dat_cb$Time %in% time_ind_train, ]
dat_cb_test <- dat_cb[dat_cb$Time %in% time_ind_test, ]

# ------------------- Prediction ---------------------
set.seed(1)
panel_TVP_model <- Panel_TVP(dat = dat_train, dat_cb = dat_cb_train, ntot = 3000)

save(panel_TVP_model, file = paste0("~/GitHub/GVAR/result/", file_date, "/panel_TVP_model.rda"))
# --------------------- Plot -----------------------
## Plot the predicted results
var_names <- c("CPI", "HUR", "GDP")
country_ID <- unique(dat_all$ID)
all_names <- c(paste0(var_names, "_", rep(country_ID, each = length(var_names))), "ECB")
beta_len <- panel_TVP_model$beta_len

mse_all <- c()
png(paste0("~/GitHub/GVAR/result/", file_date, "/panel_tvp_insample", h_pred, ".png"), width = 3000, height = 3000, res = 300)
par(mfrow = c(4, 3))
for (i in 1:beta_len) {
  mse_temp <- CalMSE(panel_TVP_model$y[, i], panel_TVP_model$fit[, i])
  plot(as.Date(unique(dat_train$Time)[-1]), panel_TVP_model$y[, i], type = "l", xlab = "Year", ylab = all_names[i], 
       main = paste0("MSE: ", sprintf("%0.4f", mse_temp)))
  lines(as.Date(unique(dat_train$Time)[-1]), panel_TVP_model$fit[, i], col = "red")
  mse_all <- c(mse_all, mse_temp)
}
dev.off()

# --------------- Show the beta
(bt_dim <- dim(panel_TVP_model$bt))
png(paste0("~/GitHub/GVAR/result/", file_date, "/panel_tvp_beta0.png"), width = 3000, height = 3000, res = 300)
par(mfrow = c(4, 3))
for (i in 1:beta_len) {
  plot(as.Date(unique(dat_train$Time)[-1]), panel_TVP_model$bt[, 2, i], type = "l", xlab = "Year", ylab = all_names[i])
}
dev.off()

png(paste0("~/GitHub/GVAR/result/", file_date, "/panel_tvp_beta1.png"), width = 3000, height = 3000, res = 300)
par(mfrow = c(4, 3))
for (i in 1:beta_len) {
  plot(as.Date(unique(dat_train$Time)[-1]), panel_TVP_model$bt[, 1, i], type = "l", xlab = "Year", ylab = all_names[i])
}
dev.off()


# ----------------------------- Prediction ----------------------------------
set.seed(1)
MSE_Methods <- NULL
methods <- c("const", "VAR", "tree", "rf", "lasso", "lstm", "gru", "lstm_gru")
for (m in methods) {
  pred <- Pred_Panel_TVP(model = panel_TVP_model, method = m, p = 1, h = h_pred)
  save(pred, file = paste0("~/GitHub/GVAR/result/", file_date, "/Pred_Panel_TVP6direct_", m, ".rda"))
  mse_all <- c()
  png(paste0("~/GitHub/GVAR/result/", file_date, "/panel_tvp_pred_all", h_pred, "_", m, ".png"), width = 3000, height = 3000, res = 300)
  par(mfrow = c(4, 3))
  for (i in 1:beta_len) {
    if (i < beta_len) {
      plot_dat <- dat_all[dat_all$ID == country_ID[((i - 1) %/% (length(var_names))) + 1], ((i - 1) %% length(var_names) + 1) + 2]
    } else {
      plot_dat <- dat_cb$CB
    }
    mse_temp <- CalMSE(plot_dat[(length(plot_dat) - h_pred + 1):length(plot_dat)], pred[i, ][-1])
    plot(as.Date(time_all), plot_dat, type = "l", xlab = "Year", ylab = all_names[i], 
         main = paste0("MSE: ", sprintf("%0.4f", mse_temp)))
    lines(as.Date(time_all[(len_time - h_pred):len_time]), pred[i, ], col = "red", ylab = all_names[i])
    mse_all <- c(mse_all, mse_temp)
  }
  dev.off()
  
  MSE_Methods <- rbind(MSE_Methods, c(m, sprintf("%0.4f", mean(mse_all))))
}

colnames(MSE_Methods) <- c("Methods", "MSE")
write.csv(MSE_Methods, file = paste0("~/GitHub/GVAR/result/", file_date, "/MSE_All_Methods_pred6.csv"),
          row.names = F)

# ----------------------------- Prediction 1 ----------------------------------
set.seed(1)
h_pred <- 6
time_all <- unique(dat_all$Time)
len_time <- length(time_all)

## Plot the predicted results
var_names <- c("CPI", "HUR", "GDP")
country_ID <- unique(dat_all$ID)
all_names <- c(paste0(var_names, "_", rep(country_ID, each = length(var_names))), "ECB")

MSE_Methods1 <- NULL
methods <- c("const", "VAR", "tree", "rf", "lasso", "lstm", "gru", "lstm_gru")
for (m in methods) {
  pred_model <- Pred_Panel_TVP1(dat_all, dat_cb, method = m, p = 1, h_pred = h_pred, slide = 100, ntot = 1000)
  save(pred_model, file = paste0("~/GitHub/GVAR/result/", file_date, "/Pred_Panel_TVP1_", m, ".rda"))
  pred <- pred_model$pred
  beta_len <- pred_model$beta_len
  
  mse_all1 <- c()
  png(paste0("~/GitHub/GVAR/result/", file_date, "/pred1_outsample", h_pred, "_", m, "_v2.png"), width = 3000, height = 3000, res = 300)
  par(mfrow = c(4, 3))
  for (i in 1:beta_len) {
    if (i < beta_len) {
      plot_dat <- dat_all[dat_all$ID == country_ID[((i - 1) %/% (length(var_names))) + 1], ((i - 1) %% length(var_names) + 1) + 2]
    } else {
      plot_dat <- dat_cb$CB
    }
    mse_temp <- CalMSE(plot_dat[(length(plot_dat) - h_pred + 1):length(plot_dat)], pred[i, ][-1])
    plot(as.Date(time_all), plot_dat, type = "l", xlab = "Year", ylab = all_names[i], 
         main = paste0("MSE: ", sprintf("%0.4f", mse_temp)))
    lines(as.Date(time_all[(len_time - h_pred):len_time]), pred[i, ], col = "red", ylab = all_names[i])
    mse_all1 <- c(mse_all1, mse_temp)
  }
  dev.off()
  
  MSE_Methods1 <- rbind(MSE_Methods1, c(m, sprintf("%0.4f", mean(mse_all1))))
}

colnames(MSE_Methods1) <- c("Methods", "MSE")
write.csv(MSE_Methods1, file = paste0("~/GitHub/GVAR/result/", file_date, "/MSE_All_Methods_pred1.csv"), 
          row.names = F)


# ----------------------------- Impulse Response Function ----------------------------------
# ------------------ Time-invariant 
model_result <- GVAR_Est(dat = dat_train, dat_cb = dat_cb_train,
                         weight_matrix = weight_mat, weight_cb = cb_weight)
result_irf <- GVAR_IRF(model_result, weight_matrix = weight_mat, weight_cb = cb_weight, n = 6)

irf_ci <- IRF_CI_asym(dat = dat_train, dat_cb = dat_cb_train,
                    panel_TVP_model, result_irf,
                    n = 6, alpha = 0.05)

irf_low_ci <- irf_ci$irf_low_ci
irf_up_ci <- irf_ci$irf_up_ci
gvar_irf <- result_irf$IRF

save(gvar_irf, irf_low_ci, irf_up_ci, 
     file = paste0("~/GitHub/GVAR/result/", file_date, "/GIRF.rda"))

# Examples of results
coef(model_result$gvar_coef[[1]]$JPN_CPI)
gvar_irf[[1]]

# A simple plot for GIRF
var_names <- paste0(rep(c("JPN", "USA", "EU27"), each = 3), "-", c("CPI", "HUR", "GDP"))
var_names <- c(var_names, "Oil")

png(paste0("~/GitHub/GVAR/result/", file_date, "/GIRF_time_invariant.png"), width = 10000, height = 10000, res = 300)
par(mfrow = c(10, 10))
for (j in 1:length(var_names)) {
  for (l in 1:length(var_names)) {
    result_jl <- sapply(gvar_irf, function(x) x[j, l])
    result_low_jl <- sapply(irf_low_ci, function(x) x[j, l])
    result_up_jl <- sapply(irf_up_ci, function(x) x[j, l])
    plot(result_jl, type = "l", ylim = c(min(result_low_jl, 0), max(result_up_jl, 0)), 
         xlab = "Time", ylab = "GIRF", 
         main = paste0(var_names[j], " ----- ", var_names[l]))
    abline(h = 0, col = 2)
    lines(result_up_jl, lty = 2, col = 3)
    lines(result_low_jl, lty = 2, col = 3)
  }
}
dev.off()

# ------------------ Time-variant 

sel_time <- "2020-07-01" # replace it with "2011-04-01" and "2007-12-01
result_irf_tv <- TV_IRF(sel_time, panel_TVP_model, result_irf, OIRF = T, n = 6)

irf_ci_tv <- IRF_TV_CI_asym(
  sel_time = sel_time,
  dat = dat_train, dat_cb = dat_cb_train,
  panel_TVP_model, result_irf, result_irf_tv,
  n = 6, alpha = 0.05
)

irf_low_ci <- irf_ci_tv$irf_low_ci
irf_up_ci <- irf_ci_tv$irf_up_ci
gvar_irf <- result_irf_tv$IRF

# Plot
png(paste0("~/GitHub/GVAR/result/", file_date, "/OIRF_time_variant_", sel_time, ".png"), width = 10000, height = 10000, res = 300)
par(mfrow = c(10, 10))
for (j in 1:length(var_names)) {
  for (l in 1:length(var_names)) {
    result_jl <- sapply(gvar_irf, function(x) x[j, l])
    result_low_jl <- sapply(irf_low_ci, function(x) x[j, l])
    result_up_jl <- sapply(irf_up_ci, function(x) x[j, l])
    plot(result_jl, type = "l", xlab = "Time", ylab = "TV-GIRF", 
         ylim = c(min(0, result_low_jl, result_jl),
                  max(0, result_up_jl, result_jl)),
         main = paste0(var_names[j], " ----- ", var_names[l]))
    abline(h = 0, col = 2)
    lines(result_up_jl, lty = 2, col = 3)
    lines(result_low_jl, lty = 2, col = 3)
    }
  }
}
dev.off()
