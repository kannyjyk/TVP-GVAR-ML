library(data.table)
library(dplyr)
library(zoo)
library(mfbvar)

setwd("~/GitHub/GVAR/data/PaperData")

ReadDataMain <- function() {
  areas <- c("JPN", "USA", "EU27_2020")
  
  dat_CPI <- fread("DP_CPI.csv", data.table = FALSE)
  dat_HUR <- fread("DP_Harmonised unemployment rate.csv", data.table = FALSE)
  dat_GDP <- fread("DP_GDP.csv", data.table = FALSE)
  
  dat_all <- data.frame()
  for (i in 1:length(areas)) {
    dat_CPI_sel <- dat_CPI %>% 
      filter(LOCATION == areas[i] & (SUBJECT == "TOT") & (FREQUENCY == "M") & (MEASURE == "AGRWTH") & (TIME >= "2000-01") & (TIME <= "2021-06"))
    dat_CPI_ts <- ts(dat_CPI_sel$Value, freq = 12, start = c(2000, 1))
    
    dat_HUR_sel <- dat_HUR %>%
      filter(LOCATION == areas[i] & (SUBJECT == "TOT") & (FREQUENCY == "M") & (TIME >= "2000-01") & (TIME <= "2021-06"))
    dat_HUR_sel$Value <- c(0, diff(dat_HUR_sel$Value) / dat_HUR_sel$Value[- nrow(dat_HUR_sel)])
    
    dat_HUR_ts <- ts(dat_HUR_sel$Value, freq = 12, start = c(2000, 1))
    
    dat_GDP_sel <- dat_GDP %>% 
      filter(LOCATION == areas[i] & (SUBJECT == "TOT") & (FREQUENCY == "Q") & (MEASURE == "PC_CHGPY") & (TIME >= "2000-Q1") & (TIME <= "2021-Q2"))
    dat_GDP_ts <- ts(dat_GDP_sel$Value, freq = 4, start = c(2000, 1))
    
    ts_all <- list(CPI = dat_CPI_ts / 100, HUR = dat_HUR_ts, GDP = dat_GDP_ts / 100)
    
    prior_obj <- set_prior(Y = ts_all, n_lags = 3, n_reps = 1000)
    mod_minn <- estimate_mfbvar(prior_obj, prior = "minn")
    
    dat_temp <- data.frame(ID = areas[i], Time = rownames(mod_minn$Z[, , 1000]), mod_minn$Z[, , 1000])
    dat_all <- rbind(dat_all, dat_temp)
    
  }
  return(dat_all)
}

# ----- Crude Oil Prices
ReadCB <- function() {
  dat_oil <- fread("Crude Oil Prices West Texas Intermediate (WTI) - Cushing, Oklahoma.csv", data.table = FALSE)
  oil_price <- dat_oil$MCOILWTICO
  oil_price_rate <- diff(oil_price) / oil_price[-length(oil_price)]
  dat_oil <- data.frame(ID = "CB", Time = dat_oil$DATE[-1], CB = oil_price_rate)
  dat_cb <- dat_oil %>% 
    filter((Time >= "2000-01-01") & (Time <= "2021-06-01"))
  dat_cb$Time <- as.character(dat_cb$Time)
  
  return(dat_cb)
}

# ---- Weight Matrix 
ReadWeights <- function(rda_file = "DOT_2020Q3.rda") {
  load(rda_file)
  
  # China, P.R.: Mainland: 924 --- United States: 111 --- European Union: 998 
  # Japan: 158 --- United States: 111 --- European Union: 998 
  weight_names <- c("JPN", "US", "EU")
  
  us_china <- dat_result[dat_result$`Country Code` == 111 & dat_result$`Counterpart Country Code` == 158, ]
  eu_china <- dat_result[dat_result$`Country Code` == 998 & dat_result$`Counterpart Country Code` == 158, ]
  eu_us <- dat_result[dat_result$`Country Code` == 998 & dat_result$`Counterpart Country Code` == 111, ]
  
  weight_mat <- matrix(0, nrow = 3, ncol = 3)
  weight_mat[2, 1] <- as.numeric(us_china$`2020Q3`[us_china$`Indicator Code` == "TXG_FOB_USD" & us_china$`2020Q3` != "" & us_china$`2020Q3` != "r"])
  weight_mat[1, 2] <- as.numeric(us_china$`2020Q3`[us_china$`Indicator Code` == "TMG_CIF_USD" & us_china$`2020Q3` != "" & us_china$`2020Q3` != "r"])
  weight_mat[3, 1] <- as.numeric(eu_china$`2020Q3`[eu_china$`Indicator Code` == "TXG_FOB_USD" & eu_china$`2020Q3` != ""])
  weight_mat[1, 3] <- as.numeric(eu_china$`2020Q3`[eu_china$`Indicator Code` == "TMG_CIF_USD" & eu_china$`2020Q3` != ""])
  weight_mat[3, 2] <- as.numeric(eu_us$`2020Q3`[eu_us$`Indicator Code` == "TXG_FOB_USD" & eu_us$`2020Q3` != ""])
  weight_mat[2, 3] <- as.numeric(eu_us$`2020Q3`[eu_us$`Indicator Code` == "TMG_CIF_USD" & eu_us$`2020Q3` != ""])
  
  cb_weight <- rowSums(weight_mat)
  cb_weight <- cb_weight / sum(cb_weight)
  weight_mat <- t(apply(weight_mat, 1, function(x) x / sum(x)))
  
  return(list(weight_mat = weight_mat, cb_weight = cb_weight))
}

# ------------------------- Load Data ------------------------------
dat_all <- ReadDataMain()
dat_cb <- ReadCB()
weights_list <- ReadWeights("DOT_2020Q3.rda")
weight_mat <- weights_list$weight_mat
cb_weight <- weights_list$cb_weight

