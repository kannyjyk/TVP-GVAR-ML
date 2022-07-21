library(expm)
library(Rcpp)
library(reticulate)

sourceCpp('~/GitHub/GVAR/code/PaperCode/threshold_functions.cpp')
source_python("~/GitHub/GVAR/code/PaperCode/lstm_function.py")

# ---------------------------------- Functions -----------------------------------------
Panel_TVP <- function(dat = dat, dat_cb = dat_cb, ntot = 3000, check.time = TRUE, b_prior = NULL, datamat = NULL) {
  # ---------------------- Data Procession ----------------------
  if (is.null(datamat)) {
    p <- ncol(dat) - 2
    p_cb <- ncol(dat_cb) - 2
    all_time <- unique(dat$Time)
    n_time <- length(all_time)
    all_country <- unique(dat$ID)
    n_country <- length(all_country)
    beta_len <- p * n_country + p_cb
    
    datamat <- matrix(0, nrow = n_time, ncol = beta_len)
    for (i in 1:n_country) {
      datamat[, ((i - 1) * p + 1):(i * p)] <- as.matrix(dat[dat$ID == all_country[i], 3:(2 + p)])
    }
    datamat[, (p * n_country + 1):beta_len] <- as.matrix(dat_cb[, 3:(2 + p_cb)])
  }
  
  N <- beta_len
  T <- nrow(datamat) - 1
  K <- 2
  y <- datamat[-1, ]
  X <- datamat[1:T, ]
  Xraw <- cbind(apply(X, 1, mean, na.rm = T), 1)
  
  X_all <- array(1, dim = c(T, N, K))
  X_all[, , 1] <- X
  
  # ------------------ Parameter Initialization -----------------
  # error specification
  epsilon <- eta <- matrix(0, T, N)
  S <- log(matrix(0.1, T, N))
  
  a_sig <- 0.01
  b_sig <- 0.01
  
  d0 <- d1 <- 0.6
  com.var <- 10
  
  b_draw <- matrix(0,K,N)
  bt_draw <- array(0,dim=c(T,K,N))
  omega_draw <- matrix(1,K,N)
  
  # prior-scaling (and upper bound for TVP) 
  XXinv <- diag(solve(crossprod(Xraw)))
  sd.shock <- apply(Xraw, 2, sd)[[1]]
  XXinv[1] <- XXinv[1]*sd.shock
  
  sc.fac <- c(100, 0.01)
  names(sc.fac) <- c("cons", "tvp")
  sc.ub <- XXinv*sc.fac["cons"]
  sc.ub.tvp <- c(sc.ub, XXinv*sc.fac["tvp"])
  sc.ub.tvp[2*K] <- sc.ub.tvp[2*K]*1e-5
  sc.ub.rho <- sc.fac["tvp"]
  
  if (is.null(b_prior)) b_prior <- matrix(0,2*K,N)
  V0 <- diag(sc.ub.tvp)
  
  # -------------------- Model Training ------------------------
  # sampling
  if(check.time){
    pb <- txtProgressBar(min = 1, max = ntot, style = 3)
  }
  
  t.start <- Sys.time()
  for(irep in 1:ntot){
    y_ <- Sy <- y
    
    # --------------------------
    # sample mean coefficients
    
    for(i in 1:N){
      normalizer <- exp(-S[, i] / 2)
      YY <- (y_[, i] - X_all[, i, ] %*% b_draw[, i]) * normalizer
      XX <- X_all[,i,]%*%diag(omega_draw[,i]) * normalizer
      b_it <- KF_fast(t(as.matrix(YY)),XX,as.matrix(exp(S[,i])),t(matrix(1,K,T)),K,1,T,matrix(0,K,1),diag(K)/10^15)
      
      Xtilde <- X_all[, i, ] * t(b_it)
      XX <- cbind(X_all[, i, ], Xtilde) * normalizer
      YY <- y_[,i] * normalizer
      
      V0inv <- diag(1/diag(V0))
      V_post <- try(solve((crossprod(XX)) + V0inv),silent=T)
      if (is(V_post,"try-error")) V_post <- ginv((crossprod(XX) + V0inv))
      b_post <- V_post %*% (crossprod(XX, YY)+ V0inv%*%b_prior[,i])
      b_i <- try(b_post+t(chol(V_post))%*%rnorm(2*K,0,1),silent=T)
      if (is(b_i,"try-error")) b_i <- mvrnorm(1,b_post,V_post)
      
      b_draw[,i] <- b_i[1:K]
      sqrtomega <- b_i[(K+1):(2*K)]
      
      # non-identified element
      uu <- runif(1, 0, 1)
      if (uu > 0.5){
        sqrtomega <- -sqrtomega
        b_it <- -b_it
      }
      omega_draw[,i] <- sqrtomega
      for(tt in 1:T){
        eta[tt,i] <- y_[tt,i]-X_all[tt,i,]%*%(as.numeric(b_draw[,i])+(sqrtomega*b_it))[,tt]
      }
      bt_draw[,,i] <- t((as.numeric(b_draw[,i])+sqrtomega*b_it))
    }
    
    fit <- matrix(0,T,N)
    for(i in 1:N){
      for(tt in 1:T){
        fit[tt,i] <- X_all[tt,i,] %*% bt_draw[tt,,i]
      }
    }
    
    # --------------------------
    # sample variances
    for(i in 1:N){
      eps <- eta[,i]
      S[,i] <- log(1/rgamma(1,a_sig + T/2,b_sig + crossprod(eps)/2))
    }
    
    if(irep %% 2000 == 1){
      par(mfrow = c(1, 2))
      ts.plot(bt_draw[,2,],main= paste0("Draw #",irep, ": beta_0 (interc.)"))
      ts.plot(bt_draw[,1,],main= paste0("Draw #",irep, ": beta_1"))
    }
    
    if(check.time == TRUE){
      setTxtProgressBar(pb, irep)
    }
  }
  return(list(fit = fit, bt = bt_draw, y = y, beta_len = beta_len, omega = omega_draw, 
              dat = dat, dat_cb = dat_cb, datamat = datamat))
}

Pred_Panel_TVP <- function(model, method = "const", p = 1, h = 1, slide = 100) {
  dim_bt <- dim(model$bt)
  n_time <- dim_bt[1]
  b <- model$bt[n_time, , ]
  x_ori <- model$y[n_time, ]
  y_t_pred <- array(1, dim = c(dim_bt[3], h + 1, 2))
  y_t_pred[, 1, 1] <- x_ori
  
  if (method == "const") {
    
    for (i in 1:h) {
      y_t_pred[, i + 1, 1] <- rowSums(y_t_pred[, i, ] * t(b))
    }
    
  } else if (method == "VAR") {
    
    bt_trans <- matrix(nrow = n_time, ncol = dim_bt[2] * dim_bt[3])
    for (i in 1:n_time) {
      bt_trans[i, ] <- c(model$bt[i, 1, ], model$bt[i, 2, ])
    }
    var.b <- VAR(bt_trans, p = p, type = "const")
    pred_var <- predict(var.b, n.ahead = h, ci = 0.95)
    bt_pred <- sapply(pred_var$fcst, function(x) x[, 1])
    for (i in 1:h) {
      if (h == 1) {
        b <- rbind(bt_pred[1:dim_bt[3]], bt_pred[(dim_bt[3] + 1):(dim_bt[2] * dim_bt[3])])
      } else {
        b <- rbind(bt_pred[i, 1:dim_bt[3]], bt_pred[i, (dim_bt[3] + 1):(dim_bt[2] * dim_bt[3])])
      }
      y_t_pred[, i + 1, 1] <- rowSums(y_t_pred[, i, ] * t(b))
    }
    
  } else if (method == "tree" | method == "rf" | method == "lasso" | method == "ridge") {
    
    for (i in 1:h) {
      n_vars <- dim_bt[2] * dim_bt[3]
      if (i == 1) {
        bt_trans <- matrix(nrow = n_time, ncol = n_vars)
        for (j in 1:n_time) {
          bt_trans[j, ] <- c(model$bt[j, 1, ], model$bt[j, 2, ])
        }
      } else {
        bt_trans <- rbind(bt_trans, pred)
        n_time <- n_time + 1
      }
      
      pred <- vector(length = n_vars)
      for (j in 1:n_vars) {
        X <- bt_trans[(n_time - slide):(n_time - 1), ]
        Y <- bt_trans[(n_time - slide + 1):n_time, j]
        colnames(X) <- paste0("X", 1:n_vars)
        dat_j <- data.frame(cbind(X, Y))
        
        if (method == "tree") {
          fit_j <- rpart(Y ~ ., data = dat_j)
          pred[j] <- predict(fit_j, data.frame(t(bt_trans[n_time, ])))
        } else if (method == "rf") {
          fit_j <- ranger(Y ~ ., data = dat_j)
          pred[j] <- predict(fit_j, data.frame(t(bt_trans[n_time, ])))$predictions
        } else if (method == "lasso") {
          fit_cv_j <- cv.glmnet(X, Y, type.measure = "mse", alpha = 1)
          pred[j] <- predict(fit_cv_j, t(bt_trans[n_time, ]))
          # pred[j] <- predict(fit_cv_j, t(bt_trans[n_time, ]), lambda = fit_cv$lambda.min)
        } else if (method == "ridge") {
          fit_cv_j <- cv.glmnet(X, Y, type.measure = "mse", alpha = 0)
          pred[j] <- predict(fit_cv_j, t(bt_trans[n_time, ]))
          # pred[j] <- predict(fit_cv_j, t(bt_trans[n_time, ]), lambda = fit_cv$lambda.min)
        }
      }
      b <- rbind(pred[1:dim_bt[3]], pred[(dim_bt[3] + 1):n_vars])
      y_t_pred[, i + 1, 1] <- rowSums(y_t_pred[, i, ] * t(b))
    }
  } else if (method == 'lstm' | method == "gru" | method == "lstm_gru") {
    n_vars <- dim_bt[2] * dim_bt[3]
    bt_trans <- matrix(nrow = n_time, ncol = n_vars)
    for (j in 1:n_time) {
      bt_trans[j, ] <- c(model$bt[j, 1, ], model$bt[j, 2, ])
    }
    if (method == "lstm") {
      pred_data <- lstm_gvar(bt_trans, h_pred = h)
    } else if (method == "gru") {
      pred_data <- gru_gvar(bt_trans, h_pred = h)
    } else if (method == "lstm_gru") {
      pred_data <- lstm_gru_gvar(bt_trans, h_pred = h)
    }
    for (i in 1:h){
      pred <- pred_data[i, ]
      pred <- t(pred)
      b <- rbind(pred[1:dim_bt[3]], pred[(dim_bt[3] + 1):n_vars])
      y_t_pred[, i + 1, 1] <- rowSums(y_t_pred[, i, ] * t(b))
    }
  }
  
  
  return(y_t_pred[, , 1])
}

Pred_Panel_TVP1 <- function(dat_all, dat_cb, method = "const", p = 1, h_pred = 6, slide = 100, ntot = 1000) {
  
  for (h_ in h_pred:1) {
    time_all1 <- unique(dat_all$Time)
    len_time1 <- length(time_all1)
    time_ind_train <- time_all1[1:(len_time1 - h_)]
    time_ind_test <- time_all1[(len_time1 - h_ + 1):len_time1]
    dat_train <- dat_all[dat_all$Time %in% time_ind_train, ]
    dat_test <- dat_all[dat_all$Time %in% time_ind_test, ]
    dat_cb_train <- dat_cb[dat_cb$Time %in% time_ind_train, ]
    dat_cb_test <- dat_cb[dat_cb$Time %in% time_ind_test, ]
    
    panel_TVP_model <- Panel_TVP(dat = dat_train, dat_cb = dat_cb_train, ntot = ntot)
    pred_temp <- Pred_Panel_TVP(model = panel_TVP_model, method = method, p = p, h = 1)
    if (h_ == h_pred) {
      pred <- pred_temp
    } else {
      pred <- cbind(pred, pred_temp[, 2])
    }
  }
  return(list(pred = pred, beta_len = panel_TVP_model$beta_len))
}


Pred_Panel_TVP_Slide <- function(model, h = 10, ntot_pred = 500, check.time = TRUE) {
  
  cat("Predicting:\n")
  
  if(check.time){
    pb <- txtProgressBar(min = 1, max = h, style = 3)
  }
  
  for (h_temp in 1:h) {
    if (h_temp == 1) {
      pred1 <- Pred_Panel_TVP(model)
    } else {
      model <- Panel_TVP(dat = model$dat, dat_cb = model$dat_cb, datamat = datamat, 
                         b_prior = b_prior, ntot = ntot_pred, check.time = FALSE)
      pred1 <- Pred_Panel_TVP(model)
    }
    datamat <- rbind(model$datamat, pred1[, 2])
    b_prior <- rbind(model$bt[dim(model$bt)[1], , ], model$omega)
    
    if(check.time == TRUE){
      setTxtProgressBar(pb, h_temp)
    }
  }
  
  return(datamat)
}

CalMSE <- function(X, Y, num = NULL) {
  if (is.null(num)) {
    return(mean((X - Y) ^ 2))
  } else {
    return(mean((X[1:num] - Y[1:num]) ^ 2))
  }
}

GVAR0 <- function(y, y_cb, exogen) {
  
  y <- as.matrix(y)
  y_names <- colnames(y)
  p <- ncol(y)
  ylags <- y[-nrow(y), ]
  colnames(ylags) <- paste0(y_names, "_lag1")
  
  y_cb_names <- names(y_cb)[!is.na(names(y_cb))]
  y_cb <- as.matrix(y_cb)
  ylags_cb <- embed(y_cb, 2)
  colnames(ylags_cb) <- paste0(y_cb_names, "_lag", rep(c(0, 1), each = length(y_cb_names)))
  
  exogen <- as.matrix(exogen)
  datamat <- data.frame(cbind(ylags, exogen, ylags_cb)) # combine the global effect, local effect and center bank effect
  
  varresult <- list()
  for (i in 1:p) {
    y.temp <- y[-1, i]
    varresult[[y_names[i]]] <- lm(y.temp ~ ., data = datamat)
  }
  return(varresult)
}

GVAR0_CB <- function(y, exogen) {
  
  y_names <- names(y)[!is.na(names(y))]
  y <- as.matrix(y)
  ylags <- y[-nrow(y), ]
  
  exogen <- as.matrix(exogen)
  datamat_ <- data.frame(cbind(ylags, exogen))
  names(datamat_) <- c(y_names, colnames(exogen))
  
  varresult <- list()
  for (i in 1:ncol(y)) {
    if (ncol(y) == 1) {
      datamat <- datamat_
    } else {
      datamat <- datamat_[, grepl(paste0("CB", i), names(datamat_))]
    }
    
    # datamat <- datamat_
    y.temp <- y[-1, i]
    varresult[[y_names[i]]] <- lm(y.temp ~ ., data = datamat)
  }
  
  return(varresult)
}

GVAR_Ft <- function(data, weight = NULL, variates = variates, CB = FALSE,
                    Year = Year, NAME = NAME, N = N, timeID = timeID) {
  
  dat <- data[, variates]
  Ft <- list()
  for (i in 1:N) {
    exo <- NULL  # Compute Exogenous Foreign Variables
    for (j in 1:length(variates)) {
      dat_matrix <- matrix(as.numeric(dat[, j]), ncol = N) # The size of matrix is (timeID * N) 
      
      if (is.vector(weight)) {
        if (CB == FALSE) {
          F.tmp <- dat_matrix %*% as.matrix(weight[, i])
        } else {
          p_cb <- 1
          F.tmp <- dat_matrix %*% as.matrix(weight)
        }
        
      } else {
        if (CB == FALSE) {
          F.tmp <- dat_matrix %*% as.matrix(weight[, i])
        } else {
          p_cb <- nrow(weight)
          F.tmp <- dat_matrix %*% t(as.matrix(weight))
        }
      }
      
      exo <- cbind(exo, F.tmp) # Foreign variables
    }
    if (CB == FALSE) {
      colnames(exo) <- paste0(NAME[i], "_", variates)
    } else {
      colnames(exo) <- paste0(NAME[i], "_cb", 1:p_cb, "_", rep(variates, each = p_cb))
    }
    Ft[[i]] <- exo
    
    if (CB == TRUE) break
  }
  return(Ft)
}

GVAR_Est <- function(dat, dat_cb, weight_matrix, weight_cb) {
  
  variates <- colnames(dat[, -(1:2)]) # Names of column variables
  Year <- unique(as.character(year(dat$Time)))
  NAME <- as.character(unique(dat$ID)) # Names of countries
  N <- length(NAME) # Number of countries
  timeID <- as.character(dat$Time[dat$ID == NAME[1]])
  p_cb <- ncol(dat_cb) - 2
  
  results <- results_cb <- list()
  Ft <- GVAR_Ft(data = dat, weight = weight_matrix, CB = F, variates, Year, NAME, N, timeID)
  Ft_cb <- GVAR_Ft(data = dat, weight = weight_cb, CB = T, variates, Year, NAME, N, timeID)[[1]]
  
  ytmp_cb <- dat_cb[, 3:(2 + p_cb)]
  if (p_cb == 1) {
    names(ytmp_cb) <- names(dat_cb)[3]
  }
  exogen_cb <- embed(Ft_cb, 2)
  colnames(exogen_cb) <- paste0(paste0("F_CB", 1:p_cb, "_", variates), "_lag", 
                                rep(c(0, 1), each = length(variates) * p_cb))
  results_cb <- GVAR0_CB(y = ytmp_cb, exogen = exogen_cb)
  
  for (i in 1:N) {
    ytmp <- dat[dat$ID == NAME[i], -(1:2)]
    exogen <- embed(Ft[[i]], 2)
    colnames(ytmp) <- paste0(NAME[i], "_", variates)
    colnames(exogen) <- paste0(paste0(NAME[i], "_F_", variates), "_lag", 
                               rep(c(0, 1), each = length(variates)))
    results[[i]] <- GVAR0(y = ytmp, y_cb = ytmp_cb, exogen = exogen)
  }
  return(list(gvar_coef = results, cb_coef = results_cb, n_var = variates))
}

GVAR_IRF <- function(model_result, weight_matrix, weight_cb, n = 10, 
                     OIRF = FALSE) {
  
  res_u <- NULL
  res_name <- NULL
  at <- NULL
  G0 <- NULL
  G1 <- NULL
  N <- length(model_result$gvar_coef)
  n_var <- length(model_result$n_var)  # numbers of column variables
  
  ## Translate the weight of matrix
  weight_matrix <- as.matrix(weight_matrix)
  
  if (is.vector(weight_cb)) {
    p_cb <- 1
  } else {
    p_cb <- nrow(weight_cb)
  }
  weight_cb <- unlist(weight_cb)
  
  for (i in 1:N) {
    model_temp <- model_result$gvar_coef[[i]]
    
    ## Compute residuals about u
    for (j in 1:n_var) {
      res_u <- cbind(res_u, resid(model_temp[[j]]))
    }
    res_name <- c(res_name, names(model_temp))
    
    ## Compute the coefficients of functions
    exo_lag0_w <- NULL
    exo_lag1_w <- NULL
    exo_cb_lag0_w <- NULL
    exo_cb_lag1_w <- NULL
    
    a_i <- NULL
    phi_i1 <- NULL
    gamma_i0 <- NULL
    gamma_i1 <- NULL
    gamma_cb_i0 <- NULL
    gamma_cb_i1 <- NULL
    
    for (j in 1:n_var) {
      coef_all <- coef(model_temp[[j]])
      coef_const <- coef_all[1]
      coef_lag1 <- coef_all[2:(n_var + 1)] # the jth row of phi_i1
      coef_exo_lag0 <- coef_all[(n_var + 2):(2 * n_var + 1)] # the jth row of gamma_i0
      coef_exo_lag1 <- coef_all[(2 * n_var + 2):(3 * n_var + 1)] # the jth row of gamma_i1
      coef_exo_cb_lag0 <- coef_all[(3 * n_var + 2):(3 * n_var + 1 + p_cb)]
      coef_exo_cb_lag1 <- coef_all[(3 * n_var + 2 + p_cb):(3 * n_var + 1 + 2 * p_cb)]
      
      a_i <- c(a_i, coef_const)
      phi_i1 <- rbind(phi_i1, coef_lag1)
      gamma_i0 <- rbind(gamma_i0, coef_exo_lag0)
      gamma_i1 <- rbind(gamma_i1, coef_exo_lag1)
      gamma_cb_i0 <- rbind(gamma_cb_i0, coef_exo_cb_lag0)
      gamma_cb_i1 <- rbind(gamma_cb_i1, coef_exo_cb_lag1)
    }
    
    ## Compute G0
    A_i0 <- cbind(diag(n_var), - gamma_i0, - gamma_cb_i0)
    A_i1 <- cbind(phi_i1, gamma_i1, gamma_cb_i1)
    
    W_i <- matrix(0, nrow = n_var, ncol = n_var * N + p_cb)
    W_i[, (n_var * (i - 1) + 1):(n_var * i)] <- diag(n_var)
    W_i_exo <- NULL
    for (j in 1:N) {
      W_i_exo <- cbind(W_i_exo, weight_matrix[j, i] * diag(n_var))
    }
    W_i_exo <- cbind(W_i_exo, matrix(0, n_var, p_cb))
    W_i <- rbind(W_i, W_i_exo)
    W_i <- rbind(W_i, cbind(matrix(0, p_cb, n_var * N), diag(p_cb)))
    
    at <- c(at, a_i)
    G0 <- rbind(G0, A_i0 %*% W_i)
    G1 <- rbind(G1, A_i1 %*% W_i)
  }
  
  for (i in 1:p_cb) {
    model_temp_cb <- model_result$cb_coef[[i]]
    coef_all_cb <- coef(model_temp_cb)
    a_cb_i <- coef_all_cb[1]
    phi2_cb_i1 <- coef_all_cb[2]
    gamma2_cb_i0 <- coef_all_cb[3:(n_var + 2)]
    gamma2_cb_i1 <- coef_all_cb[(n_var + 3):(2 * n_var + 2)]
    
    A_i0_cb <- c(1, - gamma2_cb_i0)
    A_i1_cb <- c(phi2_cb_i1, gamma2_cb_i1)
    
    W_i_cb <- matrix(0, nrow = n_var + 1, ncol = n_var * N + p_cb)
    W_i_cb[1, n_var * N + i] <- 1
    for (j in 1:n_var) {
      W_i_cb[j + 1, (1:N) * n_var + j - n_var] <- weight_cb
    }
    
    at <- c(at, a_cb_i)
    G0 <- rbind(G0, A_i0_cb %*% W_i_cb)
    G1 <- rbind(G1, A_i1_cb %*% W_i_cb)
    
    res_u <- cbind(res_u, resid(model_temp_cb))
  }
  
  colnames(res_u) <- c(res_name, paste0("CB", 1:p_cb))
  m <- n_var * N + p_cb
  invGO <- solve(G0)
  F1 <- invGO %*% G1
  bt <- invGO %*% at
  eps_all <- res_u %*% t(invGO)
  
  ## Calculate the GIRF
  A0 <- diag(m)
  A1 <- F1 %*% A0
  Sigma_u <- cov(res_u)
  
  result <- list()
  if (OIRF == TRUE) {
    P <- t(chol(Sigma_u))
    for(i in 1:n) {
      result[[i]] <- (A1 %^% i) %*% invGO %*% P
    }
  } else {
    for(i in 1:n) {
      result_temp <- (A1 %^% i) %*% invGO %*% Sigma_u
      result[[i]] <- t(apply(result_temp, 1, function(x) x / sqrt(diag(Sigma_u))))
    }
  }
  
  return(list(IRF = result, F1 = F1, bt = bt, G0 = G0, G1 = G1, 
              res_u = res_u, eps_all = eps_all))
}

# Elimination Matrix
Eli_Mat <- function(m) {
  I <- diag(rep(1, m ^ 2))
  L <- NULL
  for (i in 1:m) {
    L <- rbind(L, I[(m * (i - 1) + i):(m * i), ])
  }
  return(L)
}

# Communication Matrix
Commu_Mat <- function(m, n) {
  I <- diag(rep(1, m * n))
  K <- NULL
  for (i in 1:m) {
    K <- rbind(K, I[seq(from = i, by = m, length.out = n), ])
  }
  return(K)
}

# Duplication Matrix
Dup_Mat <- function(m) {
  I <- diag(rep(1, m * (m + 1) / 2))
  D <- NULL
  for (i in 1:m) {
    if (i > 1) {
      for (j in 1:(i - 1)) {
        D <- rbind(D, I[(m + m - j + 2) * (j - 1) / 2 + i - j + 1, ])
      }
    }
    D <- rbind(D, I[
      ((m + 1) * m / 2 - (m - i + 1 + 1) * (m - i + 1) / 2 + 1):((m + 1) * m / 2 - (m - i + 1) * (m - i) / 2), ])
  }
  return(D)
}

"%^%" = function(a, n) {
  for (i in 1:n) {
    if (i == 1){
      b = a}
    else{
      b = a %*% b}}
  return = b
}

IRF_CI_asym <- function(
  dat, dat_cb,
  panel_TVP_model, result_irf,
  n = 10, alpha = 0.05
) {
  ## The Asymptotic Distributions of Impulse Response
  ## Ljtkepohl, H. New introduction to multiple time series analysis. Berlin:
  ## Springer, 2006. ISBN 978-3-540-26239-8.
  
  n_var <- ncol(dat) - 2
  p_cb <- ncol(dat_cb) - 2
  all_time <- unique(dat$Time)
  first_time <- all_time[1]
  all_ID <- dat[dat$Time == first_time, "ID"]
  p_all <- length(all_ID) * n_var + p_cb ## K
  
  F1_t <- result_irf$F1 ## A1 = A
  
  dim_bt <- dim(panel_TVP_model$bt)
  n_time <- dim_bt[1] ## T
  
  datamat <- panel_TVP_model$datamat
  Z <- rbind(rep(1, n_time), t(datamat[-(n_time + 1), ])) ## Z
  Y <- t(datamat[-1, ]) ## Y
  sigma_eps <- (Y %*% (diag(n_time) - t(Z) %*% solve(Z %*% t(Z)) %*% Z) %*% t(Y)) / (n_time - p_all - 1) ## Sigma_u
  P_eps <- t(chol(sigma_eps)) ## P
  
  H <- t(Eli_Mat(p_all)) %*% solve(
    Eli_Mat(p_all) %*% (diag(p_all ^ 2) + Commu_Mat(p_all, p_all))
    %*% kronecker(P_eps, diag(p_all)) %*% t(Eli_Mat(p_all))
  ) ## H
  
  D_plus_K <- solve(t(Dup_Mat(p_all)) %*% Dup_Mat(p_all)) %*% t(Dup_Mat(p_all))
  sigma_sig <- 2 * D_plus_K %*% kronecker(sigma_eps, sigma_eps) %*% t(D_plus_K) ## Sigma_sig
  
  gamma_inv <- solve(Z %*% t(Z) / n_time)
  gamma_inv <- gamma_inv[-1, -1] ## Gamma_Y_inv
  sigma_beta <- kronecker(gamma_inv, sigma_eps) ## Sigma_alpha
  
  irf_low_ci <- irf_up_ci <- list()
  for (i in 1:n) {
    G_i <- 0
    for (m in 0:(i - 1)) {
      G_i <- G_i + kronecker(t(F1_t) %^% (i - 1 - m), F1_t %^% m)
    } ## G_i
    
    C_i <- kronecker(t(P_eps), diag(p_all)) %*% G_i ## C_i
    
    C_bar_i <- kronecker(diag(p_all), F1_t %^% i) %*% H ## C_bar_i
    
    X <- C_i %*% sigma_beta %*% t(C_i) + C_bar_i %*% sigma_sig %*% t(C_bar_i)
    
    half <- matrix(nrow = p_all, ncol = p_all)
    for (I in 1:p_all) {
      for (J in 1:p_all) {
        half[I, J] <- qnorm(1 - alpha / 2) * 1 / sqrt(n_time) * sqrt(X[p_all * (J - 1) + I, p_all * (J - 1) + I])
      }
    }
    
    irf_low_ci[[i]] <- result_irf$IRF[[i]] - half
    irf_up_ci[[i]] <- result_irf$IRF[[i]] + half
  }
  
  return(list(irf_low_ci = irf_low_ci, irf_up_ci = irf_up_ci, X = X))
}

GVAR_pred <- function(model, dat, dat_cb, F1, bt, h = 10) {
  
  n_var <- ncol(dat) - 2
  p_cb <- ncol(dat_cb) - 2
  all_time <- unique(dat$Time)
  last_time <- all_time[length(all_time)]
  
  all_ID <- dat[dat$Time == last_time, "ID"]
  all_variables <- dat[dat$Time == last_time, 3:(n_var + 2)]
  all_cb <- dat_cb[dat_cb$Time == last_time, 3:(p_cb + 2)]
  
  xt_all <- c(as.vector(t(as.matrix(all_variables))), unlist(all_cb))
  
  xt_pred <- NULL
  for (i in 1:h) {
    xt_all <- F1 %*% xt_all + bt
    xt_pred <- cbind(xt_pred, xt_all)
  }
  
  pred <- list()
  for (i in 1:n_var) {
    pred[[i]] <- xt_pred[(1:length(all_ID)) * n_var - n_var + i, ]
    colnames(pred[[i]]) <- paste0("time", 1:h)
    rownames(pred[[i]]) <- all_ID
  }
  pred_cb <- xt_pred[(length(all_ID) * n_var + 1):length(xt_all), ]
  if (is.vector(pred_cb)) {
    pred_cb <- matrix(pred_cb, nrow = 1)
  }
  colnames(pred_cb) <- paste0("time", 1:h)
  rownames(pred_cb) <- names(all_cb)
  
  return(list(pred = pred, pred_cb = pred_cb))
}

### TVP-IRF

TV_IRF <- function(sel_time, panel_TVP_model, result_irf, OIRF = FALSE, n = 10) {
  
  dat_cb <- panel_TVP_model$dat_cb
  ind_time <- which(dat_cb$Time == sel_time) - 1
  N <- dim(panel_TVP_model$bt)[3]
  b0_t <- as.matrix(panel_TVP_model$bt[ind_time, 2, ])
  F1_t <- diag(panel_TVP_model$bt[ind_time, 1, ])
  F1_t_inv <- diag(1 / panel_TVP_model$bt[ind_time, 1, ])
  G0_t <- result_irf$G1 %*% F1_t_inv
  G0_t_inv <- F1_t %*% solve(result_irf$G1)
  
  res_eps <- panel_TVP_model$datamat[-1, ] - panel_TVP_model$fit
  res_u <- res_eps %*% t(G0_t)
  Sigma_u <- G0_t_inv %*% cov(res_eps) %*% t(G0_t_inv)
  # Sigma_u <- cov(result_irf$res_u)
  
  A0 <- diag(N)
  A1 <- F1_t %*% A0
  
  result <- list()
  if (OIRF == TRUE) {
    P <- t(chol(Sigma_u))
    for(i in 1:n) {
      result[[i]] <- (A1 %^% i) %*% G0_t_inv %*% P
    }
  } else {
    for(i in 1:n) {
      result_temp <- (A1 %^% i) %*% G0_t_inv %*% Sigma_u
      result[[i]] <- t(apply(result_temp, 1, function(x) x / sqrt(diag(Sigma_u))))
    }
  }
  return(list(IRF = result, G0_t = G0_t, F1_t = F1_t, b0_t = b0_t, res_u = res_u,
              sel_time = sel_time))
}

IRF_TV_CI_asym <- function(
  sel_time,
  dat, dat_cb,
  panel_TVP_model, result_irf, result_irf_tv,
  n = 10, alpha = 0.05
) {
  ## The Asymptotic Distributions of Impulse Response
  ## Ljtkepohl, H. New introduction to multiple time series analysis. Berlin:
  ## Springer, 2006. ISBN 978-3-540-26239-8.
  
  n_var <- ncol(dat) - 2
  p_cb <- ncol(dat_cb) - 2
  all_time <- unique(dat$Time)
  first_time <- all_time[1]
  all_ID <- dat[dat$Time == first_time, "ID"]
  p_all <- length(all_ID) * n_var + p_cb ## K
  
  ind_time <- which(dat_cb$Time == sel_time) - 1
  F1_t <- diag(panel_TVP_model$bt[ind_time, 1, ]) ## A1 = A
  
  F1_t_inv <- diag(1 / panel_TVP_model$bt[ind_time, 1, ])
  G0_t <- result_irf$G1 %*% F1_t_inv
  G0_t_inv <- F1_t %*% solve(result_irf$G1)
  
  dim_bt <- dim(panel_TVP_model$bt)
  n_time <- dim_bt[1] ## T
  
  datamat <- panel_TVP_model$datamat
  Z <- rbind(rep(1, n_time), t(datamat[-(n_time + 1), ])) ## Z
  Y <- t(datamat[-1, ]) ## Y
  sigma_eps <- (Y %*% (diag(n_time) - t(Z) %*% solve(Z %*% t(Z)) %*% Z) %*% t(Y)) / (n_time - p_all - 1) ## Sigma_u
  P_eps <- t(chol(sigma_eps)) ## P
  
  H <- t(Eli_Mat(p_all)) %*% solve(
    Eli_Mat(p_all) %*% (diag(p_all ^ 2) + Commu_Mat(p_all, p_all))
    %*% kronecker(P_eps, diag(p_all)) %*% t(Eli_Mat(p_all))
  ) ## H

  D_plus_K <- solve(t(Dup_Mat(p_all)) %*% Dup_Mat(p_all)) %*% t(Dup_Mat(p_all))
  sigma_sig <- 2 * D_plus_K %*% kronecker(sigma_eps, sigma_eps) %*% t(D_plus_K) ## Sigma_sig
  
  gamma_inv <- solve(Z %*% t(Z) / n_time)
  gamma_inv <- gamma_inv[-1, -1] ## Gamma_Y_inv
  sigma_beta <- kronecker(gamma_inv, sigma_eps) ## Sigma_alpha
  
  irf_low_ci <- irf_up_ci <- list()
  for (i in 1:n) {
    G_i <- 0
    for (m in 0:(i - 1)) {
      G_i <- G_i + kronecker(t(F1_t) %^% (i - 1 - m), F1_t %^% m)
    } ## G_i
    
    C_i <- kronecker(t(P_eps), diag(p_all)) %*% G_i ## C_i
    
    C_bar_i <- kronecker(diag(p_all), F1_t %^% i) %*% H ## C_bar_i
    
    X <- C_i %*% sigma_beta %*% t(C_i) + C_bar_i %*% sigma_sig %*% t(C_bar_i)
    
    half <- matrix(nrow = p_all, ncol = p_all)
    for (I in 1:p_all) {
      for (J in 1:p_all) {
        half[I, J] <- qnorm(1 - alpha / 2) * 1 / sqrt(n_time) * sqrt(X[p_all * (J - 1) + I, p_all * (J - 1) + I])
      }
    }
    
    irf_low_ci[[i]] <- result_irf_tv$IRF[[i]] - half
    irf_up_ci[[i]] <- result_irf_tv$IRF[[i]] + half
  }
  
  return(list(irf_low_ci = irf_low_ci, irf_up_ci = irf_up_ci, X = X))
}
