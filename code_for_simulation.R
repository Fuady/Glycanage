#################################################################################################################### 
##                                             Code for simulation                                                ##
##       Objectives: (1) Assess the robustness of the regression calibration for the two mapping procedures       ##
##                                      when the assumptions are violated                                         ##
##                   (2) compare the performance of the glycan structure mapping approach and O2PLS mapping       ##
##                                                30 -06- 2021                                                    ##
#################################################################################################################### 

# Clean worksheet
# --------------------------------------------------------------------------------------------
rm(list = ls())


# Libraries --------------------------------------------------------------------------------------------

library(readr)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gplots)
library(stringr)
library(preprocessCore)
library(OmicsPLS)
library(msm)
library(systemfit)
library(xtable)
library(mvtnorm)
library(mice)
library(PO2PLS)
library(reshape)

# Function --------------------------------------------------------------------------------------------

logistic.regression <- function(yy, ww) {
  # Function to fit logistic regression 
  # Input: 
  #   yy: outcome 
  #   ww: predictor 
  # Output: 
  #   results: list of regression coefficient, covariance matrix, and data frame of regression coefficient
  
  my.data.main <- data.frame(y = factor(yy), ww)
  my.glm <- glm(y ~ ., data = my.data.main, family = binomial, control = glm.control(epsilon = 1e-08, maxit = 30))
  
  sumres <- summary(my.glm)
  estalp <- my.glm$coefficients
  covalp <- sumres$dispersion * sumres$cov.unscaled
  oddalp <- exp(estalp)
  sealp <- sqrt(diag(covalp))
  cloddalp <- cbind(exp(estalp - qnorm(1 - 0.025) * sealp), exp(estalp + qnorm(1 - 0.025) * sealp))
  pvalalp <- sumres$coefficient[, 4]  # pvalue of estimates
  
  resesta <- cbind(estalp, sealp, oddalp, cloddalp, pvalalp)
  dimnames(resesta) <- list(names(estalp), c("Estimate", "S.E", "Odds", "95% LL", "95%UL", "pval"))
  
  results <- list(alpha = estalp, covalpha = covalp, df.alpha = resesta)
  
  return(results)
}


linear.regression <- function(xx, ww) {
  # Function to fit linear regression 
  # Input: 
  #   xx: outcome 
  #   ww: predictor 
  # Output: 
  #   results: list of regression coefficient, covariance matrix, and data frame of regression coefficient
  
  n <- length(xx)
  results <- list()
  
  for (i in 1:n) {
    my.data.main <- data.frame(y = xx[[i]], ww[[i]])
    my.lm <- lm(y ~ ., data = my.data.main, na.action = na.omit)
    
    corx <- cor(xx[[i]], my.lm$fitted.values)
    varxwz <- summary(my.lm)$sigma^2
    n2 <- my.lm$df.residual + my.lm$rank
    
    sumres <- summary(my.lm)
    estgam <- my.lm$coefficients  # estimates
    covgam <- vcov(my.lm)
    segam <- sqrt(diag(covgam))
    pvalgam <- sumres$coefficient[, 4]  # pvalue of estimates
    
    resesta <- cbind(estgam, segam, pvalgam)
    dimnames(resesta) <- list(names(estgam), c("Estimate", "S.E", "pval"))
    
    results1 <- list(gamma = estgam, covgamma = covgam, df.gamma = resesta, n2 = n2, varxwz = varxwz)
    
    names(results1) <- paste0(names(results1), "_", i, sep = "")
    results <- c(results, results1)
  }
  
  return(results)
}


GP_prop <- function(logmodel, memodel, xx) {
  # Function to compute corrected logistic regression coefficient (beta_x) 
  # Input: 
  #   logmodel: results of the naive logistic regression 
  #   memodel: results of the measurement error model (linear regression) 
  # Output:
  # results: corrected logistic regression coefficient, estimate of naive logistic regression and its
  # covariance matrix, estimate of the measurement error model and its covariance matrix
  
  n <- length(xx)
  results <- list()
  for (i in 1:n) {
    selected <- grep(paste("_", i, sep = ""), names(memodel))
    select_par <- memodel[selected]
    names(select_par) <- c("gamma", "covgamma", "df.gamma", "n2", "varxwz")
    
    idx <- which(names(logmodel$alpha) %in% names(select_par$gamma))[-1]
    r <- length(idx)
    estaw <- logmodel$alpha[idx]
    estgw <- select_par$gamma[-1]
    covaw <- logmodel$covalpha[idx, idx]
    covgw <- select_par$covgamma[-1, -1]
    bdiag_ <- as.matrix(bdiag(covaw, covgw))
    betaw <- estaw/estgw
    
    if (length(estgw) == 1) {
      derbetaw <- rbind(1/estgw, -estaw/estgw^2)
    } else {
      derbetaw <- rbind(diag(1/estgw), diag(-estaw/estgw^2))
    }
    
    varxwz <- select_par$varxwz
    gamma0 <- select_par$gamma[1]
    
    res <- list(idx = idx, r = r, estaw = estaw, estgw = estgw, covaw = covaw, covgw = covgw, b_diag = bdiag_, 
                betaww = betaw, derbetaw = derbetaw, varxwz = varxwz, gamma0 = gamma0)
    names(res) <- paste0(names(res), "_", i, sep = "")
    results <- c(results, res)
  }
  return(results)
}


GP_mat <- function(xx, memodel, betaopt, varxwz_all, optwts, all_par) {
  # Function to compute corrected logistic regression coefficient 
  # Input: 
  #   logmodel: results of the naive logistic regression
  #   memodel: results of the measurement error model (linear regression)
  # Output: 
  #   results:
  #     corrected logistic regression coefficient, estimate of naive logistic regression and its covariance
  #     matrix, estimate of the measurement error model and its covariance matrix
  
  n <- length(xx)
  results <- list()
  estaw <- Reduce(c, all_par[grep("estaw_", names(all_par))])
  estgw <- Reduce(c, all_par[grep("estgw_", names(all_par))])
  
  for (t in 1:n) {
    selected <- grep(paste("_", t, sep = ""), names(memodel))
    select_par <- memodel[selected]
    names(select_par) <- c("gamma", "covgamma", "df.gamma", "n2", "varxwz")
    
    A0 <- optwts[t, ]
    optwts1 <- A0[A0 != 0]
    C <- -betaopt[t]
    D <- as.vector(as.vector(select_par$gamma[1] + betaopt[t] * varxwz_all[t, t]) * optwts1 * estaw[t]/estgw[t]^2)
    E <- as.vector(-optwts1 * estaw[t]/estgw[t]^2)
    
    res <- list(C = C, D = D, E = c(0, E), CandD = c(C, D))
    names(res) <- paste0(names(res), "_", t, sep = "")
    results <- c(results, res)
  }
  return(results)
}



adjfun <- function(mainy, mainw, validx, validw, wt1 = 1, wt2 = 1, digits = 3) {
  # Function to compute regression calibration in the main study 
  # Input: 
  #   mainy: the outcome in the main study
  #   mainw: the surrogate variable in the main study 
  #   validx: the true covariate of interest in the calibration study
  #   validw: the surrogate variable in the calibration study 
  # Output: 
  #   results: list of results for the
  #         measurement error in the calibration study, the naive logistic regression in the main study, the
  #         corrected logistic regression in the main study, inverse variance weighting, estimate of
  #         measurement error variance.
  
  # naive logistic model in main study between Y and W
  logmodel <- logistic.regression(yy = mainy, ww = do.call("cbind", mainw))  # ,zz=mainz,includeZ = F
  
  # measurement error model in validation study between X and W
  memodel <- linear.regression(xx = validx, ww = validw)  # , zz=validz, includeZ = F
  
  all_par <- GP_prop(logmodel = logmodel, memodel = memodel, xx = validx)
  
  if (length(validw) == 1 & length(validx) == 1) {
    covag <- do.call(rbind, all_par[grep("b_diag_", names(all_par))])
    derbetaw <- do.call(rbind, all_par[grep("derbetaw_", names(all_par))])
  } else {
    covag <- as.matrix(do.call(bdiag, all_par[grep("b_diag_", names(all_par))]))
    derbetaw <- as.matrix(do.call(bdiag, all_par[grep("derbetaw_", names(all_par))]))
  }
  
  varbetaw <- t(derbetaw) %*% covag %*% derbetaw
  stdbetaw <- sqrt(diag(as.matrix(varbetaw)))
  invvar <- solve(varbetaw, tol = 1e-25)
  
  if (sum(unlist(lapply(validw, ncol))) == 1 & length(validx) == 1) {
    onevec <- 1
    varxwz_all <- do.call(rbind, all_par[grep("varxwz_", names(all_par))])
    gamma1_all <- t(do.call(rbind, all_par[grep("estgw_", names(all_par))]))
    covall_tes <- as.matrix(bdiag(logmodel$covalpha, as.matrix(do.call(rbind, memodel[grep("covgamma", 
                                                                                           names(memodel))])), varxwz_all))
  } else if (length(validw) == 1) {
    onevec <- rep(1, ncol(do.call(rbind, validw)))
    varxwz_all <- do.call(rbind, all_par[grep("varxwz_", names(all_par))])
    gamma1_all <- (do.call(rbind, all_par[grep("estgw_", names(all_par))]))
    covall_tes <- as.matrix(bdiag(logmodel$covalpha, as.matrix(do.call(rbind, memodel[grep("covgamma", 
                                                                                           names(memodel))])), varxwz_all))
  } else {
    onevec <- as.matrix(do.call(bdiag, all_par[grep("betaww_", names(all_par))]))
    onevec[onevec != 0] <- 1
    varxwz_all <- as.matrix(do.call(bdiag, all_par[grep("varxwz_", names(all_par))]))
    gamma1_all <- t(as.matrix(do.call(bdiag, all_par[grep("estgw_", names(all_par))])))
    covall_tes <- as.matrix(bdiag(logmodel$covalpha, as.matrix(do.call(bdiag, memodel[grep("covgamma", 
                                                                                           names(memodel))])), varxwz_all))
  }
  
  optwts <- as.vector(1/(sum(invvar))) * (t(onevec) %*% invvar)
  betaw_all <- Reduce(c, all_par[grep("betaww_", names(all_par))])
  betaopt <- as.matrix(optwts %*% betaw_all)
  
  gamma0_all <- Reduce(c, all_par[grep("gamma0_", names(all_par))])
  beta_0_red <- t(betaopt) %*% varxwz_all %*% betaopt/2
  int_adj <- logmodel$alpha[1] - as.numeric(t(betaopt) %*% gamma0_all) - beta_0_red
  betaadj <- c(int_adj, betaopt)
  
  A1 <- t(as.matrix(optwts/gamma1_all))
  A1[is.na(A1)] <- 0
  A1_ready <- -A1 %*% (gamma0_all + (varxwz_all %*% betaopt))
  B1_ready <- A1
  
  mat_per_GP <- GP_mat(xx = validx, memodel = memodel, betaopt = betaopt, varxwz_all = varxwz_all, optwts = optwts, 
                       all_par = all_par)
  
  vec_CD <- Reduce(c, mat_per_GP[grep("CandD", names(mat_per_GP))])
  
  if (sum(unlist(lapply(validw, ncol))) == 1 & length(validx) == 1) {
    vect_E <- as.matrix(do.call(rbind, mat_per_GP[grep("E", names(mat_per_GP))]))
    derbeta_tes_1 <- rbind(cbind(A1_ready, B1_ready), t(rbind(vec_CD, vect_E)))
  } else if (length(validw) == 1) {
    vect_E <- as.matrix(do.call(rbind, mat_per_GP[grep("E", names(mat_per_GP))]))
    derbeta_tes_1 <- rbind(cbind(A1_ready, B1_ready), t(rbind(vec_CD, vect_E)))
  } else {
    vect_E <- as.matrix(do.call(bdiag, mat_per_GP[grep("E", names(mat_per_GP))]))
    derbeta_tes_1 <- rbind(cbind(A1_ready, B1_ready), cbind(vec_CD, vect_E))
  }
  
  derbeta_tes <- rbind(t(c(1, rep(0, (ncol(derbeta_tes_1) - 1)))), derbeta_tes_1, cbind(-betaopt^2/2, matrix(0, 
                                                                                                             ncol = (ncol(derbeta_tes_1) - 1), nrow = length(betaopt))))
  
  varbeta <- t(derbeta_tes) %*% covall_tes %*% derbeta_tes
  colnames(varbeta) <- NULL
  rownames(varbeta) <- NULL
  
  stdbeta <- sqrt(diag(varbeta))
  
  chiscores <- (betaadj/stdbeta)^2
  odds <- exp(betaadj)
  
  clodds <- cbind(exp((betaadj - qnorm(1 - 0.025) * stdbeta)), exp((betaadj + qnorm(1 - 0.025) * stdbeta)))
  pval <- 1 - pchisq(chiscores, 1)
  
  adj.coeff <- data.frame(beta = betaadj, se = stdbeta, chiscores, odds, clodds, pvalue = pval)
  
  dimnames(adj.coeff) <- list(c("(Intercept)", "GlycanAge"), c("Estimate", "S.E", "Statistic", "Odds", "95% LL", 
                                                               "95%UL", "pval"))
  
  results <- list(alpha_uncorrected = logmodel$df.alpha, alpha_corrected = adj.coeff, betaopt = betaopt, 
                  optwts = optwts, me_results = memodel$df.gamma_1, varxwz_all = varxwz_all)
  return(results)
}

compute_bias_mse <- function(df, true_val) {
  # Function to compute bias and MSE of beta1 estimate 
  # Input: 
  #   df: data frame of simulation results 
  #   true_val: true value of simulation 
  # Output: 
  #   results: data frame bias and MSE for mapping and O2PLS approach
  dataf_res <- list()
  dataf_results <- list()
  
  p = unique(df$n_main)
  q = unique(df$sigmae2)
  for (i in 1:length(p)) {
    
    for (j in 1:length(q)) {
      
      df_o2pls <- df[which(df$n_main == p[i] & df$sigmae2 == q[j]), "beta1_o2pls"]
      df_rc <- df[which(df$n_main == p[i] & df$sigmae2 == q[j]), "beta1_rc"]
      
      bias_o2pls = mean(df_o2pls) - true_val
      mse_o2pls = mse(x = df_o2pls, y = rep(true_val, length(df_o2pls)))
      
      bias_rc = mean(df_rc) - true_val
      mse_rc = mse(x = df_rc, y = rep(true_val, length(df_rc)))
      
      df_res <- data.frame(n_main = p[i], sigmae2 = q[j], bias_o2pls = bias_o2pls, bias_rc = bias_rc, 
                           mse_o2pls = mse_o2pls, mse_rc = mse_rc)
      dataf_res <- rbind(dataf_res, df_res)
    }
    dataf_results[[i]] <- dataf_res
    dataf_res <- c()
  }
  
  dataf_results_fin <- do.call(rbind, dataf_results)
  return(dataf_results_fin)
}

#### 

### Simulation

load("sim_param.RData")

n_main = c(500, 1000)
n_valid_perc = c(0.8)
sigmae2_xwz = c(1, 9.6, 20)
n_sim = 1000
df_results <- df_primary <- df_final <- c()
beta0 = c(1.18, -4)
beta1 = c(0.34, 1.43)

for (m in 1:length(n_main)) {
  for (t in 1:length(n_valid_perc)) {
    for (n in 1:length(sigmae2_xwz)) {
      nsample_main <- n_main[m]
      nsample_valid <- n_valid_perc[t] * n_main[m]
      sigmae2 = sigmae2_xwz[n]
      nsim = n_sim
      
      beta1.est_res <- matrix(0, nsim, ncol = 1)
      beta1.est_adj <- matrix(0, nsim, ncol = 1)
      beta_plugin_mapping <- matrix(0, nsim, ncol = 1)
      beta_plugin_o2pls <- matrix(0, nsim, ncol = 1)
      std.rc.beta1_est <- matrix(0, nsim, ncol = 1)
      gamma1.est_adj <- matrix(0, nsim, ncol = 10)
      std.rc.gamma1_est <- matrix(0, nsim, ncol = 10)
      varxwz_est <- matrix(0, nsim, ncol = 1)
      o2pls_beta1_est <- matrix(0, nsim, ncol = 1)
      ppls_beta1_est <- matrix(0, nsim, ncol = 1)
      varxwz_ppls <- matrix(0, nsim, ncol = 1)
      varxwz_o2pls <- matrix(0, nsim, ncol = 1)
      
      dat_sim_main <- generate_data(N = nsample_main, params = param_from_dat)
      dat_sim_main_w <- dat_sim_main$X %>% scale()  # kor lcms
      
      dat_sim_valid <- generate_data(N = nsample_valid, params = param_from_dat)
      dat_sim_valid_w <- dat_sim_valid$X %>% scale()  # vis lcms
      dat_sim_valid_x <- dat_sim_valid$Y %>% scale()  # vis uplc
      
      gamma = c(56.9787248, 4.0298086, 3.9750948, 0.8623402, -3.6172211, -2.0286535, -0.5873948, -0.675818, 
                -1.1819272, -0.4198655)
      
      for (h in 1:length(beta0)) {
        for (g in 1:length(beta1)) {
          beta = c(beta0[h], beta1[g])
          cat("n_main=", n_main[m], "sigma2=", sigmae2_xwz[n], "beta0=", beta0[h], "beta1=", beta1[g], 
              "\n")
          for (j in 1:nsim) {
            set.seed(j)
            
            # generate data sets
            W_main = list(subset(dat_sim_main_w, select = c("IgG1_G0FN", "IgG2_G0FN", "IgG4_G0FN", 
                                                            "IgG1_G2F", "IgG2_G2F", "IgG4_G2F", "IgG1_G2FN", "IgG2_G2FN", "IgG4_G2FN")))
            W_valid = list(subset(dat_sim_valid_w, select = c("IgG1_G0FN", "IgG2_G0FN", "IgG4_G0FN", 
                                                              "IgG1_G2F", "IgG2_G2F", "IgG4_G2F", "IgG1_G2FN", "IgG2_G2FN", "IgG4_G2FN")))
            
            X_main <- list(W_main[[1]] %*% gamma[-1] + rnorm(nsample_main, mean = 0, sd = sqrt(sigmae2)))
            X_valid <- list(W_valid[[1]] %*% gamma[-1] + rnorm(nsample_valid, mean = 0, sd = sqrt(sigmae2)))
            logitp <- as.matrix(cbind(1, X_main[[1]])) %*% beta
            
            probp <- exp(logitp)/(1 + exp(logitp))
            Y_main <- rbinom(nsample_main, 1, probp)
            
            # approximation regression calibration
            adjfun_contoh <- adjfun(mainy = Y_main, mainw = W_main, validx = X_valid, validw = W_valid, 
                                    wt1 = 1, wt2 = 1, digits = 3)  # ,includeZ=F
            
            # o2pls
            fit_o2pls <- o2m(X = dat_sim_valid_w %>% as.matrix, Y = dat_sim_valid_x %>% as.matrix, 
                             n = 2, nx = 1, ny = 0)
            
            # o2pls calibrated
            Vhat_o2pls <- predict(fit_o2pls, dat_sim_main_w %>% as.matrix, "X")
            T_kor <- dat_sim_main_w %*% fit_o2pls$W. %*% fit_o2pls$B_T.
            adjfun_o2pls <- adjfun(mainy = Y_main, mainw = list(T_kor), validx = X_valid, validw = list(fit_o2pls$Tt), 
                                   wt1 = 1, wt2 = 1, digits = 3)  # ,includeZ=F
            
            beta1.est_adj[j] <- adjfun_contoh$alpha_corrected[2, 1]
            gamma1.est_adj[j, ] <- adjfun_contoh$me_results[, 1]
            varxwz_est[j] <- as.numeric(adjfun_contoh$varxwz_all)
            
            o2pls_beta1_est[j] <- adjfun_o2pls$alpha_corrected[2, 1]
            varxwz_o2pls[j] <- as.numeric(adjfun_o2pls$varxwz_all)
          }
          
          df <- data.frame(n_main = n_main[m], n_valid = n_valid_perc[t] * n_main[m], sigmae2 = sigmae2_xwz[n], 
                           beta0 = beta0[h], beta1 = beta1[g], beta1_rc = beta1.est_adj, beta1_std_rc = std.rc.beta1_est, 
                           beta1_o2pls = o2pls_beta1_est, varxwz_rc = varxwz_est, varxwz_o2pls = varxwz_o2pls)
          
          df_final <- rbind(df_final, df)
        }
      }
    }
  }
}

beta0_true <- unique(df_final$beta0)
beta1_true <- unique(df_final$beta1)
df_final_11 <- df_final[which(df_final$beta0 == beta0_true[1] & df_final$beta1 == beta1_true[1]), ]
df_final_12 <- df_final[which(df_final$beta0 == beta0_true[1] & df_final$beta1 == beta1_true[2]), ]
df_final_21 <- df_final[which(df_final$beta0 == beta0_true[2] & df_final$beta1 == beta1_true[1]), ]
df_final_22 <- df_final[which(df_final$beta0 == beta0_true[2] & df_final$beta1 == beta1_true[2]), ]

df_violation_500 <- rbind(data.frame(beta0 = beta0_true[1], beta1 = beta1_true[1], compute_bias_mse(df = df_final_11, 
                                                                                                    true_val = beta1_true[1])[c(1, 2, 3), c(2, 4, 6, 3, 5)]), data.frame(beta0 = beta0_true[1], beta1 = beta1_true[2], 
                                                                                                                                                                         compute_bias_mse(df = df_final_12, true_val = beta1_true[2])[c(1, 2, 3), c(2, 4, 6, 3, 5)]), data.frame(beta0 = beta0_true[2], 
                                                                                                                                                                                                                                                                                 beta1 = beta1_true[1], compute_bias_mse(df = df_final_21, true_val = beta1_true[1])[c(1, 2, 3), c(2, 4, 
                                                                                                                                                                                                                                                                                                                                                                                   6, 3, 5)]), data.frame(beta0 = beta0_true[2], beta1 = beta1_true[2], compute_bias_mse(df = df_final_22, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                         true_val = beta1_true[2])[c(1, 2, 3), c(2, 4, 6, 3, 5)]))

df_violation_1000 <- rbind(data.frame(beta0 = beta0_true[1], beta1 = beta1_true[1], compute_bias_mse(df = df_final_11, 
                                                                                                     true_val = beta1_true[1])[c(4, 5, 6), c(2, 4, 6, 3, 5)]), data.frame(beta0 = beta0_true[1], beta1 = beta1_true[2], 
                                                                                                                                                                          compute_bias_mse(df = df_final_12, true_val = beta1_true[2])[c(4, 5, 6), c(2, 4, 6, 3, 5)]), data.frame(beta0 = beta0_true[2], 
                                                                                                                                                                                                                                                                                  beta1 = beta1_true[1], compute_bias_mse(df = df_final_21, true_val = beta1_true[1])[c(4, 5, 6), c(2, 4, 
                                                                                                                                                                                                                                                                                                                                                                                    6, 3, 5)]), data.frame(beta0 = beta0_true[2], beta1 = beta1_true[2], compute_bias_mse(df = df_final_22, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                          true_val = beta1_true[2])[c(4, 5, 6), c(2, 4, 6, 3, 5)]))



# full results
xtable(rbind(df_violation_500[1:9, ], df_violation_1000[1:9, ]), digits = 2)

