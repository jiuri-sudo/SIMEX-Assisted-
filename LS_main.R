log_dir <- "F:/MCcode/real_data/simulation_jinzhan"  
if(!dir.exists(log_dir)) dir.create(log_dir)  
log_file <- file.path(log_dir, paste0("MC_run_", Sys.Date(), ".txt"))
library(parallel)



datagenerator_linear <- function(d,K,sd_Q,sd_P,nQ,sample_prop,B_max_value,Omega_diag_value,thetaQ){
  
  nP <- nQ * sample_prop
  N <- nP+nQ
  
  SigmaXQ <- matrix(0, nrow = d, ncol = d)
  for (i in 1:d) {
    for (j in 1:d) {
      SigmaXQ[i, j] <- 0.5^abs(i - j) 
    }
  }
  SigmaXP <- SigmaXQ  
  epsilonQ <- rnorm(nQ, mean = 0, sd = sd_Q)
  epsilonP <- rnorm(nP, mean = 0, sd = sd_P)
  B <- as.matrix(round(runif(d, 0, B_max_value), 2))
  B_mat <- matrix(B, nrow = d, ncol = nP)
  
  rho <- 0.8
  Omega <- outer(1:d,1:d,function(i, j) {Omega_diag_value * rho^abs(i - j)})
  Z_mat <- matrix(rnorm(d * (K*nP), mean = 0, sd = 0.25), nrow = d, ncol = (K*nP))  
  U_mat <- matrix(rnorm(d * nP, mean = 0, sd = 0.25), nrow = d, ncol = nP)
  
  XQ <- t(MASS::mvrnorm(n = nQ, mu = rep(0, d), Sigma = SigmaXQ))
  XP <- t(MASS::mvrnorm(n = nP, mu = rep(0, d), Sigma = SigmaXP))

  yQ <- as.matrix(t(thetaQ) %*% XQ) + epsilonQ
  yP <- as.matrix(t(thetaQ) %*% XP) + epsilonP
  WP <- XP + B_mat + sqrt(Omega) %*% U_mat
  
  data<-list(
    SigmaXQ = SigmaXQ, SigmaXP = SigmaXP, epsilonQ = epsilonQ, epsilonP = epsilonP,
    B = B, B_mat= B_mat, Omega = Omega, Z_mat = Z_mat, U_mat = U_mat,
    XQ = XQ, XP = XP, yQ = yQ, yP = yP, WP = WP, N = N, nP=nP, nQ=nQ
  )
  return(data)
  
}




MC_Method<-function(d,K,N,w1_seq,nQ,XQ,yQ,nP,WP,yP,Z_mat,B,thetaQ_true,Omega){
  
  error<-c()
  error_min <-NA
  best_w1 <-NA
  MC_LS_matrix <- matrix(0, nrow = d, ncol = length(w1_seq))
  MCE1 <- matrix(0, nrow = d, ncol = 1)
  B_mat <- matrix(B, nrow = d, ncol = nP)
  
  for (i in seq_along(w1_seq)) {
    w1 <- w1_seq[i]
    
    
    sum_matrix1 <- Reduce("+", lapply(1:nQ, function(j) {
      xj <- XQ[, j, drop = FALSE]  
      xj %*% t(xj)                
    }))
    sum_matrix2 <- Reduce("+", lapply(1:nP, function(j) {
      WPj_B <- WP[, j, drop = FALSE]-B  
      WPj_B %*% t(WPj_B)                
    }))
  
    sum_matrix3 <- matrix(0, nrow = d, ncol = d) 
    for (j in 1:ncol(Z_mat)) {
      col_vec <- Z_mat[, j, drop = FALSE] 
      sum_matrix3 <- sum_matrix3 + col_vec %*% t(col_vec)
    }
    hatM <- 1/(K*N) * sum_matrix3
    Omega_half <- expm::sqrtm(Omega)
    A_part <- (w1/N) * sum_matrix1 + ((1-w1)/N) * sum_matrix2 - (1-w1)* Omega_half %*% hatM %*% Omega_half
    
    sum_matrix4 <- matrix(0, nrow = d, ncol = 1) 
    for (j in 1:nQ) {
      XjQyjQ <-  XQ[, j, drop = FALSE] * yQ[j]  
      sum_matrix4 <- sum_matrix4 + XjQyjQ
    }
    sum_matrix5 <- matrix(0, nrow = d, ncol = 1) 
    for (j in 1:nP) {
      WjpByjQ <-  as.matrix(WP[,j]-B_mat[,j]) * yP[j]  
      sum_matrix5 <- sum_matrix5 + WjpByjQ
    }
    B_part <- (w1/N) * sum_matrix4 + ((1-w1)/N) * sum_matrix5 
    thetaQ_estimator<- solve(A_part) %*% B_part 
    
    MC_LS_matrix[,i]<-thetaQ_estimator
    error[i] <- sum((thetaQ_estimator - thetaQ_true)^2)
  }
  error_min <- min(error)
  best_w1 <- w1_seq[which.min(error)]
  MCE1 <-  MC_LS_matrix[,which.min(error)]
  
  return(list(error_min = error_min, best_w1 = best_w1, MCE1 = MCE1))
}



para_weight_Method<-function(d,K,N,nQ,XQ,yQ,nP,WP,yP,w1_seq,B,thetaQ_true,Omega){
  
 
  betaQ_OLS <- solve(XQ %*% t(XQ)) %*% XQ %*% t(yQ)
  betaP_OLS <- solve(WP %*% t(WP)) %*% WP %*% t(yP)
  
  error <- c()
  for (i in seq_along(w1_seq)) {
    
    w1 <- w1_seq[i]
    weight_beta <- w1*betaQ_OLS + (1-w1)*betaP_OLS
    error[i] <- sum((weight_beta - thetaQ_true)^2)
   
  }
  weight_error <- min(error)
  best_w1 <- w1_seq[which.min(error)]
  return(list(weight_error = weight_error, best_w1 = best_w1))
}


Lossfun_weight_Method <- function(thetaQ_true,XQ, yQ, WP, yP, w1_seq){
  
  error <- c()
  for (i in seq_along(w1_seq)) {
    
    w1 <- w1_seq[i]
    A <- w1 *XQ %*% t(XQ) + (1-w1) * WP %*% t(WP)
    b <-  w1 * XQ %*% t(yQ) + (1-w1) * WP %*% t(yP)
    theta_hat <- solve(A) %*% b
    
    error[i] <- sum((theta_hat - thetaQ_true)^2)
  }
  error_min <- min(error)
  best_w1 <- w1_seq[which.min(error)]
  
  return(list(theta_hat = theta_hat, error_min = error_min,best_w1=best_w1))
}


alldata_Method <- function(XQ,WP,yQ,yP,thetaQ_true){
  
  all_X <- cbind(XQ,WP)
  all_y <- cbind(yQ,yP)
  all_beta <- solve(all_X %*% t(all_X)) %*% all_X %*% t(all_y)
  error <- sum((all_beta - thetaQ_true)^2)
  return(list(error=error))
}



modelQ_Method <- function(XQ,yQ,thetaQ_true){
  betaQ_OLS <- solve(XQ %*% t(XQ)) %*% XQ %*% t(yQ)
  error <- sum((betaQ_OLS - thetaQ_true)^2)
  return(list(error=error))
}


# study_seq_name = "Differ Omega values"
# Omega_diag_value_seq <- seq(0.01,0.05,0.01)
# study_seq <- Omega_diag_value_seq


# study_seq_name = "Differ K values"
# K_seq <-  seq(20,120,20)
# K_seq <-  c(1,20,100,500)
# study_seq <- K_seq


study_seq_name = "Differ B values"
B_max_value_seq <- seq(0.2,0.3,0.02)
study_seq <- B_max_value_seq


# study_seq_name = "sample_prop"
# sample_prop_seq <- seq(2,7,1)
# study_seq <- sample_prop_seq



sd_Q <- 0.5
sd_P <- 0.5
d <- 100  
non_zero_num <- 50    
thetaQ <- c(rep(0.4, non_zero_num), rep(0, d-non_zero_num))   
thetaP <- thetaQ
reputation <- 300



MCE_all  <- para_weight_all <- Lossfun_weight_all <- alldata_all <- modelQ_all <- matrix(NA, nrow = reputation, ncol = length(study_seq))
Omega_trace_all <- B_L1_all <- Omega_L1_all <- matrix(NA, nrow = reputation, ncol = length(study_seq))
MCE_w1 <- para_weight_w1 <- Lossfun_weight_w1 <- numeric(length(study_seq))



for (q_idx in seq_along(study_seq)) {

  # Omega_diag_value <- study_seq[q_idx]
  # sample_prop <- 2
  # nQ <- 120
  # K <- 20
  # B_max_value <- 0.2

  
  # K <- study_seq[q_idx]
  # sample_prop <- 5
  # nQ <- 110
  # B_max_value <- 0.4
  # Omega_diag_value <- 0.5
  
  
  B_max_value <- study_seq[q_idx]
  sample_prop <- 2
  nQ <- 120
  K <- 20
  Omega_diag_value <- 0.01
  
  # sample_prop <- study_seq[q_idx]
  # K <- 20
  # B_max_value <- 0.5
  # Omega_diag_value <- 0.5
  
  
  q_idx_now <- q_idx
  
  cat(paste0(
    Sys.time(), 
    " | START For study_seq = Omega_diag_value: ", Omega_diag_value,
    " | K: ", K,
    " | sample_prop: ", sample_prop,
    " | B_max_value: ", B_max_value,
    " | Omega_diag_value: ", Omega_diag_value,
    " | nQ: ", nQ,
    " (", q_idx, "/", length(study_seq), ")\n"
  ),
  file = log_file, append = TRUE)
  
  
  n_core <- 14
  cl <- makeCluster(n_core)
  clusterEvalQ(cl, {
    library(MASS)
    library(Matrix)
    library(glmnet)
    library(dplyr)
    library(expm)
  })
  clusterExport(cl, c("datagenerator_linear","MC_Method","para_weight_Method","Lossfun_weight_Method","alldata_Method","modelQ_Method",
                      "log_file", "log_dir",
                      "sd_Q","sd_P","d","thetaQ","thetaP",
                      "nQ" ,"K","sample_prop",
                      "B_max_value"
                      ,"Omega_diag_value"
                      ,"q_idx_now"
  ))
  clusterEvalQ(cl, { options(stringsAsFactors = FALSE, encoding="UTF-8") })
  
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(12345)
  clusterSetRNGStream(cl, 12345)
  
  rep_list <- 1:reputation
  

  res_par <- parLapply(cl, rep_list, function(rep) {
    
    rep_id <- rep

    q_idx_node <- q_idx_now
    
    MCE_error <-  para_weight_error <-  Lossfun_weight_error <-  alldata_error <- modelQ_error <- NA
    MCE_w1 <- para_weight_w1 <- Lossfun_weight_w1 <- NA
    
    
    data_linear<-datagenerator_linear(d,K,sd_Q,sd_P,nQ,sample_prop,B_max_value,Omega_diag_value,thetaQ)
    
    
    SigmaXQ <- data_linear$SigmaXQ
    SigmaXP <- data_linear$SigmaXP
    
    B <- data_linear$B
    B_mat <- data_linear$B_mat
    Omega <- data_linear$Omega
    
    Z_mat <- data_linear$Z_mat
    U_mat <- data_linear$U_mat
    
    XQ <- data_linear$XQ
    XP <- data_linear$XP
    yQ <- data_linear$yQ
    yP <- data_linear$yP
    WP <- data_linear$WP
    
    N  <- data_linear$N
    nP  <- data_linear$nP
    nQ  <- data_linear$nQ
    
    
    w1_seq_part <- seq(1/(1+sample_prop)*0.5,min(5 / (1 + sample_prop), 1),0.03)
    w1_seq <- seq(0.05,0.95,0.05) 
    
    modelQ_results <- modelQ_Method(XQ,yQ,thetaQ_true=thetaQ)
    modelQ_error <- modelQ_results$error
    
    MCE_results <- MC_Method(d,K,N,w1_seq=w1_seq_part,nQ,XQ,yQ,nP,WP,yP,Z_mat,B,thetaQ_true=thetaQ,Omega)
    MCE_error <- MCE_results$error_min
    MCE_w1 <- MCE_results$best_w1
    
    para_weight_results <- para_weight_Method(d,K,N,nQ,XQ,yQ,nP,WP,yP,w1_seq,B,thetaQ_true=thetaQ,Omega)
    para_weight_error <- para_weight_results$weight_error
    para_weight_w1 <- para_weight_results$best_w1
    
    Lossfun_weight_results<- Lossfun_weight_Method(thetaQ_true=thetaQ,XQ, yQ, WP, yP,w1_seq)
    Lossfun_weight_error <- Lossfun_weight_results$error_min
    Lossfun_weight_w1 <-  Lossfun_weight_results$best_w1
    
    alldata_results <- alldata_Method(XQ,WP,yQ,yP,thetaQ_true=thetaQ)
    alldata_error <- alldata_results$error
    
    cat(paste0(Sys.time(), " | END   | rep = ", rep_id, "\n"), file = log_file,append = TRUE)
    
    list(
      MCE_error = MCE_error,
      para_weight_error = para_weight_error,
      Lossfun_weight_error = Lossfun_weight_error,
      alldata_error = alldata_error,
      modelQ_error = modelQ_error,
      MCE_w1 = MCE_w1,
      para_weight_w1 = para_weight_w1,
      Lossfun_weight_w1 = Lossfun_weight_w1
      ,q_idx_node = q_idx_node  
    )
  })
  
  stopCluster(cl)
  
  for (r in 1:reputation) {
    res <- res_par[[r]]
    
    idx <- res$q_idx_node
    
    cat(
      paste0(
        Sys.time(), 
        " | Storing results | rep = ", r, 
        "\n"
      ),
      file = log_file, append = TRUE
    )
    
    
    MCE_all[r, idx]            <- res$MCE_error
    para_weight_all[r, idx]    <- res$para_weight_error
    Lossfun_weight_all[r, idx] <- res$Lossfun_weight_error
    alldata_all[r, idx]        <- res$alldata_error
    modelQ_all[r, idx]         <- res$modelQ_error
    
  }
  
  # 记录 w1
  MCE_w1[q_idx]          <- mean(sapply(res_par, function(x) x$MCE_w1), na.rm = TRUE)
  para_weight_w1[q_idx]  <- mean(sapply(res_par, function(x) x$para_weight_w1), na.rm = TRUE)
  Lossfun_weight_w1[q_idx] <- mean(sapply(res_par, function(x) x$Lossfun_weight_w1), na.rm = TRUE)
}




library(reshape2)
library(ggplot2)

make_long <- function(mat, label) {
  df <- as.data.frame(mat)
  data.frame(
    Column = rep(1:ncol(df), each = nrow(df)),
    Value = as.vector(as.matrix(df)),
    Method = label
  )
}


df_all <- rbind(
  make_long(MCE_all, "SIMEX"),
  make_long(para_weight_all, "Param-W"),
  make_long(Lossfun_weight_all, "Loss-W"),
  make_long(alldata_all, "Full"),
  make_long(modelQ_all, "Naive")
)


method_colors <- c(
  "SIMEX" = "red3",
  "Param-W" = "green4",
  "Loss-W" = "blue3",
  "Full" = "orange2",
  "Naive" = "purple2"
)

df_all$Method <- factor(df_all$Method, levels = c("Naive", "Full","Param-W", "Loss-W",  "SIMEX"))

df_all$study_seq_value <- rep(study_seq, each = reputation, times = 5)


if (study_seq_name == "Differ Omega values") {
  x_label <- expression(omega)
} else if (study_seq_name == "Differ K values") {
  x_label <- expression(k)
} else if (study_seq_name == "Differ B values") {
  x_label <- expression(b)
} else if (study_seq_name == "sample_prop") {
  x_label <- expression("Different Sample Proportion (" * n^P * " vs. " * n^Q * ")")
} else if (study_seq_name == "split_rate") {
  x_label <- expression("Different Split Ratios (" * n[1]^Q * " vs. " * n^Q * ")")
} else if (study_seq_name == "Differ d values") {
  x_label <- expression("Different " * d * " Values")
} 



# 绘图
ggplot(df_all, aes(x = as.factor(study_seq_value), y = Value, fill = Method)) +
  geom_boxplot(alpha = 0.55, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  stat_summary(aes(group = Method, color = Method),
               fun = mean, geom = "point", shape = 21, size = 2,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  labs(
    x = x_label,
    y = expression("Squared " ~ L[2] ~ " error")
  )+
  coord_cartesian(ylim = c(0, 8.20798109734609726662461))+
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.spacing.y = unit(0, "cm"),
    legend.box.margin = margin(5, 0, -100, 0),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )
