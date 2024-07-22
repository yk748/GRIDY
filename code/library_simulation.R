#--------------------------------------------------#
# DGP function
#--------------------------------------------------#
sim_model <- function(Joint_rank, Gr_Ind_rank, Num_var, Num_sample, Num_subj_per_group, Control_param){
  
  r_J <- Joint_rank
  r_G <- Gr_Ind_rank
  dd <- Num_var
  TT <- Num_sample
  KK <- Num_subj_per_group
  
  
  if (r_J == 1){
    if (Control_param$type==1){
      Psi_comm <- diag(0.8944,r_J) 
    }else if(Control_param$type==2){
      Psi_comm <- diag(0.8367,r_J)  
    }
  }else if (r_J == 2){
    if (Control_param$type==1){
      Psi_comm <- matrix(c(0.8213,-0.1142,
                           -0.1142,0.8213),ncol=r_J,byrow=TRUE) 
    }else if(Control_param$type==2){
      Psi_comm <- matrix(c(0.8367,0,
                           0,0.8367),ncol=r_J,byrow=TRUE) 
    }
  }else if (r_J == 3){
    if (Control_param$type==1){
      Psi_comm <- matrix(c(0.8154,-0.1377,-0.0298,
                           -0.1377,0.7168,-0.1377,
                           -0.0298,-0.1377,0.8154),ncol=r_J,byrow=TRUE)
    }else if(Control_param$type==2){
      Psi_comm <- matrix(c(0.8367,0,0,
                           0,0.8367,0,
                           0,0,0.8367),ncol=r_J,byrow=TRUE)
    }
  }else if (r_J == 4){
    if (Control_param$type==1){
      Psi_comm <- matrix(c(0.8148,-0.1372,-0.0241,0.0084,
                           -0.1372, 0.7118,-0.1575,-0.0241,
                           -0.0241,-0.1575,0.7118,-0.1372,
                           0.0084,-0.0241,-0.1372,0.8148),ncol=r_J,byrow=TRUE)
    }else if(Control_param$type==2){
      Psi_comm <- matrix(c(0.8367,0,0,0,
                           0,0.8367,0,0,
                           0,0,0.8367,0,
                           0,0,0,0.8367),ncol=r_J,byrow=TRUE)
    }
  }
  
  if (r_G == 1){
    if (Control_param$type==1){
      Psi_ind <- diag(0.8944,r_G) 
    }else if(Control_param$type==2){
      Psi_ind <- diag(0.8367,r_G)  
    }
  }else if (r_G == 2){
    if (Control_param$type==1){
      Psi_ind <- matrix(c(0.8213,-0.1142,
                          -0.1142,0.8213),ncol=r_G,byrow=TRUE) 
    }else if(Control_param$type==2){
      Psi_ind <- matrix(c(0.8367,0,
                          0,0.8367),ncol=r_G,byrow=TRUE) 
    }
  }else if (r_G == 3){
    if (Control_param$type==1){
      Psi_ind <- matrix(c(0.8154,-0.1377,-0.0298,
                          -0.1377,0.7168,-0.1377,
                          -0.0298,-0.1377,0.8154),ncol=r_G,byrow=TRUE)
    }else if(Control_param$type==2){
      Psi_ind <- matrix(c(0.8367,0,0,
                          0,0.8367,0,
                          0,0,0.8367),ncol=r_G,byrow=TRUE)
    }
  }else if (r_G == 4){
    if (Control_param$type==1){
      Psi_ind <- matrix(c(0.8148,-0.1372,-0.0241,0.0084,
                          -0.1372, 0.7118,-0.1575,-0.0241,
                          -0.0241,-0.1575,0.7118,-0.1372,
                          0.0084,-0.0241,-0.1372,0.8148),ncol=r_G,byrow=TRUE)
    }else if(Control_param$type==2){
      Psi_ind <- matrix(c(0.8367,0,0,0,
                          0,0.8367,0,0,
                          0,0,0.8367,0,
                          0,0,0,0.8367),ncol=r_G,byrow=TRUE)
    }
  }
  
  A_comm <- list()
  A_G_ind <- list()
  C_comm <- list()
  C_G_ind <- list()
  F_comm <- list()
  F_G_ind <- list()
  
  Burn_period <- 500 + TT
  for (kk in 1:(2*KK)){
    
    tmp_A_comm <- array(0,dim=c(r_J,(Burn_period+TT)))
    tmp_A_G_ind <- array(0,dim=c(r_G,(Burn_period+TT)))
    noise_comm <- t(rmvnorm((Burn_period+TT),rep(0,r_J),diag(Control_param$sigma_xi,r_J)))
    noise_ind <- t(rmvnorm((Burn_period+TT),rep(0,r_G),diag(Control_param$sigma_xi,r_G)))
    
    for (tt in 1:(Burn_period+TT)){
      if (tt == 1){
        tmp_A_comm[,1] <- noise_comm[,1]
        tmp_A_G_ind[,1] <- noise_ind[,1]
      }else{
        tmp_A_comm[,tt] <- Psi_comm %*% tmp_A_comm[,(tt-1)] + noise_comm[,tt]
        tmp_A_G_ind[,tt] <- Psi_ind %*% tmp_A_G_ind[,(tt-1)] + noise_ind[,tt]
      }
    }
    A_comm[[kk]] <- t(tmp_A_comm[,(Burn_period+1):(Burn_period+TT)])
    A_G_ind[[kk]] <- t(tmp_A_G_ind[,(Burn_period+1):(Burn_period+TT)])
    
    C_comm[[kk]] <- diag(sqrt(Control_param$sigma_c)*c((1+4):(r_J+4)),r_J)
    C_G_ind[[kk]] <-diag(sqrt(Control_param$sigma_c)*c((1+4):(r_G+4)),r_G)
    
    if (r_J == 1){
      F_comm[[kk]] <- C_comm[[kk]] %*% A_comm[[kk]]
    }else{
      F_comm[[kk]] <- A_comm[[kk]] %*% C_comm[[kk]]
    }
    
    if (r_G == 1){
      F_G_ind[[kk]] <- C_G_ind[[kk]] %*% A_G_ind[[kk]]
    }else{
      F_G_ind[[kk]] <- A_G_ind[[kk]] %*% C_G_ind[[kk]]
    }
  }
  
  
  idx_comm <- sort(sample(c(1:dd),round(Control_param$J_per*dd),replace=FALSE),decreasing=FALSE)
  idx_group1 <- sort(sample(c(1:dd)[-idx_comm],round((1-Control_param$J_per)/2*dd),replace=FALSE),decreasing=FALSE)
  idx_group2 <- c(1:dd)[-c(idx_group1,idx_comm)]
  
  B_comm <- matrix(runif((dd*r_J),Control_param$b_bar_min,Control_param$b_bar_max),nrow=dd,ncol=r_J)
  B_comm[-idx_comm,] <- 0
  
  B_group1 <- matrix(runif((dd*r_G),Control_param$b_tilde_min,Control_param$b_tilde_max),nrow=dd,ncol=r_G)
  B_group1[-idx_group1,] <- 0
  B_group2 <- matrix(runif((dd*r_G),Control_param$b_tilde_min,Control_param$b_tilde_max),nrow=dd,ncol=r_G)
  B_group2[-idx_group2,] <- 0
  
  X <- list()
  X_comm <- list()
  X_G_ind <- list()
  E <- list()
  
  scaled_X <- list()
  scaled_X_comm <- list()
  scaled_X_G_ind <- list()
  
  for (kk in 1:(2*KK)){
    noise_var <- Control_param$sigma_eps
    if (noise_var >0 ){
      E[[kk]] <- rmvnorm(TT,rep(0,dd),noise_var*diag(1,dd))
    }else{
      E[[kk]] <- array(0,dim=c(TT,dd))
    }
    
    if (r_J == 1){
      X_comm[[kk]] <- t(F_comm[[kk]]) %*% t(B_comm)
    }else{
      X_comm[[kk]] <- F_comm[[kk]] %*% t(B_comm)
    }
    
    if (kk <= KK){
      if (r_G == 1){
        X_G_ind[[kk]] <- t(F_G_ind[[kk]]) %*% t(B_group1)
      }else{
        X_G_ind[[kk]] <- F_G_ind[[kk]] %*% t(B_group1)
      }
    }else{
      if (r_G == 1){
        X_G_ind[[kk]] <- t(F_G_ind[[kk]]) %*% t(B_group2)
      }else{
        X_G_ind[[kk]] <- F_G_ind[[kk]] %*% t(B_group2)
      }
    }
    
    X[[kk]] <- X_comm[[kk]] + X_G_ind[[kk]] + E[[kk]]
    scaled_X[[kk]] <-  t( scale(X[[kk]], center = TRUE, scale = TRUE) )
    scaled_X[[kk]][is.nan(scaled_X[[kk]])] <- 0
    scaled_X_comm[[kk]] <-  t( scale(X_comm[[kk]], center = TRUE, scale = TRUE) )
    scaled_X_comm[[kk]][is.nan(scaled_X_comm[[kk]])] <- 0
    scaled_X_G_ind[[kk]] <-  t (scale(X_G_ind[[kk]], center = TRUE, scale = TRUE) )
    scaled_X_G_ind[[kk]][is.nan(scaled_X_G_ind[[kk]])] <- 0
  }
  
  output <- list(Xk = X, X_bar = X_comm, X_tilde = X_G_ind, Ek = E,
                 scaled_X = scaled_X, scaled_X_comm = scaled_X_comm, 
                 scaled_X_G_ind = scaled_X_G_ind, 
                 F_bar = F_comm, F_tilde = F_G_ind,
                 A_bar = A_comm, A_tilde = A_G_ind,
                 C_bar = C_comm, C_tilde = C_G_ind,
                 Psi_bar = Psi_comm, Psi_tilde = Psi_ind,
                 B_bar= B_comm, B_tilde1 = B_group1, B_tilde2 = B_group2,
                 comm_idx = idx_comm, group1_idx = idx_group1, group2_idx = idx_group2)
  return(output)
}


#--------------------------------------------------#
# Factor_regression function
#--------------------------------------------------#
factor_regression <- function(KK,r_J,r_G,DGP,Series_joint,Series_group1,Series_group2){
  
  factor_regression_joint <- list()
  factor_regression_group1 <- list()
  factor_regression_group2 <- list()
  
  for (kk in 1:(2*KK)){
    if (kk <= KK){
      B_aug <- cbind(Series_joint$B,Series_group1$B)      
    }else{
      B_aug <- cbind(Series_joint$B,Series_group2$B)
    }
    factor_tmp <- ( t(DGP$scaled_X[[kk]]) %*% B_aug ) %*% solve (t(B_aug) %*% B_aug)
    factor_regression_joint[[kk]] <- as.matrix(factor_tmp[,1:r_J])
    
    if (kk <= KK){
      factor_regression_group1[[kk]] <- as.matrix(factor_tmp[,(r_J+1):(r_J+r_G)])
    }else{
      factor_regression_group2[[(kk-KK)]] <- as.matrix(factor_tmp[,(r_J+1):(r_J+r_G)])
    }
  }
  output <- list(factor_joint = factor_regression_joint,
                 factor_group1 = factor_regression_group1,
                 factor_group2 = factor_regression_group2)
  return(output)
}


#--------------------------------------------------#
# Yule-Walker function
#--------------------------------------------------#
YK_compute <- function(dd,TT,KK,r_J,r_G,Control_param,factor_list){
  
  Psi_joint_hat <- array(NA,dim=c(r_J,r_J,(2*KK)))
  Eta_joint_hat <- array(NA,dim=c(r_J,r_J,(2*KK)))
  
  Psi_indiv_hat <- array(NA,dim=c(r_G,r_G,(2*KK)))
  Eta_G_indiv_hat <- array(NA,dim=c(r_G,r_G,(2*KK)))
  
  for (kk in 1:(2*KK)){
    Cov_joint_F0 <- t(scale(factor_list$factor_joint[[kk]],center=TRUE,scale=FALSE)) %*% scale(factor_list$factor_joint[[kk]],center=TRUE,scale=FALSE) / TT
    Cov_joint_F1 <- t(scale(factor_list$factor_joint[[kk]][-1,],center=TRUE,scale=FALSE)) %*% scale(factor_list$factor_joint[[kk]][-TT,],center=TRUE,scale=FALSE) / TT
    
    Psi_joint_hat[,,kk] <- t(solve(Cov_joint_F0) %*% t(Cov_joint_F1))
    Eta_joint_hat[,,kk] <- Cov_joint_F0 - Psi_joint_hat[,,kk]%*%Cov_joint_F1
    
    if (kk <= KK){
      Cov_indiv_F0 <- t(scale(factor_list$factor_group1[[kk]],center=TRUE,scale=FALSE)) %*% scale(factor_list$factor_group1[[kk]],center=TRUE,scale=FALSE) / TT
      Cov_indiv_F1 <- t(scale(factor_list$factor_group1[[kk]][-1,],center=TRUE,scale=FALSE)) %*% scale(factor_list$factor_group1[[kk]][-TT,],center=TRUE,scale=FALSE) / TT
    }else{
      Cov_indiv_F0 <- t(scale(factor_list$factor_group2[[(kk-KK)]],center=TRUE,scale=FALSE)) %*% scale(factor_list$factor_group2[[(kk-KK)]],center=TRUE,scale=FALSE) / TT
      Cov_indiv_F1 <- t(scale(factor_list$factor_group2[[(kk-KK)]][-1,],center=TRUE,scale=FALSE)) %*% scale(factor_list$factor_group2[[(kk-KK)]][-TT,],center=TRUE,scale=FALSE) / TT
    }
    
    Psi_indiv_hat[,,kk] <- t(solve(Cov_indiv_F0) %*% t(Cov_indiv_F1))
    Eta_G_indiv_hat[,,kk] <- Cov_indiv_F0 - Psi_indiv_hat[,,kk]%*%Cov_indiv_F1
  }
  output <- list(Psi_joint_hat = Psi_joint_hat, Eta_joint_hat = Eta_joint_hat,
                 Psi_indiv_hat = Psi_indiv_hat, Eta_G_indiv_hat = Eta_G_indiv_hat)
  return(output)
}


