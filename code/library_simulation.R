#--------------------------------------------------#
# DGP function
#--------------------------------------------------#
sim_model <- function(Joint_rank, Ind_rank, Num_var, Num_sample, Num_subj_per_group, Control_param){
  
  r_J <- Joint_rank
  r_I <- Ind_rank
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
  
  if (r_I == 1){
    if (Control_param$type==1){
      Psi_ind <- diag(0.8944,r_I) 
    }else if(Control_param$type==2){
      Psi_ind <- diag(0.8367,r_I)  
    }
  }else if (r_I == 2){
    if (Control_param$type==1){
      Psi_ind <- matrix(c(0.8213,-0.1142,
                          -0.1142,0.8213),ncol=r_I,byrow=TRUE) 
    }else if(Control_param$type==2){
      Psi_ind <- matrix(c(0.8367,0,
                          0,0.8367),ncol=r_I,byrow=TRUE) 
    }
  }else if (r_I == 3){
    if (Control_param$type==1){
      Psi_ind <- matrix(c(0.8154,-0.1377,-0.0298,
                          -0.1377,0.7168,-0.1377,
                          -0.0298,-0.1377,0.8154),ncol=r_I,byrow=TRUE)
    }else if(Control_param$type==2){
      Psi_ind <- matrix(c(0.8367,0,0,
                          0,0.8367,0,
                          0,0,0.8367),ncol=r_I,byrow=TRUE)
    }
  }else if (r_I == 4){
    if (Control_param$type==1){
      Psi_ind <- matrix(c(0.8148,-0.1372,-0.0241,0.0084,
                          -0.1372, 0.7118,-0.1575,-0.0241,
                          -0.0241,-0.1575,0.7118,-0.1372,
                          0.0084,-0.0241,-0.1372,0.8148),ncol=r_I,byrow=TRUE)
    }else if(Control_param$type==2){
      Psi_ind <- matrix(c(0.8367,0,0,0,
                          0,0.8367,0,0,
                          0,0,0.8367,0,
                          0,0,0,0.8367),ncol=r_I,byrow=TRUE)
    }
  }
  
  A_comm <- list()
  A_ind <- list()
  C_comm <- list()
  C_ind <- list()
  F_comm <- list()
  F_ind <- list()
  
  Burn_period <- 500 + TT
  for (kk in 1:(2*KK)){
    
    tmp_A_comm <- array(0,dim=c(r_J,(Burn_period+TT)))
    tmp_A_ind <- array(0,dim=c(r_I,(Burn_period+TT)))
    noise_comm <- t(rmvnorm((Burn_period+TT),rep(0,r_J),diag(Control_param$sigma_xi,r_J)))
    noise_ind <- t(rmvnorm((Burn_period+TT),rep(0,r_I),diag(Control_param$sigma_xi,r_I)))
    
    for (tt in 1:(Burn_period+TT)){
      if (tt == 1){
        tmp_A_comm[,1] <- noise_comm[,1]
        tmp_A_ind[,1] <- noise_ind[,1]
      }else{
        tmp_A_comm[,tt] <- Psi_comm %*% tmp_A_comm[,(tt-1)] + noise_comm[,tt]
        tmp_A_ind[,tt] <- Psi_ind %*% tmp_A_ind[,(tt-1)] + noise_ind[,tt]
      }
    }
    A_comm[[kk]] <- t(tmp_A_comm[,(Burn_period+1):(Burn_period+TT)])
    A_ind[[kk]] <- t(tmp_A_ind[,(Burn_period+1):(Burn_period+TT)])
    
    C_comm[[kk]] <- diag(runif(r_J,min=Control_param$c_min,
                               max=Control_param$c_max), r_J)
    C_ind[[kk]] <- diag(runif(r_I,min=Control_param$c_min,
                              max=Control_param$c_max), r_I)
    
    if (r_J == 1){
      F_comm[[kk]] <- C_comm[[kk]] %*% A_comm[[kk]]
    }else{
      F_comm[[kk]] <- A_comm[[kk]] %*% C_comm[[kk]]
    }
    
    if (r_I == 1){
      F_ind[[kk]] <- C_ind[[kk]] %*% A_ind[[kk]]
    }else{
      F_ind[[kk]] <- A_ind[[kk]] %*% C_ind[[kk]]
    }
  }
  
  
  idx_comm <- sort(sample(c(1:dd),round(0.5*dd),replace=FALSE),decreasing=FALSE)
  idx_group1 <- sort(sample(c(1:dd)[-idx_comm],round(0.25*dd),replace=FALSE),decreasing=FALSE)
  idx_group2 <- c(1:dd)[-c(idx_group1,idx_comm)]
  
  B_comm <- matrix(runif((dd*r_J),Control_param$b_min,Control_param$b_max),nrow=dd,ncol=r_J)
  B_comm[-idx_comm,] <- 0
  
  B_group1 <- matrix(runif((dd*r_I),Control_param$b_min,Control_param$b_max),nrow=dd,ncol=r_I)
  B_group1[-idx_group1,] <- 0
  B_group2 <- matrix(runif((dd*r_I),Control_param$b_min,Control_param$b_max),nrow=dd,ncol=r_I)
  B_group2[-idx_group2,] <- 0
  
  X <- list()
  X_comm <- list()
  X_ind <- list()
  E <- list()
  
  scaled_X <- list()
  scaled_X_comm <- list()
  scaled_X_ind <- list()
  
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
      if (r_I == 1){
        X_ind[[kk]] <- t(F_ind[[kk]]) %*% t(B_group1)
      }else{
        X_ind[[kk]] <- F_ind[[kk]] %*% t(B_group1)
      }
    }else{
      if (r_I == 1){
        X_ind[[kk]] <- t(F_ind[[kk]]) %*% t(B_group2)
      }else{
        X_ind[[kk]] <- F_ind[[kk]] %*% t(B_group2)
      }
    }
    
    X[[kk]] <- X_comm[[kk]] + X_ind[[kk]] + E[[kk]]
    scaled_X[[kk]] <-  t( scale(X[[kk]], center = TRUE, scale = TRUE) )
    scaled_X[[kk]][is.nan(scaled_X[[kk]])] <- 0
    scaled_X_comm[[kk]] <-  t( scale(X_comm[[kk]], center = TRUE, scale = TRUE) )
    scaled_X_comm[[kk]][is.nan(scaled_X_comm[[kk]])] <- 0
    scaled_X_ind[[kk]] <-  t (scale(X_ind[[kk]], center = TRUE, scale = TRUE) )
    scaled_X_ind[[kk]][is.nan(scaled_X_ind[[kk]])] <- 0
  }
  
  output <- list(Xk = X, X_bar = X_comm, X_tilde = X_ind, Ek = E,
                 scaled_X = scaled_X, scaled_X_comm = scaled_X_comm, scaled_X_ind = scaled_X_ind,
                 F_bar = F_comm, F_tilde = F_ind,
                 A_bar = A_comm, A_tilde = A_ind,
                 C_bar = C_comm, C_tilde = C_ind,
                 B_bar= B_comm, B_tilde1 = B_group1, B_tilde2 = B_group2,
                 comm_idx = idx_comm, group1_idx = idx_group1, group2_idx = idx_group2)
  return(output)
}


#--------------------------------------------------#
# Factor_regression function
#--------------------------------------------------#
factor_regression <- function(KK,r_J,r_I,DGP,Series_joint,Series_group1,Series_group2){
  
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
      factor_regression_group1[[kk]] <- as.matrix(factor_tmp[,(r_J+1):(r_J+r_I)])
    }else{
      factor_regression_group2[[(kk-KK)]] <- as.matrix(factor_tmp[,(r_J+1):(r_J+r_I)])
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
YK_compute <- function(dd,TT,KK,r_J,r_I,Control_param,DGP,factor_list){
  
  
  if (r_J == 1){
    if (Control_param$type==1){
      Psi_joint <- diag(0.8944,r_J) 
    }else if(Control_param$type==2){
      Psi_joint <- diag(0.8367,r_J)  
    }
  }else if (r_J == 2){
    if (Control_param$type==1){
      Psi_joint <- matrix(c(0.8213,-0.1142,
                            -0.1142,0.8213),ncol=r_J,byrow=TRUE) 
    }else if(Control_param$type==2){
      Psi_joint <- matrix(c(0.8367,0,
                            0,0.8367),ncol=r_J,byrow=TRUE) 
    }
  }else if (r_J == 3){
    if (Control_param$type==1){
      Psi_joint <- matrix(c(0.8154,-0.1377,-0.0298,
                            -0.1377,0.7168,-0.1377,
                            -0.0298,-0.1377,0.8154),ncol=r_J,byrow=TRUE)
    }else if(Control_param$type==2){
      Psi_joint <- matrix(c(0.8367,0,0,
                            0,0.8367,0,
                            0,0,0.8367),ncol=r_J,byrow=TRUE)
    }
  }else if (r_J == 4){
    if (Control_param$type==1){
      Psi_joint <- matrix(c(0.8148,-0.1372,-0.0241,0.0084,
                            -0.1372, 0.7118,-0.1575,-0.0241,
                            -0.0241,-0.1575,0.7118,-0.1372,
                            0.0084,-0.0241,-0.1372,0.8148),ncol=r_J,byrow=TRUE)
    }else if(Control_param$type==2){
      Psi_joint <- matrix(c(0.8367,0,0,0,
                            0,0.8367,0,0,
                            0,0,0.8367,0,
                            0,0,0,0.8367),ncol=r_J,byrow=TRUE)
    }
  }
  
  if (r_I == 1){
    if (Control_param$type==1){
      Psi_indiv <- diag(0.8944,r_I) 
    }else if(Control_param$type==2){
      Psi_indiv <- diag(0.8367,r_I)  
    }
  }else if (r_I == 2){
    if (Control_param$type==1){
      Psi_indiv <- matrix(c(0.8213,-0.1142,
                            -0.1142,0.8213),ncol=r_I,byrow=TRUE) 
    }else if(Control_param$type==2){
      Psi_indiv <- matrix(c(0.8367,0,
                            0,0.8367),ncol=r_I,byrow=TRUE) 
    }
  }else if (r_I == 3){
    if (Control_param$type==1){
      Psi_indiv <- matrix(c(0.8154,-0.1377,-0.0298,
                            -0.1377,0.7168,-0.1377,
                            -0.0298,-0.1377,0.8154),ncol=r_I,byrow=TRUE)
    }else if(Control_param$type==2){
      Psi_indiv <- matrix(c(0.8367,0,0,
                            0,0.8367,0,
                            0,0,0.8367),ncol=r_I,byrow=TRUE)
    }
  }else if (r_I == 4){
    if (Control_param$type==1){
      Psi_indiv <- matrix(c(0.8148,-0.1372,-0.0241,0.0084,
                            -0.1372, 0.7118,-0.1575,-0.0241,
                            -0.0241,-0.1575,0.7118,-0.1372,
                            0.0084,-0.0241,-0.1372,0.8148),ncol=r_I,byrow=TRUE)
    }else if(Control_param$type==2){
      Psi_indiv <- matrix(c(0.8367,0,0,0,
                            0,0.8367,0,0,
                            0,0,0.8367,0,
                            0,0,0,0.8367),ncol=r_I,byrow=TRUE)
    }
  }
  
  
  Psi_joint_k <- array(NA,dim=c(r_J,r_J,(2*KK)))
  Psi_indiv_k <- array(NA,dim=c(r_I,r_I,(2*KK)))
  
  Psi_joint_hat <- array(NA,dim=c(r_J,r_J,(2*KK)))
  Eta_joint_hat <- array(NA,dim=c(r_J,r_J,(2*KK)))
  
  Psi_indiv_hat <- array(NA,dim=c(r_I,r_I,(2*KK)))
  Eta_indiv_hat <- array(NA,dim=c(r_I,r_I,(2*KK)))
  
  for (kk in 1:(2*KK)){
    
    Psi_joint_k[,,kk] <- DGP$C_bar[[kk]] %*% Psi_joint %*% solve(DGP$C_bar[[kk]])
    
    Cov_joint_F0 <- t(factor_list$factor_joint[[kk]]) %*% factor_list$factor_joint[[kk]] / TT
    Cov_joint_F1 <- t(factor_list$factor_joint[[kk]][-1,]) %*% factor_list$factor_joint[[kk]][-TT,] / TT
    
    Psi_joint_hat[,,kk] <- t(solve(Cov_joint_F0) %*% t(Cov_joint_F1))
    Eta_joint_hat[,,kk] <- Cov_joint_F0 - Psi_joint_hat[,,kk]%*%Cov_joint_F1
    
    if (kk <= KK){
      Cov_indiv_F0 <- t(factor_list$factor_group1[[kk]]) %*% factor_list$factor_group1[[kk]] / TT
      Cov_indiv_F1 <- t(factor_list$factor_group1[[kk]][-1,]) %*% factor_list$factor_group1[[kk]][-TT,] / TT
    }else{
      Cov_indiv_F0 <- t(factor_list$factor_group2[[(kk-KK)]]) %*% factor_list$factor_group2[[(kk-KK)]] / TT
      Cov_indiv_F1 <- t(factor_list$factor_group2[[(kk-KK)]][-1,]) %*% factor_list$factor_group2[[(kk-KK)]][-TT,] / TT
    }
    Psi_indiv_k[,,kk] <- DGP$C_tilde[[kk]] %*% Psi_indiv %*% solve(DGP$C_tilde[[kk]])
    Psi_indiv_hat[,,kk] <- t(solve(Cov_indiv_F0) %*% t(Cov_indiv_F1))
    Eta_indiv_hat[,,kk] <- Cov_indiv_F0 - Psi_indiv_hat[,,kk]%*%Cov_indiv_F1
  }
  output <- list(Psi_joint_hat = Psi_joint_hat, Psi_joint_k = Psi_joint_k, Eta_joint_hat = Eta_joint_hat,
                 Psi_indiv_hat = Psi_indiv_hat, Psi_indiv_k = Psi_indiv_k, Eta_indiv_hat = Eta_indiv_hat)
  return(output)
}



#--------------------------------------------------#
# Rotational bootstrapping function
#--------------------------------------------------#
Rotational_bootstrap <- function(dd,TT,DGP,cull=0.5){
  
  Num_sim <- 200
  
  X_full <- DGP
  svd_full <- svd(X_full)
  
  beta <- min(TT,dd)/max(TT,dd)
  upper <- (1+sqrt(beta))
  lower <- (1-sqrt(beta))
  lambda <- seq(lower^2,upper^2,length.out=1000)
  h <- na.omit(sqrt(((1+sqrt(beta))^2-lambda)*(lambda-(1-sqrt(beta))^2))/((2*pi)*beta*lambda))
  MP_beta <- median(h)
  
  sigma_hat <- median(svd_full$d)/median(h)
  single_hat <- svd_full$d/sigma_hat
  single_shrunken <- ifelse(single_hat >= upper,1/sqrt(2)*sqrt(single_hat^2-beta-1+sqrt((single_hat^2-beta-1)^2-4*beta)),0)
  single_shrunken <- sigma_hat*single_shrunken
  
  r_hat <- max(length(single_shrunken[which(single_shrunken>0)]),1)
  A_hat <- svd_full$u[,1:r_hat] %*% diag(single_shrunken[1:r_hat],r_hat) %*% t(svd_full$v[,1:r_hat])
  E_hat <- X_full - A_hat
  E_hat_remain <- svd_full$u[,(r_hat+1):min(dd,TT)] %*% diag(svd_full$d[(r_hat+1):min(dd,TT)],length((r_hat+1):min(dd,TT))) %*% t(svd_full$v[,(r_hat+1):min(dd,TT)])
  
  Imputed_single <- vector("numeric",length=r_hat)
  MP_cumm <- cumsum(h)/sum(h)
  for (i in 1:r_hat){
    percent <- sample(seq(1,1000,1),1,FALSE)
    MP_rand <- MP_cumm[percent]
    Imputed_single[i] <- sqrt(MP_rand)
  }
  E_hat_imputed <- E_hat_remain + svd_full$u[,1:r_hat] %*% diag(Imputed_single*sigma_hat,length(1:r_hat)) %*% t(svd_full$v[,1:r_hat])
  E_hat_imputed[is.nan(E_hat_imputed)] <- 0
  
  angle_horizontal <- vector("numeric",length=1000)
  angle_vertical <- vector("numeric",length=1000)
  for (iter in 1:1000){
    vec_horizontal <- rnorm(dd,0,1)
    vec_vertical <- rnorm(TT,0,1)
    angle_horizontal[iter] <- acos(norm(vec_horizontal[1:r_hat],"2")/norm(vec_horizontal,"2"))*180/pi
    angle_vertical[iter] <- acos(norm(vec_vertical[1:r_hat],"2")/norm(vec_vertical,"2"))*180/pi
  }
  rangle_horizontal <- quantile(angle_horizontal,0.05)
  rangle_vertical <- quantile(angle_vertical,0.05)
  
  
  PC_angle_U <- array(NA,dim=c(Num_sim,r_hat))
  PC_angle_V <- array(NA,dim=c(Num_sim,r_hat))
  for (i in 1:Num_sim){
    
    Rand_U <- matrix(rnorm((TT*r_hat),0,1),nrow=TT)
    Rand_U <- (diag(1,TT) - 1/TT*matrix(rep(1,TT^2),nrow=TT)) %*% Rand_U
    Rand_U <- pracma::orth(Rand_U)
    
    Rand_V <- matrix(rnorm((dd*r_hat),0,1),nrow=dd)
    # # Our case is column-centered
    Rand_v <- (diag(1,dd) - 1/dd*matrix(rep(1,dd^2),nrow=dd)) %*% Rand_V
    Rand_V <- pracma::orth(Rand_V)
    
    Rand_X <- Rand_U %*% diag(single_shrunken[1:r_hat],r_hat) %*% t(Rand_V) + E_hat_imputed
    svd_rand <- svd(Rand_X,nu=r_hat,nv=r_hat)
    for (j in 1:r_hat){
      PC_angle_U[i,j] <- acos(min(svd(t(Rand_U) %*% svd_rand$u[,1:j])$d))*180/pi
      PC_angle_V[i,j] <- acos(min(svd(t(Rand_V) %*% svd_rand$v[,1:j])$d))*180/pi
    }
  }
  PC_angle_U[is.nan(PC_angle_U)] <- 0
  PC_angle_V[is.nan(PC_angle_V)] <- 0
  
  for (j in 1:r_hat){
    PC_angle_U[,j] <- sort(PC_angle_U[,j])
    PC_angle_V[,j] <- sort(PC_angle_V[,j])
  }
  
  r_check <- min(sum(PC_angle_U[Num_sim*0.95,] < cull*rangle_horizontal),
                 sum(PC_angle_V[Num_sim*0.95,] < cull*rangle_vertical))
  return(r_check)
}
