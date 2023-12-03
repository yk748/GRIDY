#--------------------------------------------------#
# Factor_regression function
#--------------------------------------------------#
factor_regression <- function(K1_idx,K2_idx,r_J,r_I_1,r_I_2,
                              Series_observation,Series_joint,Series_group1,Series_group2){
  
  factor_regression_joint <- list()
  noise_refitted <- list()
  if (r_I_1 != 0){
    factor_regression_group1 <- list()
  }else{
    factor_regression_group1 <- 0
  }
  if (r_I_2 != 0){
    factor_regression_group2 <- list()
  }else{
    factor_regression_group2 <- 0
  }
  
  
  for (k in 1:(length(K1_idx) + length(K2_idx))){
    if (k <= length(K1_idx)){
      if (r_I_1 != 0){
        B_aug <- cbind(Series_joint$B,Series_group1$B) 
      }else{
        B_aug <- cbind(Series_joint$B)
      }
    }else{
      if (r_I_2 != 0){
        B_aug <- cbind(Series_joint$B,Series_group2$B) 
      }else{
        B_aug <- cbind(Series_joint$B)
      }
    }
    factor_tmp <- ( t(Series_observation[[k]]) %*% B_aug ) %*% solve (t(B_aug) %*% B_aug)
    factor_regression_joint[[k]] <- as.matrix(factor_tmp[,1:r_J])
    
    if (k <= length(K1_idx)){
      if (r_I_1 != 0){
        factor_regression_group1[[k]] <- as.matrix(factor_tmp[,(r_J+1):(r_J+r_I_1)])
      }
    }else{
      if (r_I_2 != 0){
        factor_regression_group2[[(k-length(K1_idx))]] <- as.matrix(factor_tmp[,(r_J+1):(r_J+r_I_2)])
      }
    }
    
    noise_refitted[[k]] <- t(Series_observation[[k]]) - as.matrix(factor_tmp) %*% t(B_aug)

  }
  output <- list(factor_joint = factor_regression_joint,
                 factor_group1 = factor_regression_group1,
                 factor_group2 = factor_regression_group2,
                 noise_refitted = noise_refitted)
  return(output)
}

#--------------------------------------------------#
# Yule-Walker function
#--------------------------------------------------#
YK_compute <- function(K1_idx,K2_idx,r_J,r_I_1,r_I_2,Series_observation,factor_list){
  
  Psi_joint_hat <- array(NA,dim=c(r_J,r_J,(length(K1_idx) + length(K2_idx))))
  Eta_joint_hat <- array(NA,dim=c(r_J,r_J,(length(K1_idx) + length(K2_idx))))
  
  if (r_I_1 != 0){
    Psi_indiv1_hat <- array(NA,dim=c(r_I_1,r_I_1,length(K1_idx)))
    Eta_indiv1_hat <- array(NA,dim=c(r_I_1,r_I_1,length(K1_idx)))
  }else{
    Psi_indiv1_hat <- 0
    Eta_indiv1_hat <- 0
  }

  if (r_I_2 != 0){
    Psi_indiv2_hat <- array(NA,dim=c(r_I_2,r_I_2,length(K2_idx)))
    Eta_indiv2_hat <- array(NA,dim=c(r_I_2,r_I_2,length(K2_idx)))
  }else{
    Psi_indiv2_hat <- 0
    Eta_indiv2_hat <- 0
  }
  
  for (k in 1:(length(K1_idx) + length(K2_idx))){
    
    T_k <- dim(X_scale_final[[k]])[2]
    
    Cov_joint_F0 <- t(factor_list$factor_joint[[k]]) %*% factor_list$factor_joint[[k]] / T_k
    Cov_joint_F1 <- t(factor_list$factor_joint[[k]][-1,]) %*% factor_list$factor_joint[[k]][-T_k,] / T_k
    
    Psi_joint_hat[,,k] <- t(solve(Cov_joint_F0) %*% t(Cov_joint_F1))
    Eta_joint_hat[,,k] <- Cov_joint_F0 - Psi_joint_hat[,,k]%*%t(Cov_joint_F1)
    
    if (k <= length(K1_idx)){
      if (r_I_1 != 0){
        Cov_indiv1_F0 <- t(factor_list$factor_group1[[k]]) %*% factor_list$factor_group1[[k]] / T_k
        Cov_indiv1_F1 <- t(factor_list$factor_group1[[k]][-1,]) %*% factor_list$factor_group1[[k]][-T_k,] / T_k
        
        Psi_indiv1_hat[,,k] <- t(solve(Cov_indiv1_F0) %*% t(Cov_indiv1_F1))
        Eta_indiv1_hat[,,k] <- Cov_indiv1_F0 - Psi_indiv1_hat[,,k]%*%t(Cov_indiv1_F1)
      }
    }else{
      if (r_I_2 != 0){
        Cov_indiv2_F0 <- t(factor_list$factor_group2[[(k-length(K1_idx))]]) %*% factor_list$factor_group2[[(k-length(K1_idx))]] / T_k
        Cov_indiv2_F1 <- t(factor_list$factor_group2[[(k-length(K1_idx))]][-1,]) %*% factor_list$factor_group2[[(k-length(K1_idx))]][-T_k,] / T_k
        
        Psi_indiv2_hat[,,(k-length(K1_idx))] <- t(solve(Cov_indiv2_F0) %*% t(Cov_indiv2_F1))
        Eta_indiv2_hat[,,(k-length(K1_idx))] <- Cov_indiv2_F0 - Psi_indiv2_hat[,,(k-length(K1_idx))]%*%t(Cov_indiv2_F1)
      }
    }
  }
  output <- list(Psi_joint_hat = Psi_joint_hat, Eta_joint_hat = Eta_joint_hat,
                 Psi_indiv1_hat = Psi_indiv1_hat, Eta_indiv1_hat = Eta_indiv1_hat,
                 Psi_indiv2_hat = Psi_indiv2_hat, Eta_indiv2_hat = Eta_indiv2_hat)
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
    # Rand_U <- (diag(1,TT) - 1/TT*matrix(rep(1,TT^2),nrow=TT)) %*% Rand_U
    Rand_U <- pracma::orth(Rand_U)
    
    Rand_V <- matrix(rnorm((dd*r_hat),0,1),nrow=dd)
    # # Our case is column-centered
    # Rand_v <- (diag(1,dd) - 1/dd*matrix(rep(1,dd^2),nrow=dd)) %*% Rand_V
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
