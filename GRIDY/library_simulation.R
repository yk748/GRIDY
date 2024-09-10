# -----------------------------------------------------------------------------#
# Running simulation
# -----------------------------------------------------------------------------#
# Packages required
library(mvtnorm)
library(Matrix)
library(combinat)
library(multiway) # For SCA and GICA
library(ica) # For SCA and GICA
library(readxl)
library(ggplot2)
library(latex2exp)
library(ggpubr)
library(RColorBrewer)
library(reshape2)
library(gridExtra)

#--------------------------------------------------#
#   Function name : sim_model    												  
#
#   Project : "Group integrative dynamic factor models 
#             with application to multiple subject brain connectivity"
#
#   Maintainer : Younghoon Kim                     
#
#   Date : Sep. 1st, 2024
#
#   Purpose : Data generating function used in Section 4. Illustrative examples
#
#   R version 4.0.5 (2021-03-31)                                    
#
#   Input:    
#    @ Joint_rank : number of joint factors.
#    @ Gr_Ind_rank : number of group individual factors.
#    @ Num_var : number of variables.
#    @ Num_sample : number of sample lengths.
#    @ Num_subj_per_group : number of subjects in each group.
#    @ Control_param : a set of parameters that determines the model. 
#                      It includes type (1 and 2 provided. 1 stands for correlated factor series
#                      and 2 stands for independent factors), 
#                      sigma_xi (standard deviation of noise of factor models), 
#                      J_per (percentage of rows belonging to the common components), 
#                      b_bar_min (minimal value of joint loadings), 
#                      b_bar_max (maximal value of joint loadings), 
#                      b_tilde_min (minimal value of group individual loadings), 
#                      b_tilde_max (maximal value of group individual loadings), 
#                      sigma_c (standard deviation of factor series), 
#                      sigma_eps (standard deviation of noise of observations).
#                      Note that to make the factor series stable, 
#                      a pre-calculated combinations of type and sigma_xi are required and 
#                      this pair matches with the pre-calculated VAR transition matrices 
#                      of the factor series. See Section 4.1 for the detail.
# 
#   Output:
#    @ Xk : list of Num_sample x Num_sample generated observational multi-view blocks.
#    @ X_bar : list of Num_sample x Num_sample generated joint components.
#    @ X_tilde : list of Num_sample x Num_sample generated group individual components.
#    @ Ek : list of Num_sample x Num_sample observation noises.
#    @ scaled_X : centered and scaled Xk.
#    @ scaled_X_comm : centered and scaled X_bar.
#    @ scaled_X_G_ind : centered and scaled X_tilde.
#    @ F_bar : joint factor series.
#    @ F_tilde : group individual factor series.
#    @ A_bar : scaled joint factor series.
#    @ A_tilde : scaled group individual factor series.
#    @ C_bar : scalers of joint factor series.
#    @ C_tilde : scalers of group individual factor series.
#    @ Psi_bar : VAR transition matrices of joint factor series.
#    @ Psi_tilde : VAR transition matrices of group individual factor series.
#    @ B_bar : joint loadings matrix
#    @ B_tilde1 : group individual loadings matrix belonging to the first group.
#    @ B_tilde2 : group individual loadings matrix belonging to the second group.
#    @ comm_idx : index of variables belonging to joint components.
#    @ group1_idx : index of variables belonging to first group individual components.
#    @ group2_idx : index of variables belonging to second group individual components.
#
#   Required R packages : mvtnorm_1.2-5 and Matrix_1.5-1.
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
#   Function name : factor_regression    												  
#
#   Project : "Group integrative dynamic factor models 
#             with application to multiple subject brain connectivity"
#
#   Maintainer : Younghoon Kim                     
#
#   Date : Sep. 1st, 2024
#
#   Purpose : subject-wise regression used in estimation of factor series. 
#             See (18) in Section 3.3 for the detail.
#
#   R version 4.0.5 (2021-03-31)                                     
#
#   Input:    
#    @ KK : number of subjects.
#    @ r_J : number of joint factor series.
#    @ r_G : number of group individual factor series.
#    @ DGP : generated data. Observation from sim_model() is called.
#    @ Series_joint : joint components from previous step. 
#                     Estimates of joint loadings are called. 
#    @ Series_group1 : group individual components of the first group from previous step.
#                      Estimates of group individual loadings of the group 1 are called.
#    @ Series_group2 : group individual components of the second group from previous step.
#                      Estimates of group individual loadings of the group 2 are called.
# 
#   Output:
#    @ factor_joint : estimated joint factor series.
#    @ factor_group1 : estimated group individual factor series of group 1.
#    @ factor_group2 : estimated group individual factor series of group 2.
#
#   Required R packages : Matrix_1.5-1.
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
#   Function name : YW_compute   												  
#
#   Project : "Group integrative dynamic factor models 
#             with application to multiple subject brain connectivity"
#
#   Maintainer : Younghoon Kim                     
#
#   Date : Sep. 1st, 2024
#
#   Purpose : compute Yule-Walker equation for estimating 
#             vector autoregressive (VAR) transition matrices and 
#             covariance matrices of noise for all subjects 
#             by taking the output of factor_regression. 
#
#   R version 4.0.5 (2021-03-31)                                     
#
#   Input:    
#    @ dd : number of variables.
#    @ TT : sample lengths.
#    @ KK : number of subjects.
#    @ r_J : number of joint factors.
#    @ r_G : number of group individual factors.
#    @ factor_list: list of estimated factor series.
# 
#   Output:
#    @ Psi_joint_hat : VAR transition matrices of estimated joint factor series.
#    @ Psi_indiv_hat: VAR transition matrices of estimated group individual factor series.
#    @ Eta_joint_hat : covariance matrices of noises of estimated joint factor series.
#    @ Eta_G_indiv_hat : covariance matrices of noises of estimated 
#                        group individual factor series.
#
#   Required R packages : Matrix_1.5-1.
#--------------------------------------------------#
YW_compute <- function(dd,TT,KK,r_J,r_G,factor_list){
  
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


