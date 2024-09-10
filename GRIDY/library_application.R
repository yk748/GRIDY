# -----------------------------------------------------------------------------#
# Running data application
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
#             The input of this function is slightly different from that of the function in 
#             the same name under library_simulation since this function can handle 
#             unequal numbers of subjects in groups and contains the exceptional handling 
#             that can be encountered in real data applications. But their roles are the same.
#
#   R version 4.0.5 (2021-03-31)                                    
#
#   Input:    
#    @ K1_idx : subject indices belonging to group 1.
#    @ K2_idx : subject indices belonging to group 2.
#    @ r_J : number of joint factor series.
#    @ r_G1 : number of group individual factor series in group 1.
#    @ r_G2 : number of group individual factor series in group 2.
#    @ Series_observation : Observed data. Observation X is called.
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
#    @ noise_refitted : estimated residuals after fitting factor series.
#
#   Required R packages : Matrix_1.5-1.
#--------------------------------------------------#
factor_regression <- function(K1_idx,K2_idx,r_J,r_G1,r_G2,
                              Series_observation,Series_joint,Series_group1,Series_group2){
  
  factor_regression_joint <- list()
  noise_refitted <- list()
  if (r_G1 != 0){
    factor_regression_group1 <- list()
  }else{
    factor_regression_group1 <- 0
  }
  if (r_G2 != 0){
    factor_regression_group2 <- list()
  }else{
    factor_regression_group2 <- 0
  }
  
  
  for (kk in 1:(length(K1_idx) + length(K2_idx))){
    if (kk <= length(K1_idx)){
      if (r_G1 != 0){
        B_aug <- cbind(Series_joint$B,Series_group1$B) 
      }else{
        B_aug <- cbind(Series_joint$B)
      }
    }else{
      if (r_G2 != 0){
        B_aug <- cbind(Series_joint$B,Series_group2$B) 
      }else{
        B_aug <- cbind(Series_joint$B)
      }
    }
    factor_tmp <- ( t(Series_observation[[kk]]) %*% B_aug ) %*% solve (t(B_aug) %*% B_aug)
    factor_regression_joint[[kk]] <- as.matrix(factor_tmp[,1:r_J])
    
    if (kk <= length(K1_idx)){
      if (r_G1 != 0){
        factor_regression_group1[[kk]] <- as.matrix(factor_tmp[,(r_J+1):(r_J+r_G1)])
      }
    }else{
      if (r_G2 != 0){
        factor_regression_group2[[(kk-length(K1_idx))]] <- as.matrix(factor_tmp[,(r_J+1):(r_J+r_G2)])
      }
    }
    
    noise_refitted[[kk]] <- t(Series_observation[[kk]]) - as.matrix(factor_tmp) %*% t(B_aug)
    
  }
  output <- list(factor_joint = factor_regression_joint,
                 factor_group1 = factor_regression_group1,
                 factor_group2 = factor_regression_group2,
                 noise_refitted = noise_refitted)
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
#             The input of this function is slightly different from that of the function in 
#             the same name under library_simulation since this function can handle 
#             unequal numbers of subjects in groups and contains the exceptional handling 
#             that can be encountered in real data applications. But their roles are the same.
#
#   R version 4.0.5 (2021-03-31)                                    
#
#   Input:    
#    @ K1_idx : subject indices belonging to group 1.
#    @ K2_idx : subject indices belonging to group 2.
#    @ r_J : number of joint factors.
#    @ r_G1 : number of group individual factor series in group 1.
#    @ r_G2 : number of group individual factor series in group 2.
#    @ Series_observation : Observed data. Observation X is called.
#    @ factor_list: list of estimated factor series.
# 
#   Output:
#    @ Psi_joint_hat : VAR transition matrices of estimated joint factor series.
#    @ Eta_joint_hat : covariance matrices of noises of estimated joint factor series
#    @ Psi_indiv1_hat : VAR transition matrices of estimated group individual 
#                       factor series in group 1.
#    @ Eta_indiv1_hat : covariance matrices of noises of estimated 
#                        group individual factor series in group 1.
#    @ Psi_indiv2_hat : VAR transition matrices of estimated group individual 
#                       factor series in group 2.
#    @ Eta_indiv2_hat : covariance matrices of noises of estimated 
#                        group individual factor series in group 2.
# 
#   Required R packages : Matrix_1.5-1.
#--------------------------------------------------#
YW_compute <- function(K1_idx,K2_idx,r_J,r_G1,r_G2,Series_observation,factor_list){
  
  Psi_joint_hat <- array(NA,dim=c(r_J,r_J,(length(K1_idx) + length(K2_idx))))
  Eta_joint_hat <- array(NA,dim=c(r_J,r_J,(length(K1_idx) + length(K2_idx))))
  
  if (r_G1 != 0){
    Psi_indiv1_hat <- array(NA,dim=c(r_G1,r_G1,length(K1_idx)))
    Eta_indiv1_hat <- array(NA,dim=c(r_G1,r_G1,length(K1_idx)))
  }else{
    Psi_indiv1_hat <- 0
    Eta_indiv1_hat <- 0
  }
  
  if (r_G2 != 0){
    Psi_indiv2_hat <- array(NA,dim=c(r_G2,r_G2,length(K2_idx)))
    Eta_indiv2_hat <- array(NA,dim=c(r_G2,r_G2,length(K2_idx)))
  }else{
    Psi_indiv2_hat <- 0
    Eta_indiv2_hat <- 0
  }
  
  for (kk in 1:(length(K1_idx) + length(K2_idx))){
    
    T_k <- dim(Series_observation[[kk]])[2]
    
    Cov_joint_F0 <- t(scale(factor_list$factor_joint[[kk]],center=TRUE,scale=FALSE)) %*% scale(factor_list$factor_joint[[kk]],center=TRUE,scale=FALSE) / T_k
    Cov_joint_F1 <- t(scale(factor_list$factor_joint[[kk]][-1,],center=TRUE,scale=FALSE)) %*% scale(factor_list$factor_joint[[kk]][-T_k,],center=TRUE,scale=FALSE) / T_k
    
    Psi_joint_hat[,,kk] <- t(solve(Cov_joint_F0) %*% t(Cov_joint_F1))
    Eta_joint_hat[,,kk] <- Cov_joint_F0 - Psi_joint_hat[,,kk]%*%t(Cov_joint_F1)
    
    if (kk <= length(K1_idx)){
      if (r_G1 != 0){
        Cov_indiv1_F0 <- t(scale(factor_list$factor_group1[[kk]],center=TRUE,scale=FALSE)) %*% scale(factor_list$factor_group1[[kk]],center=TRUE,scale=FALSE) / T_k
        Cov_indiv1_F1 <- t(scale(factor_list$factor_group1[[kk]][-1,],center=TRUE,scale=FALSE)) %*% scale(factor_list$factor_group1[[kk]][-T_k,],center=TRUE,scale=FALSE) / T_k
        
        Psi_indiv1_hat[,,kk] <- t(solve(Cov_indiv1_F0) %*% t(Cov_indiv1_F1))
        Eta_indiv1_hat[,,kk] <- Cov_indiv1_F0 - Psi_indiv1_hat[,,kk]%*%t(Cov_indiv1_F1)
      }
    }else{
      if (r_G2 != 0){
        Cov_indiv2_F0 <- t(scale(factor_list$factor_group2[[(kk-length(K1_idx))]],center=TRUE,scale=FALSE)) %*% scale(factor_list$factor_group2[[(kk-length(K1_idx))]],center=TRUE,scale=FALSE) / T_k
        Cov_indiv2_F1 <- t(scale(factor_list$factor_group2[[(kk-length(K1_idx))]][-1,],center=TRUE,scale=FALSE)) %*% scale(factor_list$factor_group2[[(kk-length(K1_idx))]][-T_k,],center=TRUE,scale=FALSE) / T_k
        
        Psi_indiv2_hat[,,(kk-length(K1_idx))] <- t(solve(Cov_indiv2_F0) %*% t(Cov_indiv2_F1))
        Eta_indiv2_hat[,,(kk-length(K1_idx))] <- Cov_indiv2_F0 - Psi_indiv2_hat[,,(kk-length(K1_idx))]%*%t(Cov_indiv2_F1)
      }
    }
  }
  output <- list(Psi_joint_hat = Psi_joint_hat, Eta_joint_hat = Eta_joint_hat,
                 Psi_indiv1_hat = Psi_indiv1_hat, Eta_indiv1_hat = Eta_indiv1_hat,
                 Psi_indiv2_hat = Psi_indiv2_hat, Eta_indiv2_hat = Eta_indiv2_hat)
  return(output)
}


