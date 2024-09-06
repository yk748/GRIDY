#--------------------------------------------------#
#   Function name : rotational_boostrap    												  
#
#   Originally developed by Prothero et al. (2024).                    
#   
#   citation: Prothero, J., Jiang, M., Hannig, J., Tran-Dinh, Q., 
#             Ackerman, A., and Marron, J.S. (2024) 
#             Data integration via analysis of subspaces (DIVAS). TEST, 14, 1-42.
#
#   Last visit date : Sep. 1st, 2024
#
#   Purpose : this function is created to implement 
#             principal angle-based estimation methods for initial rank 
#             of multi-view data by using bootstrap 
#             The reference of this R implementation is 
#             Prothero, J., Jiang, M., Hannig, J., Tran-Dinh, Q., 
#             Ackerman, A., and Marron, J.S. (2024) 
#             Data integration via analysis of subspaces (DIVAS). TEST, 14, 1-42.
#
#   Note: to the best knowledge, the public accessible code of this function does not exist. 
#         Therefore, Younghoon Kim, the author of "Group integrative dynamic factor models 
#         with application to multiple subject brain connectivity"
#         wrote this function for his research.
#         This function is a part of what the authors of DIVAS has been implemented; 
#         the segmentation of multi-view data is not included in this version. 
#
#   R version 4.0.5 (2021-03-31)                                      
#
#   Input:    
#    @ dd : number of variables in the multi-view data.
#    @ TT : sample length of the multi-view data.
#    @ DGP : TT x dd - dimensional block observations.
#    @ cull : multiplied constant to the threshold. Default is 0.5. 
#
#   Output:
#    @ r check: estimated initial rank of a single multi-view data.
#
#   Required R packages : mvtnorm_1.2-5 and Matrix_1.5-1.
#--------------------------------------------------#
Rotational_bootstrap <- function(dd,TT,DGP,cull=0.5){
  
  Num_sim <- 200
  
  #--------------------------------------------------#
  # Obtain maximal initial rank of the multi-view data by using optimal hard threshold.
  #--------------------------------------------------#
  X_full <- DGP
  svd_full <- svd(X_full)
  
  # Compute Marchenko–Pastur distribution
  beta <- min(TT,dd)/max(TT,dd)
  upper <- (1+sqrt(beta))
  lower <- (1-sqrt(beta))
  lambda <- seq(lower^2,upper^2,length.out=1000)
  h <- na.omit(sqrt(((1+sqrt(beta))^2-lambda)*(lambda-(1-sqrt(beta))^2))/((2*pi)*beta*lambda))
  MP_beta <- median(h)
  
  # Compute singular values up to the maximal numbers
  sigma_hat <- median(svd_full$d)/median(h)
  single_hat <- svd_full$d/sigma_hat
  single_shrunken <- vector("numeric",length(single_hat))
  for (j in 1:length(single_hat)){
    if (single_hat[j] >= upper){
      single_shrunken[j] <- sigma_hat*1/sqrt(2)*sqrt(single_hat[j]^2-beta-1+sqrt((single_hat[j]^2-beta-1)^2-4*beta))
    }else{
      next
    }
  }
  
  #--------------------------------------------------#
  # Compute imputed noise matrix for boostrapping
  #--------------------------------------------------#
  # Copmute residual
  r_hat <- max(length(single_shrunken[which(single_shrunken>0)]),1)
  A_hat <- svd_full$u[,1:r_hat] %*% diag(single_shrunken[1:r_hat],r_hat) %*% t(svd_full$v[,1:r_hat])
  E_hat <- X_full - A_hat
  E_hat_remain <- svd_full$u[,(r_hat+1):min(dd,TT)] %*% diag(svd_full$d[(r_hat+1):min(dd,TT)],length((r_hat+1):min(dd,TT))) %*% t(svd_full$v[,(r_hat+1):min(dd,TT)])
  
  # Compute imputed energy
  Imputed_single <- vector("numeric",length=r_hat)
  MP_cumm <- cumsum(h)/sum(h)
  for (i in 1:r_hat){
    percent <- sample(seq(1,1000,1),1,FALSE)
    MP_rand <- MP_cumm[percent]
    Imputed_single[i] <- sqrt(MP_rand)
  }
  E_hat_imputed <- E_hat_remain + svd_full$u[,1:r_hat] %*% diag(Imputed_single*sigma_hat,length(1:r_hat)) %*% t(svd_full$v[,1:r_hat])
  E_hat_imputed[is.nan(E_hat_imputed)] <- 0
  
  #--------------------------------------------------#
  # Compute random directional bounds for row and column spaces as thresholds
  #--------------------------------------------------#
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
  
  #--------------------------------------------------#
  # Bootstrap the block signals by generating random orthogonal row and column matrices.
  #--------------------------------------------------#
  PC_angle_U <- array(NA,dim=c(Num_sim,r_hat))
  PC_angle_V <- array(NA,dim=c(Num_sim,r_hat))
  for (i in 1:Num_sim){
    
    # row-centered
    Rand_U <- matrix(rnorm((TT*r_hat),0,1),nrow=TT)
    Rand_U <- (diag(1,TT) - 1/TT*matrix(rep(1,TT^2),nrow=TT)) %*% Rand_U
    Rand_U <- pracma::orth(Rand_U)
    
    # column-centered
    Rand_V <- matrix(rnorm((dd*r_hat),0,1),nrow=dd)
    Rand_v <- (diag(1,dd) - 1/dd*matrix(rep(1,dd^2),nrow=dd)) %*% Rand_V
    Rand_V <- pracma::orth(Rand_V)
    
    #--------------------------------------------------#
    # Compute principal angles between generated block signals and their approximations
    #--------------------------------------------------#
    Rand_X <- Rand_U %*% diag(single_shrunken[1:r_hat],r_hat) %*% t(Rand_V) + E_hat_imputed
    svd_rand <- svd(Rand_X,nu=r_hat,nv=r_hat)
    for (j in 1:r_hat){
      PC_angle_U[i,j] <- acos(min(svd(t(Rand_U) %*% svd_rand$u[,1:j])$d))*180/pi
      PC_angle_V[i,j] <- acos(min(svd(t(Rand_V) %*% svd_rand$v[,1:j])$d))*180/pi
    }
  }
  PC_angle_U[is.nan(PC_angle_U)] <- 0
  PC_angle_V[is.nan(PC_angle_V)] <- 0
  
  # sort computed angles
  for (j in 1:r_hat){
    PC_angle_U[,j] <- sort(PC_angle_U[,j])
    PC_angle_V[,j] <- sort(PC_angle_V[,j])
  }
  
  #--------------------------------------------------#
  # Return the count of singular values measured by the angles above.
  #--------------------------------------------------#
  r_check <- min(sum(PC_angle_U[Num_sim*0.95,] < cull*rangle_horizontal),
                 sum(PC_angle_V[Num_sim*0.95,] < cull*rangle_vertical))
  return(r_check)
}