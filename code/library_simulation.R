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
# Auxiliary functions:
#--------------------------------------------------#
# get_wedin_bound_samples
get_wedin_bound_samples <- function(X, SVD, signal_rank, num_samples=1000){
  
  # resample for U and V
  U_perp <- SVD[['u']][ , -(1:signal_rank)]
  U_sampled_norms <- wedin_bound_resampling(X=X,
                                            perp_basis=U_perp,
                                            right_vectors=FALSE,
                                            num_samples=num_samples)
  
  V_perp <- SVD[['v']][ , -(1:signal_rank)]
  V_sampled_norms <- wedin_bound_resampling(X=X,
                                            perp_basis=V_perp,
                                            right_vectors=TRUE,
                                            num_samples=num_samples)
  
  sigma_min <- SVD[['d']][signal_rank]
  wedin_bound_samples <- mapply(function(u, v)  min(max(u, v)/sigma_min, 1)^2, U_sampled_norms, V_sampled_norms)
  
  wedin_bound_samples
}

# wedin_bound_resampling
wedin_bound_resampling <- function(X, perp_basis, right_vectors, num_samples=1000){
  
  rank <- dim(perp_basis)[2]
  resampled_norms <- rep(0, num_samples)
  
  for(s in 1:num_samples){
    
    sampled_col_index <- sample.int(n=dim(perp_basis)[2],
                                    size=rank,
                                    replace=TRUE)
    
    
    perp_resampled <- perp_basis[ , sampled_col_index]
    
    if(right_vectors){
      resampled_projection <- X %*% perp_resampled
    } else{
      resampled_projection <- t(perp_resampled) %*% X
    }
    
    # operator L2 norm
    resampled_norms[s] <- norm(resampled_projection,
                               type='2')
  }
  
  resampled_norms
}

# get_random_direction_bound
get_random_direction_bound <- function(n_obs, dims, num_samples=1000){
  
  n_blocks <- length(dims)
  rand_dir_samples <- rep(0, num_samples)
  for(s in 1:num_samples){
    rand_subspaces <- list()
    for(b in 1:n_blocks){
      X <- matrix(rnorm(n_obs * dims[b], mean=0,sd=1), n_obs, dims[b])
      U <- get_svd(X)[['u']]
      
      rand_subspaces[[b]] <- U
      
    }
    M <- do.call(cbind, rand_subspaces)
    M_svd <- get_svd(M, rank=min(dims))
    
    rand_dir_samples[s] <- M_svd[['d']][1]^2
    
  }
  
  rand_dir_samples
}

# get_svd
get_svd <- function(X, rank=NULL){
  # SVD <- get_svd(X, rank=2)
  
  if(is.null(rank)){
    svd(X)
  } else if(rank == 0){
    # TODO: what to do
    decomposition <- list()
    decomposition[['u']] <- matrix(0, ncol=1, nrow=dim(X)[1])
    decomposition[['d']] <- 0
    decomposition[['v']] <- matrix(0, ncol=1, nrow=dim(X)[2])
    decomposition
    
  } else{
    decomposition <- svd(X, nu=rank, nv=rank)
    decomposition[['d']] <- decomposition[['d']][1:rank]
    decomposition
  }
  
}

# get_sv_threshold
get_sv_threshold <- function(singular_values, rank){
  0.5*(singular_values[rank] + singular_values[rank + 1])
}

# get_joint_scores
get_joint_scores <- function(blocks, block_svd, initial_signal_ranks, sv_thresholds,
                             n_wedin_samples=1000, n_rand_dir_samples=1000,
                             joint_rank=NA){
  
  
  if(is.na(n_wedin_samples) & is.na(n_rand_dir_samples) & is.na(joint_rank)){
    stop('at least one of n_wedin_samples, n_rand_dir_samples, or joint_rank must not be NA',
         call.=FALSE)
  }
  
  K <- length(blocks)
  n_obs <- dim(blocks[[1]])[1]
  
  # SVD of the signal scores matrix -----------------------------------------
  signal_scores <- list()
  for(k in 1:K){
    signal_scores[[k]] <- block_svd[[k]][['u']][, 1:initial_signal_ranks[k]]
  }
  
  M <- do.call(cbind, signal_scores)
  M_svd <- get_svd(M, rank=min(initial_signal_ranks))
  
  
  # estimate joint rank with wedin bound and random direction bound-------------------------------------------------------------
  
  rank_sel_results  <- list()
  rank_sel_results[['obs_svals']] <- M_svd[['d']]
  
  if(is.na(joint_rank)){
    
    # maybe comptue wedin bound
    if(!is.na(n_wedin_samples)){
      
      block_wedin_samples <- matrix(NA, K, n_wedin_samples)
      
      for(k in 1:K){
        block_wedin_samples[k, ] <- get_wedin_bound_samples(X=blocks[[k]],
                                                            SVD=block_svd[[k]],
                                                            signal_rank=initial_signal_ranks[k],
                                                            num_samples=n_wedin_samples)
      }
      
      wedin_samples <-  K - colSums(block_wedin_samples)
      wedin_svsq_threshold <- quantile(wedin_samples, .05)
      
      rank_sel_results[['wedin']] <- list(block_wedin_samples=block_wedin_samples,
                                          wedin_samples=wedin_samples,
                                          wedin_svsq_threshold=wedin_svsq_threshold)
    } else{
      wedin_svsq_threshold <- NA
    }
    
    # maybe compute random direction bound
    if(!is.na(n_rand_dir_samples)){
      
      rand_dir_samples <- get_random_direction_bound(n_obs=n_obs, dims=initial_signal_ranks, num_samples=n_rand_dir_samples)
      rand_dir_svsq_threshold <- quantile(rand_dir_samples, .95)
      
      rank_sel_results[['rand_dir']] <- list(rand_dir_samples=rand_dir_samples,
                                             rand_dir_svsq_threshold=rand_dir_svsq_threshold)
      
    } else {
      rand_dir_svsq_threshold <- NA
    }
    
    overall_sv_sq_threshold <- max(wedin_svsq_threshold, rand_dir_svsq_threshold, na.rm=TRUE)
    joint_rank_estimate <- sum(M_svd[['d']]^2 > overall_sv_sq_threshold)
    
    rank_sel_results[['overall_sv_sq_threshold']] <- overall_sv_sq_threshold
    rank_sel_results[['joint_rank_estimate']] <- joint_rank_estimate
    
    
  } else { # user provided joint rank
    joint_rank_estimate <- joint_rank
    rank_sel_results[['joint_rank_estimate']] <- joint_rank
  }
  
  
  # estimate joint score space ------------------------------------
  
  if(joint_rank_estimate >= 1){
    joint_scores <- M_svd[['u']][ , 1:joint_rank_estimate, drop=FALSE]
    
    # reconsider joint score space ------------------------------------
    # remove columns of joint_scores that have a
    # trivial projection from one of the data matrices
    
    to_remove <- c()
    for(k in 1:K){
      for(j in 1:joint_rank_estimate){
        
        score <- t(blocks[[k]]) %*% joint_scores[ , j]
        sv <- norm(score)
        
        if(sv < sv_thresholds[[k]]){
          print(paste('removing column', j))
          to_remove <- c(to_remove, j)
          break
        }
      }
      
    }
    to_keep <- setdiff(1:joint_rank_estimate, to_remove)
    joint_rank <- length(to_keep)
    joint_scores <- joint_scores[ , to_keep, drop=FALSE]
  } else {
    joint_scores <- NA
  }
  
  
  list(joint_scores=joint_scores, rank_sel_results=rank_sel_results)
}

# truncate_svd
truncate_svd <- function(decomposition, rank){
  
  if(rank==0){
    n <- dim(decomposition[['u']])[1]
    d <- dim(decomposition[['v']])[1]
    decomposition[['u']] <- matrix(0, ncol=1, nrow=n)
    decomposition[['d']] <- 0
    decomposition[['v']] <- matrix(0, ncol=1, nrow=d)
  }else{
    decomposition[['u']] <- decomposition[['u']][, 1:rank, drop=FALSE]
    decomposition[['d']] <- decomposition[['d']][1:rank]
    decomposition[['v']] <- decomposition[['v']][, 1:rank, drop=FALSE]
  }
  
  decomposition
}

# svd_reconstruction
svd_reconstruction <- function(decomposition){
  
  # decomposition rank -- need to truncated singluar values
  r <- dim(decomposition[['u']])[2]
  
  decomposition[['u']]  %*%
    diag(decomposition[['d']][1:r], nrow=r, ncol=r) %*%
    t(decomposition[['v']])
  
}

# get_block_full
get_block_full <- function(ajive_output, k, type){
  
  if(! type  %in% c('joint', 'individual', 'noise')){
    stop('type must be: joint, individual, or noise')
  }
  
  if(type  == 'noise'){
    ajive_output$block_decomps[[k]][[type]]
  } else{
    ajive_output$block_decomps[[k]][[type]][['full']]
  }
  
}

# get_common_normalized_scores
get_common_normalized_scores <- function(ajive_output){
  ajive_output[['joint_scores']]
}

# get_block_scores
get_block_scores <- function(ajive_output, k, type, normalized){
  if(! type  %in% c('joint', 'individual')){
    stop('type must be: joint or individual')
  }
  
  scores <- ajive_output$block_decomps[[k]][[type]][['u']]
  
  if(!normalized){
    D <- diag(ajive_output$block_decomps[[k]][[type]][['d']],
              ncol=dim(scores)[2])
    
    scores <- scores %*% D
  }
  
  scores
}

# get_block_loadings
get_block_loadings <- function(ajive_output, k, type){
  if(! type  %in% c('joint', 'individual')){
    stop('type must be: joint or individual')
  }
  
  ajive_output$block_decomps[[k]][[type]][['v']]
}

# get_joint_rank
get_joint_rank <- function(ajive_output){
  ajive_output$block_decomps[[1]][['joint']][['rank']]
}

# get_individual_rank
get_individual_rank <- function(ajive_output, k){
  ajive_output$block_decomps[[k]][['individual']][['rank']]
}

# get_final_decomposition
get_final_decomposition <- function(X, joint_scores, sv_threshold, full=TRUE){
  
  jive_decomposition <- list()
  jive_decomposition[['individual']] <- get_individual_decomposition(X, joint_scores, sv_threshold, full)
  jive_decomposition[['joint']] <- get_joint_decomposition(X, joint_scores, full)
  
  
  if(full){
    jive_decomposition[['noise']] <- X - (jive_decomposition[['joint']][['full']] +
                                            jive_decomposition[['individual']][['full']])
  } else{
    jive_decomposition[['noise']] <- NA
  }
  
  jive_decomposition
}

# get_individual_decomposition
get_individual_decomposition <- function(X, joint_scores, sv_threshold, full=TRUE){
  
  if(is.null(joint_scores)){
    indiv_decomposition <- get_svd(X)
  } else{
    X_orthog <- (diag(dim(X)[1]) - joint_scores %*% t(joint_scores)) %*% X
    indiv_decomposition <- get_svd(X_orthog)
  }
  
  
  indiv_rank <- sum(indiv_decomposition[['d']] > sv_threshold)
  
  indiv_decomposition <- truncate_svd(decomposition=indiv_decomposition, rank=indiv_rank)
  
  if(full){
    indiv_decomposition[['full']] <- svd_reconstruction(indiv_decomposition)
  } else{
    indiv_decomposition[['full']] <- NA
  }
  
  indiv_decomposition[['rank']] <- indiv_rank
  indiv_decomposition
}

# get_joint_decomposition
get_joint_decomposition <- function(X, joint_scores, full=TRUE){
  
  if(is.null(joint_scores)){
    joint_decomposition <- list(full= NA, rank=0, u=NA, d=NA, v=NA)
    return(joint_decomposition)
  }
  joint_rank <- dim(joint_scores)[2]
  J <-  joint_scores %*% t(joint_scores) %*% X
  
  joint_decomposition <- get_svd(J, joint_rank)
  
  if(full){
    joint_decomposition[['full']] <- J
  } else{
    joint_decomposition[['full']] <- NA
  }
  
  joint_decomposition[['rank']] <- joint_rank
  joint_decomposition
  
}

#--------------------------------------------------#
# AJIVE
#--------------------------------------------------#
ajive <- function(blocks, initial_signal_ranks, full=TRUE, n_wedin_samples=1000, n_rand_dir_samples=1000, joint_rank=NA){
  
  K <- length(blocks)
  
  if(K < 2){
    stop('ajive expects at least two data matrices.')
  }
  
  if(sum(sapply(blocks,function(X) any(is.na(X)))) > 0){
    stop('Some of the blocks has missing data -- ajive expects full data matrices.')
  }
  
  # ----------------------------------------------------- #
  # step 1: initial signal space extraction 
  # ----------------------------------------------------- #
  # initial estimate of signal space with SVD
  block_svd <- list()
  sv_thresholds <- rep(0, K)
  for(k in 1:K){
    block_svd[[k]] <- get_svd(blocks[[k]])
    sv_thresholds[k] <- get_sv_threshold(singular_values = block_svd[[k]][['d']],
                                         rank=initial_signal_ranks[k])
  }
  
  # ----------------------------------------------------- #
  # step 2: joint space estimation
  # ----------------------------------------------------- #
  out <- get_joint_scores(blocks, block_svd, initial_signal_ranks, sv_thresholds,
                          n_wedin_samples=n_wedin_samples,
                          n_rand_dir_samples=n_rand_dir_samples,
                          joint_rank=joint_rank)
  joint_rank_sel_results <- out$rank_sel_results
  joint_scores <- out$joint_scores
  
  joint_rank <- out[['rank_sel_results']][['joint_rank_estimate']] # dim(joint_scores)[2]
  
  # ----------------------------------------------------- #
  # step 3: final decomposition 
  # ----------------------------------------------------- #
  block_decomps <- list()
  for(k in 1:K){
    block_decomps[[k]] <- get_final_decomposition(X=blocks[[k]],
                                                  joint_scores=joint_scores,
                                                  sv_threshold=sv_thresholds[k])
  }
  
  jive_decomposition <- list(block_decomps=block_decomps)
  jive_decomposition[['joint_scores']] <- joint_scores
  jive_decomposition[['joint_rank']] <- joint_rank
  
  jive_decomposition[['joint_rank_sel']] <- joint_rank_sel_results
  jive_decomposition
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
YK_compute <- function(dd,TT,KK,r_J,r_G,Control_param,DGP,factor_list){
  
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
  
  if (r_G == 1){
    if (Control_param$type==1){
      Psi_indiv <- diag(0.8944,r_G) 
    }else if(Control_param$type==2){
      Psi_indiv <- diag(0.8367,r_G)  
    }
  }else if (r_G == 2){
    if (Control_param$type==1){
      Psi_indiv <- matrix(c(0.8213,-0.1142,
                            -0.1142,0.8213),ncol=r_G,byrow=TRUE) 
    }else if(Control_param$type==2){
      Psi_indiv <- matrix(c(0.8367,0,
                            0,0.8367),ncol=r_G,byrow=TRUE) 
    }
  }else if (r_G == 3){
    if (Control_param$type==1){
      Psi_indiv <- matrix(c(0.8154,-0.1377,-0.0298,
                            -0.1377,0.7168,-0.1377,
                            -0.0298,-0.1377,0.8154),ncol=r_G,byrow=TRUE)
    }else if(Control_param$type==2){
      Psi_indiv <- matrix(c(0.8367,0,0,
                            0,0.8367,0,
                            0,0,0.8367),ncol=r_G,byrow=TRUE)
    }
  }else if (r_G == 4){
    if (Control_param$type==1){
      Psi_indiv <- matrix(c(0.8148,-0.1372,-0.0241,0.0084,
                            -0.1372, 0.7118,-0.1575,-0.0241,
                            -0.0241,-0.1575,0.7118,-0.1372,
                            0.0084,-0.0241,-0.1372,0.8148),ncol=r_G,byrow=TRUE)
    }else if(Control_param$type==2){
      Psi_indiv <- matrix(c(0.8367,0,0,0,
                            0,0.8367,0,0,
                            0,0,0.8367,0,
                            0,0,0,0.8367),ncol=r_G,byrow=TRUE)
    }
  }
  
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
