#--------------------------------------------------#
# AJIVE
#--------------------------------------------------#
ajive <- function(blocks, initial_signal_ranks, full=TRUE, n_wedin_samples=1000, n_rand_dir_samples=1000, joint_rank=NA){
  
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