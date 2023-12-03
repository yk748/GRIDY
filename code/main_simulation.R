rm(list = ls())

# Packages required
library(mvtnorm)
library(Matrix)
library(combinat)

# For AJIVE:
# devtools::install_github("idc9/r_jive")
library(ajive)

# For SCA:
library(multiway)

# For miscellaneous functions
source("library_simulation.R")


#-----------------------------------------------------------------------------#
# Basic settings
#-----------------------------------------------------------------------------#
# # Control parameters 1: PF2, standard
type <- 1
b_min <- 0
b_max <- 1
c_min <- 5
c_max <- 10
sigma_xi <- 0.2
sigma_eps <- 1
Control_param_1 <- data.frame(type=type,sigma_xi=sigma_xi,b_min=b_min,b_max=b_max,
                              c_min=c_min, c_max=c_max,sigma_eps=sigma_eps)
#-----------------------------------------------------------------------------#
# Control parameters 2: IND
type <- 2
b_min <- 0
b_max <- 1
c_min <- 5
c_max <- 10
sigma_xi <- 0.3
sigma_eps <- 1
Control_param_2 <- data.frame(type=type,sigma_xi=sigma_xi,b_min=b_min,b_max=b_max,
                              c_min=c_min, c_max=c_max,sigma_eps=sigma_eps)
#-----------------------------------------------------------------------------#
# # Control parameters 3: PF2, larger loading
type <- 1
b_min <- 0
b_max <- 2
c_min <- 5
c_max <- 10
sigma_xi <- 0.2
sigma_eps <- 1
Control_param_3 <- data.frame(type=type,sigma_xi=sigma_xi,b_min=b_min,b_max=b_max,
                              c_min=c_min, c_max=c_max,sigma_eps=sigma_eps)
#-----------------------------------------------------------------------------#
# # Control parameters 4: PF2, larger standard dev of factor series
type <- 1
b_min <- 0
b_max <- 1
c_min <- 10
c_max <- 20
sigma_xi <- 0.2
sigma_eps <- 1
Control_param_4 <- data.frame(type=type,sigma_xi=sigma_xi,b_min=b_min,b_max=b_max,
                              c_min=c_min, c_max=c_max,sigma_eps=sigma_eps)
#-----------------------------------------------------------------------------#
# # Control parameters 5: PF2, larger noise
type <- 1
b_min <- 0
b_max <- 1
c_min <- 5
c_max <- 10
sigma_xi <- 0.2
sigma_eps <- 4
Control_param_5 <- data.frame(type=type,sigma_xi=sigma_xi,b_min=b_min,b_max=b_max,
                              c_min=c_min, c_max=c_max,sigma_eps=sigma_eps)

# -----------------------------------------------------------------------------#
# Main section
# -----------------------------------------------------------------------------#
## Main setting: Problem size
d_list <- c(200,400,800,100,100,100,100,100,100)
T_list <- c(100,100,100,200,400,800,200,200,200)
K_list <- c(50,50,50,50,50,50,100,200,400)

min_iter <- 1
max_iter <- 100
H <- 10

r_J_list <- c(2,2,3)
r_I_list <- c(2,3,3)

for (case in 5:5){
  
  dd <- d_list[case]
  TT <- T_list[case]
  KK <- K_list[case]
  
  for (set in 1:1){
    
    if (set == 1){
      Control_param <- Control_param_1
    }else if(set == 2){
      Control_param <- Control_param_2
    }else if(set ==3){
      Control_param <- Control_param_3
    }else if(set ==4){
      Control_param <- Control_param_4
    }
    cat("Currently, set",set,"is set.\n")
    
    for (r in 1:1){
      
      r_J <- r_J_list[r]
      r_I <- r_I_list[r]
      
      # ---------------------------------------------------------------- #
      # Preparing containers:
      # ---------------------------------------------------------------- #
      list_GRIDY <- list()
      list_DSCA <- list ()
      list_PCA <- list()
      list_Unfitted <- list()
      
      for (iter in 1:100){
        
        time <- Sys.time()
        #-----------------------------------------------------------------#
        # 0-1. Data generation:
        model_dgp_ext <- sim_model(r_J, r_I, dd, (TT+H), KK, Control_param)
        # 0-2. Use this for estimation:
        model_dgp <- list()
        for (kk in 1:(2*KK)){
          model_dgp$Xk[[kk]] <- model_dgp_ext$Xk[[kk]][1:TT,]
          model_dgp$X_bar[[kk]] <- model_dgp_ext$X_bar[[kk]][1:TT,]
          model_dgp$X_tilde[[kk]] <- model_dgp_ext$X_tilde[[kk]][1:TT,]
          model_dgp$scaled_X[[kk]] <- model_dgp_ext$scaled_X[[kk]][,1:TT]
          model_dgp$scaled_X_comm[[kk]] <- model_dgp_ext$scaled_X_comm[[kk]][,1:TT]
          model_dgp$scaled_X_ind[[kk]] <- model_dgp_ext$scaled_X_ind[[kk]][,1:TT]
          model_dgp$Ek[[kk]] <- model_dgp_ext$Ek[[kk]][1:TT,]
          
          if (r_J == 1){
            model_dgp$F_bar[[kk]] <- model_dgp_ext$F_bar[[kk]][1:TT]
            model_dgp$A_bar[[kk]] <- model_dgp_ext$A_bar[[kk]][1:TT]
          }else{
            model_dgp$F_bar[[kk]] <- model_dgp_ext$F_bar[[kk]][1:TT,]
            model_dgp$A_bar[[kk]] <- model_dgp_ext$A_bar[[kk]][1:TT,]
          }
          
          if (r_I == 1){
            model_dgp$F_tilde[[kk]] <- model_dgp_ext$F_tilde[[kk]][1:TT]
            model_dgp$A_tilde[[kk]] <- model_dgp_ext$A_tilde[[kk]][1:TT]
          }else{
            model_dgp$F_tilde[[kk]] <- model_dgp_ext$F_tilde[[kk]][1:TT,]
            model_dgp$A_tilde[[kk]] <- model_dgp_ext$A_tilde[[kk]][1:TT,]
          }
        }
        model_dgp$C_bar <- model_dgp_ext$C_bar
        model_dgp$C_tilde <- model_dgp_ext$C_tilde
        model_dgp$B_bar <- model_dgp_ext$B_bar
        model_dgp$B_tilde1 <- model_dgp_ext$B_tilde1
        model_dgp$B_tilde2 <- model_dgp_ext$B_tilde2
        model_dgp$comm_idx <- model_dgp_ext$comm_idx
        model_dgp$group1_idx <- model_dgp_ext$group1_idx
        model_dgp$group2_idx <-model_dgp_ext$group2_idx
        
        #-----------------------------------------------------------------#
        # 1. Run AJIVE (GRIDY,PCA,Unfitted):
        #-----------------------------------------------------------------#
        possibleError <- tryCatch(
          GRIDY_AJIVE <- ajive(model_dgp$scaled_X,
                               initial_signal_ranks = rep((r_J+r_I),(2*KK)),
                               n_wedin_samples = 1000,
                               n_rand_dir_samples = 1000,
                               full = TRUE, joint_rank = r_J),
          error = function(iter){
            print(iter,"th iteration is failed to run AJIVE")
          }
        )
        
        Joint_block <- list()
        Group1_block <- list()
        Group2_block <- list()
        
        for (kk in 1:(2*KK)){
          Joint_block[[kk]] <- t(GRIDY_AJIVE[["block_decomps"]][[kk]]$joint$full)
          if (kk <= KK){
            Group1_block[[kk]] <- t(GRIDY_AJIVE[["block_decomps"]][[kk]]$individual$full)
          }else{
            Group2_block[[(kk-KK)]] <- t(GRIDY_AJIVE[["block_decomps"]][[kk]]$individual$full)
          }
        }
        
        #-----------------------------------------------------------------#
        # 2. Run SCA (GRIDY,Unfitted,Unfitted):
        #-----------------------------------------------------------------#
        # Run SCA for joint structure:
        possibleError <- tryCatch(
          GRIDY_PF2_Joint <- sca(Joint_block, nfac = r_J,type="sca-pf2", verbose = FALSE),
          error = function(iter){
            print(iter,"th iteration is failed to estimate Joint")
          }
        )
        # Run SCA for 1st group structure:
        possibleError <- tryCatch(
          GRIDY_PF2_Group1 <- sca(Group1_block, nfac = r_I,type="sca-pf2", verbose = FALSE),
          error = function(iter){
            print(iter,"th iteration is failed to estimate Group1")
          }
        )
        # Run SCA for 2nd group structure:
        possibleError <- tryCatch(
          GRIDY_PF2_Group2 <- sca(Group2_block, nfac = r_I,type="sca-pf2", verbose = FALSE),
          error = function(iter){
            print(iter,"th iteration is failed to estimate Group2")
          }
        )
        
        #-----------------------------------------------------------------#
        # 3. Run Yule-Walker (GRIDY,DSCA,PCA):
        #-----------------------------------------------------------------#
        # Refit factor series
        GRIDY_refitted <- factor_regression(KK,r_J,r_I,model_dgp,
                                            GRIDY_PF2_Joint,
                                            GRIDY_PF2_Group1,
                                            GRIDY_PF2_Group2)
        # Run Yule-Walker
        GRIDY_YK <- YK_compute(dd,TT,KK,r_J,r_I,Control_param,model_dgp,GRIDY_refitted)
        
        #-----------------------------------------------------------------#
        # Double SCA:
        #-----------------------------------------------------------------#
        Whole_block <- list()
        
        for (kk in 1:(2*KK)){
          Whole_block[[kk]] <- t(model_dgp$scaled_X[[kk]])
        }
        
        possibleError <- tryCatch(
          DSCA_Joint <- sca(Whole_block, nfac = r_J,type="sca-pf2", verbose = FALSE),
          error = function(iter){
            print(iter,"th iteration is failed to estimate first stage of DSCA")
          }
        )
        
        Remaining_block_group1 <- list()
        Remaining_block_group2 <- list()
        for (kk in 1:(2*KK)){
          kth_Joint_block <- DSCA_Joint$D[[kk]] %*% t(DSCA_Joint$B)
          if (kk <= KK){
            Remaining_block_group1[[kk]] <- Whole_block[[kk]] - kth_Joint_block
          }else{
            Remaining_block_group2[[(kk-KK)]] <- Whole_block[[kk]] - kth_Joint_block
          }
        }
        
        possibleError <- tryCatch(
          DSCA_Group1 <- sca(Remaining_block_group1, nfac = r_I,type="sca-pf2", verbose = FALSE),
          error = function(iter){
            print(iter,"th iteration is failed to estimate Group 1 at second stage of DSCA")
          }
        )
        
        possibleError <- tryCatch(
          DSCA_Group2 <- sca(Remaining_block_group2, nfac = r_I,type="sca-pf2", verbose = FALSE),
          error = function(iter){
            print(iter,"th iteration is failed to estimate Group 2 at second stage of DSCA")
          }
        )
        
        
        DSCA_refitted <- factor_regression(KK,r_J,r_I,model_dgp,
                                           DSCA_Joint,
                                           DSCA_Group1,
                                           DSCA_Group2)
        
        DSCA_YK <- YK_compute(dd,TT,KK,r_J,r_I,Control_param,model_dgp,DSCA_refitted)
        
        #-----------------------------------------------------------------#
        # Individual PCA:
        #-----------------------------------------------------------------#
        PCA_Joint <- list()
        PCA_Group1 <- list()
        PCA_Group2 <- list()
        
        possibleError <- tryCatch(
          for (kk in 1:(2*KK)){
            
            svd_joint <- svd(Joint_block[[kk]])
            B_k_joint <- svd_joint$v[,1:r_J]
            D_k_joint <- svd_joint$u[,1:r_J] %*% diag(svd_joint$d[1:r_J],r_J)
            Cov_F_joint <- diag(diag(t(D_k_joint) %*% D_k_joint/TT),r_J)
            
            if (kk == 1){
              PCA_Joint$D[[1]] <- D_k_joint %*% solve(sqrt(Cov_F_joint))
              PCA_Joint$B <- B_k_joint %*% sqrt(Cov_F_joint)
            }else{
              D_k_joint <- D_k_joint %*% solve(sqrt(Cov_F_joint))
              B_k_joint <- B_k_joint %*% sqrt(Cov_F_joint)
              Q_rotate_joint <- solve(t(B_k_joint) %*% B_k_joint) %*% t(B_k_joint) %*% PCA_Joint$B
              PCA_Joint$D[[kk]] <- D_k_joint %*% Q_rotate_joint
            }
            
            if (kk <= KK){
              svd_group1 <- svd(Group1_block[[kk]])
              B_k_group1 <- svd_group1$v[,1:r_I]
              D_k_group1 <- svd_group1$u[,1:r_I] %*% diag(svd_group1$d[1:r_I],r_I)
              Cov_F_group1 <- diag(diag(t(D_k_group1) %*% D_k_group1/TT),r_I)
              
              if (kk == 1){
                PCA_Group1$D[[1]] <- D_k_group1 %*% solve(sqrt(Cov_F_group1))
                PCA_Group1$B <- B_k_group1 %*% sqrt(Cov_F_group1)
                
              }else{
                D_k_group1 <- D_k_group1 %*% solve(sqrt(Cov_F_group1))
                B_k_group1 <- B_k_group1 %*% sqrt(Cov_F_group1)
                Q_rotate_group1 <- solve(t(B_k_group1) %*% B_k_group1) %*% t(B_k_group1) %*% PCA_Group1$B
                PCA_Group1$D[[kk]] <- D_k_group1 %*% Q_rotate_group1
              }
            }else{
              svd_group2 <- svd(Group2_block[[(kk-KK)]])
              B_k_group2 <- svd_group2$v[,1:r_I]
              D_k_group2 <- svd_group2$u[,1:r_I] %*% diag(svd_group2$d[1:r_I],r_I)
              Cov_F_group2 <- diag(diag(t(D_k_group2) %*% D_k_group2/TT),r_I)
              
              if (kk == (KK+1)){
                PCA_Group2$D[[1]] <- D_k_group2 %*% solve(sqrt(Cov_F_group2))
                PCA_Group2$B <- B_k_group2 %*% sqrt(Cov_F_group2)
              }else{
                D_k_group2 <- D_k_group2 %*% solve(sqrt(Cov_F_group2))
                B_k_group2 <- B_k_group2 %*% sqrt(Cov_F_group2)
                Q_rotate_group2 <- solve(t(B_k_group2) %*% B_k_group2) %*% t(B_k_group2) %*% PCA_Group2$B
                PCA_Group2$D[[(kk-KK)]] <- D_k_group2 %*% Q_rotate_group2
              }
            }
          },
          error = function(iter){
            print(iter,"th iteration is failed to estimate PCA")
          }
        )
        
        PCA_refitted <- factor_regression(KK,r_J,r_I,model_dgp,
                                          PCA_Joint,
                                          PCA_Group1,
                                          PCA_Group2)
        PCA_YK <- YK_compute(dd,TT,KK,r_J,r_I,Control_param,model_dgp,PCA_refitted)
        
        #-----------------------------------------------------------------#
        # Unfitted:
        #-----------------------------------------------------------------#
        Factor_unifitted <- list()
        Factor_unifitted$factor_joint <- GRIDY_PF2_Joint$D
        Factor_unifitted$factor_group1 <- GRIDY_PF2_Group1$D
        Factor_unifitted$factor_group2 <- GRIDY_PF2_Group2$D
        
        Unfitted_YK <- YK_compute(dd,TT,KK,r_J,r_I,Control_param,model_dgp,Factor_unifitted)
        
        
        #-----------------------------------------------------------------#
        # 4. Evaluation (GRIDY,DSCA,PCA,Unfitted):
        #-----------------------------------------------------------------#
        
        # R2: GRIDY, DSCA, PCA, Unfitted
        total_variation <- 0
        var_GRIDY <- 0
        var_DSCA <- 0
        var_PCA <- 0
        var_Unfitted <- 0
        for (kk in 1:(2*KK)){
          total_variation <- total_variation + norm(t(model_dgp$scaled_X[[kk]]),"F")^2
          if (kk <= KK){
            var_GRIDY <- var_GRIDY + norm(t(model_dgp$scaled_X[[kk]])
                                          - GRIDY_refitted$factor_joint[[kk]] %*% t(GRIDY_PF2_Joint$B)
                                          - GRIDY_refitted$factor_group1[[kk]] %*% t(GRIDY_PF2_Group1$B) ,"F")^2
            
            var_DSCA <- var_DSCA + norm(t(model_dgp$scaled_X[[kk]])
                                        - DSCA_refitted$factor_joint[[kk]] %*% t(DSCA_Joint$B)
                                        - DSCA_refitted$factor_group1[[kk]] %*% t(DSCA_Group1$B) ,"F")^2
            
            var_PCA <- var_PCA + norm(t(model_dgp$scaled_X[[kk]])
                                      - PCA_refitted$factor_joint[[kk]] %*% t(PCA_Joint$B)
                                      - PCA_refitted$factor_group1[[kk]] %*% t(PCA_Group1$B) ,"F")^2
            
            var_Unfitted <- var_Unfitted + norm(t(model_dgp$scaled_X[[kk]])
                                                - GRIDY_PF2_Joint$D[[kk]] %*% t(GRIDY_PF2_Joint$B)
                                                - GRIDY_PF2_Group1$D[[kk]] %*% t(GRIDY_PF2_Group1$B) ,"F")^2
            
          }else{
            var_GRIDY <- var_GRIDY + norm(t(model_dgp$scaled_X[[kk]])
                                          - GRIDY_refitted$factor_joint[[kk]] %*% t(GRIDY_PF2_Joint$B)
                                          - GRIDY_refitted$factor_group2[[(kk-KK)]] %*% t(GRIDY_PF2_Group2$B) ,"F")^2
            
            var_DSCA <- var_DSCA + norm(t(model_dgp$scaled_X[[kk]])
                                        - DSCA_refitted$factor_joint[[kk]] %*% t(DSCA_Joint$B)
                                        - DSCA_refitted$factor_group2[[(kk-KK)]] %*% t(DSCA_Group2$B) ,"F")^2
            
            var_PCA <- var_PCA + norm(t(model_dgp$scaled_X[[kk]])
                                      - PCA_refitted$factor_joint[[kk]] %*% t(PCA_Joint$B)
                                      - PCA_refitted$factor_group2[[(kk-KK)]] %*% t(PCA_Group2$B) ,"F")^2
            
            var_Unfitted <- var_Unfitted + norm(t(model_dgp$scaled_X[[kk]])
                                                - GRIDY_PF2_Joint$D[[kk]] %*% t(GRIDY_PF2_Joint$B)
                                                - GRIDY_PF2_Group2$D[[(kk-KK)]] %*% t(GRIDY_PF2_Group2$B) ,"F")^2
          }
        }
        list_GRIDY[["R2"]][[iter]] <- 1 - var_GRIDY / total_variation
        list_DSCA[["R2"]][[iter]] <- 1 - var_DSCA / total_variation
        list_PCA[["R2"]][[iter]] <- 1 - var_PCA / total_variation
        list_Unfitted[["R2"]][[iter]] <- 1 - var_Unfitted / total_variation
        
        
        # RMSE: GRIDY, DSCA, PCA, Unfitted
        RMSE_GRIDY_joint <- 0
        RMSE_DSCA_joint <- 0
        RMSE_PCA_joint <- 0
        RMSE_Unfitted_joint <- 0
        
        RMSE_GRIDY_group1 <- 0
        RMSE_DSCA_group1 <- 0
        RMSE_PCA_group1 <- 0
        RMSE_Unfitted_group1 <- 0
        
        RMSE_GRIDY_group2 <- 0
        RMSE_DSCA_group2 <- 0
        RMSE_PCA_group2 <- 0
        RMSE_Unfitted_group2 <- 0
        
        for (kk in 1:(2*KK)){
          RMSE_GRIDY_joint <- RMSE_GRIDY_joint + norm(t(model_dgp$scaled_X_comm[[kk]])
                                                      - GRIDY_refitted$factor_joint[[kk]] %*% t(GRIDY_PF2_Joint$B), "F")^2 / (TT*dd*(2*KK))
          
          RMSE_DSCA_joint <- RMSE_DSCA_joint + norm(t(model_dgp$scaled_X_comm[[kk]])
                                                    - DSCA_refitted$factor_joint[[kk]] %*% t(DSCA_Joint$B), "F")^2 / (TT*dd*(2*KK))
          
          RMSE_PCA_joint <- RMSE_PCA_joint + norm(t(model_dgp$scaled_X_comm[[kk]])
                                                  - PCA_refitted$factor_joint[[kk]] %*% t(PCA_Joint$B), "F")^2 / (TT*dd*(2*KK))
          
          RMSE_Unfitted_joint <- RMSE_Unfitted_joint + norm(t(model_dgp$scaled_X_comm[[kk]])
                                                            - GRIDY_PF2_Joint$D[[kk]] %*% t(GRIDY_PF2_Joint$B), "F")^2 / (TT*dd*(2*KK))
          
          if (kk <= KK){
            RMSE_GRIDY_group1 <- RMSE_GRIDY_group1 + norm(t(model_dgp$scaled_X_ind[[kk]])
                                                          - GRIDY_refitted$factor_group1[[kk]] %*% t(GRIDY_PF2_Group1$B), "F")^2 / (TT*dd*KK)
            
            RMSE_DSCA_group1 <- RMSE_DSCA_group1 + norm(t(model_dgp$scaled_X_ind[[kk]])
                                                        - DSCA_refitted$factor_group1[[kk]] %*% t(DSCA_Group1$B), "F")^2 / (TT*dd*KK)
            
            RMSE_PCA_group1 <- RMSE_PCA_group1 + norm(t(model_dgp$scaled_X_ind[[kk]])
                                                      - PCA_refitted$factor_group1[[kk]] %*% t(PCA_Group1$B), "F")^2 / (TT*dd*KK)
            
            RMSE_Unfitted_group1 <- RMSE_Unfitted_group1 + norm(t(model_dgp$scaled_X_ind[[kk]])
                                                                - GRIDY_PF2_Group1$D[[kk]] %*% t(GRIDY_PF2_Group1$B), "F")^2 / (TT*dd*KK)
            
          }else{
            RMSE_GRIDY_group2 <- RMSE_GRIDY_group2 + norm(t(model_dgp$scaled_X_ind[[kk]])
                                                          - GRIDY_refitted$factor_group2[[(kk-KK)]] %*% t(GRIDY_PF2_Group2$B), "F")^2 / (TT*dd*KK)
            
            RMSE_DSCA_group2 <- RMSE_DSCA_group2 + norm(t(model_dgp$scaled_X_ind[[kk]])
                                                        - DSCA_refitted$factor_group2[[(kk-KK)]] %*% t(DSCA_Group2$B), "F")^2 / (TT*dd*KK)
            
            RMSE_PCA_group2 <- RMSE_PCA_group2 + norm(t(model_dgp$scaled_X_ind[[kk]])
                                                      - PCA_refitted$factor_group2[[(kk-KK)]] %*% t(PCA_Group2$B), "F")^2 / (TT*dd*KK)
            
            RMSE_Unfitted_group2 <- RMSE_Unfitted_group2 + norm(t(model_dgp$scaled_X_ind[[kk]])
                                                                - GRIDY_PF2_Group2$D[[(kk-KK)]] %*% t(GRIDY_PF2_Group2$B), "F")^2 / (TT*dd*KK)
            
          }
        }
        list_GRIDY[["RMSE"]][["Joint"]][[iter]] <- sqrt(RMSE_GRIDY_joint)
        list_GRIDY[["RMSE"]][["Group1"]][[iter]] <- sqrt(RMSE_GRIDY_group1)
        list_GRIDY[["RMSE"]][["Group2"]][[iter]] <- sqrt(RMSE_GRIDY_group2)
        list_DSCA[["RMSE"]][["Joint"]][[iter]] <- sqrt(RMSE_DSCA_joint)
        list_DSCA[["RMSE"]][["Group1"]][[iter]] <- sqrt(RMSE_DSCA_group1)
        list_DSCA[["RMSE"]][["Group2"]][[iter]] <- sqrt(RMSE_DSCA_group2)
        list_PCA[["RMSE"]][["Joint"]][[iter]] <- sqrt(RMSE_PCA_joint)
        list_PCA[["RMSE"]][["Group1"]][[iter]] <- sqrt(RMSE_PCA_group1)
        list_PCA[["RMSE"]][["Group2"]][[iter]] <- sqrt(RMSE_PCA_group2)
        list_Unfitted[["RMSE"]][["Joint"]][[iter]] <- sqrt(RMSE_Unfitted_joint)
        list_Unfitted[["RMSE"]][["Group1"]][[iter]] <- sqrt(RMSE_Unfitted_group1)
        list_Unfitted[["RMSE"]][["Group2"]][[iter]] <- sqrt(RMSE_Unfitted_group2)
        
        # CC_B: GRIDY, DSCA, PCA
        W_perm_joint <- matrix(unlist(permn(r_J)),nrow=factorial(r_J),byrow=TRUE)
        Con_B_joint_GRIDY <- rep(NA, nrow(W_perm_joint))
        Con_B_joint_DSCA <- rep(NA, nrow(W_perm_joint))
        Con_B_joint_PCA <- rep(NA, nrow(W_perm_joint))
        Con_B_joint_Unfitted <- rep(NA, nrow(W_perm_joint))
        
        for(i in 1:nrow(W_perm_joint)){
          Con_B_joint_GRIDY[i] <- sum(diag(abs(as.matrix(congru(GRIDY_PF2_Joint$B[,W_perm_joint[i,]],model_dgp$B_bar)))))/r_J
          Con_B_joint_DSCA[i] <- sum(diag(abs(as.matrix(congru(DSCA_Joint$B[,W_perm_joint[i,]],model_dgp$B_bar)))))/r_J
          Con_B_joint_PCA[i] <- sum(diag(abs(as.matrix(congru(PCA_Joint$B[,W_perm_joint[i,]],model_dgp$B_bar)))))/r_J
          Con_B_joint_Unfitted[i] <- Con_B_joint_GRIDY[i]
        }
        idx_joint_GRIDY <- which.max(Con_B_joint_GRIDY)
        Con_B_joint_GRIDY <- Con_B_joint_GRIDY[idx_joint_GRIDY]
        
        idx_joint_DSCA <- which.max(Con_B_joint_DSCA)
        Con_B_joint_DSCA <- Con_B_joint_DSCA[idx_joint_DSCA]
        
        idx_joint_PCA <- which.max(Con_B_joint_PCA)
        Con_B_joint_PCA <- Con_B_joint_PCA[idx_joint_PCA]
        
        Con_B_joint_Unfitted <- Con_B_joint_GRIDY
        
        list_GRIDY[["CC_B"]][["Joint"]][[iter]] <- Con_B_joint_GRIDY
        list_DSCA[["CC_B"]][["Joint"]][[iter]] <- Con_B_joint_DSCA
        list_PCA[["CC_B"]][["Joint"]][[iter]] <- Con_B_joint_PCA
        list_Unfitted[["CC_B"]][["Joint"]][[iter]] <- Con_B_joint_Unfitted
        
        
        W_perm_group1 <- matrix(unlist(permn(r_I)),nrow=factorial(r_I),byrow=TRUE)
        Con_B_group1_GRIDY <- rep(NA, nrow(W_perm_group1))
        Con_B_group1_DSCA <- rep(NA, nrow(W_perm_group1))
        Con_B_group1_PCA <- rep(NA, nrow(W_perm_group1))
        Con_B_group1_Unfitted <- rep(NA, nrow(W_perm_group1))
        
        for(i in 1:nrow(W_perm_group1)){
          Con_B_group1_GRIDY[i] <- sum(diag(abs(as.matrix(congru(GRIDY_PF2_Group1$B[,W_perm_group1[i,]],model_dgp$B_tilde1)))))/r_I
          Con_B_group1_DSCA[i] <- sum(diag(abs(as.matrix(congru(DSCA_Group1$B[,W_perm_group1[i,]],model_dgp$B_tilde1)))))/r_I
          Con_B_group1_PCA[i] <- sum(diag(abs(as.matrix(congru(PCA_Group1$B[,W_perm_group1[i,]],model_dgp$B_tilde1)))))/r_I
          Con_B_group1_Unfitted[i] <- Con_B_group1_GRIDY[i]
        }
        idx_group1_GRIDY <- which.max(Con_B_group1_GRIDY)
        Con_B_group1_GRIDY <- Con_B_group1_GRIDY[idx_group1_GRIDY]
        
        idx_group1_DSCA <- which.max(Con_B_group1_DSCA)
        Con_B_group1_DSCA <- Con_B_group1_DSCA[idx_group1_DSCA]
        
        idx_group1_PCA <- which.max(Con_B_group1_PCA)
        Con_B_group1_PCA <- Con_B_group1_PCA[idx_group1_PCA]
        
        Con_B_group1_Unfitted <- Con_B_group1_GRIDY
        
        list_GRIDY[["CC_B"]][["Group1"]][[iter]] <- Con_B_group1_GRIDY
        list_DSCA[["CC_B"]][["Group1"]][[iter]] <- Con_B_group1_DSCA
        list_PCA[["CC_B"]][["Group1"]][[iter]] <- Con_B_group1_PCA
        list_Unfitted[["CC_B"]][["Group1"]][[iter]] <- Con_B_group1_Unfitted
        
        
        W_perm_group2 <- matrix(unlist(permn(r_I)),nrow=factorial(r_I),byrow=TRUE)
        Con_B_group2_GRIDY <- rep(NA, nrow(W_perm_group2))
        Con_B_group2_DSCA <- rep(NA, nrow(W_perm_group2))
        Con_B_group2_PCA <- rep(NA, nrow(W_perm_group2))
        Con_B_group2_Unfitted <- rep(NA, nrow(W_perm_group2))
        
        for(i in 1:nrow(W_perm_group2)){
          Con_B_group2_GRIDY[i] <- sum(diag(abs(as.matrix(congru(GRIDY_PF2_Group2$B[,W_perm_group2[i,]],model_dgp$B_tilde2)))))/r_I
          Con_B_group2_DSCA[i] <- sum(diag(abs(as.matrix(congru(DSCA_Group2$B[,W_perm_group2[i,]],model_dgp$B_tilde2)))))/r_I
          Con_B_group2_PCA[i] <- sum(diag(abs(as.matrix(congru(PCA_Group2$B[,W_perm_group2[i,]],model_dgp$B_tilde2)))))/r_I
          Con_B_group2_Unfitted[i] <- Con_B_group2_GRIDY[i]
        }
        idx_group2_GRIDY <- which.max(Con_B_group2_GRIDY)
        Con_B_group2_GRIDY <- Con_B_group2_GRIDY[idx_group2_GRIDY]
        
        idx_group2_DSCA <- which.max(Con_B_group2_DSCA)
        Con_B_group2_DSCA <- Con_B_group2_DSCA[idx_group2_DSCA]
        
        idx_group2_PCA <- which.max(Con_B_group2_PCA)
        Con_B_group2_PCA <- Con_B_group2_PCA[idx_group2_PCA]
        
        Con_B_group2_Unfitted <- Con_B_group2_GRIDY
        
        list_GRIDY[["CC_B"]][["Group2"]][[iter]] <- Con_B_group2_GRIDY
        list_DSCA[["CC_B"]][["Group2"]][[iter]] <- Con_B_group2_DSCA
        list_PCA[["CC_B"]][["Group2"]][[iter]] <- Con_B_group2_PCA
        list_Unfitted[["CC_B"]][["Group2"]][[iter]] <- Con_B_group2_Unfitted
        
        
        # CC_F: GRIDY, DSCA, PCA, Unfitted
        
        Reorder_joint_GRIDY <- lapply(X=c(1:(2*KK)),FUN=function(X) GRIDY_refitted$factor_joint[[X]][,W_perm_joint[idx_joint_GRIDY,]])
        Reorder_joint_DSCA <- lapply(X=c(1:(2*KK)),FUN=function(X) DSCA_refitted$factor_joint[[X]][,W_perm_joint[idx_joint_DSCA,]])
        Reorder_joint_PCA <- lapply(X=c(1:(2*KK)),FUN=function(X) PCA_refitted$factor_joint[[X]][,W_perm_joint[idx_joint_PCA,]])
        Reorder_joint_Unfitted <- lapply(X=c(1:(2*KK)),FUN=function(X) GRIDY_PF2_Joint$D[[X]][,W_perm_joint[idx_joint_GRIDY,]])
        
        Reorder_group1_GRIDY <- lapply(X=c(1:KK),FUN=function(X) GRIDY_refitted$factor_group1[[X]][,W_perm_group1[idx_group1_GRIDY,]])
        Reorder_group1_DSCA <- lapply(X=c(1:KK),FUN=function(X) DSCA_refitted$factor_group1[[X]][,W_perm_group1[idx_group1_DSCA,]])
        Reorder_group1_PCA <- lapply(X=c(1:KK),FUN=function(X) PCA_refitted$factor_group1[[X]][,W_perm_group1[idx_group1_PCA,]])
        Reorder_group1_Unfitted <- lapply(X=c(1:KK),FUN=function(X) GRIDY_PF2_Group1$D[[X]][,W_perm_group1[idx_group1_GRIDY,]])
        
        Reorder_group2_GRIDY <- lapply(X=c(1:KK),FUN=function(X) GRIDY_refitted$factor_group2[[X]][,W_perm_group2[idx_group2_GRIDY,]])
        Reorder_group2_DSCA <- lapply(X=c(1:KK),FUN=function(X) DSCA_refitted$factor_group2[[X]][,W_perm_group2[idx_group2_DSCA,]])
        Reorder_group2_PCA <- lapply(X=c(1:KK),FUN=function(X) PCA_refitted$factor_group2[[X]][,W_perm_group2[idx_group2_PCA,]])
        Reorder_group2_Unfitted <- lapply(X=c(1:KK),FUN=function(X) GRIDY_PF2_Group2$D[[X]][,W_perm_group2[idx_group2_GRIDY,]])
        
        Con_F_joint_GRIDY <- 0
        Con_F_joint_DSCA <- 0
        Con_F_joint_PCA <- 0
        Con_F_joint_Unfitted <- 0
        
        Con_F_group1_GRIDY <- 0
        Con_F_group1_DSCA <- 0
        Con_F_group1_PCA <- 0
        Con_F_group1_Unfitted <- 0
        
        Con_F_group2_GRIDY <- 0
        Con_F_group2_DSCA <- 0
        Con_F_group2_PCA <- 0
        Con_F_group2_Unfitted <- 0
        for(kk in 1:(2*KK)) {
          Con_F_joint_GRIDY <- Con_F_joint_GRIDY + sum(abs(diag(matrix(congru(Reorder_joint_GRIDY[[kk]],model_dgp$F_bar[[kk]]),nrow=r_J))))/ (r_J*(2*KK))
          
          Con_F_joint_DSCA <- Con_F_joint_DSCA + sum(abs(diag(matrix(congru(Reorder_joint_DSCA[[kk]],model_dgp$F_bar[[kk]]),nrow=r_J))))/ (r_J*(2*KK))
          
          Con_F_joint_PCA <- Con_F_joint_PCA + sum(abs(diag(matrix(congru(Reorder_joint_PCA[[kk]],model_dgp$F_bar[[kk]]),nrow=r_J))))/ (r_J*(2*KK))
          
          Con_F_joint_Unfitted <- Con_F_joint_Unfitted + sum(abs(diag(matrix(congru(Reorder_joint_Unfitted[[kk]],model_dgp$F_bar[[kk]]),nrow=r_J))))/ (r_J*(2*KK))
          if (kk <= KK){
            Con_F_group1_GRIDY <- Con_F_group1_GRIDY + sum(abs(diag(matrix(congru(Reorder_group1_GRIDY[[kk]],model_dgp$F_tilde[[kk]]),nrow=r_I))))/ (r_I*KK)
            
            Con_F_group1_DSCA <- Con_F_group1_DSCA +  sum(abs(diag(matrix(congru(Reorder_group1_DSCA[[kk]],model_dgp$F_tilde[[kk]]),nrow=r_I))))/ (r_I*KK)
            
            Con_F_group1_PCA <- Con_F_group1_PCA +  sum(abs(diag(matrix(congru(Reorder_group1_PCA[[kk]],model_dgp$F_tilde[[kk]]),nrow=r_I))))/ (r_I*KK)
            
            Con_F_group1_Unfitted <- Con_F_group1_Unfitted +  sum(abs(diag(matrix(congru(Reorder_group1_Unfitted[[kk]],model_dgp$F_tilde[[kk]]),nrow=r_I))))/ (r_I*KK)
          }else{
            Con_F_group2_GRIDY <- Con_F_group2_GRIDY + sum(abs(diag(matrix(congru(Reorder_group2_GRIDY[[(kk-KK)]],model_dgp$F_tilde[[kk]]),nrow=r_I))))/ (r_I*KK)
            
            Con_F_group2_DSCA <- Con_F_group2_DSCA +  sum(abs(diag(matrix(congru(Reorder_group2_DSCA[[(kk-KK)]],model_dgp$F_tilde[[kk]]),nrow=r_I))))/ (r_I*KK)
            
            Con_F_group2_PCA <- Con_F_group2_PCA +  sum(abs(diag(matrix(congru(Reorder_group2_PCA[[(kk-KK)]],model_dgp$F_tilde[[kk]]),nrow=r_I))))/ (r_I*KK)
            
            Con_F_group2_Unfitted <- Con_F_group2_Unfitted +  sum(abs(diag(matrix(congru(Reorder_group2_Unfitted[[(kk-KK)]],model_dgp$F_tilde[[kk]]),nrow=r_I))))/ (r_I*KK)
          }
        }
        list_GRIDY[["CC_F"]][["Joint"]][[iter]] <- Con_F_joint_GRIDY
        list_DSCA[["CC_F"]][["Joint"]][[iter]] <- Con_F_joint_DSCA
        list_PCA[["CC_F"]][["Joint"]][[iter]] <- Con_F_joint_PCA
        list_Unfitted[["CC_F"]][["Joint"]][[iter]] <- Con_F_joint_Unfitted
        
        list_GRIDY[["CC_F"]][["Group1"]][[iter]] <- Con_F_group1_GRIDY
        list_DSCA[["CC_F"]][["Group1"]][[iter]] <- Con_F_group1_DSCA
        list_PCA[["CC_F"]][["Group1"]][[iter]] <- Con_F_group1_PCA
        list_Unfitted[["CC_F"]][["Group1"]][[iter]] <- Con_F_group1_Unfitted
        
        list_GRIDY[["CC_F"]][["Group2"]][[iter]] <- Con_F_group2_GRIDY
        list_DSCA[["CC_F"]][["Group2"]][[iter]] <- Con_F_group2_DSCA
        list_PCA[["CC_F"]][["Group2"]][[iter]] <- Con_F_group2_PCA
        list_Unfitted[["CC_F"]][["Group2"]][[iter]] <- Con_F_group2_Unfitted
        
        
        # CC_Psi: GRIDY, DSCA, PCA, Unfitted
        Con_Psi_joint_GRIDY <- 0
        Con_Psi_joint_DSCA <- 0
        Con_Psi_joint_PCA <- 0
        Con_Psi_joint_Unfitted <- 0
        
        Con_Psi_group1_GRIDY <- 0
        Con_Psi_group1_DSCA <- 0
        Con_Psi_group1_PCA <- 0
        Con_Psi_group1_Unfitted <- 0
        
        Con_Psi_group2_GRIDY <- 0
        Con_Psi_group2_DSCA <- 0
        Con_Psi_group2_PCA <- 0
        Con_Psi_group2_Unfitted <- 0
        for (kk in 1:(2*KK)){
          Con_Psi_joint_GRIDY <- Con_Psi_joint_GRIDY + sum(abs(diag(
            congru(GRIDY_YK$Psi_joint_hat[W_perm_joint[idx_joint_GRIDY,],W_perm_joint[idx_joint_GRIDY,],kk],
                   GRIDY_YK$Psi_joint_k[W_perm_joint[idx_joint_GRIDY,],W_perm_joint[idx_joint_GRIDY,],kk])))) / (r_J*(2*KK))
          
          Con_Psi_joint_DSCA <- Con_Psi_joint_DSCA + sum(abs(diag(
            congru(DSCA_YK$Psi_joint_hat[W_perm_joint[idx_joint_DSCA,],W_perm_joint[idx_joint_DSCA,],kk],
                   DSCA_YK$Psi_joint_k[W_perm_joint[idx_joint_DSCA,],W_perm_joint[idx_joint_DSCA,],kk])))) / (r_J*(2*KK))
          
          Con_Psi_joint_PCA <- Con_Psi_joint_PCA + sum(abs(diag(
            congru(PCA_YK$Psi_joint_hat[W_perm_joint[idx_joint_PCA,],W_perm_joint[idx_joint_PCA,],kk],
                   PCA_YK$Psi_joint_k[W_perm_joint[idx_joint_PCA,],W_perm_joint[idx_joint_PCA,],kk])))) / (r_J*(2*KK))
          
          Con_Psi_joint_Unfitted <- Con_Psi_joint_Unfitted + sum(abs(diag(
            congru(Unfitted_YK$Psi_joint_hat[W_perm_joint[idx_joint_GRIDY,],W_perm_joint[idx_joint_GRIDY,],kk],
                   Unfitted_YK$Psi_joint_k[W_perm_joint[idx_joint_GRIDY,],W_perm_joint[idx_joint_GRIDY,],kk])))) / (r_J*(2*KK))
          
          if (kk <= KK){
            Con_Psi_group1_GRIDY <- Con_Psi_group1_GRIDY + sum(abs(diag(
              congru(GRIDY_YK$Psi_indiv_hat[W_perm_group1[idx_group1_GRIDY,],W_perm_group1[idx_group1_GRIDY,],kk],
                     GRIDY_YK$Psi_indiv_k[W_perm_group1[idx_group1_GRIDY,],W_perm_group1[idx_group1_GRIDY,],kk])))) / (r_I*KK)
            
            Con_Psi_group1_DSCA <- Con_Psi_group1_DSCA + sum(abs(diag(
              congru(DSCA_YK$Psi_indiv_hat[W_perm_group1[idx_group1_DSCA,],W_perm_group1[idx_group1_DSCA,],kk],
                     DSCA_YK$Psi_indiv_k[W_perm_group1[idx_group1_DSCA,],W_perm_group1[idx_group1_DSCA,],kk])))) / (r_I*KK)
            
            Con_Psi_group1_PCA <- Con_Psi_group1_PCA + sum(abs(diag(
              congru(PCA_YK$Psi_indiv_hat[W_perm_group1[idx_group1_PCA,],W_perm_group1[idx_group1_PCA,],kk],
                     PCA_YK$Psi_indiv_k[W_perm_group1[idx_group1_PCA,],W_perm_group1[idx_group1_PCA,],kk])))) / (r_I*KK)
            
            Con_Psi_group1_Unfitted <- Con_Psi_group1_Unfitted + sum(abs(diag(
              congru(Unfitted_YK$Psi_indiv_hat[W_perm_group1[idx_group1_GRIDY,],W_perm_group1[idx_group1_GRIDY,],kk],
                     Unfitted_YK$Psi_indiv_k[W_perm_group1[idx_group1_GRIDY,],W_perm_group1[idx_group1_GRIDY,],kk])))) / (r_I*KK)
          }else{
            Con_Psi_group2_GRIDY <- Con_Psi_group2_GRIDY + sum(abs(diag(
              congru(GRIDY_YK$Psi_indiv_hat[W_perm_group2[idx_group2_GRIDY,],W_perm_group2[idx_group2_GRIDY,],kk],
                     GRIDY_YK$Psi_indiv_k[W_perm_group2[idx_group2_GRIDY,],W_perm_group2[idx_group2_GRIDY,],kk])))) / (r_I*KK)
            
            Con_Psi_group2_DSCA <- Con_Psi_group2_DSCA + sum(abs(diag(
              congru(DSCA_YK$Psi_indiv_hat[W_perm_group2[idx_group2_DSCA,],W_perm_group2[idx_group2_DSCA,],kk],
                     DSCA_YK$Psi_indiv_k[W_perm_group2[idx_group2_DSCA,],W_perm_group2[idx_group2_DSCA,],kk])))) / (r_I*KK)
            
            Con_Psi_group2_PCA <- Con_Psi_group2_PCA + sum(abs(diag(
              congru(PCA_YK$Psi_indiv_hat[W_perm_group2[idx_group2_PCA,],W_perm_group2[idx_group2_PCA,],kk],
                     PCA_YK$Psi_indiv_k[W_perm_group2[idx_group2_PCA,],W_perm_group2[idx_group2_PCA,],kk])))) / (r_I*KK)
            
            Con_Psi_group2_Unfitted <- Con_Psi_group2_Unfitted + sum(abs(diag(
              congru(Unfitted_YK$Psi_indiv_hat[W_perm_group2[idx_group2_GRIDY,],W_perm_group2[idx_group2_GRIDY,],kk],
                     Unfitted_YK$Psi_indiv_k[W_perm_group2[idx_group2_GRIDY,],W_perm_group2[idx_group2_GRIDY,],kk])))) / (r_I*KK)
          }
        }
        list_GRIDY[["CC_Psi"]][["Joint"]][[iter]] <- Con_Psi_joint_GRIDY
        list_DSCA[["CC_Psi"]][["Joint"]][[iter]] <- Con_Psi_joint_DSCA
        list_PCA[["CC_Psi"]][["Joint"]][[iter]] <- Con_Psi_joint_PCA
        list_Unfitted[["CC_Psi"]][["Joint"]][[iter]] <- Con_Psi_joint_Unfitted
        
        list_GRIDY[["CC_Psi"]][["Group1"]][[iter]] <- Con_Psi_group1_GRIDY
        list_DSCA[["CC_Psi"]][["Group1"]][[iter]] <- Con_Psi_group1_DSCA
        list_PCA[["CC_Psi"]][["Group1"]][[iter]] <- Con_Psi_group1_PCA
        list_Unfitted[["CC_Psi"]][["Group1"]][[iter]] <- Con_Psi_group1_Unfitted
        
        list_GRIDY[["CC_Psi"]][["Group2"]][[iter]] <- Con_Psi_group2_GRIDY
        list_DSCA[["CC_Psi"]][["Group2"]][[iter]] <- Con_Psi_group2_DSCA
        list_PCA[["CC_Psi"]][["Group2"]][[iter]] <- Con_Psi_group2_PCA
        list_Unfitted[["CC_Psi"]][["Group2"]][[iter]] <- Con_Psi_group2_Unfitted
        
        
        # MSFE: GRIDY, DSCA, PCA, Unfitted
        X_hat_GRIDY <- array(NA,dim=c(H,dd,(2*KK)))
        F_hat_GRIDY <- array(NA,dim=c(H,(r_J+r_I),(2*KK)))
        
        X_hat_DSCA <- array(NA,dim=c(H,dd,(2*KK)))
        F_hat_DSCA <- array(NA,dim=c(H,(r_J+r_I),(2*KK)))
        
        X_hat_PCA <- array(NA,dim=c(H,dd,(2*KK)))
        F_hat_PCA <- array(NA,dim=c(H,(r_J+r_I),(2*KK)))
        
        X_hat_Unfitted <- array(NA,dim=c(H,dd,(2*KK)))
        F_hat_Unfitted <- array(NA,dim=c(H,(r_J+r_I),(2*KK)))
        
        for (kk in 1:(2*KK)){
          if (kk <= KK){
            Lambda_hat_GRIDY <- cbind(GRIDY_PF2_Joint$B,GRIDY_PF2_Group1$B)
            Lambda_hat_DSCA <- cbind(DSCA_Joint$B,DSCA_Group1$B)
            Lambda_hat_PCA <- cbind(PCA_Joint$B,PCA_Group1$B)
            Lambda_hat_Unfitted <- cbind(GRIDY_PF2_Joint$B,GRIDY_PF2_Group1$B)
            
            if (r_J == 1){
              F_T_GRIDY <- matrix( c(GRIDY_refitted$factor_joint[[kk]][TT],GRIDY_refitted$factor_group1[[kk]][TT,]), ncol=1)
              F_T_DSCA <- matrix( c(DSCA_refitted$factor_joint[[kk]][TT],DSCA_refitted$factor_group1[[kk]][TT,]), ncol=1)
              F_T_PCA <- matrix( c(PCA_refitted$factor_joint[[kk]][TT],PCA_refitted$factor_group1[[kk]][TT,]), ncol=1)
              F_T_Unfitted <- matrix( c(GRIDY_PF2_Joint$D[[kk]][TT],GRIDY_PF2_Group1$D[[kk]][TT,]), ncol=1)
            }else{
              F_T_GRIDY <- matrix( c(GRIDY_refitted$factor_joint[[kk]][TT,],GRIDY_refitted$factor_group1[[kk]][TT,]), ncol=1)
              F_T_DSCA <- matrix( c(DSCA_refitted$factor_joint[[kk]][TT,],DSCA_refitted$factor_group1[[kk]][TT,]), ncol=1)
              F_T_PCA <- matrix( c(PCA_refitted$factor_joint[[kk]][TT,],PCA_refitted$factor_group1[[kk]][TT,]), ncol=1)
              F_T_Unfitted <- matrix( c(GRIDY_PF2_Joint$D[[kk]][TT,],GRIDY_PF2_Group1$D[[kk]][TT,]), ncol=1)
            }
            
          }else{
            Lambda_hat_GRIDY <- cbind(GRIDY_PF2_Joint$B,GRIDY_PF2_Group2$B)
            Lambda_hat_DSCA <- cbind(DSCA_Joint$B,DSCA_Group2$B)
            Lambda_hat_PCA <- cbind(PCA_Joint$B,PCA_Group2$B)
            Lambda_hat_Unfitted <- cbind(GRIDY_PF2_Joint$B,GRIDY_PF2_Group2$B)
            
            if (r_J == 1){
              F_T_GRIDY <- matrix( c(GRIDY_refitted$factor_joint[[kk]][TT],GRIDY_refitted$factor_group2[[(kk-KK)]][TT,]), ncol=1)
              F_T_DSCA <- matrix( c(DSCA_refitted$factor_joint[[kk]][TT],DSCA_refitted$factor_group2[[(kk-KK)]][TT,]), ncol=1)
              F_T_PCA <- matrix( c(PCA_refitted$factor_joint[[kk]][TT],PCA_refitted$factor_group2[[(kk-KK)]][TT,]), ncol=1)
              F_T_Unfitted <- matrix( c(GRIDY_PF2_Joint$D[[kk]][TT],GRIDY_PF2_Group2$D[[(kk-KK)]][TT,]), ncol=1)
            }else{
              F_T_GRIDY <- matrix( c(GRIDY_refitted$factor_joint[[kk]][TT,],GRIDY_refitted$factor_group2[[(kk-KK)]][TT,]), ncol=1)
              F_T_DSCA <- matrix( c(DSCA_refitted$factor_joint[[kk]][TT,],DSCA_refitted$factor_group2[[(kk-KK)]][TT,]), ncol=1)
              F_T_PCA <- matrix( c(PCA_refitted$factor_joint[[kk]][TT,],PCA_refitted$factor_group2[[(kk-KK)]][TT,]), ncol=1)
              F_T_Unfitted <- matrix( c(GRIDY_PF2_Joint$D[[kk]][TT,],GRIDY_PF2_Group2$D[[(kk-KK)]][TT,]), ncol=1)
            }
            
          }
          Psi_hat_aug_GRIDY <- as.matrix(bdiag(GRIDY_YK$Psi_joint_hat[,,kk],GRIDY_YK$Psi_indiv_hat[,,kk]))
          Psi_hat_aug_DSCA <- as.matrix(bdiag(DSCA_YK$Psi_joint_hat[,,kk],DSCA_YK$Psi_indiv_hat[,,kk]))
          Psi_hat_aug_PCA <- as.matrix(bdiag(PCA_YK$Psi_joint_hat[,,kk],PCA_YK$Psi_indiv_hat[,,kk]))
          Psi_hat_aug_Unfitted <- as.matrix(bdiag(Unfitted_YK$Psi_joint_hat[,,kk],Unfitted_YK$Psi_indiv_hat[,,kk]))
          
          for (h in 1:H){
            if (h == 1){
              F_hat_GRIDY[h,,kk] <- t(Psi_hat_aug_GRIDY %*% F_T_GRIDY)
              F_hat_DSCA[h,,kk] <- t(Psi_hat_aug_DSCA %*% F_T_DSCA)
              F_hat_PCA[h,,kk] <- t(Psi_hat_aug_PCA %*% F_T_PCA)
              F_hat_Unfitted[h,,kk] <- t(Psi_hat_aug_Unfitted %*% F_T_Unfitted)
              
            }else{
              F_hat_GRIDY[h,,kk] <- t(Psi_hat_aug_GRIDY %*% as.matrix(F_hat_GRIDY[(h-1),,kk],ncol=1))
              F_hat_DSCA[h,,kk] <- t(Psi_hat_aug_DSCA %*% as.matrix(F_hat_DSCA[(h-1),,kk],ncol=1))
              F_hat_PCA[h,,kk] <- t(Psi_hat_aug_PCA %*% as.matrix(F_hat_PCA[(h-1),,kk],ncol=1))
              F_hat_Unfitted[h,,kk] <- t(Psi_hat_aug_Unfitted %*% as.matrix(F_hat_Unfitted[(h-1),,kk],ncol=1))
              
            }
            X_hat_GRIDY[h,,kk] <- t(Lambda_hat_GRIDY %*% F_hat_GRIDY[h,,kk])
            X_hat_DSCA[h,,kk] <- t(Lambda_hat_DSCA %*% F_hat_DSCA[h,,kk])
            X_hat_PCA[h,,kk] <- t(Lambda_hat_PCA %*% F_hat_PCA[h,,kk])
            X_hat_Unfitted[h,,kk] <- t(Lambda_hat_Unfitted %*% F_hat_Unfitted[h,,kk])
          }
        }
        
        MSFE_GRIDY <- vector("numeric",length=H)
        MSFE_DSCA <- vector("numeric",length=H)
        MSFE_PCA <- vector("numeric",length=H)
        MSFE_Unfitted <- vector("numeric",length=H)
        for (h in 1:H){
          for (kk in 1:(2*KK)){
            MSFE_GRIDY[h] <- norm(model_dgp_ext$scaled_X[[kk]][,(TT+1):(TT+h)] - t(X_hat_GRIDY[(1:h),,kk]),"2")^2/(2*dd*KK)
            MSFE_DSCA[h] <- norm(model_dgp_ext$scaled_X[[kk]][,(TT+1):(TT+h)] - t(X_hat_DSCA[(1:h),,kk]),"2")^2/(2*dd*KK)
            MSFE_PCA[h] <- norm(model_dgp_ext$scaled_X[[kk]][,(TT+1):(TT+h)] - t(X_hat_PCA[(1:h),,kk]),"2")^2/(2*dd*KK)
            MSFE_Unfitted[h] <- norm(model_dgp_ext$scaled_X[[kk]][,(TT+1):(TT+h)] - t(X_hat_Unfitted[(1:h),,kk]),"2")^2/(2*dd*KK)
          }
        }
        
        list_GRIDY[["MSFE"]][[iter]] <- MSFE_GRIDY
        list_DSCA[["MSFE"]][[iter]] <- MSFE_DSCA
        list_PCA[["MSFE"]][[iter]] <- MSFE_PCA
        list_Unfitted[["MSFE"]][[iter]] <- MSFE_Unfitted
        
        cat(iter,"th iteration","time taken is",Sys.time() - time,
            "and parameter set is",set,"r_J is",r_J,"r_I is",r_I,
            "d is",dd,"T is",TT," 2K is",(2*KK),"\n")
      }
      
      file_name <- paste0("main_param",set,
                          "_rJ",r_J,"rI",r_I,
                          "_d",dd,"T",TT,"K",KK,".RData")
      save( list_GRIDY,list_DSCA,list_PCA,list_Unfitted,
            file=file_name)
      
    }
  }
}




-----------------------------------------------------------------------------#
# Initial rank selection
-----------------------------------------------------------------------------#
d_list <- c(100,200,400,
            100,100)
T_list <- c(100,100,100,
            200,400)
K_list <- c(50,50,50,
            50,50)
joint_rank <- c(1,1,1,2,1,2,2,3,3,4)
indiv_rank <- c(1,2,3,2,4,3,4,3,4,4)

for (case in 1:5){

  dd <- d_list[case]
  TT <- T_list[case]
  KK <- K_list[case]

  for (set in 1:5){

    set.seed(1234)
    if (set == 1){
      Control_param <- Control_param_1
    }else if(set == 2){
      Control_param <- Control_param_2
    }else if(set ==3){
      Control_param <- Control_param_3
    }else if(set ==4){
      Control_param <- Control_param_4
    }else if(set ==5){
      Control_param <- Control_param_5
    }
    cat("Currently, set",set,"is set.\n")

    r_hat <- array(NA,dim=c((2*KK),length(joint_rank)))

    for (r in 1:10){
      cat("Currently, rank combination",r,"th has been performed.\n")
      time <- Sys.time()

      r_J <- joint_rank[r]
      r_I <- indiv_rank[r]
      model_dgp <- sim_model(r_J, r_I, dd, TT, KK, Control_param)

      for (kk in 1:(2*KK)){
        possibleError <- tryCatch(
          r_hat[kk,r] <- Rotational_bootstrap(dd,TT,model_dgp$Xk[[kk]],cull=0.5),

          error = function(kk){
            print(kk,"th iteration is failed to do Boostrap")
          }
        )
      }
      cat(r,"th iteration completed.","Time taken is",Sys.time() - time,
          "and parameter set is type",set,"r_J is",r_J,"r_I is",r_I,
          "d is",dd,"T is",TT," 2K is",(2*KK),"\n")
    }
    file_name <- paste0("bootstrap_param",set,"_d",dd,"T",TT,"K",KK,".RData")
    save(r_hat, file=file_name)
  }
}




