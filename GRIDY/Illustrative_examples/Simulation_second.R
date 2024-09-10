#-----------------------------------------------------------------------------#
#   File name : Simulation_second.R    												  
#
#   Project : "Group integrative dynamic factor models 
#             with application to multiple subject brain connectivity"
#
#   Maintainer : Younghoon Kim                     
#
#   Date : Sep. 1st, 2024
#
#   Purpose : code for producing results for the second simulation setting in 
#             Illustrative example section. Codes for plotting Figures 3 and 4 
#             are provided separately
#
#   R version 4.0.5 (2021-03-31)                                    
#
#   Input data file : ---- 
# 
#   Output data file : /Illustrative_examples/result/sim2_snr-_type-.RData
#
#   Required R packages : combinat_0.0-8, mvtnorm_1.2-5, and Matrix_1.5-1.
#-----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#
# Running simulation
# -----------------------------------------------------------------------------#
# Load rotational bootstrap from the parent directory
source(paste0(dirname(getwd()),"/","gica.R")) # For SCA and GICA
source(paste0(dirname(getwd()),"/","AJIVE_retrieve.R"))
source(paste0(dirname(getwd()),"/","library_simulation.R"))


# Problem size
dd <- 100
TT <- 200
KK <- 50

# Control parameters
sigma_eps <- 1
r_J <- 2
r_G <- 2

sigma_c_list <- c(0.25,0.5,0.75,1.25,2,4)

type_list <- c(1,2)
sigma_xi_list <- c(0.2,0.3)

# Simulation begins
max_iter <- 100
for (case in 1:6){
  
  # ---------------------------------------------------------------- #
  # Fix the strength of the signals:
  # ---------------------------------------------------------------- #
  sigma_c <- sigma_c_list[case]
  
  for (set in 1:2){
    
    type <- type_list[set]
    sigma_xi <- sigma_xi_list[set]
    
    Control_param <- data.frame(type=type,sigma_xi=sigma_xi,J_per=0.5,
                                b_bar_min=0,b_bar_max=1,
                                b_tilde_min=0,b_tilde_max=1,
                                sigma_c=sigma_c,sigma_eps=sigma_eps)
    
    # ---------------------------------------------------------------- #
    # Preparing containers:
    # ---------------------------------------------------------------- #
    list_GRIDY <- list(); list_SCA_P <- list () ;list_GICA <- list()
    list_DSCA <- list (); list_DGICA <- list ()
    list_Unfitted_SCA_PF2 <- list(); list_Unfitted_GICA <- list()
    
    cat("#------------------------------------#\n")
    cat("Currently case",case," set",set," is performing.\n")
    cat("#------------------------------------#\n")
    
    set.seed(1234)
    for (iter in 1:max_iter){
      start_time <- Sys.time()
      
      H <- 3
      model_dgp_ext <- sim_model(r_J, r_G, dd, (TT+H), KK, Control_param)
      
      model_dgp <- list()
      for (kk in 1:(2*KK)){
        model_dgp$Xk[[kk]] <- model_dgp_ext$Xk[[kk]][1:TT,]
        model_dgp$X_bar[[kk]] <- model_dgp_ext$X_bar[[kk]][1:TT,]
        model_dgp$X_tilde[[kk]] <- model_dgp_ext$X_tilde[[kk]][1:TT,]
        model_dgp$scaled_X[[kk]] <- model_dgp_ext$scaled_X[[kk]][,1:TT]
        model_dgp$scaled_X_comm[[kk]] <- model_dgp_ext$scaled_X_comm[[kk]][,1:TT]
        model_dgp$scaled_X_G_ind[[kk]] <- model_dgp_ext$scaled_X_G_ind[[kk]][,1:TT]
        model_dgp$Ek[[kk]] <- model_dgp_ext$Ek[[kk]][1:TT,]
        
        model_dgp$F_bar[[kk]] <- model_dgp_ext$F_bar[[kk]][1:TT,]
        model_dgp$A_bar[[kk]] <- model_dgp_ext$A_bar[[kk]][1:TT,]
        
        model_dgp$F_tilde[[kk]] <- model_dgp_ext$F_tilde[[kk]][1:TT,]
        model_dgp$A_tilde[[kk]] <- model_dgp_ext$A_tilde[[kk]][1:TT,]
      }
      model_dgp$C_bar <- model_dgp_ext$C_bar
      model_dgp$C_tilde <- model_dgp_ext$C_tilde
      model_dgp$Psi_bar <- model_dgp_ext$Psi_bar
      model_dgp$Psi_tilde <- model_dgp_ext$Psi_tilde
      model_dgp$B_bar <- model_dgp_ext$B_bar
      model_dgp$B_tilde1 <- model_dgp_ext$B_tilde1
      model_dgp$B_tilde2 <- model_dgp_ext$B_tilde2
      model_dgp$comm_idx <- model_dgp_ext$comm_idx
      model_dgp$group1_idx <- model_dgp_ext$group1_idx
      model_dgp$group2_idx <-model_dgp_ext$group2_idx
      #-----------------------------------------------------------------#
      # 1. Run AJIVE (GRIDY,SCA_P,GICA):
      #-----------------------------------------------------------------#
      AJIVE_step <- ajive(model_dgp$scaled_X,
                          initial_signal_ranks = rep((r_J+r_G),(2*KK)),
                          n_wedin_samples = 1000,
                          n_rand_dir_samples = 1000,
                          full = TRUE, joint_rank = r_J)
      
      Joint_block <- list()
      Group1_block <- list()
      Group2_block <- list()
      for (kk in 1:(2*KK)){
        Joint_block[[kk]] <- t(AJIVE_step[["block_decomps"]][[kk]]$joint$full)
        if (kk <= KK){
          Group1_block[[kk]] <- t(AJIVE_step[["block_decomps"]][[kk]]$individual$full)
        }else{
          Group2_block[[(kk-KK)]] <- t(AJIVE_step[["block_decomps"]][[kk]]$individual$full)
        }
      }
      
      #-----------------------------------------------------------------#
      # 2. Run SCA_PF2, SCA_P and GICA:
      #-----------------------------------------------------------------#
      # SCA-PF2:
      # Run SCA for joint structure:
      SCA_PF2_Joint <- sca(Joint_block, nfac = r_J,type="sca-pf2", verbose = FALSE)
      # Run SCA for 1st group structure:
      SCA_PF2_Group1 <- sca(Group1_block, nfac = r_G,type="sca-pf2", verbose = FALSE)
      # Run SCA for 2nd group structure:
      SCA_PF2_Group2 <- sca(Group2_block, nfac = r_G,type="sca-pf2", verbose = FALSE)
      
      # SCA-P:
      # Run SCA for joint structure:
      SCA_P_Joint <- sca(Joint_block, nfac = r_J,type="sca-p", rotation="varimax",verbose = FALSE)
      # Run SCA for 1st group structure:
      SCA_P_Group1 <- sca(Group1_block, nfac = r_G,type="sca-p", rotation="varimax",verbose = FALSE)
      # Run SCA for 2nd group structure:
      SCA_P_Group2 <- sca(Group2_block, nfac = r_G,type="sca-p", rotation="varimax",verbose = FALSE)
      
      # GICA:
      # Run SCA for joint structure:
      GICA_Joint <- gica(Joint_block, nc = r_J, dual.reg = FALSE, center = FALSE)
      # Run SCA for 1st group structure:
      GICA_Group1 <- gica(Group1_block, nc = r_G, dual.reg = FALSE, center = FALSE)
      # Run SCA for 2nd group structure:
      GICA_Group2 <- gica(Group2_block, nc = r_G, dual.reg = FALSE, center = FALSE)
      
      #-----------------------------------------------------------------#
      # 3. Construct dynamics (GRIDY,SCA_P,GICA):
      #-----------------------------------------------------------------#
      # GRIDY:
      # Refit factor series & Run Yule-Walker:
      GRIDY_refitted <- factor_regression(KK,r_J,r_G,model_dgp,
                                          SCA_PF2_Joint,
                                          SCA_PF2_Group1,
                                          SCA_PF2_Group2)
      GRIDY_YK <- YW_compute(dd,TT,KK,r_J,r_G,GRIDY_refitted)
      
      # SCA-P:
      # Refit factor series & Run Yule-Walker:
      SCA_P_refitted <- factor_regression(KK,r_J,r_G,model_dgp,
                                          SCA_P_Joint,
                                          SCA_P_Group1,
                                          SCA_P_Group2)
      SCA_P_YK <- YW_compute(dd,TT,KK,r_J,r_G,SCA_P_refitted)
      
      # GICA:
      # Refit factor series & Run Yule-Walker:
      GICA_refitted <- factor_regression(KK,r_J,r_G,model_dgp,
                                         GICA_Joint,
                                         GICA_Group1,
                                         GICA_Group2)
      GICA_YK <- YW_compute(dd,TT,KK,r_J,r_G,GICA_refitted)
      
      #-----------------------------------------------------------------#
      # 4. Construct dynamics (Unfitted_SCA_PF2,Unfitted_GICA):
      #-----------------------------------------------------------------#
      # Unfitted_SCA_PF2:
      # Run Yule-Walker:
      Factor_unfitted_SCA_PF2 <- list()
      Factor_unfitted_SCA_PF2$factor_joint <- SCA_PF2_Joint$D
      Factor_unfitted_SCA_PF2$factor_group1 <- SCA_PF2_Group1$D
      Factor_unfitted_SCA_PF2$factor_group2 <- SCA_PF2_Group2$D
      
      Unfitted_SCA_PF2_YK <- YW_compute(dd,TT,KK,r_J,r_G,Factor_unfitted_SCA_PF2)
      
      # Unfitted_GICA:
      # Run Yule-Walker:
      Factor_unfitted_GICA <- list()
      Factor_unfitted_GICA$factor_joint <- GICA_Joint$D
      Factor_unfitted_GICA$factor_group1 <- GICA_Group1$D
      Factor_unfitted_GICA$factor_group2 <- GICA_Group2$D
      Unfitted_GICA_YK <- YW_compute(dd,TT,KK,r_J,r_G,Factor_unfitted_GICA)
      
      #-----------------------------------------------------------------#
      # 5. Double SCA and Double GICA:
      #-----------------------------------------------------------------#
      Whole_block <- list()
      
      for (kk in 1:(2*KK)){
        Whole_block[[kk]] <- t(model_dgp$scaled_X[[kk]])
      }
      
      DSCA_Joint <- sca(Whole_block, nfac = r_J,type="sca-pf2", verbose = FALSE)
      DGICA_Joint <- gica(Whole_block, nc = r_J, dual.reg = FALSE, center = FALSE)
      
      Rem_DSCA_block_group1 <- list(); Rem_DSCA_block_group2 <- list()
      Rem_DGICA_block_group1 <- list(); Rem_DGICA_block_group2 <- list()
      for (kk in 1:(2*KK)){
        kth_Joint_block_DSCA <- DSCA_Joint$D[[kk]] %*% t(DSCA_Joint$B)
        kth_Joint_block_DGICA <- DGICA_Joint$D[[kk]] %*% t(DGICA_Joint$B)
        if (kk <= KK){
          Rem_DSCA_block_group1[[kk]] <- Whole_block[[kk]] - kth_Joint_block_DSCA
          Rem_DGICA_block_group1[[kk]] <- Whole_block[[kk]] - kth_Joint_block_DGICA
        }else{
          Rem_DSCA_block_group2[[(kk-KK)]] <- Whole_block[[kk]] - kth_Joint_block_DSCA
          Rem_DGICA_block_group2[[(kk-KK)]] <- Whole_block[[kk]] - kth_Joint_block_DGICA
        }
      }
      
      DSCA_Group1 <- sca(Rem_DSCA_block_group1, nfac = r_G,type="sca-pf2", verbose = FALSE)
      DSCA_Group2 <- sca(Rem_DSCA_block_group2, nfac = r_G,type="sca-pf2", verbose = FALSE)
      
      DGICA_Group1 <- gica(Rem_DGICA_block_group1, nc = r_G, dual.reg = FALSE, center = FALSE)
      DGICA_Group2 <- gica(Rem_DGICA_block_group2, nc = r_G, dual.reg = FALSE, center = FALSE)
      
      DSCA_refitted <- factor_regression(KK,r_J,r_G,model_dgp,
                                         DSCA_Joint,
                                         DSCA_Group1,
                                         DSCA_Group2)
      DSCA_YK <- YW_compute(dd,TT,KK,r_J,r_G,DSCA_refitted)
      
      DGICA_refitted <- factor_regression(KK,r_J,r_G,model_dgp,
                                          DGICA_Joint,
                                          DGICA_Group1,
                                          DGICA_Group2)
      DGICA_YK <- YW_compute(dd,TT,KK,r_J,r_G,DGICA_refitted)
      
      #-----------------------------------------------------------------#
      # 4. Evaluation:
      #-----------------------------------------------------------------#
      #-----------------------------------------------------------------#
      # R2: 
      #-----------------------------------------------------------------#
      R2_GRIDY_Joint <- 0; R2_SCA_P_Joint <- 0; R2_GICA_Joint <- 0; 
      R2_DSCA_Joint <- 0; R2_DGICA_Joint <- 0;
      R2_Unfitted_SCA_PF2_Joint <- 0; R2_Unfitted_GICA_Joint <- 0;
      
      R2_GRIDY_Group1 <- 0; R2_SCA_P_Group1 <- 0; R2_GICA_Group1 <- 0; 
      R2_DSCA_Group1 <- 0; R2_DGICA_Group1 <- 0;
      R2_Unfitted_SCA_PF2_Group1 <- 0; R2_Unfitted_GICA_Group1 <- 0;
      
      R2_GRIDY_Group2 <- 0; R2_SCA_P_Group2 <- 0; R2_GICA_Group2 <- 0; 
      R2_DSCA_Group2 <- 0; R2_DGICA_Group2 <- 0;
      R2_Unfitted_SCA_PF2_Group2 <- 0; R2_Unfitted_GICA_Group2 <- 0;
      
      for (kk in 1:(2*KK)){
        total_variation_k <- norm(t(model_dgp$scaled_X[[kk]]),"F")^2
        
        # GRIDY:
        R2_GRIDY_Joint <- R2_GRIDY_Joint + 1- norm(t(model_dgp$scaled_X[[kk]]) - GRIDY_refitted$factor_joint[[kk]] %*% t(SCA_PF2_Joint$B),"F")^2/total_variation_k
        # SCA_P:
        R2_SCA_P_Joint <- R2_SCA_P_Joint + 1- norm(t(model_dgp$scaled_X[[kk]]) - SCA_P_refitted$factor_joint[[kk]] %*% t(SCA_P_Joint$B),"F")^2/total_variation_k
        # GICA:
        R2_GICA_Joint <- R2_GICA_Joint + 1- norm(t(model_dgp$scaled_X[[kk]]) - GICA_refitted$factor_joint[[kk]] %*% t(GICA_Joint$B),"F")^2/total_variation_k
        # DSCA:
        R2_DSCA_Joint <- R2_DSCA_Joint + 1- norm(t(model_dgp$scaled_X[[kk]]) - DSCA_refitted$factor_joint[[kk]] %*% t(DSCA_Joint$B),"F")^2/total_variation_k
        # DGICA:
        R2_DGICA_Joint <- R2_DGICA_Joint + 1- norm(t(model_dgp$scaled_X[[kk]])  - DGICA_refitted$factor_joint[[kk]] %*% t(DGICA_Joint$B),"F")^2/total_variation_k
        # Unfitted_SCA_PF2:
        R2_Unfitted_SCA_PF2_Joint <- R2_Unfitted_SCA_PF2_Joint + 1- norm(t(model_dgp$scaled_X[[kk]])  - Factor_unfitted_SCA_PF2$factor_joint[[kk]] %*% t(SCA_PF2_Joint$B),"F")^2/total_variation_k
        # Unfitted_GICA:
        R2_Unfitted_GICA_Joint <- R2_Unfitted_GICA_Joint + 1- norm(t(model_dgp$scaled_X[[kk]])  - Factor_unfitted_GICA$factor_joint[[kk]] %*% t(GICA_Joint$B),"F")^2/total_variation_k
        
        if (kk <= KK){
          # GRIDY:
          R2_GRIDY_Group1 <- R2_GRIDY_Group1 + 1- norm(t(model_dgp$scaled_X[[kk]])  - GRIDY_refitted$factor_group1[[kk]] %*% t(SCA_PF2_Group1$B) ,"F")^2/total_variation_k
          # SCA_P:
          R2_SCA_P_Group1 <- R2_SCA_P_Group1 + 1- norm(t(model_dgp$scaled_X[[kk]]) - SCA_P_refitted$factor_group1[[kk]] %*% t(SCA_P_Group1$B) ,"F")^2/total_variation_k
          # GICA:
          R2_GICA_Group1 <- R2_GICA_Group1 + 1- norm(t(model_dgp$scaled_X[[kk]]) - GICA_refitted$factor_group1[[kk]] %*% t(GICA_Group1$B) ,"F")^2/total_variation_k
          # DSCA:
          R2_DSCA_Group1 <- R2_DSCA_Group1 + 1- norm(t(model_dgp$scaled_X[[kk]]) - DSCA_refitted$factor_group1[[kk]] %*% t(DSCA_Group1$B) ,"F")^2/total_variation_k
          # DGICA:
          R2_DGICA_Group1 <- R2_DGICA_Group1 + 1- norm(t(model_dgp$scaled_X[[kk]])  - DGICA_refitted$factor_group1[[kk]] %*% t(DGICA_Group1$B) ,"F")^2/total_variation_k
          # Unfitted_SCA_PF2:
          R2_Unfitted_SCA_PF2_Group1 <- R2_Unfitted_SCA_PF2_Group1 + 1- norm(t(model_dgp$scaled_X[[kk]])  - Factor_unfitted_SCA_PF2$factor_group1[[kk]] %*% t(SCA_PF2_Group1$B) ,"F")^2/total_variation_k
          # Unfitted_GICA:
          R2_Unfitted_GICA_Group1 <- R2_Unfitted_GICA_Group1 + 1- norm(t(model_dgp$scaled_X[[kk]])  - Factor_unfitted_GICA$factor_group1[[kk]] %*% t(GICA_Group1$B) ,"F")^2/total_variation_k
          
        }else{
          # GRIDY:
          R2_GRIDY_Group2 <- R2_GRIDY_Group2 + 1- norm(t(model_dgp$scaled_X[[kk]])  - GRIDY_refitted$factor_group2[[(kk-KK)]] %*% t(SCA_PF2_Group2$B) ,"F")^2/total_variation_k
          # SCA_P:
          R2_SCA_P_Group2 <- R2_SCA_P_Group2 + 1- norm(t(model_dgp$scaled_X[[kk]])  - SCA_P_refitted$factor_group2[[(kk-KK)]] %*% t(SCA_P_Group2$B) ,"F")^2/total_variation_k
          # GICA:
          R2_GICA_Group2 <- R2_GICA_Group2 + 1- norm(t(model_dgp$scaled_X[[kk]])  - GICA_refitted$factor_group2[[(kk-KK)]] %*% t(GICA_Group2$B) ,"F")^2/total_variation_k
          # DSCA:
          R2_DSCA_Group2 <- R2_DSCA_Group2 + 1- norm(t(model_dgp$scaled_X[[kk]]) - DSCA_refitted$factor_group2[[(kk-KK)]] %*% t(DSCA_Group2$B) ,"F")^2/total_variation_k
          # DGICA:
          R2_DGICA_Group2 <- R2_DGICA_Group2 + 1- norm(t(model_dgp$scaled_X[[kk]]) - DGICA_refitted$factor_group2[[(kk-KK)]] %*% t(DGICA_Group2$B) ,"F")^2/total_variation_k
          # Unfitted_SCA_PF2:
          R2_Unfitted_SCA_PF2_Group2 <- R2_Unfitted_SCA_PF2_Group2 + 1- norm(t(model_dgp$scaled_X[[kk]]) - Factor_unfitted_SCA_PF2$factor_group2[[(kk-KK)]] %*% t(SCA_PF2_Group2$B) ,"F")^2/total_variation_k
          # Unfitted_GICA:    
          R2_Unfitted_GICA_Group2 <- R2_Unfitted_GICA_Group2 + 1- norm(t(model_dgp$scaled_X[[kk]]) - Factor_unfitted_GICA$factor_group2[[(kk-KK)]] %*% t(GICA_Group2$B) ,"F")^2/total_variation_k
        }
      }
      
      list_GRIDY[["R2"]][["Joint"]][[iter]] <- R2_GRIDY_Joint / (2*KK)
      list_SCA_P[["R2"]][["Joint"]][[iter]] <- R2_SCA_P_Joint / (2*KK)
      list_GICA[["R2"]][["Joint"]][[iter]] <- R2_GICA_Joint / (2*KK)
      list_DSCA[["R2"]][["Joint"]][[iter]] <- R2_DSCA_Joint / (2*KK)
      list_DGICA[["R2"]][["Joint"]][[iter]] <- R2_DGICA_Joint / (2*KK)
      list_Unfitted_SCA_PF2[["R2"]][["Joint"]][[iter]] <- R2_Unfitted_SCA_PF2_Joint / (2*KK)
      list_Unfitted_GICA[["R2"]][["Joint"]][[iter]] <- R2_Unfitted_GICA_Joint / (2*KK)
      
      list_GRIDY[["R2"]][["Group1"]][[iter]] <- R2_GRIDY_Group1 / (KK)
      list_SCA_P[["R2"]][["Group1"]][[iter]] <- R2_SCA_P_Group1 / (KK)
      list_GICA[["R2"]][["Group1"]][[iter]] <- R2_GICA_Group1 / (KK)
      list_DSCA[["R2"]][["Group1"]][[iter]] <- R2_DSCA_Group1 / (KK)
      list_DGICA[["R2"]][["Group1"]][[iter]] <- R2_DGICA_Group1 / (KK)
      list_Unfitted_SCA_PF2[["R2"]][["Group1"]][[iter]] <- R2_Unfitted_SCA_PF2_Group1 / (KK)
      list_Unfitted_GICA[["R2"]][["Group1"]][[iter]] <- R2_Unfitted_GICA_Group1 / (KK)
      
      list_GRIDY[["R2"]][["Group2"]][[iter]] <- R2_GRIDY_Group2 / (KK)
      list_SCA_P[["R2"]][["Group2"]][[iter]] <- R2_SCA_P_Group2 / (KK)
      list_GICA[["R2"]][["Group2"]][[iter]] <- R2_GICA_Group2 / (KK)
      list_DSCA[["R2"]][["Group2"]][[iter]] <- R2_DSCA_Group2 / (KK)
      list_DGICA[["R2"]][["Group2"]][[iter]] <- R2_DGICA_Group2 / (KK)
      list_Unfitted_SCA_PF2[["R2"]][["Group2"]][[iter]] <- R2_Unfitted_SCA_PF2_Group2 / (KK)
      list_Unfitted_GICA[["R2"]][["Group2"]][[iter]] <- R2_Unfitted_GICA_Group2 / (KK)
      
      #-----------------------------------------------------------------#
      # RMSE:
      #-----------------------------------------------------------------#
      RMSE_GRIDY_Joint <- 0; RMSE_SCA_P_Joint <- 0; RMSE_GICA_Joint <- 0; 
      RMSE_DSCA_Joint <- 0; RMSE_DGICA_Joint <- 0;
      RMSE_Unfitted_SCA_PF2_Joint <- 0; RMSE_Unfitted_GICA_Joint <- 0;
      
      RMSE_GRIDY_Group1 <- 0; RMSE_SCA_P_Group1 <- 0; RMSE_GICA_Group1 <- 0; 
      RMSE_DSCA_Group1 <- 0; RMSE_DGICA_Group1 <- 0;
      RMSE_Unfitted_SCA_PF2_Group1 <- 0; RMSE_Unfitted_GICA_Group1 <- 0;
      
      RMSE_GRIDY_Group2 <- 0; RMSE_SCA_P_Group2 <- 0; RMSE_GICA_Group2 <- 0; 
      RMSE_DSCA_Group2 <- 0; RMSE_DGICA_Group2 <- 0;
      RMSE_Unfitted_SCA_PF2_Group2 <- 0; RMSE_Unfitted_GICA_Group2 <- 0;
      
      for (kk in 1:(2*KK)){
        # GRIDY:
        RMSE_GRIDY_Joint <- RMSE_GRIDY_Joint + norm(t(model_dgp$scaled_X_comm[[kk]]) - GRIDY_refitted$factor_joint[[kk]] %*% t(SCA_PF2_Joint$B), "F")^2 / (TT*dd*(2*KK))
        # SCA_P:
        RMSE_SCA_P_Joint <- RMSE_SCA_P_Joint + norm(t(model_dgp$scaled_X_comm[[kk]])  - SCA_P_refitted$factor_joint[[kk]] %*% t(SCA_P_Joint$B), "F")^2 / (TT*dd*(2*KK))
        # GICA:
        RMSE_GICA_Joint <- RMSE_GICA_Joint + norm(t(model_dgp$scaled_X_comm[[kk]]) - GICA_refitted$factor_joint[[kk]] %*% t(GICA_Joint$B), "F")^2 / (TT*dd*(2*KK))
        # DSCA:
        RMSE_DSCA_Joint <- RMSE_DSCA_Joint + norm(t(model_dgp$scaled_X_comm[[kk]]) - DSCA_refitted$factor_joint[[kk]] %*% t(DSCA_Joint$B), "F")^2 / (TT*dd*(2*KK))
        # DGICA:
        RMSE_DGICA_Joint <- RMSE_DGICA_Joint + norm(t(model_dgp$scaled_X_comm[[kk]]) - DGICA_refitted$factor_joint[[kk]] %*% t(DGICA_Joint$B), "F")^2 / (TT*dd*(2*KK))
        # Unfitted_SCA_PF2:
        RMSE_Unfitted_SCA_PF2_Joint <- RMSE_Unfitted_SCA_PF2_Joint + norm(t(model_dgp$scaled_X_comm[[kk]]) - Factor_unfitted_SCA_PF2$factor_joint[[kk]] %*% t(SCA_PF2_Joint$B), "F")^2 / (TT*dd*(2*KK))
        # Unfitted_GICA:
        RMSE_Unfitted_GICA_Joint <- RMSE_Unfitted_GICA_Joint + norm(t(model_dgp$scaled_X_comm[[kk]]) - Factor_unfitted_GICA$factor_joint[[kk]] %*% t(GICA_Joint$B), "F")^2 / (TT*dd*(2*KK))
        
        if (kk <= KK){
          # GRIDY:
          RMSE_GRIDY_Group1 <- RMSE_GRIDY_Group1 + norm(t(model_dgp$scaled_X_G_ind[[kk]]) - GRIDY_refitted$factor_group1[[kk]] %*% t(SCA_PF2_Group1$B), "F")^2 / (TT*dd*(KK))
          # SCA_P:
          RMSE_SCA_P_Group1 <- RMSE_SCA_P_Group1 + norm(t(model_dgp$scaled_X_G_ind[[kk]]) - SCA_P_refitted$factor_group1[[kk]] %*% t(SCA_P_Group1$B), "F")^2 / (TT*dd*(KK))
          # GICA:
          RMSE_GICA_Group1 <- RMSE_GICA_Group1 + norm(t(model_dgp$scaled_X_G_ind[[kk]]) - GICA_refitted$factor_group1[[kk]] %*% t(GICA_Group1$B), "F")^2 / (TT*dd*(KK))
          # DSCA:
          RMSE_DSCA_Group1 <- RMSE_DSCA_Group1 + norm(t(model_dgp$scaled_X_G_ind[[kk]]) - DSCA_refitted$factor_group1[[kk]] %*% t(DSCA_Group1$B), "F")^2 / (TT*dd*(KK))
          # DGICA:
          RMSE_DGICA_Group1 <- RMSE_DGICA_Group1 + norm(t(model_dgp$scaled_X_G_ind[[kk]]) - DGICA_refitted$factor_group1[[kk]] %*% t(DGICA_Group1$B), "F")^2 / (TT*dd*(KK))
          # Unfitted_SCA_PF2:
          RMSE_Unfitted_SCA_PF2_Group1 <- RMSE_Unfitted_SCA_PF2_Group1 + norm(t(model_dgp$scaled_X_G_ind[[kk]]) - Factor_unfitted_SCA_PF2$factor_group1[[kk]] %*% t(SCA_PF2_Group1$B), "F")^2 / (TT*dd*(KK))
          # Unfitted_GICA:
          RMSE_Unfitted_GICA_Group1 <- RMSE_Unfitted_GICA_Group1 + norm(t(model_dgp$scaled_X_G_ind[[kk]]) - Factor_unfitted_GICA$factor_group1[[kk]] %*% t(GICA_Group1$B), "F")^2 / (TT*dd*(KK))
          
        }else{
          # GRIDY:
          RMSE_GRIDY_Group2 <- RMSE_GRIDY_Group2 + norm(t(model_dgp$scaled_X_G_ind[[kk]]) - GRIDY_refitted$factor_group2[[(kk-KK)]] %*% t(SCA_PF2_Group2$B), "F")^2 / (TT*dd*(KK))
          # SCA_P:
          RMSE_SCA_P_Group2 <- RMSE_SCA_P_Group2 + norm(t(model_dgp$scaled_X_G_ind[[kk]])  - SCA_P_refitted$factor_group2[[(kk-KK)]] %*% t(SCA_P_Group2$B), "F")^2 / (TT*dd*(KK))
          # GICA:
          RMSE_GICA_Group2 <- RMSE_GICA_Group2 + norm(t(model_dgp$scaled_X_G_ind[[kk]]) - GICA_refitted$factor_group2[[(kk-KK)]] %*% t(GICA_Group2$B), "F")^2 / (TT*dd*(KK))
          # DSCA:
          RMSE_DSCA_Group2 <- RMSE_DSCA_Group2 + norm(t(model_dgp$scaled_X_G_ind[[kk]]) - DSCA_refitted$factor_group2[[(kk-KK)]] %*% t(DSCA_Group2$B), "F")^2 / (TT*dd*(KK))
          # DGICA:
          RMSE_DGICA_Group2 <- RMSE_DGICA_Group2 + norm(t(model_dgp$scaled_X_G_ind[[kk]]) - DGICA_refitted$factor_group2[[(kk-KK)]] %*% t(DGICA_Group2$B), "F")^2 / (TT*dd*(KK))
          # Unfitted_SCA_PF2:
          RMSE_Unfitted_SCA_PF2_Group2 <- RMSE_Unfitted_SCA_PF2_Group2 + norm(t(model_dgp$scaled_X_G_ind[[kk]]) - Factor_unfitted_SCA_PF2$factor_group2[[(kk-KK)]] %*% t(SCA_PF2_Group2$B), "F")^2 / (TT*dd*(KK))
          # Unfitted_GICA:
          RMSE_Unfitted_GICA_Group2 <- RMSE_Unfitted_GICA_Group2 + norm(t(model_dgp$scaled_X_G_ind[[kk]]) - Factor_unfitted_GICA$factor_group2[[(kk-KK)]] %*% t(GICA_Group2$B), "F")^2 / (TT*dd*(KK))
          
        }
      }
      
      list_GRIDY[["RMSE"]][["Joint"]][[iter]] <- sqrt(RMSE_GRIDY_Joint)
      list_SCA_P[["RMSE"]][["Joint"]][[iter]] <- sqrt(RMSE_SCA_P_Joint)
      list_GICA[["RMSE"]][["Joint"]][[iter]] <- sqrt(RMSE_GICA_Joint)
      list_DSCA[["RMSE"]][["Joint"]][[iter]] <- sqrt(RMSE_DSCA_Joint)
      list_DGICA[["RMSE"]][["Joint"]][[iter]] <- sqrt(RMSE_DGICA_Joint)
      list_Unfitted_SCA_PF2[["RMSE"]][["Joint"]][[iter]] <- sqrt(RMSE_Unfitted_SCA_PF2_Joint)
      list_Unfitted_GICA[["RMSE"]][["Joint"]][[iter]] <- sqrt(RMSE_Unfitted_GICA_Joint)
      
      list_GRIDY[["RMSE"]][["Group1"]][[iter]] <- sqrt(RMSE_GRIDY_Group1)
      list_SCA_P[["RMSE"]][["Group1"]][[iter]] <- sqrt(RMSE_SCA_P_Group1)
      list_GICA[["RMSE"]][["Group1"]][[iter]] <- sqrt(RMSE_GICA_Group1)
      list_DSCA[["RMSE"]][["Group1"]][[iter]] <- sqrt(RMSE_DSCA_Group1)
      list_DGICA[["RMSE"]][["Group1"]][[iter]] <- sqrt(RMSE_DGICA_Group1)
      list_Unfitted_SCA_PF2[["RMSE"]][["Group1"]][[iter]] <- sqrt(RMSE_Unfitted_SCA_PF2_Group1)
      list_Unfitted_GICA[["RMSE"]][["Group1"]][[iter]] <- sqrt(RMSE_Unfitted_GICA_Group1)
      
      list_GRIDY[["RMSE"]][["Group2"]][[iter]] <- sqrt(RMSE_GRIDY_Group2)
      list_SCA_P[["RMSE"]][["Group2"]][[iter]] <- sqrt(RMSE_SCA_P_Group2)
      list_GICA[["RMSE"]][["Group2"]][[iter]] <- sqrt(RMSE_GICA_Group2)
      list_DSCA[["RMSE"]][["Group2"]][[iter]] <- sqrt(RMSE_DSCA_Group2)
      list_DGICA[["RMSE"]][["Group2"]][[iter]] <- sqrt(RMSE_DGICA_Group2)
      list_Unfitted_SCA_PF2[["RMSE"]][["Group2"]][[iter]] <- sqrt(RMSE_Unfitted_SCA_PF2_Group2)
      list_Unfitted_GICA[["RMSE"]][["Group2"]][[iter]] <- sqrt(RMSE_Unfitted_GICA_Group2)
      
      #-----------------------------------------------------------------#
      # CC_B: 
      #-----------------------------------------------------------------#
      # Joint loadings:
      W_perm_joint <- matrix(unlist(permn(r_J)),nrow=factorial(r_J),byrow=TRUE)
      Con_B_joint_GRIDY <- rep(NA, nrow(W_perm_joint))
      Con_B_joint_SCA_P <- rep(NA, nrow(W_perm_joint))
      Con_B_joint_GICA <- rep(NA, nrow(W_perm_joint))
      Con_B_joint_DSCA <- rep(NA, nrow(W_perm_joint))
      Con_B_joint_DGICA <- rep(NA, nrow(W_perm_joint))
      
      for(i in 1:nrow(W_perm_joint)){
        Con_B_joint_GRIDY[i] <- sum(abs(diag(as.matrix(congru(SCA_PF2_Joint$B[,W_perm_joint[i,]],
                                                              model_dgp$B_bar)))))/r_J
        Con_B_joint_SCA_P[i] <- sum(abs(diag(as.matrix(congru(SCA_P_Joint$B[,W_perm_joint[i,]],
                                                              model_dgp$B_bar)))))/r_J
        Con_B_joint_GICA[i] <- sum(abs(diag(as.matrix(congru(GICA_Joint$B[,W_perm_joint[i,]],
                                                             model_dgp$B_bar)))))/r_J
        Con_B_joint_DSCA[i] <- sum(abs(diag(as.matrix(congru(DSCA_Joint$B[,W_perm_joint[i,]],
                                                             model_dgp$B_bar)))))/r_J
        Con_B_joint_DGICA[i] <- sum(abs(diag(as.matrix(congru(DGICA_Joint$B[,W_perm_joint[i,]],
                                                              model_dgp$B_bar)))))/r_J
      }
      idx_joint_GRIDY <- which.min(Con_B_joint_GRIDY)
      Con_B_joint_GRIDY <- Con_B_joint_GRIDY[idx_joint_GRIDY]
      idx_joint_SCA_P <- which.min(Con_B_joint_SCA_P)
      Con_B_joint_SCA_P <- Con_B_joint_SCA_P[idx_joint_SCA_P]
      idx_joint_GICA <- which.min(Con_B_joint_GICA)
      Con_B_joint_GICA <- Con_B_joint_GICA[idx_joint_GICA]
      idx_joint_DSCA <- which.min(Con_B_joint_DSCA)
      Con_B_joint_DSCA <- Con_B_joint_DSCA[idx_joint_DSCA]
      idx_joint_DGICA <- which.min(Con_B_joint_DGICA)
      Con_B_joint_DGICA <- Con_B_joint_DGICA[idx_joint_DGICA]
      Con_B_joint_Unfitted_SCA_PF2 <- Con_B_joint_GRIDY
      Con_B_joint_Unfitted_GICA <- Con_B_joint_GICA
      
      
      # Group 1 individual loadings:
      W_perm_group1 <- matrix(unlist(permn(r_G)),nrow=factorial(r_G),byrow=TRUE)
      Con_B_group1_GRIDY <- rep(NA, nrow(W_perm_group1))
      Con_B_group1_SCA_P <- rep(NA, nrow(W_perm_group1))
      Con_B_group1_GICA <- rep(NA, nrow(W_perm_group1))
      Con_B_group1_DSCA <- rep(NA, nrow(W_perm_group1))
      Con_B_group1_DGICA <- rep(NA, nrow(W_perm_group1))
      
      for(i in 1:nrow(W_perm_group1)){
        Con_B_group1_GRIDY[i] <- sum(abs(diag(as.matrix(congru(SCA_PF2_Group1$B[,W_perm_group1[i,]],
                                                               model_dgp$B_tilde1)))))/r_G
        Con_B_group1_SCA_P[i] <- sum(abs(diag(as.matrix(congru(SCA_P_Group1$B[,W_perm_group1[i,]],
                                                               model_dgp$B_tilde1)))))/r_G
        Con_B_group1_GICA[i] <- sum(abs(diag(as.matrix(congru(GICA_Group1$B[,W_perm_group1[i,]],
                                                              model_dgp$B_tilde1)))))/r_G
        Con_B_group1_DSCA[i] <- sum(abs(diag(as.matrix(congru(DSCA_Group1$B[,W_perm_group1[i,]],
                                                              model_dgp$B_tilde1)))))/r_G
        Con_B_group1_DGICA[i] <- sum(abs(diag(as.matrix(congru(DGICA_Group1$B[,W_perm_group1[i,]],
                                                               model_dgp$B_tilde1)))))/r_G
      }
      idx_group1_GRIDY <- which.max(Con_B_group1_GRIDY)
      Con_B_group1_GRIDY <- Con_B_group1_GRIDY[idx_group1_GRIDY]
      idx_group1_SCA_P <- which.max(Con_B_group1_SCA_P)
      Con_B_group1_SCA_P <- Con_B_group1_SCA_P[idx_group1_SCA_P]
      idx_group1_GICA <- which.max(Con_B_group1_GICA)
      Con_B_group1_GICA <- Con_B_group1_GICA[idx_group1_GICA]
      idx_group1_DSCA <- which.max(Con_B_group1_DSCA)
      Con_B_group1_DSCA <- Con_B_group1_DSCA[idx_group1_DSCA]
      idx_group1_DGICA <- which.max(Con_B_group1_DGICA)
      Con_B_group1_DGICA <- Con_B_group1_DGICA[idx_group1_DGICA]
      Con_B_group1_Unfitted_SCA_PF2 <- Con_B_group1_GRIDY
      Con_B_group1_Unfitted_GICA <- Con_B_group1_GICA
      
      # Group 2 individual loadings:
      W_perm_group2 <- matrix(unlist(permn(r_G)),nrow=factorial(r_G),byrow=TRUE)
      Con_B_group2_GRIDY <- rep(NA, nrow(W_perm_group2))
      Con_B_group2_SCA_P <- rep(NA, nrow(W_perm_group2))
      Con_B_group2_GICA <- rep(NA, nrow(W_perm_group2))
      Con_B_group2_DSCA <- rep(NA, nrow(W_perm_group2))
      Con_B_group2_DGICA <- rep(NA, nrow(W_perm_group2))
      
      for(i in 1:nrow(W_perm_group2)){
        Con_B_group2_GRIDY[i] <- sum(abs(diag(as.matrix(congru(SCA_PF2_Group2$B[,W_perm_group2[i,]],
                                                               model_dgp$B_tilde2)))))/r_G
        Con_B_group2_SCA_P[i] <- sum(abs(diag(as.matrix(congru(SCA_P_Group2$B[,W_perm_group2[i,]],
                                                               model_dgp$B_tilde2)))))/r_G
        Con_B_group2_GICA[i] <- sum(abs(diag(as.matrix(congru(GICA_Group2$B[,W_perm_group2[i,]],
                                                              model_dgp$B_tilde2)))))/r_G
        Con_B_group2_DSCA[i] <- sum(abs(diag(as.matrix(congru(DSCA_Group2$B[,W_perm_group2[i,]],
                                                              model_dgp$B_tilde2)))))/r_G
        Con_B_group2_DGICA[i] <- sum(abs(diag(as.matrix(congru(DGICA_Group2$B[,W_perm_group2[i,]],
                                                               model_dgp$B_tilde2)))))/r_G
      }
      idx_group2_GRIDY <- which.max(Con_B_group2_GRIDY)
      Con_B_group2_GRIDY <- Con_B_group2_GRIDY[idx_group2_GRIDY]
      idx_group2_SCA_P <- which.max(Con_B_group2_SCA_P)
      Con_B_group2_SCA_P <- Con_B_group2_SCA_P[idx_group2_SCA_P]
      idx_group2_GICA <- which.max(Con_B_group2_GICA)
      Con_B_group2_GICA <- Con_B_group2_GICA[idx_group2_GICA]
      idx_group2_DSCA <- which.max(Con_B_group2_DSCA)
      Con_B_group2_DSCA <- Con_B_group2_DSCA[idx_group2_DSCA]
      idx_group2_DGICA <- which.max(Con_B_group2_DGICA)
      Con_B_group2_DGICA <- Con_B_group2_DGICA[idx_group2_DGICA]
      Con_B_group2_Unfitted_SCA_PF2 <- Con_B_group2_GRIDY
      Con_B_group2_Unfitted_GICA <- Con_B_group2_GICA
      
      
      list_GRIDY[["CC_B"]][["Joint"]][[iter]] <- Con_B_joint_GRIDY
      list_SCA_P[["CC_B"]][["Joint"]][[iter]] <- Con_B_joint_SCA_P
      list_GICA[["CC_B"]][["Joint"]][[iter]] <- Con_B_joint_GICA
      list_DSCA[["CC_B"]][["Joint"]][[iter]] <- Con_B_joint_DSCA
      list_DGICA[["CC_B"]][["Joint"]][[iter]] <- Con_B_joint_DGICA
      list_Unfitted_SCA_PF2[["CC_B"]][["Joint"]][[iter]] <- Con_B_joint_GRIDY
      list_Unfitted_GICA[["CC_B"]][["Joint"]][[iter]] <- Con_B_joint_GICA
      
      list_GRIDY[["CC_B"]][["Group1"]][[iter]] <- Con_B_group1_GRIDY
      list_SCA_P[["CC_B"]][["Group1"]][[iter]] <- Con_B_group1_SCA_P
      list_GICA[["CC_B"]][["Group1"]][[iter]] <- Con_B_group1_GICA
      list_DSCA[["CC_B"]][["Group1"]][[iter]] <- Con_B_group1_DSCA
      list_DGICA[["CC_B"]][["Group1"]][[iter]] <- Con_B_group1_DGICA
      list_Unfitted_SCA_PF2[["CC_B"]][["Group1"]][[iter]] <- Con_B_group1_GRIDY
      list_Unfitted_GICA[["CC_B"]][["Group1"]][[iter]] <- Con_B_group1_GICA
      
      list_GRIDY[["CC_B"]][["Group2"]][[iter]] <- Con_B_group2_GRIDY
      list_SCA_P[["CC_B"]][["Group2"]][[iter]] <- Con_B_group2_SCA_P
      list_GICA[["CC_B"]][["Group2"]][[iter]] <- Con_B_group2_GICA
      list_DSCA[["CC_B"]][["Group2"]][[iter]] <- Con_B_group2_DSCA
      list_DGICA[["CC_B"]][["Group2"]][[iter]] <- Con_B_group2_DGICA
      list_Unfitted_SCA_PF2[["CC_B"]][["Group2"]][[iter]] <- Con_B_group2_GRIDY
      list_Unfitted_GICA[["CC_B"]][["Group2"]][[iter]] <- Con_B_group2_GICA
      
      #-----------------------------------------------------------------#
      # CC_F: 
      #-----------------------------------------------------------------#
      Reorder_joint_GRIDY <- lapply(X=c(1:(2*KK)),FUN=function(X) 
        GRIDY_refitted$factor_joint[[X]][,W_perm_joint[idx_joint_GRIDY,]])
      Reorder_joint_SCA_P <- lapply(X=c(1:(2*KK)),FUN=function(X) 
        SCA_P_refitted$factor_joint[[X]][,W_perm_joint[idx_joint_SCA_P,]])
      Reorder_joint_GICA <- lapply(X=c(1:(2*KK)),FUN=function(X) 
        GICA_refitted$factor_joint[[X]][,W_perm_joint[idx_joint_GICA,]])
      Reorder_joint_DSCA <- lapply(X=c(1:(2*KK)),FUN=function(X) 
        DSCA_refitted$factor_joint[[X]][,W_perm_joint[idx_joint_DSCA,]])
      Reorder_joint_DGICA <- lapply(X=c(1:(2*KK)),FUN=function(X) 
        DGICA_refitted$factor_joint[[X]][,W_perm_joint[idx_joint_DGICA,]])
      Reorder_joint_Unfitted_SCA_PF2 <- lapply(X=c(1:(2*KK)),FUN=function(X) 
        Factor_unfitted_SCA_PF2$factor_joint[[X]][,W_perm_joint[idx_joint_GRIDY,]])
      Reorder_joint_Unfitted_GICA <- lapply(X=c(1:(2*KK)),FUN=function(X) 
        Factor_unfitted_GICA$factor_joint[[X]][,W_perm_joint[idx_joint_GICA,]])
      
      Reorder_group1_GRIDY <- lapply(X=c(1:KK),FUN=function(X) 
        GRIDY_refitted$factor_group1[[X]][,W_perm_group1[idx_group1_GRIDY,]])
      Reorder_group1_SCA_P <- lapply(X=c(1:KK),FUN=function(X) 
        SCA_P_refitted$factor_group1[[X]][,W_perm_group1[idx_group1_SCA_P,]])
      Reorder_group1_GICA <- lapply(X=c(1:KK),FUN=function(X) 
        GICA_refitted$factor_group1[[X]][,W_perm_group1[idx_group1_GICA,]])
      Reorder_group1_DSCA <- lapply(X=c(1:KK),FUN=function(X) 
        DSCA_refitted$factor_group1[[X]][,W_perm_group1[idx_group1_DSCA,]])
      Reorder_group1_DGICA <- lapply(X=c(1:KK),FUN=function(X) 
        DGICA_refitted$factor_group1[[X]][,W_perm_group1[idx_group1_DGICA,]])
      Reorder_group1_Unfitted_SCA_PF2 <- lapply(X=c(1:KK),FUN=function(X) 
        Factor_unfitted_SCA_PF2$factor_group1[[X]][,W_perm_group1[idx_group1_GRIDY,]])
      Reorder_group1_Unfitted_GICA <- lapply(X=c(1:KK),FUN=function(X) 
        Factor_unfitted_GICA$factor_group1[[X]][,W_perm_group1[idx_group1_GICA,]])
      
      Reorder_group2_GRIDY <- lapply(X=c(1:KK),FUN=function(X) 
        GRIDY_refitted$factor_group2[[X]][,W_perm_group2[idx_group2_GRIDY,]])
      Reorder_group2_SCA_P <- lapply(X=c(1:KK),FUN=function(X) 
        SCA_P_refitted$factor_group2[[X]][,W_perm_group2[idx_group2_SCA_P,]])
      Reorder_group2_GICA <- lapply(X=c(1:KK),FUN=function(X) 
        GICA_refitted$factor_group2[[X]][,W_perm_group2[idx_group2_GICA,]])
      Reorder_group2_DSCA <- lapply(X=c(1:KK),FUN=function(X) 
        DSCA_refitted$factor_group2[[X]][,W_perm_group2[idx_group2_DSCA,]])
      Reorder_group2_DGICA <- lapply(X=c(1:KK),FUN=function(X) 
        DGICA_refitted$factor_group2[[X]][,W_perm_group2[idx_group2_DGICA,]])
      Reorder_group2_Unfitted_SCA_PF2 <- lapply(X=c(1:KK),FUN=function(X) 
        Factor_unfitted_SCA_PF2$factor_group2[[X]][,W_perm_group2[idx_group2_GRIDY,]])
      Reorder_group2_Unfitted_GICA <- lapply(X=c(1:KK),FUN=function(X) 
        Factor_unfitted_GICA$factor_group2[[X]][,W_perm_group2[idx_group2_GICA,]])
      
      
      Con_F_joint_GRIDY <- 0; Con_F_joint_SCA_P <- 0; Con_F_joint_GICA <- 0; 
      Con_F_joint_DSCA <- 0; Con_F_joint_DGICA <- 0;
      Con_F_joint_Unfitted_SCA_PF2 <- 0; Con_F_joint_Unfitted_GICA <- 0;
      
      Con_F_group1_GRIDY <- 0; Con_F_group1_SCA_P <- 0; Con_F_group1_GICA <- 0; 
      Con_F_group1_DSCA <- 0; Con_F_group1_DGICA <- 0;
      Con_F_group1_Unfitted_SCA_PF2 <- 0; Con_F_group1_Unfitted_GICA <- 0;
      
      Con_F_group2_GRIDY <- 0; Con_F_group2_SCA_P <- 0; Con_F_group2_GICA <- 0; 
      Con_F_group2_DSCA <- 0; Con_F_group2_DGICA <- 0;
      Con_F_group2_Unfitted_SCA_PF2 <- 0; Con_F_group2_Unfitted_GICA <- 0;
      
      for (kk in 1:(2*KK)) {
        # GRIDY:
        Con_F_joint_GRIDY <- Con_F_joint_GRIDY + sum(abs(diag(matrix(congru(Reorder_joint_GRIDY[[kk]],
                                                                            model_dgp$F_bar[[kk]]),nrow=r_J))))/ (r_J*(2*KK))
        # SCA_P:
        Con_F_joint_SCA_P <- Con_F_joint_SCA_P + sum(abs(diag(matrix(congru(Reorder_joint_SCA_P[[kk]],
                                                                            model_dgp$F_bar[[kk]]),nrow=r_J))))/ (r_J*(2*KK))
        # GICA:
        Con_F_joint_GICA <- Con_F_joint_GICA + sum(abs(diag(matrix(congru(Reorder_joint_GICA[[kk]],
                                                                          model_dgp$F_bar[[kk]]),nrow=r_J))))/ (r_J*(2*KK))
        # DSCA:
        Con_F_joint_DSCA <- Con_F_joint_DSCA + sum(abs(diag(matrix(congru(Reorder_joint_DSCA[[kk]],
                                                                          model_dgp$F_bar[[kk]]),nrow=r_J))))/ (r_J*(2*KK))
        # DGICA:
        Con_F_joint_DGICA <- Con_F_joint_DGICA + sum(abs(diag(matrix(congru(Reorder_joint_DGICA[[kk]],
                                                                            model_dgp$F_bar[[kk]]),nrow=r_J))))/ (r_J*(2*KK))
        # Unfitted_SCA_PF2:
        Con_F_joint_Unfitted_SCA_PF2 <- Con_F_joint_Unfitted_SCA_PF2 + sum(abs(diag(matrix(congru(Reorder_joint_Unfitted_SCA_PF2[[kk]],
                                                                                                  model_dgp$F_bar[[kk]]),nrow=r_J))))/ (r_J*(2*KK))
        # Unfitted_GICA:
        Con_F_joint_Unfitted_GICA <- Con_F_joint_Unfitted_GICA + sum(abs(diag(matrix(congru(Reorder_joint_Unfitted_GICA[[kk]],
                                                                                            model_dgp$F_bar[[kk]]),nrow=r_J))))/ (r_J*(2*KK))
        
        if (kk <= KK){
          # GRIDY:
          Con_F_group1_GRIDY <- Con_F_group1_GRIDY + sum(abs(diag(matrix(congru(Reorder_group1_GRIDY[[kk]],
                                                                                model_dgp$F_tilde[[kk]]),nrow=r_G))))/ (r_G*KK)
          # SCA_P:
          Con_F_group1_SCA_P <- Con_F_group1_SCA_P + sum(abs(diag(matrix(congru(Reorder_group1_SCA_P[[kk]],
                                                                                model_dgp$F_tilde[[kk]]),nrow=r_G))))/ (r_G*KK)
          # GICA:
          Con_F_group1_GICA <- Con_F_group1_GICA + sum(abs(diag(matrix(congru(Reorder_group1_GICA[[kk]],
                                                                              model_dgp$F_tilde[[kk]]),nrow=r_G))))/ (r_G*KK)
          # DSCA:
          Con_F_group1_DSCA <- Con_F_group1_DSCA + sum(abs(diag(matrix(congru(Reorder_group1_DSCA[[kk]],
                                                                              model_dgp$F_tilde[[kk]]),nrow=r_G))))/ (r_G*KK)
          # DGICA:
          Con_F_group1_DGICA <- Con_F_group1_DGICA + sum(abs(diag(matrix(congru(Reorder_group1_DGICA[[kk]],
                                                                                model_dgp$F_tilde[[kk]]),nrow=r_G))))/ (r_G*KK)
          # Unfitted_SCA_PF2:
          Con_F_group1_Unfitted_SCA_PF2 <- Con_F_group1_Unfitted_SCA_PF2 + sum(abs(diag(matrix(congru(Reorder_group1_Unfitted_SCA_PF2[[kk]],
                                                                                                      model_dgp$F_tilde[[kk]]),nrow=r_G))))/ (r_G*KK)
          # Unfitted_GICA:
          Con_F_group1_Unfitted_GICA <- Con_F_group1_Unfitted_GICA + sum(abs(diag(matrix(congru(Reorder_group1_Unfitted_GICA[[kk]],
                                                                                                model_dgp$F_tilde[[kk]]),nrow=r_G))))/ (r_G*KK)
          
        }else{
          # GRIDY:
          Con_F_group2_GRIDY <- Con_F_group2_GRIDY + sum(abs(diag(matrix(congru(Reorder_group2_GRIDY[[(kk-KK)]],
                                                                                model_dgp$F_tilde[[kk]]),nrow=r_G))))/ (r_G*KK)
          # SCA_P:
          Con_F_group2_SCA_P <- Con_F_group2_SCA_P + sum(abs(diag(matrix(congru(Reorder_group2_SCA_P[[(kk-KK)]],
                                                                                model_dgp$F_tilde[[kk]]),nrow=r_G))))/ (r_G*KK)
          # GICA:
          Con_F_group2_GICA <- Con_F_group2_GICA + sum(abs(diag(matrix(congru(Reorder_group2_GICA[[(kk-KK)]],
                                                                              model_dgp$F_tilde[[kk]]),nrow=r_G))))/ (r_G*KK)
          # DSCA:
          Con_F_group2_DSCA <- Con_F_group2_DSCA + sum(abs(diag(matrix(congru(Reorder_group2_DSCA[[(kk-KK)]],
                                                                              model_dgp$F_tilde[[kk]]),nrow=r_G))))/ (r_G*KK)
          # DGICA:
          Con_F_group2_DGICA <- Con_F_group2_DGICA + sum(abs(diag(matrix(congru(Reorder_group2_DGICA[[(kk-KK)]],
                                                                                model_dgp$F_tilde[[kk]]),nrow=r_G))))/ (r_G*KK)
          # Unfitted_SCA_PF2:
          Con_F_group2_Unfitted_SCA_PF2 <- Con_F_group2_Unfitted_SCA_PF2 + sum(abs(diag(matrix(congru(Reorder_group2_Unfitted_SCA_PF2[[(kk-KK)]],
                                                                                                      model_dgp$F_tilde[[kk]]),nrow=r_G))))/ (r_G*KK)
          # Unfitted_GICA:
          Con_F_group2_Unfitted_GICA <- Con_F_group2_Unfitted_GICA + sum(abs(diag(matrix(congru(Reorder_group2_Unfitted_GICA[[(kk-KK)]],
                                                                                                model_dgp$F_tilde[[kk]]),nrow=r_G))))/ (r_G*KK)
          
        }
      }
      
      list_GRIDY[["CC_F"]][["Joint"]][[iter]] <- Con_F_joint_GRIDY
      list_SCA_P[["CC_F"]][["Joint"]][[iter]] <- Con_F_joint_SCA_P
      list_GICA[["CC_F"]][["Joint"]][[iter]] <- Con_F_joint_GICA
      list_DSCA[["CC_F"]][["Joint"]][[iter]] <- Con_F_joint_DSCA
      list_DGICA[["CC_F"]][["Joint"]][[iter]] <- Con_F_joint_DGICA
      list_Unfitted_SCA_PF2[["CC_F"]][["Joint"]][[iter]] <- Con_F_joint_Unfitted_SCA_PF2
      list_Unfitted_GICA[["CC_F"]][["Joint"]][[iter]] <- Con_F_joint_Unfitted_GICA
      
      list_GRIDY[["CC_F"]][["Group1"]][[iter]] <- Con_F_group1_GRIDY
      list_SCA_P[["CC_F"]][["Group1"]][[iter]] <- Con_F_group1_SCA_P
      list_GICA[["CC_F"]][["Group1"]][[iter]] <- Con_F_group1_GICA
      list_DSCA[["CC_F"]][["Group1"]][[iter]] <- Con_F_group1_DSCA
      list_DGICA[["CC_F"]][["Group1"]][[iter]] <- Con_F_group1_DGICA
      list_Unfitted_SCA_PF2[["CC_F"]][["Group1"]][[iter]] <- Con_F_group1_Unfitted_SCA_PF2
      list_Unfitted_GICA[["CC_F"]][["Group1"]][[iter]] <- Con_F_group1_Unfitted_GICA
      
      list_GRIDY[["CC_F"]][["Group2"]][[iter]] <- Con_F_group2_GRIDY
      list_SCA_P[["CC_F"]][["Group2"]][[iter]] <- Con_F_group2_SCA_P
      list_GICA[["CC_F"]][["Group2"]][[iter]] <- Con_F_group2_GICA
      list_DSCA[["CC_F"]][["Group2"]][[iter]] <- Con_F_group2_DSCA
      list_DGICA[["CC_F"]][["Group2"]][[iter]] <- Con_F_group2_DGICA
      list_Unfitted_SCA_PF2[["CC_F"]][["Group2"]][[iter]] <- Con_F_group2_Unfitted_SCA_PF2
      list_Unfitted_GICA[["CC_F"]][["Group2"]][[iter]] <- Con_F_group2_Unfitted_GICA
      
      #-----------------------------------------------------------------#
      # CC_Psi: 
      #-----------------------------------------------------------------#
      Con_Psi_joint_GRIDY <- 0; Con_Psi_joint_SCA_P <- 0; Con_Psi_joint_GICA <- 0; 
      Con_Psi_joint_DSCA <- 0; Con_Psi_joint_DGICA <- 0;
      Con_Psi_joint_Unfitted_SCA_PF2 <- 0; Con_Psi_joint_Unfitted_GICA <- 0;
      
      Con_Psi_group1_GRIDY <- 0; Con_Psi_group1_SCA_P <- 0; Con_Psi_group1_GICA <- 0; 
      Con_Psi_group1_DSCA <- 0; Con_Psi_group1_DGICA <- 0;
      Con_Psi_group1_Unfitted_SCA_PF2 <- 0; Con_Psi_group1_Unfitted_GICA <- 0;
      
      Con_Psi_group2_GRIDY <- 0; Con_Psi_group2_SCA_P <- 0; Con_Psi_group2_GICA <- 0; 
      Con_Psi_group2_DSCA <- 0; Con_Psi_group2_DGICA <- 0;
      Con_Psi_group2_Unfitted_SCA_PF2 <- 0; Con_Psi_group2_Unfitted_GICA <- 0;
      
      for (kk in 1:(2*KK)){
        # GRIDY:
        Con_Psi_joint_GRIDY <- Con_Psi_joint_GRIDY + abs(congru(as.vector(model_dgp$Psi_bar),
                                                                as.vector(GRIDY_YK$Psi_joint_hat[W_perm_joint[idx_joint_GRIDY,],
                                                                                                 W_perm_joint[idx_joint_GRIDY,],kk]))) / (2*KK)
        # SCA_P:
        Con_Psi_joint_SCA_P <- Con_Psi_joint_SCA_P + abs(congru(as.vector(model_dgp$Psi_bar),
                                                                as.vector(SCA_P_YK$Psi_joint_hat[W_perm_joint[idx_joint_SCA_P,],
                                                                                                 W_perm_joint[idx_joint_SCA_P,],kk]))) / (2*KK)
        # GICA:
        Con_Psi_joint_GICA <- Con_Psi_joint_GICA + abs(congru(as.vector(model_dgp$Psi_bar),
                                                              as.vector(GICA_YK$Psi_joint_hat[W_perm_joint[idx_joint_GICA,],
                                                                                              W_perm_joint[idx_joint_GICA,],kk]))) / (2*KK)
        # DSCA:
        Con_Psi_joint_DSCA <- Con_Psi_joint_DSCA + abs(congru(as.vector(model_dgp$Psi_bar),
                                                              as.vector(DSCA_YK$Psi_joint_hat[W_perm_joint[idx_joint_DSCA,],
                                                                                              W_perm_joint[idx_joint_DSCA,],kk]))) / (2*KK)
        # DGICA:
        Con_Psi_joint_DGICA <- Con_Psi_joint_DGICA + abs(congru(as.vector(model_dgp$Psi_bar),
                                                                as.vector(DGICA_YK$Psi_joint_hat[W_perm_joint[idx_joint_DGICA,],
                                                                                                 W_perm_joint[idx_joint_DGICA,],kk]))) / (2*KK)
        # Unfitted_SCA_PF2:
        Con_Psi_joint_Unfitted_SCA_PF2 <- Con_Psi_joint_Unfitted_SCA_PF2 + abs(congru(as.vector(model_dgp$Psi_bar),
                                                                                      as.vector(Unfitted_SCA_PF2_YK$Psi_joint_hat[W_perm_joint[idx_joint_GRIDY,],
                                                                                                                                  W_perm_joint[idx_joint_GRIDY,],kk]))) / (2*KK)
        # Unfitted_GICA:
        Con_Psi_joint_Unfitted_GICA <- Con_Psi_joint_Unfitted_GICA + abs(congru(as.vector(model_dgp$Psi_bar),
                                                                                as.vector(Unfitted_GICA_YK$Psi_joint_hat[W_perm_joint[idx_joint_GICA,],
                                                                                                                         W_perm_joint[idx_joint_GICA,],kk]))) / (2*KK)
        
        if (kk <= KK){
          # GRIDY:
          Con_Psi_group1_GRIDY <- Con_Psi_group1_GRIDY + abs(congru(as.vector(model_dgp$Psi_tilde),
                                                                    as.vector(GRIDY_YK$Psi_indiv_hat[W_perm_joint[idx_group1_GRIDY,],
                                                                                                     W_perm_joint[idx_group1_GRIDY,],kk]))) / KK
          
          # SCA_P:
          Con_Psi_group1_SCA_P <- Con_Psi_group1_SCA_P + abs(congru(as.vector(model_dgp$Psi_tilde),
                                                                    as.vector(SCA_P_YK$Psi_indiv_hat[W_perm_joint[idx_group1_SCA_P,],
                                                                                                     W_perm_joint[idx_group1_SCA_P,],kk]))) / KK
          # GICA:
          Con_Psi_group1_GICA <- Con_Psi_group1_GICA + abs(congru(as.vector(model_dgp$Psi_tilde),
                                                                  as.vector(GICA_YK$Psi_indiv_hat[W_perm_joint[idx_group1_GICA,],
                                                                                                  W_perm_joint[idx_group1_GICA,],kk]))) / KK
          # DSCA:
          Con_Psi_group1_DSCA <- Con_Psi_group1_DSCA + abs(congru(as.vector(model_dgp$Psi_tilde),
                                                                  as.vector(DSCA_YK$Psi_indiv_hat[W_perm_joint[idx_group1_DSCA,],
                                                                                                  W_perm_joint[idx_group1_DSCA,],kk]))) / KK
          # DGICA:
          Con_Psi_group1_DGICA <- Con_Psi_group1_DGICA + abs(congru(as.vector(model_dgp$Psi_tilde),
                                                                    as.vector(DGICA_YK$Psi_indiv_hat[W_perm_joint[idx_group1_DGICA,],
                                                                                                     W_perm_joint[idx_group1_DGICA,],kk]))) / KK
          # Unfitted_SCA_PF2:
          Con_Psi_group1_Unfitted_SCA_PF2 <- Con_Psi_group1_Unfitted_SCA_PF2 + abs(congru(as.vector(model_dgp$Psi_tilde),
                                                                                          as.vector(Unfitted_SCA_PF2_YK$Psi_indiv_hat[W_perm_joint[idx_group1_GRIDY,],
                                                                                                                                      W_perm_joint[idx_group1_GRIDY,],kk]))) / KK
          # Unfitted_GICA:
          Con_Psi_group1_Unfitted_GICA <- Con_Psi_group1_Unfitted_GICA + abs(congru(as.vector(model_dgp$Psi_tilde),
                                                                                    as.vector(Unfitted_GICA_YK$Psi_indiv_hat[W_perm_joint[idx_group1_GICA,],
                                                                                                                             W_perm_joint[idx_group1_GICA,],kk]))) / KK
          
        }else{
          # GRIDY:
          Con_Psi_group2_GRIDY <- Con_Psi_group2_GRIDY + abs(congru(as.vector(model_dgp$Psi_tilde),
                                                                    as.vector(GRIDY_YK$Psi_indiv_hat[W_perm_joint[idx_group2_GRIDY,],
                                                                                                     W_perm_joint[idx_group2_GRIDY,],kk]))) / KK
          
          # SCA_P:
          Con_Psi_group2_SCA_P <- Con_Psi_group2_SCA_P + abs(congru(as.vector(model_dgp$Psi_tilde),
                                                                    as.vector(SCA_P_YK$Psi_indiv_hat[W_perm_joint[idx_group2_SCA_P,],
                                                                                                     W_perm_joint[idx_group2_SCA_P,],kk]))) / KK
          # GICA:
          Con_Psi_group2_GICA <- Con_Psi_group2_GICA + abs(congru(as.vector(model_dgp$Psi_tilde),
                                                                  as.vector(GICA_YK$Psi_indiv_hat[W_perm_joint[idx_group2_GICA,],
                                                                                                  W_perm_joint[idx_group2_GICA,],kk]))) / KK
          # DSCA:
          Con_Psi_group2_DSCA <- Con_Psi_group2_DSCA + abs(congru(as.vector(model_dgp$Psi_tilde),
                                                                  as.vector(DSCA_YK$Psi_indiv_hat[W_perm_joint[idx_group2_DSCA,],
                                                                                                  W_perm_joint[idx_group2_DSCA,],kk]))) / KK
          # DGICA:
          Con_Psi_group2_DGICA <- Con_Psi_group2_DGICA + abs(congru(as.vector(model_dgp$Psi_tilde),
                                                                    as.vector(DGICA_YK$Psi_indiv_hat[W_perm_joint[idx_group2_DGICA,],
                                                                                                     W_perm_joint[idx_group2_DGICA,],kk]))) / KK
          # Unfitted_SCA_PF2:
          Con_Psi_group2_Unfitted_SCA_PF2 <- Con_Psi_group2_Unfitted_SCA_PF2 + abs(congru(as.vector(model_dgp$Psi_tilde),
                                                                                          as.vector(Unfitted_SCA_PF2_YK$Psi_indiv_hat[W_perm_joint[idx_group2_GRIDY,],
                                                                                                                                      W_perm_joint[idx_group2_GRIDY,],kk]))) / KK
          # Unfitted_GICA:
          Con_Psi_group2_Unfitted_GICA <- Con_Psi_group2_Unfitted_GICA + abs(congru(as.vector(model_dgp$Psi_tilde),
                                                                                    as.vector(Unfitted_GICA_YK$Psi_indiv_hat[W_perm_joint[idx_group2_GICA,],
                                                                                                                             W_perm_joint[idx_group2_GICA,],kk]))) / KK
        }
      }
      
      list_GRIDY[["CC_PSI"]][["Joint"]][[iter]] <- Con_Psi_joint_GRIDY
      list_SCA_P[["CC_PSI"]][["Joint"]][[iter]] <- Con_Psi_joint_SCA_P
      list_GICA[["CC_PSI"]][["Joint"]][[iter]] <- Con_Psi_joint_GICA
      list_DSCA[["CC_PSI"]][["Joint"]][[iter]] <- Con_Psi_joint_DSCA
      list_DGICA[["CC_PSI"]][["Joint"]][[iter]] <- Con_Psi_joint_DGICA
      list_Unfitted_SCA_PF2[["CC_PSI"]][["Joint"]][[iter]] <- Con_Psi_joint_Unfitted_SCA_PF2
      list_Unfitted_GICA[["CC_PSI"]][["Joint"]][[iter]] <- Con_Psi_joint_Unfitted_GICA
      
      list_GRIDY[["CC_PSI"]][["Group1"]][[iter]] <- Con_Psi_group1_GRIDY
      list_SCA_P[["CC_PSI"]][["Group1"]][[iter]] <- Con_Psi_group1_SCA_P
      list_GICA[["CC_PSI"]][["Group1"]][[iter]] <- Con_Psi_group1_GICA
      list_DSCA[["CC_PSI"]][["Group1"]][[iter]] <- Con_Psi_group1_DSCA
      list_DGICA[["CC_PSI"]][["Group1"]][[iter]] <- Con_Psi_group1_DGICA
      list_Unfitted_SCA_PF2[["CC_PSI"]][["Group1"]][[iter]] <- Con_Psi_group1_Unfitted_SCA_PF2
      list_Unfitted_GICA[["CC_PSI"]][["Group1"]][[iter]] <- Con_Psi_group1_Unfitted_GICA
      
      list_GRIDY[["CC_PSI"]][["Group2"]][[iter]] <- Con_Psi_group2_GRIDY
      list_SCA_P[["CC_PSI"]][["Group2"]][[iter]] <- Con_Psi_group2_SCA_P
      list_GICA[["CC_PSI"]][["Group2"]][[iter]] <- Con_Psi_group2_GICA
      list_DSCA[["CC_PSI"]][["Group2"]][[iter]] <- Con_Psi_group2_DSCA
      list_DGICA[["CC_PSI"]][["Group2"]][[iter]] <- Con_Psi_group2_DGICA
      list_Unfitted_SCA_PF2[["CC_PSI"]][["Group2"]][[iter]] <- Con_Psi_group2_Unfitted_SCA_PF2
      list_Unfitted_GICA[["CC_PSI"]][["Group2"]][[iter]] <- Con_Psi_group2_Unfitted_GICA
      
      #-----------------------------------------------------------------#
      # RMSFE: 
      #-----------------------------------------------------------------#
      X_hat_GRIDY <- array(NA,dim=c(H,dd,(2*KK))); F_hat_GRIDY <- array(NA,dim=c(H,(r_J+r_G),(2*KK)))
      X_hat_SCA_P <- array(NA,dim=c(H,dd,(2*KK))); F_hat_SCA_P <- array(NA,dim=c(H,(r_J+r_G),(2*KK)))
      X_hat_GICA <- array(NA,dim=c(H,dd,(2*KK))); F_hat_GICA <- array(NA,dim=c(H,(r_J+r_G),(2*KK)))
      
      X_hat_DSCA <- array(NA,dim=c(H,dd,(2*KK))); F_hat_DSCA <- array(NA,dim=c(H,(r_J+r_G),(2*KK)))
      X_hat_DGICA <- array(NA,dim=c(H,dd,(2*KK))); F_hat_DGICA <- array(NA,dim=c(H,(r_J+r_G),(2*KK)))
      
      X_hat_Unfitted_SCA_PF2 <- array(NA,dim=c(H,dd,(2*KK))); F_hat_Unfitted_SCA_PF2 <- array(NA,dim=c(H,(r_J+r_G),(2*KK)))
      X_hat_Unfitted_GICA <- array(NA,dim=c(H,dd,(2*KK))); F_hat_Unfitted_GICA <- array(NA,dim=c(H,(r_J+r_G),(2*KK)))
      
      for (kk in 1:(2*KK)){
        if (kk <= KK){
          Lambda_hat_GRIDY <- cbind(SCA_PF2_Joint$B,SCA_PF2_Group1$B)
          Lambda_hat_SCA_P <- cbind(SCA_P_Joint$B,SCA_P_Group1$B)
          Lambda_hat_GICA <- cbind(GICA_Joint$B,GICA_Group1$B)
          Lambda_hat_DSCA <- cbind(DSCA_Joint$B,DSCA_Group1$B)
          Lambda_hat_DGICA <- cbind(DGICA_Joint$B,DGICA_Group1$B)
          Lambda_hat_Unfitted_SCA_PF2 <- Lambda_hat_GRIDY
          Lambda_hat_Unfitted_GICA <- Lambda_hat_GICA 
          
          F_T_GRIDY <- matrix( c(GRIDY_refitted$factor_joint[[kk]][TT,],GRIDY_refitted$factor_group1[[kk]][TT,]), ncol=1)
          F_T_SCA_P <- matrix( c(SCA_P_refitted$factor_joint[[kk]][TT,],SCA_P_refitted$factor_group1[[kk]][TT,]), ncol=1)
          F_T_GICA <- matrix( c(GICA_refitted$factor_joint[[kk]][TT,],GICA_refitted$factor_group1[[kk]][TT,]), ncol=1)
          F_T_DSCA <- matrix( c(DSCA_refitted$factor_joint[[kk]][TT,],DSCA_refitted$factor_group1[[kk]][TT,]), ncol=1)
          F_T_DGICA <- matrix( c(DGICA_refitted$factor_joint[[kk]][TT,],DGICA_refitted$factor_group1[[kk]][TT,]), ncol=1)
          F_T_Unfitted_SCA_PF2 <- matrix( c(Factor_unfitted_SCA_PF2$factor_joint[[kk]][TT,],
                                            Factor_unfitted_SCA_PF2$factor_group1[[kk]][TT,]), ncol=1)
          F_T_Unfitted_GICA <- matrix( c(Factor_unfitted_GICA$factor_joint[[kk]][TT,],
                                         Factor_unfitted_GICA$factor_group1[[kk]][TT,]), ncol=1)
          
        }else{
          Lambda_hat_GRIDY <- cbind(SCA_PF2_Joint$B,SCA_PF2_Group2$B)
          Lambda_hat_SCA_P <- cbind(SCA_P_Joint$B,SCA_P_Group2$B)
          Lambda_hat_GICA <- cbind(GICA_Joint$B,GICA_Group2$B)
          Lambda_hat_DSCA <- cbind(DSCA_Joint$B,DSCA_Group2$B)
          Lambda_hat_DGICA <- cbind(DGICA_Joint$B,DGICA_Group2$B)
          Lambda_hat_Unfitted_SCA_PF2 <- Lambda_hat_GRIDY
          Lambda_hat_Unfitted_GICA <- Lambda_hat_GICA 
          
          F_T_GRIDY <- matrix( c(GRIDY_refitted$factor_joint[[kk]][TT,],GRIDY_refitted$factor_group2[[(kk-KK)]][TT,]), ncol=1)
          F_T_SCA_P <- matrix( c(SCA_P_refitted$factor_joint[[kk]][TT,],SCA_P_refitted$factor_group2[[(kk-KK)]][TT,]), ncol=1)
          F_T_GICA <- matrix( c(GICA_refitted$factor_joint[[kk]][TT,],GICA_refitted$factor_group2[[(kk-KK)]][TT,]), ncol=1)
          F_T_DSCA <- matrix( c(DSCA_refitted$factor_joint[[kk]][TT,],DSCA_refitted$factor_group2[[(kk-KK)]][TT,]), ncol=1)
          F_T_DGICA <- matrix( c(DGICA_refitted$factor_joint[[kk]][TT,],DGICA_refitted$factor_group2[[(kk-KK)]][TT,]), ncol=1)
          F_T_Unfitted_SCA_PF2 <- matrix( c(Factor_unfitted_SCA_PF2$factor_joint[[kk]][TT,],
                                            Factor_unfitted_SCA_PF2$factor_group2[[(kk-KK)]][TT,]), ncol=1)
          F_T_Unfitted_GICA <- matrix( c(Factor_unfitted_GICA$factor_joint[[kk]][TT,],
                                         Factor_unfitted_GICA$factor_group2[[(kk-KK)]][TT,]), ncol=1)
        }
        
        Psi_hat_aug_GRIDY <- as.matrix(bdiag(GRIDY_YK$Psi_joint_hat[,,kk],GRIDY_YK$Psi_indiv_hat[,,kk]))
        Psi_hat_aug_SCA_P <- as.matrix(bdiag(SCA_P_YK$Psi_joint_hat[,,kk],SCA_P_YK$Psi_indiv_hat[,,kk]))
        Psi_hat_aug_GICA <- as.matrix(bdiag(GICA_YK$Psi_joint_hat[,,kk],GICA_YK$Psi_indiv_hat[,,kk]))
        Psi_hat_aug_DSCA <- as.matrix(bdiag(DSCA_YK$Psi_joint_hat[,,kk],DSCA_YK$Psi_indiv_hat[,,kk]))
        Psi_hat_aug_DGICA <- as.matrix(bdiag(DGICA_YK$Psi_joint_hat[,,kk],DGICA_YK$Psi_indiv_hat[,,kk]))
        Psi_hat_aug_Unfitted_SCA_PF2 <- as.matrix(bdiag(Unfitted_SCA_PF2_YK$Psi_joint_hat[,,kk],
                                                        Unfitted_SCA_PF2_YK$Psi_indiv_hat[,,kk]))
        Psi_hat_aug_Unfitted_GICA <- as.matrix(bdiag(Unfitted_GICA_YK$Psi_joint_hat[,,kk],
                                                     Unfitted_GICA_YK$Psi_indiv_hat[,,kk]))
        
        for (h in 1:H){
          if (h == 1){
            F_hat_GRIDY[h,,kk] <- t(Psi_hat_aug_GRIDY %*% F_T_GRIDY)
            F_hat_SCA_P[h,,kk] <- t(Psi_hat_aug_SCA_P %*% F_T_SCA_P)
            F_hat_GICA[h,,kk] <- t(Psi_hat_aug_GICA %*% F_T_GICA)
            F_hat_DSCA[h,,kk] <- t(Psi_hat_aug_DSCA %*% F_T_DSCA)
            F_hat_DGICA[h,,kk] <- t(Psi_hat_aug_DGICA %*% F_T_DGICA)
            F_hat_Unfitted_SCA_PF2[h,,kk] <- t(Psi_hat_aug_Unfitted_SCA_PF2 %*% F_T_Unfitted_SCA_PF2)
            F_hat_Unfitted_GICA[h,,kk] <- t(Psi_hat_aug_Unfitted_GICA %*% F_T_Unfitted_GICA)
          }else{
            F_hat_GRIDY[h,,kk] <- t(Psi_hat_aug_GRIDY %*% as.matrix(F_hat_GRIDY[(h-1),,kk],ncol=1))
            F_hat_SCA_P[h,,kk] <- t(Psi_hat_aug_SCA_P %*% as.matrix(F_hat_SCA_P[(h-1),,kk],ncol=1))
            F_hat_GICA[h,,kk] <- t(Psi_hat_aug_GICA %*% as.matrix(F_hat_GICA[(h-1),,kk],ncol=1))
            F_hat_DSCA[h,,kk] <- t(Psi_hat_aug_DSCA %*% as.matrix(F_hat_DSCA[(h-1),,kk],ncol=1))
            F_hat_DGICA[h,,kk] <- t(Psi_hat_aug_DGICA %*% as.matrix(F_hat_DGICA[(h-1),,kk],ncol=1))
            F_hat_Unfitted_SCA_PF2[h,,kk] <- t(Psi_hat_aug_Unfitted_SCA_PF2 %*% as.matrix(F_hat_Unfitted_SCA_PF2[(h-1),,kk],ncol=1))
            F_hat_Unfitted_GICA[h,,kk] <- t(Psi_hat_aug_Unfitted_GICA %*% as.matrix(F_hat_Unfitted_GICA[(h-1),,kk],ncol=1))
          }
          X_hat_GRIDY[h,,kk] <- t(Lambda_hat_GRIDY %*% F_hat_GRIDY[h,,kk])
          X_hat_SCA_P[h,,kk] <- t(Lambda_hat_SCA_P %*% F_hat_SCA_P[h,,kk])
          X_hat_GICA[h,,kk] <- t(Lambda_hat_GICA %*% F_hat_GICA[h,,kk])
          X_hat_DSCA[h,,kk] <- t(Lambda_hat_DSCA %*% F_hat_DSCA[h,,kk])
          X_hat_DGICA[h,,kk] <- t(Lambda_hat_DGICA %*% F_hat_DGICA[h,,kk])
          X_hat_Unfitted_SCA_PF2[h,,kk] <- t(Lambda_hat_Unfitted_SCA_PF2 %*% F_hat_Unfitted_SCA_PF2[h,,kk])
          X_hat_Unfitted_GICA[h,,kk] <- t(Lambda_hat_Unfitted_GICA %*% F_hat_Unfitted_GICA[h,,kk])
        }
      }
      
      MSFE_GRIDY <- vector("numeric",length=H)
      MSFE_SCA_P <- vector("numeric",length=H)
      MSFE_GICA <- vector("numeric",length=H)
      MSFE_DSCA <- vector("numeric",length=H)
      MSFE_DGICA <- vector("numeric",length=H)
      MSFE_Unfitted_SCA_PF2 <- vector("numeric",length=H)
      MSFE_Unfitted_GICA <- vector("numeric",length=H)
      for (h in 1:H){
        for (kk in 1:(2*KK)){
          MSFE_GRIDY[h] <- MSFE_GRIDY[h] + norm(model_dgp_ext$scaled_X[[kk]][,(TT+1):(TT+h)] - t(X_hat_GRIDY[(1:h),,kk]),"2")^2/(2*dd*KK)
          MSFE_SCA_P[h] <- MSFE_SCA_P[h] + norm(model_dgp_ext$scaled_X[[kk]][,(TT+1):(TT+h)] - t(X_hat_SCA_P[(1:h),,kk]),"2")^2/(2*dd*KK)
          MSFE_GICA[h] <- MSFE_GICA[h] + norm(model_dgp_ext$scaled_X[[kk]][,(TT+1):(TT+h)] - t(X_hat_GICA[(1:h),,kk]),"2")^2/(2*dd*KK)
          MSFE_DSCA[h] <- MSFE_DSCA[h] + norm(model_dgp_ext$scaled_X[[kk]][,(TT+1):(TT+h)] - t(X_hat_DSCA[(1:h),,kk]),"2")^2/(2*dd*KK)
          MSFE_DGICA[h] <- MSFE_DGICA[h] + norm(model_dgp_ext$scaled_X[[kk]][,(TT+1):(TT+h)] - t(X_hat_DGICA[(1:h),,kk]),"2")^2/(2*dd*KK)
          MSFE_Unfitted_SCA_PF2[h] <- MSFE_Unfitted_SCA_PF2[h] + norm(model_dgp_ext$scaled_X[[kk]][,(TT+1):(TT+h)] 
                                                                      - t(X_hat_Unfitted_SCA_PF2[(1:h),,kk]),"2")^2/(2*dd*KK)
          MSFE_Unfitted_GICA[h] <- MSFE_Unfitted_GICA[h] + norm(model_dgp_ext$scaled_X[[kk]][,(TT+1):(TT+h)] 
                                                                - t(X_hat_Unfitted_GICA[(1:h),,kk]),"2")^2/(2*dd*KK)
        }
      }
      
      list_GRIDY[["RMSFE"]][[iter]] <- sqrt(MSFE_GRIDY)
      list_SCA_P[["RMSFE"]][[iter]] <- sqrt(MSFE_SCA_P)
      list_GICA[["RMSFE"]][[iter]] <- sqrt(MSFE_GICA)
      list_DSCA[["RMSFE"]][[iter]] <- sqrt(MSFE_DSCA)
      list_DGICA[["RMSFE"]][[iter]] <- sqrt(MSFE_DGICA)
      list_Unfitted_SCA_PF2[["RMSFE"]][[iter]] <- sqrt(MSFE_Unfitted_SCA_PF2)
      list_Unfitted_GICA[["RMSFE"]][[iter]] <- sqrt(MSFE_Unfitted_GICA)
      
      # ----------------------------------------------------------- #
      time_taken <- Sys.time() - start_time
      cat("Currently ",iter,"th iteration is done. Taken time is,",time_taken,"\n")
    }
    
    file_name <- paste0("./result/sim2_snr",case,"_type",set,".RData")
    save( list_GRIDY,list_SCA_P,list_GICA,list_DSCA,list_DGICA,list_Unfitted_SCA_PF2,list_Unfitted_GICA,
          file=file_name)
  }
}