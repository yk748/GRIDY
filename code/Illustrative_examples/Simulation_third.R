#-----------------------------------------------------------------------------#
# Code for producing results for the third simulation setting in Illustrative example section
# Codes for plotting Figure 5 is provided separately
#-----------------------------------------------------------------------------#




# -----------------------------------------------------------------------------#
# Running rank selection
# -----------------------------------------------------------------------------#
# Packages required
library(mvtnorm)
library(Matrix)
library(combinat)
library(multiway) # For SCA and GICA
library(ica) # For SCA and GICA

# Load rotational bootstrap from the parent directory
source(paste0(dirname(getwd()),"/","Rotational_bootstrap.R"))
source(paste0(dirname(getwd()),"/","library_simulation.R"))

sigma_eps <- 1

d_list <- c(100,200,300,100,100)
T_list <- c(100,100,100,200,300)
K_list <- c(50,50,50,50,50)

sigma_c_list <- c(0.25,0.5,0.75,1.25,2,4)

joint_rank <- c(1,1,1,2,1,2,2,3,3,4)
G_indiv_rank <- c(1,2,3,2,4,3,4,3,4,4)

for (case in 1:5){
  
  dd <- d_list[case]
  TT <- T_list[case]
  KK <- K_list[case]
  
  for (set in 1:6){
    
    sigma_c <- sigma_c_list[set]
    cat("Currently, case",case,",set",set,"is set.\n")
    
    r_hat <- array(NA,dim=c((2*KK),length(joint_rank)))
    
    set.seed(1234)
    for (r in 1:10){
      cat("Currently, rank combination",r,"th has been performed.\n")
      start_time <- Sys.time()
      
      r_J <- joint_rank[r]
      r_G <- G_indiv_rank[r]
      
      Control_param <- data.frame(type=1,sigma_xi=0.2,J_per=0.5,
                                  b_bar_min=0,b_bar_max=1,
                                  b_tilde_min=0,b_tilde_max=1,
                                  sigma_c=sigma_c,sigma_eps=sigma_eps)
      
      model_dgp <- sim_model(r_J, r_G, dd, TT, KK, Control_param)
      
      for (kk in 1:(2*KK)){
        possibleError <- tryCatch(
          r_hat[kk,r] <- Rotational_bootstrap(dd,TT,model_dgp$Xk[[kk]],cull=0.5),
          
          error = function(kk){
            print(kk,"th iteration is failed to do Boostrap")
          }
        )
      }
      #-----------------------------------------------------------------#
      time_taken <- Sys.time() - start_time
      cat(r,"th iteration completed.","Time taken is",time_taken,
          "and SNR is",set,"r_J is",r_J,"r_G is",r_G,"\n")
    }
    file_name <- paste0("./result/sim3_set",set,"_d",dd,"T",TT,"K",KK,".RData")
    save(r_hat, file=file_name)
  }
}
