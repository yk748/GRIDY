#-----------------------------------------------------------------------------#
# Code for producing results in Data application section
# Codes for plotting each figure are provided separately
#-----------------------------------------------------------------------------#



# -----------------------------------------------------------------------------#
# Step 0: Set up
# -----------------------------------------------------------------------------#
# Packages required
library(readxl)
library(mvtnorm)
library(Matrix)
library(combinat)
library(multiway) # For SCA and GICA
library(ica) # For SCA and GICA

# Load source code from the parent directory
source(paste0(dirname(getwd()),"/","Rotational_bootstrap.R"))
source(paste0(dirname(getwd()),"/","gica.R")) # For SCA and GICA
source(paste0(dirname(getwd()),"/","AJIVE_retrieve.R"))
source(paste0(dirname(getwd()),"/","library_application.R"))


# Load the list of the file names
file_name_table <- data.frame(read_excel("./data/file_name_list.xlsx", col_names=TRUE))
colnames(file_name_table) <- c("site_idx","file_idx","group_idx")

# Load the list of the variables of the ROIs
voxel_name_table <- data.frame(read_excel("./data/dos160_labels.xlsx", col_names=TRUE))
colnames(voxel_name_table) <- c("ROI_number","ROI_label","Network","Net","Name")


# -----------------------------------------------------------------------------#
# Step 1: Data preprocessing
# -----------------------------------------------------------------------------#
Tot_subject <- dim(file_name_table)[1]
G1_subject <- dim(file_name_table[which(file_name_table$group_idx==1),])[1]
G2_subject <- dim(file_name_table[which(file_name_table$group_idx==2),])[1]

list_scaled_X <- list()
for (kk in 1:Tot_subject){
  
  file_name <- paste0(file_name_table[kk,2],"_ccs_filt_noglobal_rois_dosenbach160")
  file_directory <- paste0("./data/",file_name,".csv")
  data_entry <- read.csv(file_directory)
  
  if (dim(data_entry)[2] == 160){
    list_scaled_X[[kk]] <- t( scale(data_entry, center = TRUE, scale = TRUE) )
    # impute 0 for the observations of the variables whose paths are zero. 
    list_scaled_X[[kk]][is.na(list_scaled_X[[kk]])] <- 0
    list_scaled_X[[kk]] <- list_scaled_X[[kk]]
  }else{
    list_scaled_X[[kk]] <- t( scale(data_entry[,-dim(data_entry)[2]], center = TRUE, scale = TRUE) )
    # impute 0 for the observations of the variables whose paths are zero.
    list_scaled_X[[kk]][is.na(list_scaled_X[[kk]])] <- 0
    list_scaled_X[[kk]] <- list_scaled_X[[kk]]
  }
}

# -----------------------------------------------------------------------------#
# Step 2: Rank selection
# -----------------------------------------------------------------------------#
# rotational bootstrap:
r_hat <- vector("numeric",length(Tot_subject))
d <- 160

set.seed(1234)
for (kk in 1:length(Tot_subject)){
  time_start <- Sys.time()
  T_k <- dim(list_scaled_X[[kk]])[2]
  
  possibleError <- tryCatch(
    r_hat[kk] <- Rotational_bootstrap(d,T_k,t(list_scaled_X[[kk]]),cull=0.5),
    
    error = function(kk){
      print(kk,"th subject is failed to do Boostrap")
    }
  )
  
  cat(kk,"th subject is completed.","Time taken is",Sys.time() - time_start,
      "and estimated rank is",r_hat[kk],"\n")
}

# index data frame:
df_index <- data.frame(subject=c(1:Tot_subject))
df_index$group <- file_name_table$group_idx
df_index$rank <- r_hat
df_index$site <- file_name_table$site_idx

# Save the intermediate result
# This data is used for drawing Figure S1
save(df_index,r_hat,Tot_subject,file_name_table,file="./result/result_intermediate.RData")

# -----------------------------------------------------------------------------#
# Step 3: Screening subjects from the rank selection result
# -----------------------------------------------------------------------------#
# Finalize subject index for the rest of analysis
df_index_final <- df_index[which(df_index$rank > 0 & df_index$rank <= 15),]
X_scale_ext <- list()

# Organize Group 1
extract_idx1 <- df_index_final[which(df_index_final$group == 1),"subject"]
cnt_group1 <- 0
for (kk in extract_idx1){
  cnt_group1 <- cnt_group1+1
  X_scale_ext[[cnt_group1]] <- list_scaled_X[[kk]]
}
K1_idx <- c(1:cnt_group1)


# Organize Group 2
extract_idx2 <- df_index_final[which(df_index_final$group == 2),"subject"]
cnt_group2 <- 0
for (kk in extract_idx2){
  cnt_group2 <- cnt_group2+1
  X_scale_ext[[cnt_group1+cnt_group2]] <- list_scaled_X[[kk]]
}
K2_idx <- cnt_group1+c(1:cnt_group2)


# -----------------------------------------------------------------------------#
# Step 4: running AJIVE
# -----------------------------------------------------------------------------#
r1_total <- 2; r2_total <- 2
initial_ranks <- c(rep(r1_total,length(K1_idx)),rep(r2_total,length(K2_idx)))

set.seed(1234)
AJIVE_decomp <- ajive(X_scale_ext,
                      initial_signal_ranks = initial_ranks,
                      n_wedin_samples = 1000,
                      n_rand_dir_samples = 1000,
                      full = TRUE, joint_rank = NA)
AJIVE_decomp$joint_rank


# Exclude subjects whose group individual rank is zero
Joint_block<- list()
Group1_block <- list()
Group2_block <- list()

exclude_idx <- which(c(mapply(x=K1_idx,function(x){AJIVE_decomp$block_decomps[[x]]$individual$rank}),
                       mapply(x=K2_idx,function(x){AJIVE_decomp$block_decomps[[x]]$individual$rank})) == 0)
effect_tot_idx <- 0 
effect_g1_idx <- 0
effect_g2_idx <- 0

tmp_X_scale_ext <- list()
for (kk in 1:(length(K1_idx)+length(K2_idx))){
  
  if (is.element(kk,exclude_idx)){
    next
  }else{
    effect_tot_idx <- effect_tot_idx +1
    tmp_X_scale_ext[[effect_tot_idx]] <- X_scale_ext[[kk]]
    Joint_block[[effect_tot_idx]] <- t(AJIVE_decomp[["block_decomps"]][[kk]]$joint$full)
    if ( kk <= length(K1_idx) ){
      effect_g1_idx <- effect_g1_idx +1
      Group1_block[[effect_g1_idx]] <- t(AJIVE_decomp[["block_decomps"]][[kk]]$individual$full)
    }else{
      effect_g2_idx <- effect_g2_idx +1
      Group2_block[[effect_g2_idx]] <- t(AJIVE_decomp[["block_decomps"]][[kk]]$individual$full)
    }
  }
}
X_scale_ext <- tmp_X_scale_ext
K1_idx <- c(1:effect_g1_idx)
K2_idx <- effect_g1_idx+c(1:effect_g2_idx)

# -----------------------------------------------------------------------------#
# Step 5: fitting SCA_PF2, SCA_P and GICA
# -----------------------------------------------------------------------------#
rJ <- AJIVE_decomp$joint_rank
rG1 <- 1; rG2 <- 1

SCA_PF2_Joint <- sca(Joint_block, nfac = rJ, type="sca-pf2", verbose = FALSE)
SCA_PF2_Group1 <- sca(Group1_block,nfac = rG1, type="sca-pf2", verbose = FALSE)
SCA_PF2_Group2 <- sca(Group2_block,nfac = rG2, type="sca-pf2", verbose = FALSE)

SCA_P_Joint <- sca(Joint_block, nfac = rJ,type="sca-p", rotation="varimax",verbose = FALSE)
SCA_P_Group1 <- sca(Group1_block, nfac = rG1,type="sca-p",verbose = FALSE)
SCA_P_Group2 <- sca(Group2_block, nfac = rG2,type="sca-p",verbose = FALSE)

GICA_Joint <- gica(Joint_block, nc = rJ, dual.reg = FALSE, center = FALSE)
GICA_Group1 <- gica(Group1_block, nc = rG1, dual.reg = FALSE, center = FALSE)
GICA_Group2 <- gica(Group2_block, nc = rG2, dual.reg = FALSE, center = FALSE)


# -----------------------------------------------------------------------------#
# Step 6: Constructing dynamics for GRIDY,SCA_P, and GICA
# -----------------------------------------------------------------------------#
GRIDY_refitted <- factor_regression(K1_idx,K2_idx,rJ,rG1,rG2,
                                    X_scale_ext,SCA_PF2_Joint,SCA_PF2_Group1,SCA_PF2_Group2)
GRIDY_YK <- YK_compute(K1_idx,K2_idx,rJ,rG1,rG2,X_scale_ext,GRIDY_refitted)

SCA_P_refitted <- factor_regression(K1_idx,K2_idx,rJ,rG1,rG2,
                                    X_scale_ext,SCA_P_Joint,SCA_P_Group1,SCA_P_Group2)
SCA_P_YK <- YK_compute(K1_idx,K2_idx,rJ,rG1,rG2,X_scale_ext,SCA_P_refitted)

GICA_refitted <- factor_regression(K1_idx,K2_idx,rJ,rG1,rG2,
                                   X_scale_ext,GICA_Joint,GICA_Group1,GICA_Group2)
GICA_YK <- YK_compute(K1_idx,K2_idx,rJ,rG1,rG2,X_scale_ext,GICA_refitted)


# Save the final result
# This data is used for drawing Figures S2-S3, 6, and 7
save(voxel_name_table,X_scale_ext,K1_idx,K2_idx,
     SCA_PF2_Joint,SCA_PF2_Group1,SCA_PF2_Group2,
     SCA_P_Joint,SCA_P_Group1,SCA_P_Group2,
     GICA_Joint,GICA_Group1,GICA_Group2,
     GRIDY_refitted,GRIDY_YK,
     SCA_P_refitted,SCA_P_YK,
     GICA_refitted,GICA_YK,
     file="./result/result_final.RData")


