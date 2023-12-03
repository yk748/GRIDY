rm(list = ls())

# Packages required
library("readxl");
library("mvtnorm");
library("Matrix");
library("combinat");
library("latex2exp");
library("ggplot2");
library("patchwork");
library("reshape2");
library("gridExtra");
library("ggpubr");
library("RColorBrewer");
library("forcats");
library("ggdendro");
library("tidyr");
library("matrixcalc");
library("lmtest");
library("cowplot");
library("grid");

# For AJIVE:
# devtools::install_github("idc9/r_jive")
library(ajive)

# For SCA:
library(multiway)

# For miscellaneous functions
source("library_application.R")

# -----------------------------------------------------------------------------#
# Step 0 : Data preprocessing
# -----------------------------------------------------------------------------#

file_name_table <- data.frame(read_excel("file_name_list.xlsx", col_names=TRUE))
colnames(file_name_table) <- c("site_idx","file_idx","group_idx")

voxel_name_table <- data.frame(read_excel("dos160_labels.xlsx", col_names=TRUE))
colnames(voxel_name_table) <- c("ROI_number","ROI_label","Network","Net","Name")

Tot_subject <- dim(file_name_table)[1]
G1_subject <- dim(file_name_table[which(file_name_table$group_idx==1),])[1]
G2_subject <- dim(file_name_table[which(file_name_table$group_idx==2),])[1]
Tot_subject - (G1_subject+G2_subject)

list_X <- list()
list_scaled_X <- list()
directory_folder <- "data"
for (k in 1:Tot_subject){
  
  file_name <- paste0(file_name_table[k,2],"_ccs_filt_noglobal_rois_dosenbach160")
  file_directory <- paste0(directory_folder,"/",file_name,".csv")
  data_entry <- read.csv(file_directory)
  
  if (dim(data_entry)[2] == 160){
    list_X[[k]] <- data_entry
    list_scaled_X[[k]] <- t( scale(list_X[[k]], center = TRUE, scale = TRUE) )
    list_scaled_X[[k]][is.na(list_scaled_X[[k]])] <- 0
  }else{
    list_X[[k]] <- data_entry[,-dim(data_entry)[2]] 
    list_scaled_X[[k]] <- t( scale(list_X[[k]], center = TRUE, scale = TRUE) )
    list_scaled_X[[k]][is.na(list_scaled_X[[k]])] <- 0
  }
}

avg_T <- 0
for (k in 1:Tot_subject){
  avg_T <- avg_T + dim(list_X[[k]])[1]/Tot_subject
}
avg_T 


# -----------------------------------------------------------------------------#
# Step 1 : Initial rank selection
# -----------------------------------------------------------------------------#
r_hat <- vector("numeric",Tot_subject)
d <- 160
for (k in 1:Tot_subject){
  
  time_start <- Sys.time()
  T_k <- dim(list_scaled_X[[k]])[2]
  
  possibleError <- tryCatch(
    r_hat[k] <- Rotational_bootstrap(d,T_k,t(list_scaled_X[[k]]),cull=0.5),
    
    error = function(k){
      print(k,"th subject is failed to do Boostrap")
    }
  )
  
  cat(k,"th subject is completed.","Time taken is",Sys.time() - time_start,
      "and estimated rank is",r_hat[k],"\n")
}

# rank display:
df_index <- data.frame(subject=c(1:Tot_subject))
df_index$group <- file_name_table$group_idx
df_index$rank <- r_hat
df_index$site <- file_name_table$site_idx


# rank table:
rank_table <- data.frame(rank = unique(df_index$rank)[order(unique(df_index$rank))],
                         frequency = tabulate(match(r_hat,unique(df_index$rank)[order(unique(df_index$rank))])))
rank_table

rank_table_G1 <- data.frame(rank = unique(df_index$rank[which(df_index$group == 1)])[order(unique(df_index$rank[which(df_index$group == 1)]))],
                            frequency = tabulate(match(r_hat[which(df_index$group == 1)],
                                                       unique(df_index$rank[which(df_index$group == 1)])[order(unique(df_index$rank[which(df_index$group == 1)]))])))
rank_table_G1 

rank_table_G2 <- data.frame(rank = unique(df_index$rank[which(df_index$group == 2)])[order(unique(df_index$rank[which(df_index$group == 2)]))],
                            frequency = tabulate(match(r_hat[which(df_index$group == 2)],
                                                       unique(df_index$rank[which(df_index$group == 2)])[order(unique(df_index$rank[which(df_index$group == 2)]))])))
rank_table_G2 

# Histogram for the paper:
ABIDE_color <- c("#1C9E77","#D95F02","#7570B3","#E72A8A","#19398A",
                 "#2B5592","#3E729A","#5480A0","#706B9F","#8D569E",
                 "#AA3C94","#C62482","#E22671","#E12456","#DB3C38",
                 "#D85B1F","#E36922","#EE7723","#F68D2C","#FAAC3F")
df_hist <- df_index
df_hist_G1 <- df_index[which(df_index$group==1),]
df_hist_G1$site <- factor(df_hist_G1$site,labels=sort(unique(df_hist_G1$site)))
df_hist_G2 <- df_index[which(df_index$group==2),]
df_hist_G2$site <- factor(df_hist_G2$site,labels=sort(unique(df_hist_G2$site)))

G_hist <- ggplot(df_hist) +
  geom_bar(aes(x=rank,fill=site),alpha=1) +
  scale_fill_manual("Sites",
                    values=ABIDE_color) + 
  scale_x_continuous(name=TeX("Estimated initial ranks $\\hat{r}$"),
                     breaks=unique(df_hist$rank)) +
  scale_y_continuous("Frequency",limits=c(0,250),
                     sec.axis = sec_axis(~./5*2, name = "Cumulative frequency", 
                                         labels = function(b){paste0(round(b,0),"%")})) +
  ggtitle("Total") +
  theme(legend.position="bottom",axis.text.x = element_text(size = 6)) + 
  guides(fill = guide_legend(nrow = 2)) +
  geom_step(data=rank_table_G1,aes(x=rank, y=cumsum(frequency/dim(df_hist_G1)[1]*250)),alpha=0.5,linetype="dashed")

G1_hist <- ggplot(df_hist_G1) +
  geom_bar(aes(x=rank,fill=site),alpha=1) +
  scale_fill_manual("Sites",
                    values=ABIDE_color) + 
  scale_x_continuous(name=TeX("Estimated initial ranks $\\hat{r}$"),
                     breaks=unique(df_hist_G1$rank)) +
  scale_y_continuous("Frequency",limits=c(0,100),
                     sec.axis = sec_axis(~., name = "Cumulative frequency", 
                                         labels = function(b){paste0(round(b,0),"%")})) +
  ggtitle("Group 1") +
  theme(legend.position="bottom",axis.text.x = element_text(size = 6)) + 
  guides(fill = guide_legend(nrow = 2)) +
  geom_step(data=rank_table_G1,aes(x=rank, y=cumsum(frequency/dim(df_hist_G1)[1]*100)),alpha=0.5,linetype="dashed")

G2_hist <- ggplot(df_hist_G2) +
  geom_bar(aes(x=rank,fill=site),alpha=1) +
  scale_fill_manual("Sites",
                    values=ABIDE_color) + 
  scale_x_continuous(name=TeX("Estimated initial ranks $\\hat{r}$"),
                     breaks=unique(df_hist_G2$rank)) +
  scale_y_continuous("Frequency",limits=c(0,150),
                     sec.axis = sec_axis(~./3*2, name = "Cumulative frequency", 
                                         labels = function(b){paste0(round(b,0),"%")})) +
  ggtitle("Group 2") +
  theme(legend.position="bottom",axis.text.x = element_text(size = 6)) + 
  guides(fill = guide_legend(nrow = 2)) +
  geom_step(data=rank_table_G2,aes(x=rank, y=cumsum(frequency/dim(df_hist_G2)[1]*150)),alpha=0.5,linetype="dashed")

par(mar=c(0,0,0,0))
ggarrange(G_hist,G1_hist,G2_hist,nrow=3,common.legend=TRUE,legend="bottom")


# Finalize index:
df_index_final <- df_index[which(df_index$rank > 0 & df_index$rank <= 15),]

X_final_ext <- list()
X_scale_final_ext <- list()
scaler <- list()

extract_idx1 <- df_index_final[which(df_index_final$group == 1),"subject"]
extract_idx2 <- df_index_final[which(df_index_final$group == 2),"subject"]

cnt_group1 <- 0
for (k in extract_idx1){
  
  cnt_group1 <- cnt_group1+1
  X_final_ext[[cnt_group1]] <- list_X[[k]]
  X_scale_final_ext[[cnt_group1]] <- list_scaled_X[[k]]
  scaler[[cnt_group1]] <- attr(list_scaled_X[[k]],"scaled:scale")
}

cnt_group2 <- 0
for (k in extract_idx2){
  
  cnt_group2 <- cnt_group2+1
  X_final_ext[[cnt_group1+cnt_group2]] <- list_X[[k]]
  X_scale_final_ext[[cnt_group1+cnt_group2]] <- list_scaled_X[[k]]
  scaler[[cnt_group1+cnt_group2]] <- attr(list_scaled_X[[k]],"scaled:scale")
}

K1_idx <- c(1:cnt_group1)
K2_idx <- cnt_group1+c(1:cnt_group2)

# Hold out last 10 observations
H <- 10

X_final <- list()
X_scale_final <- list()

for (k in 1:(cnt_group1+cnt_group2)){
  
  X_final[[k]] <- X_final_ext[[k]][-c((dim(X_final_ext[[k]])[1]-H+1):dim(X_final_ext[[k]])[1]),]
  X_scale_final[[k]] <- X_scale_final_ext[[k]][,-c((dim(X_scale_final_ext[[k]])[2]-H+1):dim(X_scale_final_ext[[k]])[2])]
}e


# -----------------------------------------------------------------------------#
# Step 2 : AJIVE
# -----------------------------------------------------------------------------#
r1_total <- 2
r2_total <- 2
initial_ranks <- c(rep(r1_total,length(K1_idx)),rep(r2_total,length(K2_idx)))

GRIDY_AJIVE <- ajive(X_scale_final,
                     initial_signal_ranks = initial_ranks,
                     n_wedin_samples = 1000,
                     n_rand_dir_samples = 1000,
                     full = TRUE, joint_rank = NA)

GRIDY_AJIVE$joint_rank
par(mfrow=c(1,2))
hist(mapply(x=K1_idx,function(x){GRIDY_AJIVE$block_decomps[[x]]$individual$rank}),
     main="Group individual ranks for group 1",xlab=TeX(r'($\hat{r}_{G}$)'),xlim=c(0,4),ylim=c(0,120),col="dodgerblue1",border="dodgerblue1")
hist(mapply(x=K2_idx,function(x){GRIDY_AJIVE$block_decomps[[x]]$individual$rank}),
     main="Group individual ranks for group 2",xlab=TeX(r'($\hat{r}_{G}$)'),xlim=c(0,4),ylim=c(0,120),col="firebrick1",border="firebrick1")


df_tmp <- data.frame(estim=c(mapply(x=K1_idx,function(x){GRIDY_AJIVE$block_decomps[[x]]$individual$rank}),
                   mapply(x=K2_idx,function(x){GRIDY_AJIVE$block_decomps[[x]]$individual$rank})),
           group=c(rep(1,length(K1_idx)),rep(2,length(K2_idx))))




# Exclude subjects whose individual rank is zero
Joint_block<- list()
Group1_block <- list()
Group2_block <- list()

exclude_idx <- which(c(mapply(x=K1_idx,function(x){GRIDY_AJIVE$block_decomps[[x]]$individual$rank}),
                       mapply(x=K2_idx,function(x){GRIDY_AJIVE$block_decomps[[x]]$individual$rank})) == 0)
cnt_effect_tot_idx <- 0 
cnt_effect_g1_idx <- 0
cnt_effect_g2_idx <- 0

tmp_X_final_ext <- list()
tmp_X_scale_final_ext <- list()
tmp_scaler <- list()
tmp_X_final <- list()
tmp_X_scale_final <- list()

for (k in 1:(length(K1_idx)+length(K2_idx))){
  
  if (is.element(k,exclude_idx)){
    next
    
  }else{
    
    cnt_effect_tot_idx <- cnt_effect_tot_idx +1
    
    tmp_X_final_ext[[cnt_effect_tot_idx]] <- X_final_ext[[k]]
    tmp_X_scale_final_ext[[cnt_effect_tot_idx]] <- X_scale_final_ext[[k]]
    tmp_scaler[[cnt_effect_tot_idx]] <- scaler[[k]]
    tmp_X_final[[cnt_effect_tot_idx]] <- X_final[[k]]
    tmp_X_scale_final[[cnt_effect_tot_idx]] <- X_scale_final[[k]]
    
    Joint_block[[cnt_effect_tot_idx]] <- t(GRIDY_AJIVE[["block_decomps"]][[k]]$joint$full)
    
    if ( k <= length(K1_idx) ){
      
      cnt_effect_g1_idx <- cnt_effect_g1_idx +1
      
      Group1_block[[cnt_effect_g1_idx]] <- t(GRIDY_AJIVE[["block_decomps"]][[k]]$individual$full)
      
    }else{
      
      cnt_effect_g2_idx <- cnt_effect_g2_idx +1
      
      Group2_block[[cnt_effect_g2_idx]] <- t(GRIDY_AJIVE[["block_decomps"]][[k]]$individual$full)
    }
  }
}

X_final_ext <- tmp_X_final_ext
X_scale_final_ext <- tmp_X_scale_final_ext
scaler <- tmp_scaler
X_final <- tmp_X_final
X_scale_final <- tmp_X_scale_final
K1_idx <- c(1:cnt_effect_g1_idx)
K2_idx <- cnt_effect_g1_idx+c(1:cnt_effect_g2_idx)

rm(tmp_X_final_ext,tmp_X_scale_final_ext,tmp_scaler,tmp_X_final,tmp_X_scale_final,
   cnt_effect_tot_idx,cnt_effect_g1_idx,cnt_effect_g2_idx)



#-----------------------------------------------------------------#
# Step 3: SCA
#-----------------------------------------------------------------#
rJ <- GRIDY_AJIVE$joint_rank
rG1 <- 1
rG2 <- 1
# rG1 <- r1_total - rJ
# rG2 <- r2_total - rJ

possibleError <- tryCatch(
  GRIDY_Joint <- sca(Joint_block, nfac = rJ, type="sca-pf2", verbose = FALSE),
  
  error = function(x){
    print("Failed to estimate Joint")
  }
  
)

possibleError <- tryCatch(
  GRIDY_Group1 <- sca(Group1_block,nfac = rG1, type="sca-pf2", verbose = FALSE),
  
  error = function(x){
    print("Failed to estimate Group Individual")
  }
)


possibleError <- tryCatch(
  GRIDY_Group2 <- sca(Group2_block,nfac = rG2, type="sca-pf2", verbose = FALSE),
  
  error = function(x){
    print("Failed to estimate Group Individual")
  }
)

#-----------------------------------------------------------------#
# Step 4. Refitting and Yule-Walker eqs
#-----------------------------------------------------------------#
GRIDY_refitted <- factor_regression(K1_idx,K2_idx,
                                    rJ,rG1,rG2,
                                    X_scale_final,
                                    GRIDY_Joint,
                                    GRIDY_Group1,
                                    GRIDY_Group2)

GRIDY_YK <- YK_compute(K1_idx,K2_idx,
                       rJ,rG1,rG2,
                       X_scale_final,
                       GRIDY_refitted)


#-----------------------------------------------------------------#
# Comparative methods
#-----------------------------------------------------------------#
# Double SCA:
Whole_block <- list()
for (k in 1:(length(K1_idx)+length(K2_idx))){
  Whole_block[[k]] <- t(X_scale_final[[k]])
}

possibleError <- tryCatch(
  DSCA_Joint <- sca(Whole_block, nfac = rJ,type="sca-pf2", verbose = FALSE),
  error = function(x){
    print(x,"failed to estimate first stage of DSCA")
  }
)


Remaining_block_group1 <- list()
Remaining_block_group2 <- list()
for (k in 1:(length(K1_idx)+length(K2_idx))){
  kth_Joint_block <- DSCA_Joint$D[[k]] %*% t(DSCA_Joint$B)
  if (k <= length(K1_idx)){
    Remaining_block_group1[[k]] <- Whole_block[[k]] - kth_Joint_block
  }else{
    Remaining_block_group2[[(k-length(K1_idx))]] <- Whole_block[[k]] - kth_Joint_block
  }
}


possibleError <- tryCatch(
  DSCA_Group1 <- sca(Remaining_block_group1, nfac = rG1,type="sca-pf2", verbose = FALSE),
  error = function(x){
    print(x,"failed to estimate Group 1 at second stage of DSCA")
  }
)

possibleError <- tryCatch(
  DSCA_Group2 <- sca(Remaining_block_group2, nfac = rG2,type="sca-pf2", verbose = FALSE),
  error = function(x){
    print(x,"failed to estimate Group 2 at second stage of DSCA")
  }
)


DSCA_refitted <- factor_regression(K1_idx,K2_idx,
                                   rJ,rG1,rG2,
                                   X_scale_final,
                                   DSCA_Joint,
                                   DSCA_Group1,
                                   DSCA_Group2)

DSCA_YK <- YK_compute(K1_idx,K2_idx,
                      rJ,rG1,rG2,
                      X_scale_final,
                      DSCA_refitted)


# Individual PCA:
PCA_Joint <- list()
PCA_Group1 <- list()
PCA_Group2 <- list()

possibleError <- tryCatch(
  for (k in 1:(length(K1_idx)+length(K2_idx))){
    
    T_k <- dim(X_scale_final[[k]])[2]
    
    svd_joint <- svd(Joint_block[[k]])
    B_k_joint <- svd_joint$v[,1:rJ]
    D_k_joint <- svd_joint$u[,1:rJ] %*% diag(svd_joint$d[1:rJ],rJ)
    Cov_F_joint <- diag(diag(t(D_k_joint) %*% D_k_joint/T_k),rJ)
    
    if (k == 1){
      PCA_Joint$D[[1]] <- D_k_joint %*% solve(sqrt(Cov_F_joint))
      PCA_Joint$B <- B_k_joint %*% sqrt(Cov_F_joint)
    }else{
      D_k_joint <- D_k_joint %*% solve(sqrt(Cov_F_joint))
      B_k_joint <- B_k_joint %*% sqrt(Cov_F_joint)
      Q_rotate_joint <- solve(t(B_k_joint) %*% B_k_joint) %*% t(B_k_joint) %*% PCA_Joint$B
      PCA_Joint$D[[k]] <- D_k_joint %*% Q_rotate_joint
    }
    
    if (k <= length(K1_idx)){
      svd_group1 <- svd(Group1_block[[k]])
      B_k_group1 <- svd_group1$v[,1:rG1]
      D_k_group1 <- svd_group1$u[,1:rG1] %*% diag(svd_group1$d[1:rG1],rG1)
      Cov_F_group1 <- diag(diag(t(D_k_group1) %*% D_k_group1/T_k),rG1)
      
      if (k == 1){
        PCA_Group1$D[[1]] <- D_k_group1 %*% solve(sqrt(Cov_F_group1))
        PCA_Group1$B <- B_k_group1 %*% sqrt(Cov_F_group1)
      }else{
        D_k_group1 <- D_k_group1 %*% solve(sqrt(Cov_F_group1))
        B_k_group1 <- B_k_group1 %*% sqrt(Cov_F_group1)
        Q_rotate_group1 <- solve(t(B_k_group1) %*% B_k_group1) %*% t(B_k_group1) %*% PCA_Group1$B
        PCA_Group1$D[[k]] <- D_k_group1 %*% Q_rotate_group1
        
      }
    }else{
      svd_group2 <- svd(Group2_block[[(k-length(K1_idx))]])
      B_k_group2 <- svd_group2$v[,1:rG2]
      D_k_group2 <- svd_group2$u[,1:rG2] %*% diag(svd_group2$d[1:rG2],rG2)
      Cov_F_group2 <- diag(diag(t(D_k_group2) %*% D_k_group2/T_k),rG2)
      
      if (k == (length(K1_idx)+1)){
        PCA_Group2$D[[1]] <- D_k_group2 %*% solve(sqrt(Cov_F_group2))
        PCA_Group2$B <- B_k_group2 %*% sqrt(Cov_F_group2)
      }else{
        D_k_group2 <- D_k_group2 %*% solve(sqrt(Cov_F_group2))
        B_k_group2 <- B_k_group2 %*% sqrt(Cov_F_group2)
        Q_rotate_group2 <- solve(t(B_k_group2) %*% B_k_group2) %*% t(B_k_group2) %*% PCA_Group2$B
        PCA_Group2$D[[(k-length(K1_idx))]] <- D_k_group2 %*% Q_rotate_group2
        
      }
    }
  },
  error = function(x){
    print(x,"failed to estimate PCA")
  }
)

PCA_refitted <- factor_regression(K1_idx,K2_idx,
                                  rJ,rG1,rG2,
                                  X_scale_final,
                                  PCA_Joint,
                                  PCA_Group1,
                                  PCA_Group2)

PCA_YK <- YK_compute(K1_idx,K2_idx,
                     rJ,rG1,rG2,
                     X_scale_final,
                     PCA_refitted)

#-----------------------------------------------------------------#
# Step 5. Forecasting
#-----------------------------------------------------------------#
H <- 10

X_hat_GRIDY <- array(NA,dim=c(H,160,(length(K1_idx)+length(K2_idx))))
F_group1_hat_GRIDY <- array(NA,dim=c(H,(rJ+rG1),length(K1_idx)))
F_group2_hat_GRIDY <- array(NA,dim=c(H,(rJ+rG2),length(K2_idx)))

X_hat_DSCA <- array(NA,dim=c(H,160,(length(K1_idx)+length(K2_idx))))
F_group1_hat_DSCA <- array(NA,dim=c(H,(rJ+rG1),length(K1_idx)))
F_group2_hat_DSCA <- array(NA,dim=c(H,(rJ+rG2),length(K2_idx)))

X_hat_PCA <- array(NA,dim=c(H,160,(length(K1_idx)+length(K2_idx))))
F_group1_hat_PCA <- array(NA,dim=c(H,(rJ+rG1),length(K1_idx)))
F_group2_hat_PCA <- array(NA,dim=c(H,(rJ+rG2),length(K2_idx)))

for (k in 1:(length(K1_idx)+length(K2_idx))){
  
  T_k <- dim(X_scale_final[[k]])[2]
  
  if (k <= length(K1_idx)){
    
    if (rG1 != 0){
      Lambda_hat_GRIDY <- cbind(GRIDY_Joint$B,GRIDY_Group1$B)
      F_T_GRIDY <- matrix( c(GRIDY_refitted$factor_joint[[k]][T_k,],
                             GRIDY_refitted$factor_group1[[k]][T_k,]), ncol=1)
      Psi_hat_aug_GRIDY <- as.matrix(bdiag(GRIDY_YK$Psi_joint_hat[,,k],
                                           GRIDY_YK$Psi_indiv1_hat[,,k]))
      
      Lambda_hat_DSCA <- cbind(DSCA_Joint$B,DSCA_Group1$B)
      F_T_DSCA <- matrix( c(DSCA_refitted$factor_joint[[k]][T_k,],
                            DSCA_refitted$factor_group1[[k]][T_k,]), ncol=1)
      Psi_hat_aug_DSCA <- as.matrix(bdiag(DSCA_YK$Psi_joint_hat[,,k],
                                          DSCA_YK$Psi_indiv1_hat[,,k]))
      
      Lambda_hat_PCA <- cbind(PCA_Joint$B,PCA_Group1$B)
      F_T_PCA <- matrix( c(PCA_refitted$factor_joint[[k]][T_k,],
                           PCA_refitted$factor_group1[[k]][T_k,]), ncol=1)
      Psi_hat_aug_PCA <- as.matrix(bdiag(PCA_YK$Psi_joint_hat[,,k],
                                         PCA_YK$Psi_indiv1_hat[,,k]))
    }else{
      Lambda_hat_GRIDY <- cbind(GRIDY_Joint$B)
      F_T_GRIDY <- matrix( c(GRIDY_refitted$factor_joint[[k]][T_k,]), ncol=1)
      Psi_hat_aug_GRIDY <- as.matrix(bdiag(GRIDY_YK$Psi_joint_hat[,,k]))
      
      Lambda_hat_PCA <- cbind(PCA_Joint$B)
      F_T_PCA <- matrix( c(PCA_refitted$factor_joint[[k]][T_k,]), ncol=1)
      Psi_hat_aug_PCA <- as.matrix(bdiag(PCA_YK$Psi_joint_hat[,,k]))
      
      Lambda_hat_DSCA <- cbind(DSCA_Joint$B)
      F_T_DSCA <- matrix( c(DSCA_refitted$factor_joint[[k]][T_k,]), ncol=1)
      Psi_hat_aug_DSCA <- as.matrix(bdiag(DSCA_YK$Psi_joint_hat[,,k]))
    }
  }else{
    
    if (rG2 != 0){
      Lambda_hat_GRIDY <- cbind(GRIDY_Joint$B,GRIDY_Group2$B)
      F_T_GRIDY <- matrix( c(GRIDY_refitted$factor_joint[[k]][T_k,],
                             GRIDY_refitted$factor_group2[[(k-length(K1_idx))]][T_k,]), ncol=1)
      Psi_hat_aug_GRIDY <- as.matrix(bdiag(GRIDY_YK$Psi_joint_hat[,,k],
                                           GRIDY_YK$Psi_indiv2_hat[,,(k-length(K1_idx))]))
      
      Lambda_hat_PCA <- cbind(PCA_Joint$B,PCA_Group2$B)
      F_T_PCA <- matrix( c(PCA_refitted$factor_joint[[k]][T_k,],
                           PCA_refitted$factor_group2[[(k-length(K1_idx))]][T_k,]), ncol=1)
      Psi_hat_aug_PCA <- as.matrix(bdiag(PCA_YK$Psi_joint_hat[,,k],
                                         PCA_YK$Psi_indiv2_hat[,,(k-length(K1_idx))]))
      
      Lambda_hat_DSCA <- cbind(DSCA_Joint$B,DSCA_Group2$B)
      F_T_DSCA <- matrix( c(DSCA_refitted$factor_joint[[k]][T_k,],
                            DSCA_refitted$factor_group2[[(k-length(K1_idx))]][T_k,]), ncol=1)
      Psi_hat_aug_DSCA <- as.matrix(bdiag(DSCA_YK$Psi_joint_hat[,,k],
                                          DSCA_YK$Psi_indiv2_hat[,,(k-length(K1_idx))]))
    }else{
      Lambda_hat_GRIDY <- cbind(GRIDY_Joint$B)
      F_T_GRIDY <- matrix( c(GRIDY_refitted$factor_joint[[k]][T_k,]), ncol=1)
      Psi_hat_aug_GRIDY <- as.matrix(bdiag(GRIDY_YK$Psi_joint_hat[,,k]))
      
      Lambda_hat_PCA <- cbind(PCA_Joint$B)
      F_T_PCA <- matrix( c(PCA_refitted$factor_joint[[k]][T_k,]), ncol=1)
      Psi_hat_aug_PCA <- as.matrix(bdiag(PCA_YK$Psi_joint_hat[,,k]))
      
      Lambda_hat_DSCA <- cbind(DSCA_Joint$B)
      F_T_DSCA <- matrix( c(DSCA_refitted$factor_joint[[k]][T_k,]), ncol=1)
      Psi_hat_aug_DSCA <- as.matrix(bdiag(DSCA_YK$Psi_joint_hat[,,k]))
    }
  }
  
  
  for (h in 1:H){
    if (h == 1){
      if (k <= length(K1_idx)){
        F_group1_hat_GRIDY[h,,k] <- t(Psi_hat_aug_GRIDY %*% F_T_GRIDY)  
        
        F_group1_hat_PCA[h,,k] <- t(Psi_hat_aug_PCA %*% F_T_PCA) 
        
        F_group1_hat_DSCA[h,,k] <- t(Psi_hat_aug_DSCA %*% F_T_DSCA)
      }else{
        F_group2_hat_GRIDY[h,,(k-length(K1_idx))] <- t(Psi_hat_aug_GRIDY %*% F_T_GRIDY)
        
        F_group2_hat_PCA[h,,(k-length(K1_idx))] <- t(Psi_hat_aug_PCA %*% F_T_PCA)
        
        F_group2_hat_DSCA[h,,(k-length(K1_idx))] <- t(Psi_hat_aug_DSCA %*% F_T_DSCA)
      }
    }else{
      if (k <= length(K1_idx)){
        F_group1_hat_GRIDY[h,,k] <- t(Psi_hat_aug_GRIDY %*% as.matrix(F_group1_hat_GRIDY[(h-1),,k],ncol=1)) 
        
        F_group1_hat_PCA[h,,k] <- t(Psi_hat_aug_PCA %*% as.matrix(F_group1_hat_PCA[(h-1),,k],ncol=1)) 
        
        F_group1_hat_DSCA[h,,k] <- t(Psi_hat_aug_DSCA %*% as.matrix(F_group1_hat_DSCA[(h-1),,k],ncol=1)) 
      }else{
        F_group2_hat_GRIDY[h,,(k-length(K1_idx))] <- t(Psi_hat_aug_GRIDY %*% as.matrix(F_group2_hat_GRIDY[(h-1),,(k-length(K1_idx))],ncol=1))
        
        F_group2_hat_PCA[h,,(k-length(K1_idx))] <- t(Psi_hat_aug_PCA %*% as.matrix(F_group2_hat_PCA[(h-1),,(k-length(K1_idx))],ncol=1))
        
        F_group2_hat_DSCA[h,,(k-length(K1_idx))] <- t(Psi_hat_aug_DSCA %*% as.matrix(F_group2_hat_DSCA[(h-1),,(k-length(K1_idx))],ncol=1))
        
      }
    }
    if (k <= length(K1_idx)){
      X_hat_GRIDY[h,,k] <- t(Lambda_hat_GRIDY %*% F_group1_hat_GRIDY[h,,k])
      
      X_hat_PCA[h,,k] <- t(Lambda_hat_PCA %*% F_group1_hat_PCA[h,,k])
      
      X_hat_DSCA[h,,k] <- t(Lambda_hat_DSCA %*% F_group1_hat_DSCA[h,,k])
    }else{
      X_hat_GRIDY[h,,k] <- t(Lambda_hat_GRIDY %*% F_group2_hat_GRIDY[h,,(k-length(K1_idx))])
      
      X_hat_PCA[h,,k] <- t(Lambda_hat_PCA %*% F_group2_hat_PCA[h,,(k-length(K1_idx))])
      
      X_hat_DSCA[h,,k] <- t(Lambda_hat_DSCA %*% F_group2_hat_DSCA[h,,(k-length(K1_idx))])
    }
  }
}



#-----------------------------------------------------------------#
# Step 6. Evaluation & Illustration
#-----------------------------------------------------------------#
color_table <- data.frame(color=c(rep("red",length(which(voxel_name_table$Net=="dfm"))),
                                  rep("blue",length(which(voxel_name_table$Net=="fp"))),
                                  rep("green",length(which(voxel_name_table$Net=="co"))),
                                  rep("purple",length(which(voxel_name_table$Net=="sm"))),
                                  rep("orange",length(which(voxel_name_table$Net=="cb"))),
                                  rep("black",length(which(voxel_name_table$Net=="oc")))),
                          ROI=c(which(voxel_name_table$Net=="dfm"),
                                which(voxel_name_table$Net=="fp"),
                                which(voxel_name_table$Net=="co"),
                                which(voxel_name_table$Net=="sm"),
                                which(voxel_name_table$Net=="cb"),
                                which(voxel_name_table$Net=="oc")))
color_table <- color_table[c(which(color_table$color=="red"),
                             which(color_table$color=="blue"),
                             which(color_table$color=="green"),
                             which(color_table$color=="purple"),
                             which(color_table$color=="orange"),
                             which(color_table$color=="black")),]

# loadings
df_B_joint_1 <- data.frame(value = GRIDY_Joint$B[color_table$ROI,1],
                           ROI = color_table$ROI,
                           color = color_table$color,
                           factor = rep(1,160),
                           structure=rep("Joint",160))
df_B_joint_2 <- data.frame(value = GRIDY_Joint$B[color_table$ROI,2],
                           ROI = color_table$ROI,
                           color = color_table$color,
                           factor = rep(2,160),
                           structure=rep("Joint",160))
df_B_joint <- rbind(df_B_joint_1,df_B_joint_2)
df_B_joint$ROI <- factor(df_B_joint$ROI,levels=color_table$ROI)
df_B_joint$factor <- factor(df_B_joint$factor,levels=c(1,2))


df_B_group1_1 <- data.frame(value = GRIDY_Group1$B[color_table$ROI,1],
                            ROI = color_table$ROI,
                            color = color_table$color,
                            factor = rep(1,160),
                            structure=rep("Group Individual 1",160))
df_B_group1 <- rbind(df_B_group1_1)
df_B_group1$ROI <- factor(df_B_group1$ROI,levels=color_table$ROI)
df_B_group1$factor <- factor(df_B_group1$factor,levels=1)


df_B_group2_1 <- data.frame(value = GRIDY_Group2$B[color_table$ROI,1],
                            ROI = color_table$ROI,
                            color = color_table$color,
                            factor = rep(1,160),
                            structure=rep("Group Individual 2",160))
df_B_group2 <- rbind(df_B_group2_1)
df_B_group2$ROI <- factor(df_B_group2$ROI,levels=color_table$ROI)
df_B_group2$factor <- factor(df_B_group2$factor,levels=1)

df_B <- rbind(df_B_joint,df_B_group1,df_B_group2)
df_B$structure <- factor(df_B$structure,levels=c("Joint","Group Individual 1","Group Individual 2"))
factor.labs <- c("Factor 1", "Factor 2", "Factor 3")
names(factor.labs) <- c(1,2)
structure.labs <- c("Joint", "Group Individual 1", "Group Individual 2")
names(structure.labs) <- c("Joint", "Group Individual 1", "Group Individual 2")

par(mar=c(0,0,0,0))
ggplot(df_B,aes(x=as.factor(ROI),y=value,fill=color)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="Network",
                    labels=c("dfm","fp","co","sm","cb","oc"),
                    breaks=c("red","blue","green","purple","orange","black"),
                    values=c("red","blue","green","purple","orange","black")) +
  scale_y_continuous(name="Value") + 
  theme(legend.position="bottom",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + guides(fill = guide_legend(nrow = 1)) +
  facet_wrap(~structure+factor,nrow=4,labeller=labeller(structure=structure.labs,factor=factor.labs))


# Transition matrices and covariance matrices
df_Psi_group1_joint <- data.frame(value = as.vector(GRIDY_YK$Psi_joint_hat[,,K1_idx]),
                                  index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K1_idx)),
                                  group = rep("A",length(K1_idx)),
                                  structure = rep("trans",length(K1_idx)))
df_Psi_group2_joint <- data.frame(value = as.vector(GRIDY_YK$Psi_joint_hat[,,K2_idx]),
                                  index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K2_idx)),
                                  group = rep("B",length(K2_idx)),
                                  structure = rep("trans",length(K2_idx)))



df_Psi_group1_individual <- data.frame(value = as.vector(GRIDY_YK$Psi_indiv1_hat[,,K1_idx]),
                                       index = rep(c("(1,1)"),length(K1_idx)),
                                       group = rep("C",length(K1_idx)),
                                       structure = rep("trans",length(K1_idx)))
df_Psi_group2_individual <- data.frame(value = as.vector(GRIDY_YK$Psi_indiv2_hat[,,1:length(K2_idx)]),
                                       index = rep(c("(1,1)"),length(K2_idx)),
                                       group = rep("D",length(K2_idx)),
                                       structure = rep("trans",length(K2_idx)))


df_sig_group1_joint <- data.frame(value = as.vector(GRIDY_YK$Eta_joint_hat[,,K1_idx]),
                                  index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K1_idx)),
                                  group = rep("A",length(K1_idx)),
                                  structure = rep("sigma",length(K1_idx)))
df_sig_group2_joint <- data.frame(value = as.vector(GRIDY_YK$Eta_joint_hat[,,K2_idx]),
                                  index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K2_idx)),
                                  group = rep("B",length(K2_idx)),
                                  structure = rep("sigma",length(K2_idx)))

df_sig_group1_individual <- data.frame(value = as.vector(GRIDY_YK$Eta_indiv1_hat[,,K1_idx]),
                                       index = rep(c("(1,1)"),length(K1_idx)),
                                       group = rep("C",length(K1_idx)),
                                       structure = rep("sigma",length(K1_idx)))
df_sig_group2_individual <- data.frame(value = as.vector(GRIDY_YK$Eta_indiv2_hat[,,1:length(K2_idx)]),
                                       index = rep(c("(1,1)"),length(K2_idx)),
                                       group = rep("D",length(K2_idx)),
                                       structure = rep("sigma",length(K2_idx)))

df_fs <- rbind(df_Psi_group1_joint,df_Psi_group2_joint,df_Psi_group1_individual,df_Psi_group2_individual,
               df_sig_group1_joint,df_sig_group2_joint,df_sig_group1_individual,df_sig_group2_individual)
df_fs$structure <- factor(df_fs$structure,levels=c("trans","sigma"))
mycolors_fs <- colorRampPalette(brewer.pal(8,"Spectral"))(8)
str.labs <- c("VAR transition matrices", "Covariance matrices")
names(str.labs) <- c("trans","sigma")
gr.labs <- c("Group 1 Joint", "Group 2 Joint",
             "Group 1 Group Individual", "Group 2 Group Individual")
names(gr.labs) <- c("A","B","C","D")

par(mar=c(0,0,0,0))
ggplot(df_fs,aes(y=value,x=index,fill=factor(group):factor(structure))) +
  geom_boxplot(cex=0.1) +
  scale_fill_manual(name="",values=mycolors_fs) +
  scale_y_continuous(name="Value") +
  scale_x_discrete(name="Entries") +
  theme(legend.position="none") +
  facet_grid(structure~group,
             scales = "free",
             labeller=labeller(structure=str.labs,group=gr.labs))



# Factor series
for(k in 1:(length(K1_idx)+length(K2_idx))) {
  if (k <= length(K1_idx)){
    if ( k %% 20 == 1){
      print(paste0(k,"th individual's observation in group 1 is ",dim(X_scale_final[[k]])[2]))
    }
  }else{
    if ( (k-length(K1_idx)) %% 20 == 1){
      print(paste0(k,"th individual's observation in group 2 is ",dim(X_scale_final[[k]])[2]))
    }
  }
}

# joint factors:
df_F1_group1_joint <- data.frame(value=c(GRIDY_refitted$factor_joint[[1]][,1],
                                         GRIDY_refitted$factor_joint[[21]][,1],
                                         GRIDY_refitted$factor_joint[[41]][,1],
                                         GRIDY_refitted$factor_joint[[61]][,1],
                                         GRIDY_refitted$factor_joint[[81]][,1],
                                         GRIDY_refitted$factor_joint[[101]][,1],
                                         GRIDY_refitted$factor_joint[[121]][,1]),
                                 time=c(c(1:length(GRIDY_refitted$factor_joint[[1]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[21]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[41]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[61]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[81]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[101]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[121]][,1]))),
                                 indiv=c(rep(1,length(GRIDY_refitted$factor_joint[[1]][,1])),
                                         rep(2,length(GRIDY_refitted$factor_joint[[21]][,1])),
                                         rep(3,length(GRIDY_refitted$factor_joint[[41]][,1])),
                                         rep(4,length(GRIDY_refitted$factor_joint[[61]][,1])),
                                         rep(5,length(GRIDY_refitted$factor_joint[[81]][,1])),
                                         rep(6,length(GRIDY_refitted$factor_joint[[101]][,1])),
                                         rep(7,length(GRIDY_refitted$factor_joint[[121]][,1]))))
df_F1_group1_joint$factor <- rep(1,dim(df_F1_group1_joint)[1])
df_F1_group1_joint$group <- rep(1,dim(df_F1_group1_joint)[1])

df_F2_group1_joint <- data.frame(value=c(GRIDY_refitted$factor_joint[[1]][,2],
                                         GRIDY_refitted$factor_joint[[21]][,2],
                                         GRIDY_refitted$factor_joint[[41]][,2],
                                         GRIDY_refitted$factor_joint[[61]][,2],
                                         GRIDY_refitted$factor_joint[[81]][,2],
                                         GRIDY_refitted$factor_joint[[101]][,2],
                                         GRIDY_refitted$factor_joint[[121]][,2]),
                                 time=c(c(1:length(GRIDY_refitted$factor_joint[[1]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[21]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[41]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[61]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[81]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[101]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[121]][,1]))),
                                 indiv=c(rep(1,length(GRIDY_refitted$factor_joint[[1]][,1])),
                                         rep(2,length(GRIDY_refitted$factor_joint[[21]][,1])),
                                         rep(3,length(GRIDY_refitted$factor_joint[[41]][,1])),
                                         rep(4,length(GRIDY_refitted$factor_joint[[61]][,1])),
                                         rep(5,length(GRIDY_refitted$factor_joint[[81]][,1])),
                                         rep(6,length(GRIDY_refitted$factor_joint[[101]][,1])),
                                         rep(7,length(GRIDY_refitted$factor_joint[[121]][,1]))))
df_F2_group1_joint$factor <- rep(2,dim(df_F2_group1_joint)[1])
df_F2_group1_joint$group <- rep(1,dim(df_F2_group1_joint)[1])

df_F_group1_joint <- rbind(df_F1_group1_joint,df_F2_group1_joint)
df_F_group1_joint$structure <- rep("Joint",dim(df_F_group1_joint)[1])
df_F_group1_joint$factor <- factor(df_F_group1_joint$factor,levels=c(1,2))
df_F_group1_joint$indiv <- factor(df_F_group1_joint$indiv,levels=c(1,2,3,4,5,6,7))

df_F1_group2_joint <- data.frame(value=c(GRIDY_refitted$factor_joint[[154]][,1],
                                         GRIDY_refitted$factor_joint[[174]][,1],
                                         GRIDY_refitted$factor_joint[[194]][,1],
                                         GRIDY_refitted$factor_joint[[214]][,1],
                                         GRIDY_refitted$factor_joint[[234]][,1],
                                         GRIDY_refitted$factor_joint[[254]][,1],
                                         GRIDY_refitted$factor_joint[[274]][,1]),
                                 time=c(c(1:length(GRIDY_refitted$factor_joint[[154]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[174]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[194]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[214]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[234]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[254]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[274]][,1]))),
                                 indiv=c(rep(1,length(GRIDY_refitted$factor_joint[[154]][,1])),
                                         rep(2,length(GRIDY_refitted$factor_joint[[174]][,1])),
                                         rep(3,length(GRIDY_refitted$factor_joint[[194]][,1])),
                                         rep(4,length(GRIDY_refitted$factor_joint[[214]][,1])),
                                         rep(5,length(GRIDY_refitted$factor_joint[[234]][,1])),
                                         rep(6,length(GRIDY_refitted$factor_joint[[254]][,1])),
                                         rep(7,length(GRIDY_refitted$factor_joint[[274]][,1]))))
df_F1_group2_joint$factor <- rep(1,dim(df_F1_group2_joint)[1])
df_F1_group2_joint$group <- rep(2,dim(df_F1_group2_joint)[1])

df_F2_group2_joint <- data.frame(value=c(GRIDY_refitted$factor_joint[[154]][,2],
                                         GRIDY_refitted$factor_joint[[174]][,2],
                                         GRIDY_refitted$factor_joint[[194]][,2],
                                         GRIDY_refitted$factor_joint[[214]][,2],
                                         GRIDY_refitted$factor_joint[[234]][,2],
                                         GRIDY_refitted$factor_joint[[254]][,2],
                                         GRIDY_refitted$factor_joint[[274]][,2]),
                                 time=c(c(1:length(GRIDY_refitted$factor_joint[[154]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[174]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[194]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[214]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[234]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[254]][,1])),
                                        c(1:length(GRIDY_refitted$factor_joint[[274]][,1]))),
                                 indiv=c(rep(1,length(GRIDY_refitted$factor_joint[[154]][,1])),
                                         rep(2,length(GRIDY_refitted$factor_joint[[174]][,1])),
                                         rep(3,length(GRIDY_refitted$factor_joint[[194]][,1])),
                                         rep(4,length(GRIDY_refitted$factor_joint[[214]][,1])),
                                         rep(5,length(GRIDY_refitted$factor_joint[[234]][,1])),
                                         rep(6,length(GRIDY_refitted$factor_joint[[254]][,1])),
                                         rep(7,length(GRIDY_refitted$factor_joint[[274]][,1]))))
df_F2_group2_joint$factor <- rep(2,dim(df_F2_group2_joint)[1])
df_F2_group2_joint$group <- rep(2,dim(df_F2_group2_joint)[1])

df_F_group2_joint <- rbind(df_F1_group2_joint,df_F2_group2_joint)
df_F_group2_joint$structure <- rep("Joint",dim(df_F_group2_joint)[1])
df_F_group2_joint$factor <- factor(df_F_group2_joint$factor,levels=c(1,2))
df_F_group2_joint$indiv <- factor(df_F_group2_joint$indiv,levels=c(1,2,3,4,5,6,7))

df_F_joint <- rbind(df_F_group1_joint,df_F_group2_joint)


# group individual factors:
df_F1_group1_individual <- data.frame(value=c(GRIDY_refitted$factor_group1[[1]][,1],
                                              GRIDY_refitted$factor_group1[[21]][,1],
                                              GRIDY_refitted$factor_group1[[41]][,1],
                                              GRIDY_refitted$factor_group1[[61]][,1],
                                              GRIDY_refitted$factor_group1[[81]][,1],
                                              GRIDY_refitted$factor_group1[[101]][,1],
                                              GRIDY_refitted$factor_group1[[121]][,1]),
                                      time=c(c(1:length(GRIDY_refitted$factor_group1[[1]][,1])),
                                             c(1:length(GRIDY_refitted$factor_group1[[21]][,1])),
                                             c(1:length(GRIDY_refitted$factor_group1[[41]][,1])),
                                             c(1:length(GRIDY_refitted$factor_group1[[61]][,1])),
                                             c(1:length(GRIDY_refitted$factor_group1[[81]][,1])),
                                             c(1:length(GRIDY_refitted$factor_group1[[101]][,1])),
                                             c(1:length(GRIDY_refitted$factor_group1[[121]][,1]))),
                                      indiv=c(rep(1,length(GRIDY_refitted$factor_group1[[1]][,1])),
                                              rep(2,length(GRIDY_refitted$factor_group1[[21]][,1])),
                                              rep(3,length(GRIDY_refitted$factor_group1[[41]][,1])),
                                              rep(4,length(GRIDY_refitted$factor_group1[[61]][,1])),
                                              rep(5,length(GRIDY_refitted$factor_group1[[81]][,1])),
                                              rep(6,length(GRIDY_refitted$factor_group1[[101]][,1])),
                                              rep(7,length(GRIDY_refitted$factor_group1[[121]][,1]))))
df_F1_group1_individual$factor <- rep(1,dim(df_F1_group1_individual)[1])
df_F1_group1_individual$group <- rep(1,dim(df_F1_group1_individual)[1])

df_F_group1_individual <- rbind(df_F1_group1_individual)
df_F_group1_individual$structure <- rep("Group Individual 1",dim(df_F1_group1_individual)[1])
df_F_group1_individual$factor <- factor(df_F1_group1_individual$factor,levels=1)
df_F_group1_individual$indiv <- factor(df_F1_group1_individual$indiv,levels=c(1,2,3,4,5,6,7))


df_F1_group2_individual <- data.frame(value=c(GRIDY_refitted$factor_group2[[1]][,1],
                                              GRIDY_refitted$factor_group2[[21]][,1],
                                              GRIDY_refitted$factor_group2[[41]][,1],
                                              GRIDY_refitted$factor_group2[[61]][,1],
                                              GRIDY_refitted$factor_group2[[81]][,1],
                                              GRIDY_refitted$factor_group2[[101]][,1],
                                              GRIDY_refitted$factor_group2[[121]][,1]),
                                      time=c(c(1:length(GRIDY_refitted$factor_group2[[1]][,1])),
                                             c(1:length(GRIDY_refitted$factor_group2[[21]][,1])),
                                             c(1:length(GRIDY_refitted$factor_group2[[41]][,1])),
                                             c(1:length(GRIDY_refitted$factor_group2[[61]][,1])),
                                             c(1:length(GRIDY_refitted$factor_group2[[81]][,1])),
                                             c(1:length(GRIDY_refitted$factor_group2[[101]][,1])),
                                             c(1:length(GRIDY_refitted$factor_group2[[121]][,1]))),
                                      indiv=c(rep(1,length(GRIDY_refitted$factor_group2[[1]][,1])),
                                              rep(2,length(GRIDY_refitted$factor_group2[[21]][,1])),
                                              rep(3,length(GRIDY_refitted$factor_group2[[41]][,1])),
                                              rep(4,length(GRIDY_refitted$factor_group2[[61]][,1])),
                                              rep(5,length(GRIDY_refitted$factor_group2[[81]][,1])),
                                              rep(6,length(GRIDY_refitted$factor_group2[[101]][,1])),
                                              rep(7,length(GRIDY_refitted$factor_group2[[121]][,1]))))
df_F1_group2_individual$factor <- rep(1,dim(df_F1_group2_individual)[1])
df_F1_group2_individual$group <- rep(2,dim(df_F1_group2_individual)[1])

df_F_group2_individual <- rbind(df_F1_group2_individual)
df_F_group2_individual$structure <- rep("Group Individual 2",dim(df_F_group2_individual)[1])
df_F_group2_individual$factor <- factor(df_F_group2_individual$factor,levels=1)
df_F_group2_individual$indiv <- factor(df_F_group2_individual$indiv,levels=c(1,2,3,4,5,6,7))

df_F <- rbind(df_F_joint,df_F_group1_individual,df_F_group2_individual)
df_F$structure <- factor(df_F$structure,levels=c("Joint","Group Individual 1","Group Individual 2"))

factor.labs <- c("Factor 1", "Factor 2")
names(factor.labs) <- c(1,2)
group.labs <- c("Group 1", "Group 2")
names(group.labs) <- c(1,2)
structure.labs <- c("Joint", "Group Individual 1", "Group Individual 2")
names(structure.labs) <- c("Joint","Group Individual 1", "Group Individual 2")

mycolors <- colorRampPalette(brewer.pal(7,"Spectral"))(7)
par(mar=c(0,0,0,0))
Plot_F <- ggplot(df_F,aes(y=value,x=as.factor(time),group=as.factor(indiv),color=as.factor(indiv))) +
  geom_line(cex=0.5) +
  scale_color_manual(name="",values=mycolors) +
  scale_x_discrete(name="Time",
                   breaks=c(1:320)[which(c(1:320)%%20 == 1)],
                   labels=c(1:320)[which(c(1:320)%%20 == 1)]) +
  scale_y_continuous(name="Value") +
  theme(legend.position="none") +
  facet_grid(structure+factor~group,
             labeller=labeller(structure=structure.labs,factor=factor.labs,group=group.labs))

gt <- plot_to_gtable(Plot_F)
gt$grobs[c(which(gt$layout$t == 14 & gt$layout$r == 5 & grepl('panel', gt$layout$name)),
           which(gt$layout$t == 12 & gt$layout$r == 7 & grepl('panel', gt$layout$name)))] <- NULL
gt$layout <- gt$layout[-c(which(gt$layout$t == 14 & gt$layout$r == 5 & grepl('panel', gt$layout$name)),
                          which(gt$layout$t == 12 & gt$layout$r == 7 & grepl('panel', gt$layout$name))),]
grid.newpage()
grid.draw(gt)



# R2 statistics
R2_group1 <- array(NA,dim=c(160,3,length(K1_idx)))
R2_group2 <- array(NA,dim=c(160,5,length(K2_idx)))
for (k in 1:(length(K1_idx)+length(K2_idx))){
  
  X_block <- t(X_scale_final[[k]])
  
  if (k <= length(K1_idx)){
    Components_factor1 <- as.matrix(GRIDY_refitted$factor_joint[[k]][,1]) %*% t(as.matrix(GRIDY_Joint$B[,1]))
    Components_factor2 <- as.matrix(GRIDY_refitted$factor_joint[[k]][,2]) %*% t(as.matrix(GRIDY_Joint$B[,2]))
    Components_factor3 <- as.matrix(GRIDY_refitted$factor_group1[[k]][,1]) %*% t(as.matrix(GRIDY_Group1$B[,1]))
    
    for (i in 1:160){
      R2_group1[i,1,k] <- sum((Components_factor1[,i])^2) / sum((X_block[,i])^2)
      R2_group1[i,2,k] <- sum((Components_factor2[,i])^2) / sum((X_block[,i])^2)
      R2_group1[i,3,k] <- sum((Components_factor3[,i])^2) / sum((X_block[,i])^2)
    }
    
  }else{
    
    Components_factor1 <- as.matrix(GRIDY_refitted$factor_joint[[k]][,1]) %*% t(as.matrix(GRIDY_Joint$B[,1]))
    Components_factor2 <- as.matrix(GRIDY_refitted$factor_joint[[k]][,2]) %*% t(as.matrix(GRIDY_Joint$B[,2]))
    Components_factor3 <- as.matrix(GRIDY_refitted$factor_group2[[(k-length(K1_idx))]][,1]) %*% t(as.matrix(GRIDY_Group2$B[,1]))
    for (i in 1:160){
      R2_group2[i,1,(k-length(K1_idx))] <- sum((Components_factor1[,i])^2) / sum((X_block[,i])^2)
      R2_group2[i,2,(k-length(K1_idx))] <- sum((Components_factor2[,i])^2) / sum((X_block[,i])^2)
      R2_group2[i,3,(k-length(K1_idx))] <- sum((Components_factor3[,i])^2) / sum((X_block[,i])^2)
      
    }
  }
}

reorder_idx <- color_table
R2_group1_factor1 <- data.frame(t(R2_group1[reorder_idx$ROI,1,]))
R2_group1_factor2 <- data.frame(t(R2_group1[reorder_idx$ROI,2,]))
R2_group1_factor3 <- data.frame(t(R2_group1[reorder_idx$ROI,3,]))

colnames(R2_group1_factor1) <- reorder_idx$ROI
colnames(R2_group1_factor2) <- reorder_idx$ROI
colnames(R2_group1_factor3) <- reorder_idx$ROI

R2_group2_factor1 <- data.frame(t(R2_group2[reorder_idx$ROI,1,]))
R2_group2_factor2 <- data.frame(t(R2_group2[reorder_idx$ROI,2,]))
R2_group2_factor3 <- data.frame(t(R2_group2[reorder_idx$ROI,3,]))

colnames(R2_group2_factor1) <- reorder_idx$ROI
colnames(R2_group2_factor2) <- reorder_idx$ROI
colnames(R2_group2_factor3) <- reorder_idx$ROI

df_R2_group1_factor1 <- melt(R2_group1_factor1)
df_R2_group1_factor2 <- melt(R2_group1_factor2)
df_R2_group1_factor3 <- melt(R2_group1_factor3)

df_R2_group1_factor1$color <- rep(reorder_idx$color,each=length(K1_idx))
df_R2_group1_factor2$color <- rep(reorder_idx$color,each=length(K1_idx))
df_R2_group1_factor3$color <- rep(reorder_idx$color,each=length(K1_idx))

df_R2_group2_factor1 <- melt(R2_group2_factor1)
df_R2_group2_factor2 <- melt(R2_group2_factor2)
df_R2_group2_factor3 <- melt(R2_group2_factor3)

df_R2_group2_factor1$color <- rep(reorder_idx$color,each=length(K2_idx))
df_R2_group2_factor2$color <- rep(reorder_idx$color,each=length(K2_idx))
df_R2_group2_factor3$color <- rep(reorder_idx$color,each=length(K2_idx))

df_R2_group1_factor1 <- df_R2_group1_factor1[-which(is.infinite(df_R2_group1_factor1$value)),]
df_R2_group1_factor2 <- df_R2_group1_factor2[-which(is.infinite(df_R2_group1_factor2$value)),]
df_R2_group1_factor3 <- df_R2_group1_factor3[-which(is.infinite(df_R2_group1_factor3$value)),]

df_R2_group2_factor1 <- df_R2_group2_factor1[-which(is.infinite(df_R2_group2_factor1$value)),]
df_R2_group2_factor2 <- df_R2_group2_factor2[-which(is.infinite(df_R2_group2_factor2$value)),]
df_R2_group2_factor3 <- df_R2_group2_factor3[-which(is.infinite(df_R2_group2_factor3$value)),]

df_R2_group1_factor1$factor <- rep(1,dim(df_R2_group1_factor1)[1])
df_R2_group1_factor2$factor <- rep(2,dim(df_R2_group1_factor2)[1])
df_R2_group1_factor3$factor <- rep(1,dim(df_R2_group1_factor3)[1])
df_R2_group1_factor1$structure <- rep("Joint",dim(df_R2_group1_factor1)[1])
df_R2_group1_factor2$structure <- rep("Joint",dim(df_R2_group1_factor2)[1])
df_R2_group1_factor3$structure <- rep("Group Individual 1",dim(df_R2_group1_factor3)[1])

df_R2_group2_factor1$factor <- rep(1,dim(df_R2_group2_factor1)[1])
df_R2_group2_factor2$factor <- rep(2,dim(df_R2_group2_factor2)[1])
df_R2_group2_factor3$factor <- rep(1,dim(df_R2_group2_factor3)[1])
df_R2_group2_factor1$structure <- rep("Joint",dim(df_R2_group2_factor1)[1])
df_R2_group2_factor2$structure <- rep("Joint",dim(df_R2_group2_factor2)[1])
df_R2_group2_factor3$structure <- rep("Group Individual 2",dim(df_R2_group2_factor3)[1])

df_R2_group1 <- rbind(df_R2_group1_factor1,df_R2_group1_factor2,df_R2_group1_factor3)
df_R2_group2 <- rbind(df_R2_group2_factor1,df_R2_group2_factor2,df_R2_group2_factor3)

df_R2_group1$group <- rep(1,dim(df_R2_group1)[1])
df_R2_group2$group <- rep(2,dim(df_R2_group2)[1])

df_R2 <- rbind(df_R2_group1,df_R2_group2)

df_R2$factor <- factor(df_R2$factor,levels=c(1,2))
df_R2$group <- factor(df_R2$group,levels=c(1,2))
df_R2$structure <- factor(df_R2$structure,levels=c("Joint","Group Individual 1","Group Individual 2"))

factor.labs <- c("Factor 1", "Factor 2", "Factor 2")
names(factor.labs) <- c(1,2,3)
group.labs <- c("Group 1", "Group 2")
names(group.labs) <- c(1,2)
structure.labs <- c("Joint","Group Individual 1","Group Individual 2")
names(structure.labs) <- c("Joint","Group Individual 1","Group Individual 2")

par(mar=c(0,0,0,0))
Plot_R2 <- ggplot(df_R2,aes(x=variable,y=value,color=color)) +
  geom_boxplot(size=0.1,alpha=0.5) +
  scale_color_manual(name="Network",
                     labels=c("dfm","fp","co","sm","cb","oc"),
                     breaks=c("red","blue","green","purple","orange","black"),
                     values=c("red","blue","green","purple","orange","black")) +
  scale_y_continuous(name=TeX("$\\R^2$"),limits=c(0,1)) +
  theme(legend.position="bottom",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  guides(colour = guide_legend(nrow = 1)) + 
  facet_grid(structure+factor~group,
             labeller=labeller(structure=structure.labs,factor=factor.labs,group=group.labs))

gt_R <- plot_to_gtable(Plot_R2)
gt_R$grobs[c(which(gt_R$layout$t == 14 & gt_R$layout$r == 5 & grepl('panel', gt_R$layout$name)),
             which(gt_R$layout$t == 12 & gt_R$layout$r == 7 & grepl('panel', gt_R$layout$name)))] <- NULL
gt_R$layout <- gt_R$layout[-c(which(gt_R$layout$t == 14 & gt_R$layout$r == 5 & grepl('panel', gt_R$layout$name)),
                              which(gt_R$layout$t == 12 & gt_R$layout$r == 7 & grepl('panel', gt_R$layout$name))),]
grid.newpage()
grid.draw(gt_R)


# VAR
list_tr_joint <- list()
list_tr_group1 <- list()
list_tr_group2 <- list()
list_error_joint <- list()
list_error_group1 <- list()
list_error_group2 <- list()

for (k in 1:(length(K1_idx)+length(K2_idx))){
  
  T_k <- dim(X_scale_final[[k]])[2]
  
  list_tr_joint[[k]] <- GRIDY_Joint$B %*% GRIDY_YK$Psi_joint_hat[,,k] %*% solve(t(GRIDY_Joint$B) %*% GRIDY_Joint$B) %*% t(GRIDY_Joint$B)
  
  if (k <= length(K1_idx)){
    if( rG1 !=0 ){
      list_tr_group1[[k]] <- GRIDY_Group1$B %*% GRIDY_YK$Psi_indiv1_hat[,,k] %*% solve(t(GRIDY_Group1$B) %*% GRIDY_Group1$B) %*% t(GRIDY_Group1$B)
    }else{
      next
      
    } 
  }else{
    if (rG2 != 0){
      list_tr_group2[[(k-length(K1_idx))]] <- GRIDY_Group2$B %*% GRIDY_YK$Psi_indiv2_hat[,,(k-length(K1_idx))] %*% solve(t(GRIDY_Group2$B) %*% GRIDY_Group2$B) %*% t(GRIDY_Group2$B)
    }else{
      next
      
    }
  }
  
  Cov_eps_refitted <- (t(GRIDY_refitted$noise_refitted[[k]]) %*% GRIDY_refitted$noise_refitted[[k]]) / T_k
  list_error_joint[[k]] <- list_tr_joint[[k]] %*% (diag(1,160) + Cov_eps_refitted) + GRIDY_Joint$B %*% GRIDY_YK$Eta_joint_hat[,,k] %*% t(GRIDY_Joint$B)
  
  if (k <= length(K1_idx)){
    if (rG1 != 0){
      list_error_group1[[k]] <- list_tr_group1[[k]] %*% (diag(1,160) + Cov_eps_refitted) + GRIDY_Group1$B %*% GRIDY_YK$Eta_indiv1_hat[,,k] %*% t(GRIDY_Group1$B)
      
    }else{
      next
      
    }  
  }else{
    if (rG2 != 0){
      list_error_group2[[(k-length(K1_idx))]] <- list_tr_group2[[(k-length(K1_idx))]] %*% (diag(1,160) + Cov_eps_refitted) + GRIDY_Group2$B %*% GRIDY_YK$Eta_indiv2_hat[,,(k-length(K1_idx))] %*% t(GRIDY_Group2$B)
      
    }else{
      next
      
    }  
  }
}

Avg_tr_group1 <- array(0,dim=c(160,160))
Avg_tr_group2 <- array(0,dim=c(160,160))
Avg_tr_joint_group1 <- array(0,dim=c(160,160))
Avg_tr_joint_group2 <- array(0,dim=c(160,160))
Avg_tr_indiv_group1 <- array(0,dim=c(160,160))
Avg_tr_indiv_group2 <- array(0,dim=c(160,160))

Avg_error_group1 <- array(0,dim=c(160,160))
Avg_error_group2 <- array(0,dim=c(160,160))
Avg_error_joint_group1 <- array(0,dim=c(160,160))
Avg_error_joint_group2 <- array(0,dim=c(160,160))
Avg_error_indiv_group1 <- array(0,dim=c(160,160))
Avg_error_indiv_group2 <- array(0,dim=c(160,160))

for (k in 1:(length(K1_idx)+length(K2_idx))){
  
  if (k <= length(K1_idx)){
    Avg_tr_joint_group1  <- Avg_tr_joint_group1  + list_tr_joint[[k]]/length(K1_idx)
    Avg_tr_indiv_group1  <- Avg_tr_indiv_group1  + list_tr_group1[[k]]/length(K1_idx)
    Avg_error_joint_group1 <- Avg_error_joint_group1 + list_error_joint[[k]]/length(K1_idx)
    Avg_error_indiv_group1 <- Avg_error_indiv_group1 + list_error_group1[[k]]/length(K1_idx)
    
    Avg_tr_group1 <- Avg_tr_joint_group1 + Avg_tr_indiv_group1 
    Avg_error_group1 <- Avg_error_joint_group1 + Avg_error_indiv_group1 
    
  }else{
    Avg_tr_joint_group2  <- Avg_tr_joint_group2  + list_tr_joint[[k]]/length(K2_idx)
    Avg_tr_indiv_group2  <- Avg_tr_indiv_group2  + list_tr_group2[[(k-length(K1_idx))]]/length(K2_idx)
    Avg_error_joint_group2 <- Avg_error_joint_group2 + list_error_joint[[k]]/length(K2_idx)
    Avg_error_indiv_group2 <- Avg_error_indiv_group2 + list_error_group2[[(k-length(K1_idx))]]/length(K2_idx)
    
    Avg_tr_group2 <- Avg_tr_joint_group2 + Avg_tr_indiv_group2 
    Avg_error_group2 <- Avg_error_joint_group2 + Avg_error_indiv_group2 
    
  }
}


df_tr_joint_group1 <- data.frame(value=as.vector(Avg_tr_joint_group1[color_table$ROI,color_table$ROI]),
                                 row = rep(color_table$ROI,each=1,160),
                                 column = rep(color_table$ROI,each=160))
df_tr_joint_group2 <- data.frame(value=as.vector(Avg_tr_joint_group2[color_table$ROI,color_table$ROI]),
                                 row = rep(color_table$ROI,each=1,160),
                                 column = rep(color_table$ROI,each=160))
df_tr_indiv_group1 <- data.frame(value=as.vector(Avg_tr_indiv_group1[color_table$ROI,color_table$ROI]),
                                 row = rep(color_table$ROI,each=1,160),
                                 column = rep(color_table$ROI,each=160))
df_tr_indiv_group2 <- data.frame(value=as.vector(Avg_tr_indiv_group2[color_table$ROI,color_table$ROI]),
                                 row = rep(color_table$ROI,each=1,160),
                                 column = rep(color_table$ROI,each=160))
df_tr_group1 <- data.frame(value=as.vector(Avg_tr_group1[color_table$ROI,color_table$ROI]),
                           row = rep(color_table$ROI,each=1,160),
                           column = rep(color_table$ROI,each=160))
df_tr_group2 <- data.frame(value=as.vector(Avg_tr_group2[color_table$ROI,color_table$ROI]),
                           row = rep(color_table$ROI,each=1,160),
                           column = rep(color_table$ROI,each=160))

df_error_joint_group1 <- data.frame(value=as.vector(Avg_error_joint_group1[color_table$ROI,color_table$ROI]),
                                    row = rep(color_table$ROI,each=1,160),
                                    column = rep(color_table$ROI,each=160))
df_error_joint_group2 <- data.frame(value=as.vector(Avg_error_joint_group2[color_table$ROI,color_table$ROI]),
                                    row = rep(color_table$ROI,each=1,160),
                                    column = rep(color_table$ROI,each=160))
df_error_indiv_group1 <- data.frame(value=as.vector(Avg_error_indiv_group1[color_table$ROI,color_table$ROI]),
                                    row = rep(color_table$ROI,each=1,160),
                                    column = rep(color_table$ROI,each=160))
df_error_indiv_group2 <- data.frame(value=as.vector(Avg_error_indiv_group2[color_table$ROI,color_table$ROI]),
                                    row = rep(color_table$ROI,each=1,160),
                                    column = rep(color_table$ROI,each=160))
df_error_group1 <- data.frame(value=as.vector(Avg_error_group1[color_table$ROI,color_table$ROI]),
                              row = rep(color_table$ROI,each=1,160),
                              column = rep(color_table$ROI,each=160))
df_error_group2 <- data.frame(value=as.vector(Avg_error_group2[color_table$ROI,color_table$ROI]),
                              row = rep(color_table$ROI,each=1,160),
                              column = rep(color_table$ROI,each=160))

df_tr_joint_group1$row <- factor(df_tr_joint_group1$row,levels=unique(factor(df_tr_joint_group1$row)))
df_tr_joint_group1$column <- factor(df_tr_joint_group1$column,levels=unique(factor(df_tr_joint_group1$column)))
df_tr_joint_group2$row <- factor(df_tr_joint_group2$row,levels=unique(factor(df_tr_joint_group2$row)))
df_tr_joint_group2$column <- factor(df_tr_joint_group2$column,levels=unique(factor(df_tr_joint_group2$column)))

df_tr_indiv_group1$row <- factor(df_tr_indiv_group1$row,levels=unique(factor(df_tr_indiv_group1$row)))
df_tr_indiv_group1$column <- factor(df_tr_indiv_group1$column,levels=unique(factor(df_tr_indiv_group1$column)))
df_tr_indiv_group2$row <- factor(df_tr_indiv_group2$row,levels=unique(factor(df_tr_indiv_group2$row)))
df_tr_indiv_group2$column <- factor(df_tr_indiv_group2$column,levels=unique(factor(df_tr_indiv_group2$column)))

df_tr_group1$row <- factor(df_tr_group1$row,levels=unique(factor(df_tr_group1$row)))
df_tr_group1$column <- factor(df_tr_group1$column,levels=unique(factor(df_tr_group1$column)))
df_tr_group2$row <- factor(df_tr_group2$row,levels=unique(factor(df_tr_group2$row)))
df_tr_group2$column <- factor(df_tr_group2$column,levels=unique(factor(df_tr_group2$column)))

df_error_joint_group1$row <- factor(df_error_joint_group1$row,levels=unique(factor(df_error_joint_group1$row)))
df_error_joint_group1$column <- factor(df_error_joint_group1$column,levels=unique(factor(df_error_joint_group1$column)))
df_error_joint_group2$row <- factor(df_error_joint_group2$row,levels=unique(factor(df_error_joint_group2$row)))
df_error_joint_group2$column <- factor(df_error_joint_group2$column,levels=unique(factor(df_error_joint_group2$column)))

df_error_indiv_group1$row <- factor(df_error_indiv_group1$row,levels=unique(factor(df_error_indiv_group1$row)))
df_error_indiv_group1$column <- factor(df_error_indiv_group1$column,levels=unique(factor(df_error_indiv_group1$column)))
df_error_indiv_group2$row <- factor(df_error_indiv_group2$row,levels=unique(factor(df_error_indiv_group2$row)))
df_error_indiv_group2$column <- factor(df_error_indiv_group2$column,levels=unique(factor(df_error_indiv_group2$column)))

df_error_group1$row <- factor(df_error_group1$row,levels=unique(factor(df_error_group1$row)))
df_error_group1$column <- factor(df_error_group1$column,levels=unique(factor(df_error_group1$column)))
df_error_group2$row <- factor(df_error_group2$row,levels=unique(factor(df_error_group2$row)))
df_error_group2$column <- factor(df_error_group2$column,levels=unique(factor(df_error_group2$column)))

color_x_group1_idx <- (color_table$color)[which(c(1:160)%%4==1)]
color_y_group1_idx <- (color_table$color)[which(c(1:160)%%4==1)]
color_x_group2_idx <- (color_table$color)[which(c(1:160)%%4==1)]
color_y_group2_idx <- (color_table$color)[which(c(1:160)%%4==1)]
color_x_group1_error_idx <- (color_table$color)[which(c(1:160)%%4==1)]
color_y_group1_error_idx <- (color_table$color)[which(c(1:160)%%4==1)]
color_x_group2_error_idx <- (color_table$color)[which(c(1:160)%%4==1)]
color_y_group2_error_idx <- (color_table$color)[which(c(1:160)%%4==1)]

Group1_joint_tr <- ggplot(df_tr_joint_group1,
                          aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile("directed,joint,group 1")+
  scale_fill_gradient2(name="Value",low="darkslateblue",high="firebrick4")+
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_tr_joint_group1$row)), 
                   breaks = unique(as.factor(df_tr_joint_group1$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_tr_joint_group1$column)), 
                   breaks = unique(as.factor(df_tr_joint_group1$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,
                                 color=color_x_group1_idx),
        axis.text.y=element_text(size=8,
                                 color=color_y_group1_idx))

Group2_joint_tr <- ggplot(df_tr_joint_group2,
                          aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile("directed,joint,group 2")+
  scale_fill_gradient2(name="Value",low="darkslateblue",high="firebrick4")+
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_tr_joint_group2$row)), 
                   breaks = unique(as.factor(df_tr_joint_group2$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_tr_joint_group2$column)), 
                   breaks = unique(as.factor(df_tr_joint_group2$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,
                                 color=color_x_group1_idx),
        axis.text.y=element_text(size=8,
                                 color=color_y_group1_idx))

Group1_indiv_tr <- ggplot(df_tr_indiv_group1,
                          aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile("directed,group individual,group 1")+
  scale_fill_gradient2(name="Value",low="darkslateblue",high="firebrick4")+
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_tr_indiv_group1$row)), 
                   breaks = unique(as.factor(df_tr_indiv_group1$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_tr_indiv_group1$column)), 
                   breaks = unique(as.factor(df_tr_indiv_group1$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,
                                 color=color_x_group1_idx),
        axis.text.y=element_text(size=8,
                                 color=color_y_group1_idx))

Group2_indiv_tr <- ggplot(df_tr_indiv_group2,
                          aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile("directed,group individual,group 2")+
  scale_fill_gradient2(name="Value",low="darkslateblue",high="firebrick4")+
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_tr_indiv_group2$row)), 
                   breaks = unique(as.factor(df_tr_indiv_group2$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_tr_indiv_group2$column)), 
                   breaks = unique(as.factor(df_tr_indiv_group2$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,
                                 color=color_x_group1_idx),
        axis.text.y=element_text(size=8,
                                 color=color_y_group1_idx))


Group1_joint_error <- ggplot(df_error_joint_group1,
                             aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile("contemporaneous,joint,group 1")+
  scale_fill_gradient2(name="Value",low="darkslateblue",high="firebrick4")+
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_error_joint_group1$row)), 
                   breaks = unique(as.factor(df_error_joint_group1$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_error_joint_group1$column)), 
                   breaks = unique(as.factor(df_error_joint_group1$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,
                                 color=color_x_group1_idx),
        axis.text.y=element_text(size=8,
                                 color=color_y_group1_idx))

Group2_joint_error <- ggplot(df_error_joint_group2,
                             aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile("contemporaneous,joint,group 2")+
  scale_fill_gradient2(name="Value",low="darkslateblue",high="firebrick4")+
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_error_joint_group2$row)), 
                   breaks = unique(as.factor(df_error_joint_group2$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_error_joint_group2$column)), 
                   breaks = unique(as.factor(df_error_joint_group2$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,
                                 color=color_x_group1_idx),
        axis.text.y=element_text(size=8,
                                 color=color_y_group1_idx))

Group1_indiv_error <- ggplot(df_error_indiv_group1,
                             aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile("contemporaneous,group individual,group 1")+
  scale_fill_gradient2(name="Value",low="darkslateblue",high="firebrick4")+
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_error_indiv_group1$row)), 
                   breaks = unique(as.factor(df_error_indiv_group1$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_error_indiv_group1$column)), 
                   breaks = unique(as.factor(df_error_indiv_group1$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,
                                 color=color_x_group1_idx),
        axis.text.y=element_text(size=8,
                                 color=color_y_group1_idx))

Group2_indiv_error <- ggplot(df_error_indiv_group2,
                             aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile("contemporaneous,group individual,group 2")+
  scale_fill_gradient2(name="Value",low="darkslateblue",high="firebrick4")+
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_error_indiv_group2$row)), 
                   breaks = unique(as.factor(df_error_indiv_group2$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_error_indiv_group2$column)), 
                   breaks = unique(as.factor(df_error_indiv_group2$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,
                                 color=color_x_group1_idx),
        axis.text.y=element_text(size=8,
                                 color=color_y_group1_idx))



par(mar=c(0,0,0,0))
ggarrange(Group1_joint_tr,Group2_joint_tr,Group1_indiv_tr,Group2_indiv_tr,
          nrow=2,ncol=2,common.legend=TRUE,legend="bottom")
ggarrange(Group1_joint_error,Group2_joint_error,Group1_indiv_error,Group2_indiv_error,
          nrow=2,ncol=2,common.legend=TRUE,legend="bottom")




Group1_tr <- ggplot(df_tr_group1,
                    aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile()+
  scale_fill_gradient2(name="Value",low="darkslateblue",high="firebrick4")+
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_tr_group1$row)), 
                   breaks = unique(as.factor(df_tr_group1$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_tr_group1$column)), 
                   breaks = unique(as.factor(df_tr_group1$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,
                                 color=color_x_group1_idx),
        axis.text.y=element_text(size=8,
                                 color=color_y_group1_idx))

Group2_tr <- ggplot(df_tr_group2,
                    aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile()+
  scale_fill_gradient2(name="Value",low="darkslateblue",high="firebrick4")+
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_tr_group2$row)), 
                   breaks = unique(as.factor(df_tr_group2$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_tr_group2$column)), 
                   breaks = unique(as.factor(df_tr_group2$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,
                                 color=color_x_group1_idx),
        axis.text.y=element_text(size=8,
                                 color=color_y_group1_idx))

par(mar=c(0,0,0,0))
ggarrange(Group1_tr,Group2_tr,
          nrow=1,ncol=2,common.legend=TRUE,legend="bottom")




Group1_error <- ggplot(df_error_group1,
                       aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile()+
  scale_fill_gradient2(name="Value",low="darkslateblue",high="firebrick4")+
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_error_group1$row)), 
                   breaks = unique(as.factor(df_error_group1$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_error_group1$column)), 
                   breaks = unique(as.factor(df_error_group1$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,
                                 color=color_x_group1_idx),
        axis.text.y=element_text(size=8,
                                 color=color_y_group1_idx))

Group2_error <- ggplot(df_error_group2,
                       aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile()+
  scale_fill_gradient2(name="Value",low="darkslateblue",high="firebrick4")+
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_error_group2$row)), 
                   breaks = unique(as.factor(df_error_group2$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_error_group2$column)), 
                   breaks = unique(as.factor(df_error_group2$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,
                                 color=color_x_group1_idx),
        axis.text.y=element_text(size=8,
                                 color=color_y_group1_idx))

par(mar=c(0,0,0,0))
ggarrange(Group1_error,Group2_error,
          nrow=1,ncol=2,common.legend=TRUE,legend="bottom")

# Forecast
MSFE_GRIDY <- array(0,dim=c(2,H))
MSFE_DSCA <- array(0,dim=c(2,H))
MSFE_PCA <- array(0,dim=c(2,H))

for (h in 1:H){
  for (k in 1:(length(K1_idx)+length(K2_idx))){
    
    T_k <- dim(X_scale_final[[k]])[2]
    
    if (k <= length(K1_idx)){
      MSFE_GRIDY[1,h] <- norm(X_scale_final_ext[[k]][,(T_k+1):(T_k+h)] - t(X_hat_GRIDY[(1:h),,k]),"2")^2/(2*160*length(K1_idx))
      
      MSFE_DSCA[1,h] <- norm(X_scale_final_ext[[k]][,(T_k+1):(T_k+h)] - t(X_hat_DSCA[(1:h),,k]),"2")^2/(2*160*length(K1_idx))
      
      MSFE_PCA[1,h] <- norm(X_scale_final_ext[[k]][,(T_k+1):(T_k+h)] - t(X_hat_PCA[(1:h),,k]),"2")^2/(2*160*length(K1_idx))
    }else{
      MSFE_GRIDY[2,h] <- norm(X_scale_final_ext[[k]][,(T_k+1):(T_k+h)] - t(X_hat_GRIDY[(1:h),,k]),"2")^2/(2*160*length(K2_idx))
      
      MSFE_DSCA[2,h] <- norm(X_scale_final_ext[[k]][,(T_k+1):(T_k+h)] - t(X_hat_DSCA[(1:h),,k]),"2")^2/(2*160*length(K2_idx))
      
      MSFE_PCA[2,h] <- norm(X_scale_final_ext[[k]][,(T_k+1):(T_k+h)] - t(X_hat_PCA[(1:h),,k]),"2")^2/(2*160*length(K2_idx))
    }
  }
}

sqrt(MSFE_GRIDY)
sqrt(MSFE_DSCA)
sqrt(MSFE_PCA)




