rm(list = ls())

# Packages required
library(readxl)
library(mvtnorm)
library(Matrix)
library(combinat)

# For SCA and GICA:
library(multiway)
source('gica.R')
library(ica)

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

# Change the directory folder to where the ABIDE data is stored:
directory_folder <- "D:/GRIDY/data"
for (kk in 1:Tot_subject){
  
  file_name <- paste0(file_name_table[kk,2],"_ccs_filt_noglobal_rois_dosenbach160")
  file_directory <- paste0(directory_folder,"/",file_name,".csv")
  data_entry <- read.csv(file_directory)
  
  if (dim(data_entry)[2] == 160){
    list_X[[kk]] <- data_entry
    list_scaled_X[[kk]] <- t( scale(list_X[[kk]], center = TRUE, scale = TRUE) )
    list_scaled_X[[kk]][is.na(list_scaled_X[[kk]])] <- 0
    list_scaled_X[[kk]] <- list_scaled_X[[kk]]
  }else{
    list_X[[kk]] <- data_entry[,-dim(data_entry)[2]] 
    list_scaled_X[[kk]] <- t( scale(list_X[[kk]], center = TRUE, scale = TRUE) )
    list_scaled_X[[kk]][is.na(list_scaled_X[[kk]])] <- 0
    list_scaled_X[[kk]] <- list_scaled_X[[kk]]
  }
}

avg_T <- 0
for (kk in 1:Tot_subject){
  avg_T <- avg_T + dim(list_X[[kk]])[1]/Tot_subject
}
avg_T 


# -----------------------------------------------------------------------------#
# 0. Initial rank selection
# -----------------------------------------------------------------------------#
r_hat <- vector("numeric",length(Tot_subject))
d <- 160;
for (kk in 1:Tot_subject){
  cnt <- cnt + 1
  time_start <- Sys.time()
  T_k <- dim(list_scaled_X[[kk]])[2]
  
  possibleError <- tryCatch(
    r_hat[cnt] <- Rotational_bootstrap(d,T_k,t(list_scaled_X[[kk]]),cull=0.5),
    
    error = function(k){
      print(kk,"th subject is failed to do Boostrap")
    }
  )
  
  cat(kk,"th subject is completed.","Time taken is",Sys.time() - time_start,
      "and estimated rank is",r_hat[cnt],"\n")
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

X_ext <- list()
X_scale_ext <- list()
scaler <- list()

extract_idx1 <- df_index_final[which(df_index_final$group == 1),"subject"]
extract_idx2 <- df_index_final[which(df_index_final$group == 2),"subject"]

cnt_group1 <- 0
for (kk in extract_idx1){
  cnt_group1 <- cnt_group1+1
  X_ext[[cnt_group1]] <- list_X[[kk]]
  X_scale_ext[[cnt_group1]] <- list_scaled_X[[kk]]
  scaler[[cnt_group1]] <- attr(list_scaled_X[[kk]],"scaled:scale")
}

cnt_group2 <- 0
for (kk in extract_idx2){
  
  cnt_group2 <- cnt_group2+1
  X_ext[[cnt_group1+cnt_group2]] <- list_X[[kk]]
  X_scale_ext[[cnt_group1+cnt_group2]] <- list_scaled_X[[kk]]
  scaler[[cnt_group1+cnt_group2]] <- attr(list_scaled_X[[kk]],"scaled:scale")
}

K1_idx <- c(1:cnt_group1)
K2_idx <- cnt_group1+c(1:cnt_group2)

#-----------------------------------------------------------------#
# 1. Run AJIVE (GRIDY,SCA_P,GICA):
#-----------------------------------------------------------------#
r1_total <- 2
r2_total <- 2
initial_ranks <- c(rep(r1_total,length(K1_idx)),rep(r2_total,length(K2_idx)))

AJIVE_decomp <- ajive(X_scale_ext,
                     initial_signal_ranks = initial_ranks,
                     n_wedin_samples = 1000,
                     n_rand_dir_samples = 1000,
                     full = TRUE, joint_rank = NA)


AJIVE_decomp$joint_rank
par(mfrow=c(1,2))
hist(mapply(x=K1_idx,function(x){AJIVE_decomp$block_decomps[[x]]$individual$rank}),
     main="Group individual ranks for group 1",xlab=TeX(r'($\hat{r}_{G}$)'),
     xlim=c(0,4),ylim=c(0,120),col="blue",border="blue")
hist(mapply(x=K2_idx,function(x){AJIVE_decomp$block_decomps[[x]]$individual$rank}),
     main="Group individual ranks for group 2",xlab=TeX(r'($\hat{r}_{G}$)'),
     xlim=c(0,4),ylim=c(0,120),col="red",border="red")



# Exclude subjects whose group individual rank is zero
Joint_block<- list()
Group1_block <- list()
Group2_block <- list()

exclude_idx <- which(c(mapply(x=K1_idx,function(x){AJIVE_decomp$block_decomps[[x]]$individual$rank}),
                       mapply(x=K2_idx,function(x){AJIVE_decomp$block_decomps[[x]]$individual$rank})) == 0)
cnt_effect_tot_idx <- 0 
cnt_effect_g1_idx <- 0
cnt_effect_g2_idx <- 0

tmp_X_ext <- list()
tmp_X_scale_ext <- list()
tmp_scaler <- list()
for (kk in 1:(length(K1_idx)+length(K2_idx))){
  
  if (is.element(kk,exclude_idx)){
    next
  }else{
    cnt_effect_tot_idx <- cnt_effect_tot_idx +1
    
    tmp_X_ext[[cnt_effect_tot_idx]] <- X_ext[[kk]]
    tmp_X_scale_ext[[cnt_effect_tot_idx]] <- X_scale_ext[[kk]]
    tmp_scaler[[cnt_effect_tot_idx]] <- scaler[[kk]]
    
    Joint_block[[cnt_effect_tot_idx]] <- t(AJIVE_decomp[["block_decomps"]][[kk]]$joint$full)
    
    if ( kk <= length(K1_idx) ){
      cnt_effect_g1_idx <- cnt_effect_g1_idx +1
      Group1_block[[cnt_effect_g1_idx]] <- t(AJIVE_decomp[["block_decomps"]][[kk]]$individual$full)
      
    }else{
      cnt_effect_g2_idx <- cnt_effect_g2_idx +1
      Group2_block[[cnt_effect_g2_idx]] <- t(AJIVE_decomp[["block_decomps"]][[kk]]$individual$full)
    }
  }
}

X_ext <- tmp_X_ext
X_scale_ext <- tmp_X_scale_ext
scaler <- tmp_scaler
K1_idx <- c(1:cnt_effect_g1_idx)
K2_idx <- cnt_effect_g1_idx+c(1:cnt_effect_g2_idx)
rm(tmp_X_ext,tmp_X_scale_ext,tmp_scaler,cnt_effect_tot_idx,cnt_effect_g1_idx,cnt_effect_g2_idx)


#-----------------------------------------------------------------#
# 2. Run SCA_PF2, SCA_P and GICA:
#-----------------------------------------------------------------#
rJ <- AJIVE_decomp$joint_rank
rG1 <- 1
rG2 <- 1

SCA_PF2_Joint <- sca(Joint_block, nfac = rJ, type="sca-pf2", verbose = FALSE)
SCA_PF2_Group1 <- sca(Group1_block,nfac = rG1, type="sca-pf2", verbose = FALSE)
SCA_PF2_Group2 <- sca(Group2_block,nfac = rG2, type="sca-pf2", verbose = FALSE)

SCA_P_Joint <- sca(Joint_block, nfac = rJ,type="sca-p", rotation="varimax",verbose = FALSE)
SCA_P_Group1 <- sca(Group1_block, nfac = rG1,type="sca-p",verbose = FALSE)
SCA_P_Group2 <- sca(Group2_block, nfac = rG2,type="sca-p",verbose = FALSE)

GICA_Joint <- gica(Joint_block, nc = rJ, dual.reg = FALSE, center = FALSE)
GICA_Group1 <- gica(Group1_block, nc = rG1, dual.reg = FALSE, center = FALSE)
GICA_Group2 <- gica(Group2_block, nc = rG2, dual.reg = FALSE, center = FALSE)


#-----------------------------------------------------------------#
# 3. Construct dynamics (GRIDY,SCA_P,GICA):
#-----------------------------------------------------------------#
GRIDY_refitted <- factor_regression(K1_idx,K2_idx,rJ,rG1,rG2,
                                    X_scale_ext,SCA_PF2_Joint,SCA_PF2_Group1,SCA_PF2_Group2)
GRIDY_YK <- YK_compute(K1_idx,K2_idx,rJ,rG1,rG2,X_scale_ext,GRIDY_refitted)

SCA_P_refitted <- factor_regression(K1_idx,K2_idx,rJ,rG1,rG2,
                                    X_scale_ext,SCA_P_Joint,SCA_P_Group1,SCA_P_Group2)
SCA_P_YK <- YK_compute(K1_idx,K2_idx,rJ,rG1,rG2,X_scale_ext,SCA_P_refitted)

GICA_refitted <- factor_regression(K1_idx,K2_idx,rJ,rG1,rG2,
                                    X_scale_ext,GICA_Joint,GICA_Group1,GICA_Group2)
GICA_YK <- YK_compute(K1_idx,K2_idx,rJ,rG1,rG2,X_scale_ext,GICA_refitted)


save("rev_application_new.RData")

#-----------------------------------------------------------------#
# 4. Evaluation & Illustration
#-----------------------------------------------------------------#
# # To use the application completed data:
# load("./Data_application.RData")

# Packages required
library(readxl)
library(mvtnorm)
library(Matrix)
library(combinat)
library(ggplot2)
library(patchwork)
library(reshape2)
library(latex2exp)
library(ggpubr)
library(RColorBrewer)
library(cowplot)
library(grid)

# Label setup:
color_table <- data.frame(color=c(rep("red",length(which(voxel_name_table$Net=="dfm"))),
                                  rep("blue",length(which(voxel_name_table$Net=="fp"))),
                                  rep("green",length(which(voxel_name_table$Net=="co"))),
                                  rep("purple",length(which(voxel_name_table$Net=="sm"))),
                                  rep("orange",length(which(voxel_name_table$Net=="cb"))),
                                  rep("brown",length(which(voxel_name_table$Net=="oc")))),
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
                             which(color_table$color=="brown")),]

# loadings
df_B_joint_1 <- data.frame(value = c(SCA_PF2_Joint$B[color_table$ROI,1],
                                     SCA_P_Joint$B[color_table$ROI,1],
                                     GICA_Joint$B[color_table$ROI,1]),
                           ROI = rep(color_table$ROI,3),
                           color = rep(color_table$color,3),
                           factor = rep(rep(1,160),3),
                           structure=rep(rep("Joint",160),3),
                           method=c(rep("GRIDY",160),rep("SCA-P",160),rep("GICA",160)))
df_B_joint_2 <- data.frame(value = c(SCA_PF2_Joint$B[color_table$ROI,2],
                                     SCA_P_Joint$B[color_table$ROI,2],
                                     GICA_Joint$B[color_table$ROI,2]),
                           ROI = color_table$ROI,
                           color = color_table$color,
                           factor = rep(rep(2,160),3),
                           structure=rep(rep("Joint",160),3),
                           method=c(rep("GRIDY",160),rep("SCA-P",160),rep("GICA",160)))
df_B_joint <- rbind(df_B_joint_1,df_B_joint_2)
df_B_joint$ROI <- factor(df_B_joint$ROI,levels=color_table$ROI)
df_B_joint$factor <- factor(df_B_joint$factor,levels=c(1,2),labels=c("Loadings 1","Loadings 2"))


df_B_group1_1 <- data.frame(value = c(SCA_PF2_Group1$B[color_table$ROI],
                                      SCA_P_Group1$B[color_table$ROI],
                                      GICA_Group1$B[color_table$ROI]),
                            ROI = rep(color_table$ROI,3),
                            color = rep(color_table$color,3),
                            factor = rep(rep(1,160),3),
                            structure=rep(rep("Group Individual 1",160),3),
                            method=c(rep("GRIDY",160),rep("SCA-P",160),rep("GICA",160)))
df_B_group1 <- rbind(df_B_group1_1)
df_B_group1$ROI <- factor(df_B_group1$ROI,levels=color_table$ROI)
df_B_group1$factor <- factor(df_B_group1$factor,levels=1,labels="Loadings 1")


df_B_group2_1 <- data.frame(value = c(SCA_PF2_Group2$B[color_table$ROI],
                                      SCA_P_Group2$B[color_table$ROI],
                                      GICA_Group2$B[color_table$ROI]),
                            ROI = rep(color_table$ROI,3),
                            color = rep(color_table$color,3),
                            factor = rep(rep(1,160),3),
                            structure=rep(rep("Group Individual 2",160),3),
                            method=c(rep("GRIDY",160),rep("SCA-P",160),rep("GICA",160)))
df_B_group2 <- rbind(df_B_group2_1)
df_B_group2$ROI <- factor(df_B_group2$ROI,levels=color_table$ROI)
df_B_group2$factor <- factor(df_B_group2$factor,levels=1,labels="Loadings 1")

df_B <- rbind(df_B_joint,df_B_group1,df_B_group2)
df_B$structure <- factor(df_B$structure,levels=c("Joint","Group Individual 1","Group Individual 2"))
df_B$method <- factor(df_B$method,levels=c("GRIDY","SCA-P","GICA"))
factor.labs <- c("Loadings 1","Loadings 2")
names(factor.labs) <- c("Loadings 1","Loadings 2")
structure.labs <- c("Joint", "Group Individual 1", "Group Individual 2")
names(structure.labs) <- c("Joint", "Group Individual 1", "Group Individual 2")
method.labs <- c("GRIDY", "SCA-P", "GICA")
names(method.labs) <- c("GRIDY", "SCA-P", "GICA")


# Plotting:
par(mar=c(0,0,0,0))
ggplot(df_B,aes(x=as.factor(ROI),y=value,fill=color)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="Network",
                    labels=c("dfm","fp","co","sm","cb","oc"),
                    breaks=c("red","blue","green","purple","orange","brown"),
                    values=c("red","blue","green","purple","orange","brown")) +
  scale_y_continuous(name="Value") + 
  theme(legend.position="bottom",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + guides(fill = guide_legend(nrow = 1)) +
  facet_grid(structure+factor~method,
             labeller=labeller(structure=structure.labs,factor=factor.labs,method=method.labs))



# Covariance of factors
Cov_joint_GRIDY <- Cov_joint_SCA_P <- Cov_joint_GICA <- array(NA,dim=c(2,2,(length(K1_idx) + length(K2_idx))))
Cov_group1_GRIDY <- Cov_group1_SCA_P  <- Cov_group1_GICA <- vector("numeric",length(K1_idx))
Cov_group2_GRIDY <- Cov_group2_SCA_P  <- Cov_group2_GICA <- vector("numeric",length(K2_idx))
for (kk in 1:(length(K1_idx) + length(K2_idx))){
  
  T_k <- dim(GRIDY_refitted$factor_joint[[kk]])[1]
  
  Cov_joint_GRIDY[,,kk] <- t(scale(GRIDY_refitted$factor_joint[[kk]],center=TRUE,scale=FALSE)) %*% scale(GRIDY_refitted$factor_joint[[kk]],center=TRUE,scale=FALSE) / T_k
  Cov_joint_SCA_P[,,kk] <- t(scale(SCA_P_refitted$factor_joint[[kk]],center=TRUE,scale=FALSE)) %*% scale(SCA_P_refitted$factor_joint[[kk]],center=TRUE,scale=FALSE) / T_k
  Cov_joint_GICA[,,kk] <- t(scale(GICA_refitted$factor_joint[[kk]],center=TRUE,scale=FALSE)) %*% scale(GICA_refitted$factor_joint[[kk]],center=TRUE,scale=FALSE) / T_k
  
  if (kk <= length(K1_idx)){
    Cov_group1_GRIDY[kk] <- t(scale(GRIDY_refitted$factor_group1[[kk]],center=TRUE,scale=FALSE)) %*% scale(GRIDY_refitted$factor_group1[[kk]],center=TRUE,scale=FALSE) / T_k
    Cov_group1_SCA_P[kk] <- t(scale(SCA_P_refitted$factor_group1[[kk]],center=TRUE,scale=FALSE)) %*% scale(SCA_P_refitted$factor_group1[[kk]],center=TRUE,scale=FALSE) / T_k
    Cov_group1_GICA[kk] <- t(scale(GICA_refitted$factor_group1[[kk]],center=TRUE,scale=FALSE)) %*% scale(GICA_refitted$factor_group1[[kk]],center=TRUE,scale=FALSE) / T_k
    
  }else{
    Cov_group2_GRIDY[(kk-length(K1_idx))] <- t(scale(GRIDY_refitted$factor_group2[[(kk-length(K1_idx))]],center=TRUE,scale=FALSE)) %*% scale(GRIDY_refitted$factor_group2[[(kk-length(K1_idx))]],center=TRUE,scale=FALSE) / T_k
    Cov_group2_SCA_P[(kk-length(K1_idx))] <- t(scale(SCA_P_refitted$factor_group2[[(kk-length(K1_idx))]],center=TRUE,scale=FALSE)) %*% scale(SCA_P_refitted$factor_group2[[(kk-length(K1_idx))]],center=TRUE,scale=FALSE) / T_k
    Cov_group2_GICA[(kk-length(K1_idx))] <- t(scale(GICA_refitted$factor_group2[[(kk-length(K1_idx))]],center=TRUE,scale=FALSE)) %*% scale(GICA_refitted$factor_group2[[(kk-length(K1_idx))]],center=TRUE,scale=FALSE) / T_k
    
  }
}


GRIDY_Cov_group1_joint <- data.frame(value = as.vector(Cov_joint_GRIDY[,,K1_idx]),
                                     index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K1_idx)),
                                     factor = rep(1,length(K1_idx)),
                                     method = rep("GRIDY",length(K1_idx)),
                                     structure = rep("Joint",length(K1_idx)))
GRIDY_Cov_group2_joint <- data.frame(value = as.vector(Cov_joint_GRIDY[,,K2_idx]),
                                     index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K2_idx)),
                                     factor = rep(2,length(K2_idx)),
                                     method = rep("GRIDY",length(K2_idx)),
                                     structure = rep("Joint",length(K2_idx)))
GRIDY_Cov_group1_individual <- data.frame(value = Cov_group1_GRIDY,
                                          index = rep(c("(1,1)"),length(K1_idx)),
                                          factor = rep(1,length(K1_idx)),
                                          method = rep("GRIDY",length(K1_idx)),
                                          structure = rep("Group Individual 1",length(K1_idx)))
GRIDY_Cov_group2_individual <- data.frame(value = Cov_group2_GRIDY,
                                          index = rep(c("(1,1)"),length(K2_idx)),
                                          factor = rep(1,length(K2_idx)),
                                          method = rep("GRIDY",length(K2_idx)),
                                          structure = rep("Group Individual 2",length(K2_idx)))

SCA_P_Cov_group1_joint <- data.frame(value = as.vector(Cov_joint_SCA_P[,,K1_idx]),
                                     index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K1_idx)),
                                     factor = rep(1,length(K1_idx)),
                                     method = rep("SCA-P",length(K1_idx)),
                                     structure = rep("Joint",length(K1_idx)))
SCA_P_Cov_group2_joint <- data.frame(value = as.vector(Cov_joint_SCA_P[,,K2_idx]),
                                     index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K2_idx)),
                                     factor = rep(2,length(K2_idx)),
                                     method = rep("SCA-P",length(K2_idx)),
                                     structure = rep("Joint",length(K2_idx)))
SCA_P_Cov_group1_individual <- data.frame(value = Cov_group1_SCA_P,
                                          index = rep(c("(1,1)"),length(K1_idx)),
                                          factor = rep(1,length(K1_idx)),
                                          method = rep("SCA-P",length(K1_idx)),
                                          structure = rep("Group Individual 1",length(K1_idx)))
SCA_P_Cov_group2_individual <- data.frame(value = Cov_group2_SCA_P,
                                          index = rep(c("(1,1)"),length(K2_idx)),
                                          factor = rep(1,length(K2_idx)),
                                          method = rep("SCA-P",length(K2_idx)),
                                          structure = rep("Group Individual 2",length(K2_idx)))

GICA_Cov_group1_joint <- data.frame(value = as.vector(Cov_joint_GICA[,,K1_idx]),
                                    index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K1_idx)),
                                    factor = rep(1,length(K1_idx)),
                                    method = rep("GICA",length(K1_idx)),
                                    structure = rep("Joint",length(K1_idx)))
GICA_Cov_group2_joint <- data.frame(value = as.vector(Cov_joint_GICA[,,K2_idx]),
                                    index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K2_idx)),
                                    factor = rep(2,length(K2_idx)),
                                    method = rep("GICA",length(K2_idx)),
                                    structure = rep("Joint",length(K2_idx)))
GICA_Cov_group1_individual <- data.frame(value = Cov_group1_GICA,
                                         index = rep(c("(1,1)"),length(K1_idx)),
                                         factor = rep(1,length(K1_idx)),
                                         method = rep("GICA",length(K1_idx)),
                                         structure = rep("Group Individual 1",length(K1_idx)))
GICA_Cov_group2_individual <- data.frame(value = Cov_group2_GICA,
                                         index = rep(c("(1,1)"),length(K2_idx)),
                                         factor = rep(1,length(K2_idx)),
                                         method = rep("GICA",length(K2_idx)),
                                         structure = rep("Group Individual 2",length(K2_idx)))

df_Cov <- rbind(GRIDY_Cov_group1_joint,GRIDY_Cov_group2_joint,
                GRIDY_Cov_group1_individual,GRIDY_Cov_group2_individual,
                SCA_P_Cov_group1_joint,SCA_P_Cov_group2_joint,
                SCA_P_Cov_group1_individual,SCA_P_Cov_group2_individual,
                GICA_Cov_group1_joint,GICA_Cov_group2_joint,
                GICA_Cov_group1_individual,GICA_Cov_group2_individual)
df_Cov$structure <- factor(df_Cov$structure,levels=c("Joint","Group Individual 1","Group Individual 2"))
df_Cov$method <- factor(df_Cov$method,levels=c("GRIDY","SCA-P","GICA"))
df_Cov$index <- factor(df_Cov$index,levels=c("(1,1)","(2,1)","(1,2)","(2,2)"))
df_Cov$factor <- factor(df_Cov$factor,levels=c(1,2),labels=c("Factor 1","Factor 2"))

str.labs <- c("Joint","Group Individual 1","Group Individual 2")
names(str.labs) <- c("Joint","Group Individual 1","Group Individual 2")
method.labs <- c("Joint","Group Individual 1","Group Individual 2")
names(method.labs) <- c("Joint","Group Individual 1","Group Individual 2")
factor.labs <- c("Factor 1","Factor 2")
names(factor.labs) <- c("Factor 1","Factor 2")


# Transition matrices
GRIDY_Psi_group1_joint <- data.frame(value = as.vector(GRIDY_YK$Psi_joint_hat[,,K1_idx]),
                                  index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K1_idx)),
                                  factor = rep(1,length(K1_idx)),
                                  method = rep("GRIDY",length(K1_idx)),
                                  structure = rep("Joint",length(K1_idx)))
GRIDY_Psi_group2_joint <- data.frame(value = as.vector(GRIDY_YK$Psi_joint_hat[,,K2_idx]),
                                  index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K2_idx)),
                                  factor = rep(2,length(K2_idx)),
                                  method = rep("GRIDY",length(K2_idx)),
                                  structure = rep("Joint",length(K2_idx)))
GRIDY_Psi_group1_individual <- data.frame(value = as.vector(GRIDY_YK$Psi_indiv1_hat[,,K1_idx]),
                                       index = rep(c("(1,1)"),length(K1_idx)),
                                       factor = rep(1,length(K1_idx)),
                                       method = rep("GRIDY",length(K1_idx)),
                                       structure = rep("Group Individual 1",length(K1_idx)))
GRIDY_Psi_group2_individual <- data.frame(value = as.vector(GRIDY_YK$Psi_indiv2_hat[,,1:length(K2_idx)]),
                                       index = rep(c("(1,1)"),length(K2_idx)),
                                       factor = rep(1,length(K2_idx)),
                                       method = rep("GRIDY",length(K2_idx)),
                                       structure = rep("Group Individual 2",length(K2_idx)))

SCA_P_Psi_group1_joint <- data.frame(value = as.vector(SCA_P_YK$Psi_joint_hat[,,K1_idx]),
                                     index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K1_idx)),
                                     factor = rep(1,length(K1_idx)),
                                     method = rep("SCA-P",length(K1_idx)),
                                     structure = rep("Joint",length(K1_idx)))
SCA_P_Psi_group2_joint <- data.frame(value = as.vector(SCA_P_YK$Psi_joint_hat[,,K2_idx]),
                                     index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K2_idx)),
                                     factor = rep(2,length(K2_idx)),
                                     method = rep("SCA-P",length(K2_idx)),
                                     structure = rep("Joint",length(K2_idx)))
SCA_P_Psi_group1_individual <- data.frame(value = as.vector(SCA_P_YK$Psi_indiv1_hat[,,K1_idx]),
                                          index = rep(c("(1,1)"),length(K1_idx)),
                                          factor = rep(1,length(K1_idx)),
                                          method = rep("SCA-P",length(K1_idx)),
                                          structure = rep("Group Individual 1",length(K1_idx)))
SCA_P_Psi_group2_individual <- data.frame(value = as.vector(SCA_P_YK$Psi_indiv2_hat[,,1:length(K2_idx)]),
                                          index = rep(c("(1,1)"),length(K2_idx)),
                                          factor = rep(1,length(K2_idx)),
                                          method = rep("SCA-P",length(K2_idx)),
                                          structure = rep("Group Individual 2",length(K2_idx)))

GICA_Psi_group1_joint <- data.frame(value = as.vector(GICA_YK$Psi_joint_hat[,,K1_idx]),
                                    index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K1_idx)),
                                    factor = rep(1,length(K1_idx)),
                                    method = rep("GICA",length(K1_idx)),
                                    structure = rep("Joint",length(K1_idx)))
GICA_Psi_group2_joint <- data.frame(value = as.vector(GICA_YK$Psi_joint_hat[,,K2_idx]),
                                    index = rep(c("(1,1)","(2,1)","(1,2)","(2,2)"),length(K2_idx)),
                                    factor = rep(2,length(K2_idx)),
                                    method = rep("GICA",length(K2_idx)),
                                    structure = rep("Joint",length(K2_idx)))
GICA_Psi_group1_individual <- data.frame(value = as.vector(GICA_YK$Psi_indiv1_hat[,,K1_idx]),
                                        index = rep(c("(1,1)"),length(K1_idx)),
                                        factor = rep(1,length(K1_idx)),
                                        method = rep("GICA",length(K1_idx)),
                                        structure = rep("Group Individual 1",length(K1_idx)))
GICA_Psi_group2_individual <- data.frame(value = as.vector(GICA_YK$Psi_indiv2_hat[,,1:length(K2_idx)]),
                                         index = rep(c("(1,1)"),length(K2_idx)),
                                         factor = rep(1,length(K2_idx)),
                                         method = rep("GICA",length(K2_idx)),
                                         structure = rep("Group Individual 2",length(K2_idx)))


df_Psi <- rbind(GRIDY_Psi_group1_joint,GRIDY_Psi_group2_joint,
                GRIDY_Psi_group1_individual,GRIDY_Psi_group2_individual,
                SCA_P_Psi_group1_joint,SCA_P_Psi_group2_joint,
                SCA_P_Psi_group1_individual,SCA_P_Psi_group2_individual,
                GICA_Psi_group1_joint,GICA_Psi_group2_joint,
                GICA_Psi_group1_individual,GICA_Psi_group2_individual)
df_Psi$structure <- factor(df_Psi$structure,levels=c("Joint","Group Individual 1","Group Individual 2"))
df_Psi$method <- factor(df_Psi$method,levels=c("GRIDY","SCA-P","GICA"))
df_Psi$index <- factor(df_Psi$index,levels=c("(1,1)","(2,1)","(1,2)","(2,2)"))
df_Psi$factor <- factor(df_Psi$factor,levels=c(1,2),labels=c("Factor 1","Factor 2"))

str.labs <- c("Joint","Group Individual 1","Group Individual 2")
names(str.labs) <- c("Joint","Group Individual 1","Group Individual 2")
method.labs <- c("Joint","Group Individual 1","Group Individual 2")
names(method.labs) <- c("Joint","Group Individual 1","Group Individual 2")
factor.labs <- c("Factor 1","Factor 2")
names(factor.labs) <- c("Factor 1","Factor 2")


# Plotting:
par(mar=c(0,0,0,0))
ggarrange(
ggplot(df_Cov,aes(y=value,x=index,fill=factor(index))) +
  geom_boxplot(cex=0.1) +
  scale_fill_manual(name="",values=colorRampPalette(brewer.pal(4,"Set1"))(4)) +
  scale_y_continuous(name="Value") +
  scale_x_discrete(name="Entries") +
  ggtitle("Covariance Matrices of Factors") +
  theme(legend.position="none") +
  facet_grid(structure+factor~method,
             scales = "free",
             labeller=labeller(structure=str.labs,method=method.labs,factor=factor.labs)),
ggplot(df_Psi,aes(y=value,x=index,fill=factor(index))) +
  geom_boxplot(cex=0.1) +
  scale_fill_manual(name="",values=colorRampPalette(brewer.pal(4,"Set1"))(4)) +
  scale_y_continuous(name="Value") +
  scale_x_discrete(name="Entries") +
  ggtitle("Transition Matrices of Factors") +
  theme(legend.position="none") +
  facet_grid(structure+factor~method,
             scales = "free",
             labeller=labeller(structure=str.labs,method=method.labs,factor=factor.labs)),
ncol=2,common.legend=TRUE,legend="bottom")


# R2 statistics
R2_group1_GRIDY <- R2_group1_SCA_P <- R2_group1_GICA <- array(NA,dim=c(160,3,length(K1_idx)))
R2_group2_GRIDY <- R2_group2_SCA_P <- R2_group2_GICA <- array(NA,dim=c(160,3,length(K2_idx)))
for (kk in 1:(length(K1_idx)+length(K2_idx))){
  
  X_block <- t(X_scale_ext[[kk]])
  
  GRIDY_joint_factor1 <- GRIDY_refitted$factor_joint[[kk]][,1] %*% t(SCA_PF2_Joint$B[,1])
  GRIDY_joint_factor2 <- GRIDY_refitted$factor_joint[[kk]][,2] %*% t(SCA_PF2_Joint$B[,2])
  SCA_P_joint_factor1 <- SCA_P_refitted$factor_joint[[kk]][,1] %*% t(SCA_P_Joint$B[,1])
  SCA_P_joint_factor2 <- SCA_P_refitted$factor_joint[[kk]][,2] %*% t(SCA_P_Joint$B[,2])
  GICA_joint_factor1 <- GICA_refitted$factor_joint[[kk]][,1] %*% t(GICA_Joint$B[,1])
  GICA_joint_factor2 <- GICA_refitted$factor_joint[[kk]][,2] %*% t(GICA_Joint$B[,2])
  
  if (kk <= length(K1_idx)){
    GRIDY_group_factor <- GRIDY_refitted$factor_group1[[kk]] %*% t(SCA_PF2_Group1$B)
    SCA_P_group_factor <- SCA_P_refitted$factor_group1[[kk]] %*% t(SCA_P_Group1$B)
    GICA_group_factor <- GICA_refitted$factor_group1[[kk]] %*% t(GICA_Group1$B)
    
    for (i in 1:160){
      R2_group1_GRIDY[i,1,kk] <- sum((GRIDY_joint_factor1[,i])^2) / sum((X_block[,i])^2)
      R2_group1_GRIDY[i,2,kk] <- sum((GRIDY_joint_factor2[,i])^2) / sum((X_block[,i])^2)
      R2_group1_GRIDY[i,3,kk] <- sum((GRIDY_group_factor[,i])^2) / sum((X_block[,i])^2)
      R2_group1_SCA_P[i,1,kk] <- sum((SCA_P_joint_factor1[,i])^2) / sum((X_block[,i])^2)
      R2_group1_SCA_P[i,2,kk] <- sum((SCA_P_joint_factor2[,i])^2) / sum((X_block[,i])^2)
      R2_group1_SCA_P[i,3,kk] <- sum((SCA_P_group_factor[,i])^2) / sum((X_block[,i])^2)
      R2_group1_GICA[i,1,kk] <- sum((GICA_joint_factor1[,i])^2) / sum((X_block[,i])^2)
      R2_group1_GICA[i,2,kk] <- sum((GICA_joint_factor2[,i])^2) / sum((X_block[,i])^2)
      R2_group1_GICA[i,3,kk] <- sum((GICA_group_factor[,i])^2) / sum((X_block[,i])^2)
    }
    
  }else{
    GRIDY_group_factor <- GRIDY_refitted$factor_group2[[(kk-length(K1_idx))]] %*% t(SCA_PF2_Group2$B)
    SCA_P_group_factor <- SCA_P_refitted$factor_group2[[(kk-length(K1_idx))]] %*% t(SCA_P_Group2$B)
    GICA_group_factor <- GICA_refitted$factor_group2[[(kk-length(K1_idx))]] %*% t(GICA_Group2$B)
    
    for (i in 1:160){
      R2_group2_GRIDY[i,1,(kk-length(K1_idx))] <- sum((GRIDY_joint_factor1[,i])^2) / sum((X_block[,i])^2)
      R2_group2_GRIDY[i,2,(kk-length(K1_idx))] <- sum((GRIDY_joint_factor2[,i])^2) / sum((X_block[,i])^2)
      R2_group2_GRIDY[i,3,(kk-length(K1_idx))] <- sum((GRIDY_group_factor[,i])^2) / sum((X_block[,i])^2)
      R2_group2_SCA_P[i,1,(kk-length(K1_idx))] <- sum((SCA_P_joint_factor1[,i])^2) / sum((X_block[,i])^2)
      R2_group2_SCA_P[i,2,(kk-length(K1_idx))] <- sum((SCA_P_joint_factor2[,i])^2) / sum((X_block[,i])^2)
      R2_group2_SCA_P[i,3,(kk-length(K1_idx))] <- sum((SCA_P_group_factor[,i])^2) / sum((X_block[,i])^2)
      R2_group2_GICA[i,1,(kk-length(K1_idx))] <- sum((GICA_joint_factor1[,i])^2) / sum((X_block[,i])^2)
      R2_group2_GICA[i,2,(kk-length(K1_idx))] <- sum((GICA_joint_factor2[,i])^2) / sum((X_block[,i])^2)
      R2_group2_GICA[i,3,(kk-length(K1_idx))] <- sum((GICA_group_factor[,i])^2) / sum((X_block[,i])^2)
    }
  }
}

R2_group1_factor1_GRIDY <- data.frame(t(R2_group1_GRIDY[color_table$ROI,1,]))
R2_group1_factor2_GRIDY <- data.frame(t(R2_group1_GRIDY[color_table$ROI,2,]))
R2_group1_factor3_GRIDY <- data.frame(t(R2_group1_GRIDY[color_table$ROI,3,]))
R2_group2_factor1_GRIDY <- data.frame(t(R2_group2_GRIDY[color_table$ROI,1,]))
R2_group2_factor2_GRIDY <- data.frame(t(R2_group2_GRIDY[color_table$ROI,2,]))
R2_group2_factor3_GRIDY <- data.frame(t(R2_group2_GRIDY[color_table$ROI,3,]))
R2_group1_factor1_SCA_P <- data.frame(t(R2_group1_SCA_P[color_table$ROI,1,]))
R2_group1_factor2_SCA_P <- data.frame(t(R2_group1_SCA_P[color_table$ROI,2,]))
R2_group1_factor3_SCA_P <- data.frame(t(R2_group1_SCA_P[color_table$ROI,3,]))
R2_group2_factor1_SCA_P <- data.frame(t(R2_group2_SCA_P[color_table$ROI,1,]))
R2_group2_factor2_SCA_P <- data.frame(t(R2_group2_SCA_P[color_table$ROI,2,]))
R2_group2_factor3_SCA_P <- data.frame(t(R2_group2_SCA_P[color_table$ROI,3,]))
R2_group1_factor1_GICA <- data.frame(t(R2_group1_GICA[color_table$ROI,1,]))
R2_group1_factor2_GICA <- data.frame(t(R2_group1_GICA[color_table$ROI,2,]))
R2_group1_factor3_GICA <- data.frame(t(R2_group1_GICA[color_table$ROI,3,]))
R2_group2_factor1_GICA <- data.frame(t(R2_group2_GICA[color_table$ROI,1,]))
R2_group2_factor2_GICA <- data.frame(t(R2_group2_GICA[color_table$ROI,2,]))
R2_group2_factor3_GICA <- data.frame(t(R2_group2_GICA[color_table$ROI,3,]))


colnames(R2_group1_factor1_GRIDY) <- color_table$ROI
colnames(R2_group1_factor2_GRIDY) <- color_table$ROI
colnames(R2_group1_factor3_GRIDY) <- color_table$ROI
colnames(R2_group2_factor1_GRIDY) <- color_table$ROI
colnames(R2_group2_factor2_GRIDY) <- color_table$ROI
colnames(R2_group2_factor3_GRIDY) <- color_table$ROI
colnames(R2_group1_factor1_SCA_P) <- color_table$ROI
colnames(R2_group1_factor2_SCA_P) <- color_table$ROI
colnames(R2_group1_factor3_SCA_P) <- color_table$ROI
colnames(R2_group2_factor1_SCA_P) <- color_table$ROI
colnames(R2_group2_factor2_SCA_P) <- color_table$ROI
colnames(R2_group2_factor3_SCA_P) <- color_table$ROI
colnames(R2_group1_factor1_GICA) <- color_table$ROI
colnames(R2_group1_factor2_GICA) <- color_table$ROI
colnames(R2_group1_factor3_GICA) <- color_table$ROI
colnames(R2_group2_factor1_GICA) <- color_table$ROI
colnames(R2_group2_factor2_GICA) <- color_table$ROI
colnames(R2_group2_factor3_GICA) <- color_table$ROI


df_R2_group1_factor1_GRIDY <- melt(R2_group1_factor1_GRIDY)
df_R2_group1_factor2_GRIDY <- melt(R2_group1_factor2_GRIDY)
df_R2_group1_factor3_GRIDY <- melt(R2_group1_factor3_GRIDY)
df_R2_group2_factor1_GRIDY <- melt(R2_group2_factor1_GRIDY)
df_R2_group2_factor2_GRIDY <- melt(R2_group2_factor2_GRIDY)
df_R2_group2_factor3_GRIDY <- melt(R2_group2_factor3_GRIDY)
df_R2_group1_factor1_SCA_P <- melt(R2_group1_factor1_SCA_P)
df_R2_group1_factor2_SCA_P <- melt(R2_group1_factor2_SCA_P)
df_R2_group1_factor3_SCA_P <- melt(R2_group1_factor3_SCA_P)
df_R2_group2_factor1_SCA_P <- melt(R2_group2_factor1_SCA_P)
df_R2_group2_factor2_SCA_P <- melt(R2_group2_factor2_SCA_P)
df_R2_group2_factor3_SCA_P <- melt(R2_group2_factor3_SCA_P)
df_R2_group1_factor1_GICA <- melt(R2_group1_factor1_GICA)
df_R2_group1_factor2_GICA <- melt(R2_group1_factor2_GICA)
df_R2_group1_factor3_GICA <- melt(R2_group1_factor3_GICA)
df_R2_group2_factor1_GICA <- melt(R2_group2_factor1_GICA)
df_R2_group2_factor2_GICA <- melt(R2_group2_factor2_GICA)
df_R2_group2_factor3_GICA <- melt(R2_group2_factor3_GICA)


df_R2_group1_factor1_GRIDY$color <- rep(color_table$color,each=length(K1_idx))
df_R2_group1_factor2_GRIDY$color <- rep(color_table$color,each=length(K1_idx))
df_R2_group1_factor3_GRIDY$color <- rep(color_table$color,each=length(K1_idx))
df_R2_group2_factor1_GRIDY$color <- rep(color_table$color,each=length(K2_idx))
df_R2_group2_factor2_GRIDY$color <- rep(color_table$color,each=length(K2_idx))
df_R2_group2_factor3_GRIDY$color <- rep(color_table$color,each=length(K2_idx))
df_R2_group1_factor1_SCA_P$color <- rep(color_table$color,each=length(K1_idx))
df_R2_group1_factor2_SCA_P$color <- rep(color_table$color,each=length(K1_idx))
df_R2_group1_factor3_SCA_P$color <- rep(color_table$color,each=length(K1_idx))
df_R2_group2_factor1_SCA_P$color <- rep(color_table$color,each=length(K2_idx))
df_R2_group2_factor2_SCA_P$color <- rep(color_table$color,each=length(K2_idx))
df_R2_group2_factor3_SCA_P$color <- rep(color_table$color,each=length(K2_idx))
df_R2_group1_factor1_GICA$color <- rep(color_table$color,each=length(K1_idx))
df_R2_group1_factor2_GICA$color <- rep(color_table$color,each=length(K1_idx))
df_R2_group1_factor3_GICA$color <- rep(color_table$color,each=length(K1_idx))
df_R2_group2_factor1_GICA$color <- rep(color_table$color,each=length(K2_idx))
df_R2_group2_factor2_GICA$color <- rep(color_table$color,each=length(K2_idx))
df_R2_group2_factor3_GICA$color <- rep(color_table$color,each=length(K2_idx))


df_R2_group1_factor1_GRIDY <- df_R2_group1_factor1_GRIDY[-which(is.infinite(df_R2_group1_factor1_GRIDY$value)),]
df_R2_group1_factor2_GRIDY <- df_R2_group1_factor2_GRIDY[-which(is.infinite(df_R2_group1_factor2_GRIDY$value)),]
df_R2_group1_factor3_GRIDY <- df_R2_group1_factor3_GRIDY[-which(is.infinite(df_R2_group1_factor3_GRIDY$value)),]
df_R2_group2_factor1_GRIDY <- df_R2_group2_factor1_GRIDY[-which(is.infinite(df_R2_group2_factor1_GRIDY$value)),]
df_R2_group2_factor2_GRIDY <- df_R2_group2_factor2_GRIDY[-which(is.infinite(df_R2_group2_factor2_GRIDY$value)),]
df_R2_group2_factor3_GRIDY <- df_R2_group2_factor3_GRIDY[-which(is.infinite(df_R2_group2_factor3_GRIDY$value)),]
df_R2_group1_factor1_SCA_P <- df_R2_group1_factor1_SCA_P[-which(is.infinite(df_R2_group1_factor1_SCA_P$value)),]
df_R2_group1_factor2_SCA_P <- df_R2_group1_factor2_SCA_P[-which(is.infinite(df_R2_group1_factor2_SCA_P$value)),]
df_R2_group1_factor3_SCA_P <- df_R2_group1_factor3_SCA_P[-which(is.infinite(df_R2_group1_factor3_SCA_P$value)),]
df_R2_group2_factor1_SCA_P <- df_R2_group2_factor1_SCA_P[-which(is.infinite(df_R2_group2_factor1_SCA_P$value)),]
df_R2_group2_factor2_SCA_P <- df_R2_group2_factor2_SCA_P[-which(is.infinite(df_R2_group2_factor2_SCA_P$value)),]
df_R2_group2_factor3_SCA_P <- df_R2_group2_factor3_SCA_P[-which(is.infinite(df_R2_group2_factor3_SCA_P$value)),]
df_R2_group1_factor1_GICA <- df_R2_group1_factor1_GICA[-which(is.infinite(df_R2_group1_factor1_GICA$value)),]
df_R2_group1_factor2_GICA <- df_R2_group1_factor2_GICA[-which(is.infinite(df_R2_group1_factor2_GICA$value)),]
df_R2_group1_factor3_GICA <- df_R2_group1_factor3_GICA[-which(is.infinite(df_R2_group1_factor3_GICA$value)),]
df_R2_group2_factor1_GICA <- df_R2_group2_factor1_GICA[-which(is.infinite(df_R2_group2_factor1_GICA$value)),]
df_R2_group2_factor2_GICA <- df_R2_group2_factor2_GICA[-which(is.infinite(df_R2_group2_factor2_GICA$value)),]
df_R2_group2_factor3_GICA <- df_R2_group2_factor3_GICA[-which(is.infinite(df_R2_group2_factor3_GICA$value)),]


df_R2_group1_factor1_GRIDY$factor <- rep(1,dim(df_R2_group1_factor1_GRIDY)[1])
df_R2_group1_factor2_GRIDY$factor <- rep(2,dim(df_R2_group1_factor2_GRIDY)[1])
df_R2_group1_factor3_GRIDY$factor <- rep(1,dim(df_R2_group1_factor3_GRIDY)[1])
df_R2_group1_factor1_GRIDY$structure <- rep("Joint",dim(df_R2_group1_factor1_GRIDY)[1])
df_R2_group1_factor2_GRIDY$structure <- rep("Joint",dim(df_R2_group1_factor2_GRIDY)[1])
df_R2_group1_factor3_GRIDY$structure <- rep("Group Individual 1",dim(df_R2_group1_factor3_GRIDY)[1])
df_R2_group2_factor1_GRIDY$factor <- rep(1,dim(df_R2_group2_factor1_GRIDY)[1])
df_R2_group2_factor2_GRIDY$factor <- rep(2,dim(df_R2_group2_factor2_GRIDY)[1])
df_R2_group2_factor3_GRIDY$factor <- rep(1,dim(df_R2_group2_factor3_GRIDY)[1])
df_R2_group2_factor1_GRIDY$structure <- rep("Joint",dim(df_R2_group2_factor1_GRIDY)[1])
df_R2_group2_factor2_GRIDY$structure <- rep("Joint",dim(df_R2_group2_factor2_GRIDY)[1])
df_R2_group2_factor3_GRIDY$structure <- rep("Group Individual 2",dim(df_R2_group2_factor3_GRIDY)[1])
df_R2_group1_factor1_SCA_P$factor <- rep(1,dim(df_R2_group1_factor1_SCA_P)[1])
df_R2_group1_factor2_SCA_P$factor <- rep(2,dim(df_R2_group1_factor2_SCA_P)[1])
df_R2_group1_factor3_SCA_P$factor <- rep(1,dim(df_R2_group1_factor3_SCA_P)[1])
df_R2_group1_factor1_SCA_P$structure <- rep("Joint",dim(df_R2_group1_factor1_SCA_P)[1])
df_R2_group1_factor2_SCA_P$structure <- rep("Joint",dim(df_R2_group1_factor2_SCA_P)[1])
df_R2_group1_factor3_SCA_P$structure <- rep("Group Individual 1",dim(df_R2_group1_factor3_SCA_P)[1])
df_R2_group2_factor1_SCA_P$factor <- rep(1,dim(df_R2_group2_factor1_SCA_P)[1])
df_R2_group2_factor2_SCA_P$factor <- rep(2,dim(df_R2_group2_factor2_SCA_P)[1])
df_R2_group2_factor3_SCA_P$factor <- rep(1,dim(df_R2_group2_factor3_SCA_P)[1])
df_R2_group2_factor1_SCA_P$structure <- rep("Joint",dim(df_R2_group2_factor1_SCA_P)[1])
df_R2_group2_factor2_SCA_P$structure <- rep("Joint",dim(df_R2_group2_factor2_SCA_P)[1])
df_R2_group2_factor3_SCA_P$structure <- rep("Group Individual 2",dim(df_R2_group2_factor3_SCA_P)[1])
df_R2_group1_factor1_GICA$factor <- rep(1,dim(df_R2_group1_factor1_GICA)[1])
df_R2_group1_factor2_GICA$factor <- rep(2,dim(df_R2_group1_factor2_GICA)[1])
df_R2_group1_factor3_GICA$factor <- rep(1,dim(df_R2_group1_factor3_GICA)[1])
df_R2_group1_factor1_GICA$structure <- rep("Joint",dim(df_R2_group1_factor1_GICA)[1])
df_R2_group1_factor2_GICA$structure <- rep("Joint",dim(df_R2_group1_factor2_GICA)[1])
df_R2_group1_factor3_GICA$structure <- rep("Group Individual 1",dim(df_R2_group1_factor3_GICA)[1])
df_R2_group2_factor1_GICA$factor <- rep(1,dim(df_R2_group2_factor1_GICA)[1])
df_R2_group2_factor2_GICA$factor <- rep(2,dim(df_R2_group2_factor2_GICA)[1])
df_R2_group2_factor3_GICA$factor <- rep(1,dim(df_R2_group2_factor3_GICA)[1])
df_R2_group2_factor1_GICA$structure <- rep("Joint",dim(df_R2_group2_factor1_GICA)[1])
df_R2_group2_factor2_GICA$structure <- rep("Joint",dim(df_R2_group2_factor2_GICA)[1])
df_R2_group2_factor3_GICA$structure <- rep("Group Individual 2",dim(df_R2_group2_factor3_GICA)[1])


df_R2_group1_GRIDY <- rbind(df_R2_group1_factor1_GRIDY,df_R2_group1_factor2_GRIDY,df_R2_group1_factor3_GRIDY)
df_R2_group2_GRIDY <- rbind(df_R2_group2_factor1_GRIDY,df_R2_group2_factor2_GRIDY,df_R2_group2_factor3_GRIDY)
df_R2_group1_GRIDY$group <- rep(1,dim(df_R2_group1_GRIDY)[1])
df_R2_group2_GRIDY$group <- rep(2,dim(df_R2_group2_GRIDY)[1])
df_R2_group1_SCA_P <- rbind(df_R2_group1_factor1_SCA_P,df_R2_group1_factor2_SCA_P,df_R2_group1_factor3_SCA_P)
df_R2_group2_SCA_P <- rbind(df_R2_group2_factor1_SCA_P,df_R2_group2_factor2_SCA_P,df_R2_group2_factor3_SCA_P)
df_R2_group1_SCA_P$group <- rep(1,dim(df_R2_group1_SCA_P)[1])
df_R2_group2_SCA_P$group <- rep(2,dim(df_R2_group2_SCA_P)[1])
df_R2_group1_GICA <- rbind(df_R2_group1_factor1_GICA,df_R2_group1_factor2_GICA,df_R2_group1_factor3_GICA)
df_R2_group2_GICA <- rbind(df_R2_group2_factor1_GICA,df_R2_group2_factor2_GICA,df_R2_group2_factor3_GICA)
df_R2_group1_GICA$group <- rep(1,dim(df_R2_group1_GICA)[1])
df_R2_group2_GICA$group <- rep(2,dim(df_R2_group2_GICA)[1])


df_R2_GRIDY <- rbind(df_R2_group1_GRIDY,df_R2_group2_GRIDY)
df_R2_GRIDY$factor <- factor(df_R2_GRIDY$factor,levels=c(1,2),labels=c("Factor 1","Factor 2"))
df_R2_GRIDY$group <- factor(df_R2_GRIDY$group,levels=c(1,2),labels=c("Group 1","Group 2"))
df_R2_GRIDY$structure <- factor(df_R2_GRIDY$structure,levels=c("Joint","Group Individual 1","Group Individual 2"))
df_R2_GRIDY$method <- "GRIDY"
df_R2_SCA_P <- rbind(df_R2_group1_SCA_P,df_R2_group2_SCA_P)
df_R2_SCA_P$factor <- factor(df_R2_SCA_P$factor,levels=c(1,2),labels=c("Factor 1","Factor 2"))
df_R2_SCA_P$group <- factor(df_R2_SCA_P$group,levels=c(1,2),labels=c("Group 1","Group 2"))
df_R2_SCA_P$structure <- factor(df_R2_SCA_P$structure,levels=c("Joint","Group Individual 1","Group Individual 2"))
df_R2_SCA_P$method <- "SCA-P"
df_R2_GICA <- rbind(df_R2_group1_GICA,df_R2_group2_GICA)
df_R2_GICA$factor <- factor(df_R2_GICA$factor,levels=c(1,2),labels=c("Factor 1","Factor 2"))
df_R2_GICA$group <- factor(df_R2_GICA$group,levels=c(1,2),labels=c("Group 1","Group 2"))
df_R2_GICA$structure <- factor(df_R2_GICA$structure,levels=c("Joint","Group Individual 1","Group Individual 2"))
df_R2_GICA$method <- "GICA"

df_R2 <- rbind(df_R2_GRIDY,df_R2_SCA_P,df_R2_GICA)
df_R2$method <- factor(df_R2$method,levels=c("GRIDY","SCA-P","GICA"))

factor.labs <- c("Factor 1", "Factor 2")
names(factor.labs) <- c(1,2)
group.labs <- c("Group 1", "Group 2")
names(group.labs) <- c(1,2)
structure.labs <- c("Joint","Group Individual 1","Group Individual 2")
names(structure.labs) <- c("Joint","Group Individual 1","Group Individual 2")
method.labs <- c("GRIDY","SCA-P","GICA")
names(method.labs) <- c("GRIDY","SCA-P","GICA")


# Plotting:
par(mar=c(0,0,0,0))
ggplot(df_R2,aes(x=variable,y=value,color=color,fill=color)) +
  geom_boxplot(size=0.1,alpha=0.5) +
  scale_color_manual(name="Network",
                     labels=c("dfm","fp","co","sm","cb","oc"),
                     breaks=c("red","blue","green","purple","orange","brown"),
                     values=c("red","blue","green","purple","orange","brown")) +
  scale_fill_manual(name="Network",
                     labels=c("dfm","fp","co","sm","cb","oc"),
                     breaks=c("red","blue","green","purple","orange","brown"),
                     values=c("red","blue","green","purple","orange","brown")) +
  scale_y_continuous(name=TeX("$\\R^2$")) +
  theme(legend.position="bottom",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)) + 
  facet_grid(structure+factor~method+group,
             labeller=labeller(method=method.labs,structure=structure.labs,
                               factor=factor.labs,group=group.labs))



# VAR
GRIDY_tr_joint <- GRIDY_tr_group  <-list()
GRIDY_error_joint <- GRIDY_error_group <- list()
for (kk in 1:(length(K1_idx)+length(K2_idx))){
  
  T_k <- dim(X_scale_ext[[kk]])[2]
  
  GRIDY_tr_joint[[kk]] <- SCA_PF2_Joint$B %*% GRIDY_YK$Psi_joint_hat[,,kk] %*% solve(t(SCA_PF2_Joint$B) %*% SCA_PF2_Joint$B) %*% t(SCA_PF2_Joint$B)
  
  if (kk <= length(K1_idx)){
    GRIDY_tr_group[[kk]] <- SCA_PF2_Group1$B %*% GRIDY_YK$Psi_indiv1_hat[,,kk] %*% solve(t(SCA_PF2_Group1$B) %*% SCA_PF2_Group1$B) %*% t(SCA_PF2_Group1$B)
    
  }else{
    GRIDY_tr_group[[kk]] <- SCA_PF2_Group2$B %*% GRIDY_YK$Psi_indiv2_hat[,,(kk-length(K1_idx))] %*% solve(t(SCA_PF2_Group2$B) %*% SCA_PF2_Group2$B) %*% t(SCA_PF2_Group2$B)
    
  }
  
  GRIDY_error_joint[[kk]] <- GRIDY_tr_joint[[kk]] %*% (diag(1,160) + (t(GRIDY_refitted$noise_refitted[[kk]]) %*% GRIDY_refitted$noise_refitted[[kk]]) / T_k) + SCA_PF2_Joint$B %*% GRIDY_YK$Eta_joint_hat[,,kk] %*% t(SCA_PF2_Joint$B)
  
  if (kk <= length(K1_idx)){
    GRIDY_error_group[[kk]] <- GRIDY_tr_group[[kk]] %*% (diag(1,160) + (t(GRIDY_refitted$noise_refitted[[kk]]) %*% GRIDY_refitted$noise_refitted[[kk]]) / T_k) + SCA_PF2_Group1$B %*% GRIDY_YK$Eta_indiv1_hat[,,kk] %*% t(SCA_PF2_Group1$B)
    
  }else{
    GRIDY_error_group[[kk]] <- GRIDY_tr_group[[kk]] %*% (diag(1,160) + (t(GRIDY_refitted$noise_refitted[[kk]]) %*% GRIDY_refitted$noise_refitted[[kk]]) / T_k) + SCA_PF2_Group2$B %*% GRIDY_YK$Eta_indiv2_hat[,,(kk-length(K1_idx))] %*% t(SCA_PF2_Group2$B)
    
  }
}

GRIDY_tr_group1 <- array(0,dim=c(160,160))
GRIDY_tr_group2 <- array(0,dim=c(160,160))
GRIDY_error_group1 <- array(0,dim=c(160,160))
GRIDY_error_group2 <- array(0,dim=c(160,160))

for (kk in 1:(length(K1_idx)+length(K2_idx))){
  
  if (kk <= length(K1_idx)){
    GRIDY_tr_group1 <- GRIDY_tr_group1 + GRIDY_tr_joint[[kk]]/length(K1_idx) + GRIDY_tr_group[[kk]]/length(K1_idx)
    GRIDY_error_group1 <- GRIDY_error_group1 + GRIDY_error_joint[[kk]]/length(K1_idx) + GRIDY_error_group[[kk]]/length(K1_idx)
    
  }else{
    GRIDY_tr_group2 <- GRIDY_tr_group2 + GRIDY_tr_joint[[kk]]/length(K2_idx) + GRIDY_tr_group[[kk]]/length(K2_idx)
    GRIDY_error_group2 <- GRIDY_error_group2 + GRIDY_error_joint[[kk]]/length(K2_idx) + GRIDY_error_group[[kk]]/length(K2_idx)
    
  }
}


df_GRIDY_tr_group1 <- data.frame(value=as.vector(GRIDY_tr_group1[color_table$ROI,color_table$ROI]),
                                 row = rep(color_table$ROI,each=1,160),
                                 column = rep(color_table$ROI,each=160))
df_GRIDY_tr_group2 <- data.frame(value=as.vector(GRIDY_tr_group2[color_table$ROI,color_table$ROI]),
                                 row = rep(color_table$ROI,each=1,160),
                                 column = rep(color_table$ROI,each=160))
df_GRIDY_error_group1 <- data.frame(value=as.vector(GRIDY_error_group1[color_table$ROI,color_table$ROI]),
                                    row = rep(color_table$ROI,each=1,160),
                                    column = rep(color_table$ROI,each=160))
df_GRIDY_error_group2 <- data.frame(value=as.vector(GRIDY_error_group2[color_table$ROI,color_table$ROI]),
                                    row = rep(color_table$ROI,each=1,160),
                                    column = rep(color_table$ROI,each=160))

df_GRIDY_tr_group1$row <- factor(df_GRIDY_tr_group1$row,levels=unique(df_GRIDY_tr_group1$row))
df_GRIDY_tr_group2$row <- factor(df_GRIDY_tr_group2$row,levels=unique(df_GRIDY_tr_group2$row))
df_GRIDY_error_group1$row <- factor(df_GRIDY_error_group1$row,levels=unique(df_GRIDY_error_group1$row))
df_GRIDY_error_group2$row <- factor(df_GRIDY_error_group2$row,levels=unique(df_GRIDY_error_group2$row))
df_GRIDY_tr_group1$column <- factor(df_GRIDY_tr_group1$column,levels=unique(df_GRIDY_tr_group1$column))
df_GRIDY_tr_group2$column <- factor(df_GRIDY_tr_group2$column,levels=unique(df_GRIDY_tr_group2$column))
df_GRIDY_error_group1$column <- factor(df_GRIDY_error_group1$column,levels=unique(df_GRIDY_error_group1$column))
df_GRIDY_error_group2$column <- factor(df_GRIDY_error_group2$column,levels=unique(df_GRIDY_error_group2$column))


color_idx <- (color_table$color)[which(c(1:160)%%4==1)]
GRIDY_group1_tr <- ggplot(df_GRIDY_tr_group1,aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile()+
  scale_fill_gradient2(name="Value",low=colorRampPalette(brewer.pal(3,"Set1"))(3)[2],
                       high=colorRampPalette(brewer.pal(3,"Set1"))(3)[1]) +
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_GRIDY_tr_group1$row)), 
                   breaks = unique(as.factor(df_GRIDY_tr_group1$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_GRIDY_tr_group1$column)), 
                   breaks = unique(as.factor(df_GRIDY_tr_group1$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,color=color_idx),
        axis.text.y=element_text(size=8,color=color_idx))

GRIDY_group2_tr <- ggplot(df_GRIDY_tr_group2,aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile()+
  scale_fill_gradient2(name="Value",low=colorRampPalette(brewer.pal(3,"Set1"))(3)[2],
                       high=colorRampPalette(brewer.pal(3,"Set1"))(3)[1]) +
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_GRIDY_tr_group2$row)), 
                   breaks = unique(as.factor(df_GRIDY_tr_group2$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_GRIDY_tr_group2$column)), 
                   breaks = unique(as.factor(df_GRIDY_tr_group2$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,color=color_idx),
        axis.text.y=element_text(size=8,color=color_idx))


GRIDY_group1_error <- ggplot(df_GRIDY_error_group1,aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile()+
  scale_fill_gradient2(name="Value",low=colorRampPalette(brewer.pal(3,"Set1"))(3)[2],
                       high=colorRampPalette(brewer.pal(3,"Set1"))(3)[1]) +
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_GRIDY_error_group1$row)), 
                   breaks = unique(as.factor(df_GRIDY_error_group1$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_GRIDY_error_group1$column)), 
                   breaks = unique(as.factor(df_GRIDY_error_group1$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,color=color_idx),
        axis.text.y=element_text(size=8,color=color_idx))

GRIDY_group2_error <- ggplot(df_GRIDY_error_group2,aes(x=as.factor(row),y=as.factor(column),fill=value))+
  geom_tile()+
  scale_fill_gradient2(name="Value",low=colorRampPalette(brewer.pal(3,"Set1"))(3)[2],
                       high=colorRampPalette(brewer.pal(3,"Set1"))(3)[1]) +
  scale_y_discrete(name="ROI labels",
                   limits = levels(as.factor(df_GRIDY_error_group2$row)), 
                   breaks = unique(as.factor(df_GRIDY_error_group2$row))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+
  scale_x_discrete(name="ROI labels",
                   limits = levels(as.factor(df_GRIDY_error_group2$column)), 
                   breaks = unique(as.factor(df_GRIDY_error_group2$column))[which(c(1:160)%%4==1)],
                   labels = voxel_name_table$ROI_label[color_table$ROI][which(c(1:160)%%4==1)])+  
  theme(legend.text=element_text(size=7),legend.position="bottom",
        axis.text.x=element_text(size=8,angle=45,vjust=0.5,color=color_idx),
        axis.text.y=element_text(size=8,color=color_idx))


# Plotting:
par(mar=c(0,0,0,0))
ggarrange(GRIDY_group1_tr,GRIDY_group2_tr,
          GRIDY_group1_error,GRIDY_group2_error,
          nrow=2,ncol=2,common.legend=TRUE,legend="bottom")


