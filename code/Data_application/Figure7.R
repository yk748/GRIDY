#-----------------------------------------------------------------------------#
# Code for producing Figure 7 in Data application section
#-----------------------------------------------------------------------------#


# -----------------------------------------------------------------------------#
# Preparation
# -----------------------------------------------------------------------------#
# Packages required
library(ggplot2)
library(latex2exp)
library(ggpubr)
library(RColorBrewer)
library(reshape2)

# Load data from the previous step
load("./result/result_final_X.RData")
load("./result/result_final_GRIDY.RData")
load("./result/result_final_SCA_P.RData")
load("./result/result_final_GICA.RData")

# -----------------------------------------------------------------------------#
# Creating Figure 7
# -----------------------------------------------------------------------------#
# Making color labels for the ROIs:
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


# VAR construction
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


# Plotting: This produces Figure 7
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


# -----------------------------------------------------------------------------#
# Saving Figure 7
# -----------------------------------------------------------------------------#
plot_tmp <- ggarrange(GRIDY_group1_tr,GRIDY_group2_tr,
                      GRIDY_group1_error,GRIDY_group2_error,
                      nrow=2,ncol=2,common.legend=TRUE,legend="bottom")

pdf(file = "./Figure7.pdf", width = 11.5, height = 10.5)
par(mar=c(0,0,0,0))
plot_tmp
dev.off()

