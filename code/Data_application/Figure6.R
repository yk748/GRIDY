#-----------------------------------------------------------------------------#
# Code for producing Figure 6 in Data application section
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
load("./result/result_final.RData")

# -----------------------------------------------------------------------------#
# Creating Figure 6
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


# Computing R2 statistics
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


# Plotting: This produces Figure 6
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


# -----------------------------------------------------------------------------#
# Saving Figure 6
# -----------------------------------------------------------------------------#
pdf(file = "./Figure6.pdf", width = 11.5, height = 10.5)
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
dev.off()
