#-----------------------------------------------------------------------------#
# Code for producing Figure S3 in Supplemental material
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
# Creating Figure S3
# -----------------------------------------------------------------------------#
# Computing covariance matrices
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


# Plotting: This produces Figure S3
# Covariance matrices of factors
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


# Transition matrices of factors
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


# -----------------------------------------------------------------------------#
# Saving Figure S3
# -----------------------------------------------------------------------------#
plot_tmp <- ggarrange(
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

pdf(file = "./FigureS3.pdf", width = 11.5, height = 10.5)
par(mar=c(0,0,0,0))
plot_tmp
dev.off()

