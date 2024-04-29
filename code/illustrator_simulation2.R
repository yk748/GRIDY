rm(list=ls())

# Change the directory where the data is stored:
setwd("D:/GRIDY/data_sim")

library(ggplot2)
library(patchwork)
library(reshape2)
library(latex2exp)
library(ggpubr)
library(RColorBrewer)
library(gridExtra)

#-----------------------------------------------------------------------------#
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
#-----------------------------------------------------------------------------#
# 2. Strengths of signals
#-----------------------------------------------------------------------------#
load("./sim2_type1_snr1.RData"); 
list_type1_snr1 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./sim2_type1_snr2.RData"); 
list_type1_snr2 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./sim2_type1_snr3.RData"); 
list_type1_snr3 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./sim2_type1_snr4.RData"); 
list_type1_snr4 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./sim2_type1_snr5.RData"); 
list_type1_snr5 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./sim2_type1_snr6.RData"); 
list_type1_snr6 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)

load("./sim2_type2_snr1.RData"); 
list_type2_snr1 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./sim2_type2_snr2.RData"); 
list_type2_snr2 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./sim2_type2_snr3.RData"); 
list_type2_snr3 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./sim2_type2_snr4.RData"); 
list_type2_snr4 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./sim2_type2_snr5.RData"); 
list_type2_snr5 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./sim2_type2_snr6.RData"); 
list_type2_snr6 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)


rm(list_GRIDY,list_DSCA,list_SCA_P,list_GICA)


#-----------------------------------------------------------------------------#
# type 1:
#-----------------------------------------------------------------------------#
# R2:
R2_J_snr1_type1 <- data.frame(value=c(unlist(list_type1_snr1$GRIDY$R2$Joint),
                                      unlist(list_type1_snr1$SCA_P$R2$Joint),
                                      unlist(list_type1_snr1$GICA$R2$Joint),
                                      unlist(list_type1_snr1$DSCA$R2$Joint),
                                      unlist(list_type1_snr1$DGICA$R2$Joint),
                                      unlist(list_type1_snr1$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_type1_snr1$Unfitted_GICA$R2$Joint)))
R2_J_snr1_type1$SNR <- 0.25
R2_J_snr2_type1 <- data.frame(value=c(unlist(list_type1_snr2$GRIDY$R2$Joint),
                                      unlist(list_type1_snr2$SCA_P$R2$Joint),
                                      unlist(list_type1_snr2$GICA$R2$Joint),
                                      unlist(list_type1_snr2$DSCA$R2$Joint),
                                      unlist(list_type1_snr2$DGICA$R2$Joint),
                                      unlist(list_type1_snr2$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_type1_snr2$Unfitted_GICA$R2$Joint)))
R2_J_snr2_type1$SNR <- 0.5
R2_J_snr3_type1 <- data.frame(value=c(unlist(list_type1_snr3$GRIDY$R2$Joint),
                                      unlist(list_type1_snr3$SCA_P$R2$Joint),
                                      unlist(list_type1_snr3$GICA$R2$Joint),
                                      unlist(list_type1_snr3$DSCA$R2$Joint),
                                      unlist(list_type1_snr3$DGICA$R2$Joint),
                                      unlist(list_type1_snr3$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_type1_snr3$Unfitted_GICA$R2$Joint)))
R2_J_snr3_type1$SNR <- 0.75
R2_J_snr4_type1 <- data.frame(value=c(unlist(list_type1_snr4$GRIDY$R2$Joint),
                                      unlist(list_type1_snr4$SCA_P$R2$Joint),
                                      unlist(list_type1_snr4$GICA$R2$Joint),
                                      unlist(list_type1_snr4$DSCA$R2$Joint),
                                      unlist(list_type1_snr4$DGICA$R2$Joint),
                                      unlist(list_type1_snr4$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_type1_snr4$Unfitted_GICA$R2$Joint)))
R2_J_snr4_type1$SNR <- 1.25
R2_J_snr5_type1 <- data.frame(value=c(unlist(list_type1_snr5$GRIDY$R2$Joint),
                                      unlist(list_type1_snr5$SCA_P$R2$Joint),
                                      unlist(list_type1_snr5$GICA$R2$Joint),
                                      unlist(list_type1_snr5$DSCA$R2$Joint),
                                      unlist(list_type1_snr5$DGICA$R2$Joint),
                                      unlist(list_type1_snr5$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_type1_snr5$Unfitted_GICA$R2$Joint)))
R2_J_snr5_type1$SNR <- 2
R2_J_snr6_type1 <- data.frame(value=c(unlist(list_type1_snr6$GRIDY$R2$Joint),
                                      unlist(list_type1_snr6$SCA_P$R2$Joint),
                                      unlist(list_type1_snr6$GICA$R2$Joint),
                                      unlist(list_type1_snr6$DSCA$R2$Joint),
                                      unlist(list_type1_snr6$DGICA$R2$Joint),
                                      unlist(list_type1_snr6$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_type1_snr6$Unfitted_GICA$R2$Joint)))
R2_J_snr6_type1$SNR <- 4
R2_J_type1 <- rbind(R2_J_snr1_type1,R2_J_snr2_type1,R2_J_snr3_type1,
                    R2_J_snr4_type1,R2_J_snr5_type1,R2_J_snr6_type1)
R2_J_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                               "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

R2_G1_snr1_type1 <- data.frame(value=c(unlist(list_type1_snr1$GRIDY$R2$Group1),
                                       unlist(list_type1_snr1$SCA_P$R2$Group1),
                                       unlist(list_type1_snr1$GICA$R2$Group1),
                                       unlist(list_type1_snr1$DSCA$R2$Group1),
                                       unlist(list_type1_snr1$DGICA$R2$Group1),
                                       unlist(list_type1_snr1$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_type1_snr1$Unfitted_GICA$R2$Group1)))
R2_G1_snr1_type1$SNR <- 0.25
R2_G1_snr2_type1 <- data.frame(value=c(unlist(list_type1_snr2$GRIDY$R2$Group1),
                                       unlist(list_type1_snr2$SCA_P$R2$Group1),
                                       unlist(list_type1_snr2$GICA$R2$Group1),
                                       unlist(list_type1_snr2$DSCA$R2$Group1),
                                       unlist(list_type1_snr2$DGICA$R2$Group1),
                                       unlist(list_type1_snr2$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_type1_snr2$Unfitted_GICA$R2$Group1)))
R2_G1_snr2_type1$SNR <- 0.5
R2_G1_snr3_type1 <- data.frame(value=c(unlist(list_type1_snr3$GRIDY$R2$Group1),
                                       unlist(list_type1_snr3$SCA_P$R2$Group1),
                                       unlist(list_type1_snr3$GICA$R2$Group1),
                                       unlist(list_type1_snr3$DSCA$R2$Group1),
                                       unlist(list_type1_snr3$DGICA$R2$Group1),
                                       unlist(list_type1_snr3$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_type1_snr3$Unfitted_GICA$R2$Group1)))
R2_G1_snr3_type1$SNR <- 0.75
R2_G1_snr4_type1 <- data.frame(value=c(unlist(list_type1_snr4$GRIDY$R2$Group1),
                                       unlist(list_type1_snr4$SCA_P$R2$Group1),
                                       unlist(list_type1_snr4$GICA$R2$Group1),
                                       unlist(list_type1_snr4$DSCA$R2$Group1),
                                       unlist(list_type1_snr4$DGICA$R2$Group1),
                                       unlist(list_type1_snr4$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_type1_snr4$Unfitted_GICA$R2$Group1)))
R2_G1_snr4_type1$SNR <- 1.25
R2_G1_snr5_type1 <- data.frame(value=c(unlist(list_type1_snr5$GRIDY$R2$Group1),
                                       unlist(list_type1_snr5$SCA_P$R2$Group1),
                                       unlist(list_type1_snr5$GICA$R2$Group1),
                                       unlist(list_type1_snr5$DSCA$R2$Group1),
                                       unlist(list_type1_snr5$DGICA$R2$Group1),
                                       unlist(list_type1_snr5$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_type1_snr5$Unfitted_GICA$R2$Group1)))
R2_G1_snr5_type1$SNR <- 2
R2_G1_snr6_type1 <- data.frame(value=c(unlist(list_type1_snr6$GRIDY$R2$Group1),
                                       unlist(list_type1_snr6$SCA_P$R2$Group1),
                                       unlist(list_type1_snr6$GICA$R2$Group1),
                                       unlist(list_type1_snr6$DSCA$R2$Group1),
                                       unlist(list_type1_snr6$DGICA$R2$Group1),
                                       unlist(list_type1_snr6$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_type1_snr6$Unfitted_GICA$R2$Group1)))
R2_G1_snr6_type1$SNR <- 4
R2_G1_type1 <- rbind(R2_G1_snr1_type1,R2_G1_snr2_type1,R2_G1_snr3_type1,
                     R2_G1_snr4_type1,R2_G1_snr5_type1,R2_G1_snr6_type1)
R2_G1_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


R2_G2_snr1_type1 <- data.frame(value=c(unlist(list_type1_snr1$GRIDY$R2$Group2),
                                       unlist(list_type1_snr1$SCA_P$R2$Group2),
                                       unlist(list_type1_snr1$GICA$R2$Group2),
                                       unlist(list_type1_snr1$DSCA$R2$Group2),
                                       unlist(list_type1_snr1$DGICA$R2$Group2),
                                       unlist(list_type1_snr1$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_type1_snr1$Unfitted_GICA$R2$Group2)))
R2_G2_snr1_type1$SNR <- 0.25
R2_G2_snr2_type1 <- data.frame(value=c(unlist(list_type1_snr2$GRIDY$R2$Group2),
                                       unlist(list_type1_snr2$SCA_P$R2$Group2),
                                       unlist(list_type1_snr2$GICA$R2$Group2),
                                       unlist(list_type1_snr2$DSCA$R2$Group2),
                                       unlist(list_type1_snr2$DGICA$R2$Group2),
                                       unlist(list_type1_snr2$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_type1_snr2$Unfitted_GICA$R2$Group2)))
R2_G2_snr2_type1$SNR <- 0.5
R2_G2_snr3_type1 <- data.frame(value=c(unlist(list_type1_snr3$GRIDY$R2$Group2),
                                       unlist(list_type1_snr3$SCA_P$R2$Group2),
                                       unlist(list_type1_snr3$GICA$R2$Group2),
                                       unlist(list_type1_snr3$DSCA$R2$Group2),
                                       unlist(list_type1_snr3$DGICA$R2$Group2),
                                       unlist(list_type1_snr3$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_type1_snr3$Unfitted_GICA$R2$Group2)))
R2_G2_snr3_type1$SNR <- 0.75
R2_G2_snr4_type1 <- data.frame(value=c(unlist(list_type1_snr4$GRIDY$R2$Group2),
                                       unlist(list_type1_snr4$SCA_P$R2$Group2),
                                       unlist(list_type1_snr4$GICA$R2$Group2),
                                       unlist(list_type1_snr4$DSCA$R2$Group2),
                                       unlist(list_type1_snr4$DGICA$R2$Group2),
                                       unlist(list_type1_snr4$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_type1_snr4$Unfitted_GICA$R2$Group2)))
R2_G2_snr4_type1$SNR <- 1.25
R2_G2_snr5_type1 <- data.frame(value=c(unlist(list_type1_snr5$GRIDY$R2$Group2),
                                       unlist(list_type1_snr5$SCA_P$R2$Group2),
                                       unlist(list_type1_snr5$GICA$R2$Group2),
                                       unlist(list_type1_snr5$DSCA$R2$Group2),
                                       unlist(list_type1_snr5$DGICA$R2$Group2),
                                       unlist(list_type1_snr5$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_type1_snr5$Unfitted_GICA$R2$Group2)))
R2_G2_snr5_type1$SNR <- 2
R2_G2_snr6_type1 <- data.frame(value=c(unlist(list_type1_snr6$GRIDY$R2$Group2),
                                       unlist(list_type1_snr6$SCA_P$R2$Group2),
                                       unlist(list_type1_snr6$GICA$R2$Group2),
                                       unlist(list_type1_snr6$DSCA$R2$Group2),
                                       unlist(list_type1_snr6$DGICA$R2$Group2),
                                       unlist(list_type1_snr6$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_type1_snr6$Unfitted_GICA$R2$Group2)))
R2_G2_snr6_type1$SNR <- 4
R2_G2_type1 <- rbind(R2_G2_snr1_type1,R2_G2_snr2_type1,R2_G2_snr3_type1,
                     R2_G2_snr4_type1,R2_G2_snr5_type1,R2_G2_snr6_type1)
R2_G2_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


SNR_list <- c("0.25","0.5","0.75","1.25","2","4")
R2_J_type1$SNR <- factor(R2_J_type1$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
R2_G1_type1$SNR <- factor(R2_G1_type1$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
R2_G2_type1$SNR <- factor(R2_G2_type1$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)

method_list <- c("GRIDY","SCA-P","GICA","DSCA","DGICA","Unfitted_SCA-PF2","Unfitted_GICA")
R2_J_type1$method <- factor(R2_J_type1$method,levels=method_list,labels=method_list)
R2_G1_type1$method <- factor(R2_G1_type1$method,levels=method_list,labels=method_list)
R2_G2_type1$method <- factor(R2_G2_type1$method,levels=method_list,labels=method_list)

pl_R2_J_type1 <- ggplot(R2_J_type1,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$R^2$ (Joint)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_R2_G1_type1 <- ggplot(R2_G1_type1,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$R^2$ (Group1)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_R2_G2_type1 <- ggplot(R2_G2_type1,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$R^2$ (Group2)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))


# RMSE:
RMSE_J_snr1_type1 <- data.frame(value=c(unlist(list_type1_snr1$GRIDY$RMSE$Joint),
                                        unlist(list_type1_snr1$SCA_P$RMSE$Joint),
                                        unlist(list_type1_snr1$GICA$RMSE$Joint),
                                        unlist(list_type1_snr1$DSCA$RMSE$Joint),
                                        unlist(list_type1_snr1$DGICA$RMSE$Joint),
                                        unlist(list_type1_snr1$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_type1_snr1$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr1_type1$SNR <- 0.25
RMSE_J_snr2_type1 <- data.frame(value=c(unlist(list_type1_snr2$GRIDY$RMSE$Joint),
                                        unlist(list_type1_snr2$SCA_P$RMSE$Joint),
                                        unlist(list_type1_snr2$GICA$RMSE$Joint),
                                        unlist(list_type1_snr2$DSCA$RMSE$Joint),
                                        unlist(list_type1_snr2$DGICA$RMSE$Joint),
                                        unlist(list_type1_snr2$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_type1_snr2$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr2_type1$SNR <- 0.5
RMSE_J_snr3_type1 <- data.frame(value=c(unlist(list_type1_snr3$GRIDY$RMSE$Joint),
                                        unlist(list_type1_snr3$SCA_P$RMSE$Joint),
                                        unlist(list_type1_snr3$GICA$RMSE$Joint),
                                        unlist(list_type1_snr3$DSCA$RMSE$Joint),
                                        unlist(list_type1_snr3$DGICA$RMSE$Joint),
                                        unlist(list_type1_snr3$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_type1_snr3$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr3_type1$SNR <- 0.75
RMSE_J_snr4_type1 <- data.frame(value=c(unlist(list_type1_snr4$GRIDY$RMSE$Joint),
                                        unlist(list_type1_snr4$SCA_P$RMSE$Joint),
                                        unlist(list_type1_snr4$GICA$RMSE$Joint),
                                        unlist(list_type1_snr4$DSCA$RMSE$Joint),
                                        unlist(list_type1_snr4$DGICA$RMSE$Joint),
                                        unlist(list_type1_snr4$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_type1_snr4$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr4_type1$SNR <- 1.25
RMSE_J_snr5_type1 <- data.frame(value=c(unlist(list_type1_snr5$GRIDY$RMSE$Joint),
                                        unlist(list_type1_snr5$SCA_P$RMSE$Joint),
                                        unlist(list_type1_snr5$GICA$RMSE$Joint),
                                        unlist(list_type1_snr5$DSCA$RMSE$Joint),
                                        unlist(list_type1_snr5$DGICA$RMSE$Joint),
                                        unlist(list_type1_snr5$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_type1_snr5$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr5_type1$SNR <- 2
RMSE_J_snr6_type1 <- data.frame(value=c(unlist(list_type1_snr6$GRIDY$RMSE$Joint),
                                        unlist(list_type1_snr6$SCA_P$RMSE$Joint),
                                        unlist(list_type1_snr6$GICA$RMSE$Joint),
                                        unlist(list_type1_snr6$DSCA$RMSE$Joint),
                                        unlist(list_type1_snr6$DGICA$RMSE$Joint),
                                        unlist(list_type1_snr6$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_type1_snr6$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr6_type1$SNR <- 4
RMSE_J_type1 <- rbind(RMSE_J_snr1_type1,RMSE_J_snr2_type1,RMSE_J_snr3_type1,
                      RMSE_J_snr4_type1,RMSE_J_snr5_type1,RMSE_J_snr6_type1)
RMSE_J_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                 "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

RMSE_G1_snr1_type1 <- data.frame(value=c(unlist(list_type1_snr1$GRIDY$RMSE$Group1),
                                         unlist(list_type1_snr1$SCA_P$RMSE$Group1),
                                         unlist(list_type1_snr1$GICA$RMSE$Group1),
                                         unlist(list_type1_snr1$DSCA$RMSE$Group1),
                                         unlist(list_type1_snr1$DGICA$RMSE$Group1),
                                         unlist(list_type1_snr1$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_type1_snr1$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr1_type1$SNR <- 0.25
RMSE_G1_snr2_type1 <- data.frame(value=c(unlist(list_type1_snr2$GRIDY$RMSE$Group1),
                                         unlist(list_type1_snr2$SCA_P$RMSE$Group1),
                                         unlist(list_type1_snr2$GICA$RMSE$Group1),
                                         unlist(list_type1_snr2$DSCA$RMSE$Group1),
                                         unlist(list_type1_snr2$DGICA$RMSE$Group1),
                                         unlist(list_type1_snr2$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_type1_snr2$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr2_type1$SNR <- 0.5
RMSE_G1_snr3_type1 <- data.frame(value=c(unlist(list_type1_snr3$GRIDY$RMSE$Group1),
                                         unlist(list_type1_snr3$SCA_P$RMSE$Group1),
                                         unlist(list_type1_snr3$GICA$RMSE$Group1),
                                         unlist(list_type1_snr3$DSCA$RMSE$Group1),
                                         unlist(list_type1_snr3$DGICA$RMSE$Group1),
                                         unlist(list_type1_snr3$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_type1_snr3$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr3_type1$SNR <- 0.75
RMSE_G1_snr4_type1 <- data.frame(value=c(unlist(list_type1_snr4$GRIDY$RMSE$Group1),
                                         unlist(list_type1_snr4$SCA_P$RMSE$Group1),
                                         unlist(list_type1_snr4$GICA$RMSE$Group1),
                                         unlist(list_type1_snr4$DSCA$RMSE$Group1),
                                         unlist(list_type1_snr4$DGICA$RMSE$Group1),
                                         unlist(list_type1_snr4$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_type1_snr4$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr4_type1$SNR <- 1.25
RMSE_G1_snr5_type1 <- data.frame(value=c(unlist(list_type1_snr5$GRIDY$RMSE$Group1),
                                         unlist(list_type1_snr5$SCA_P$RMSE$Group1),
                                         unlist(list_type1_snr5$GICA$RMSE$Group1),
                                         unlist(list_type1_snr5$DSCA$RMSE$Group1),
                                         unlist(list_type1_snr5$DGICA$RMSE$Group1),
                                         unlist(list_type1_snr5$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_type1_snr5$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr5_type1$SNR <- 2
RMSE_G1_snr6_type1 <- data.frame(value=c(unlist(list_type1_snr6$GRIDY$RMSE$Group1),
                                         unlist(list_type1_snr6$SCA_P$RMSE$Group1),
                                         unlist(list_type1_snr6$GICA$RMSE$Group1),
                                         unlist(list_type1_snr6$DSCA$RMSE$Group1),
                                         unlist(list_type1_snr6$DGICA$RMSE$Group1),
                                         unlist(list_type1_snr6$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_type1_snr6$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr6_type1$SNR <- 4
RMSE_G1_type1 <- rbind(RMSE_G1_snr1_type1,RMSE_G1_snr2_type1,RMSE_G1_snr3_type1,
                       RMSE_G1_snr4_type1,RMSE_G1_snr5_type1,RMSE_G1_snr6_type1)
RMSE_G1_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                  "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


RMSE_G2_snr1_type1 <- data.frame(value=c(unlist(list_type1_snr1$GRIDY$RMSE$Group2),
                                         unlist(list_type1_snr1$SCA_P$RMSE$Group2),
                                         unlist(list_type1_snr1$GICA$RMSE$Group2),
                                         unlist(list_type1_snr1$DSCA$RMSE$Group2),
                                         unlist(list_type1_snr1$DGICA$RMSE$Group2),
                                         unlist(list_type1_snr1$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_type1_snr1$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr1_type1$SNR <- 0.25
RMSE_G2_snr2_type1 <- data.frame(value=c(unlist(list_type1_snr2$GRIDY$RMSE$Group2),
                                         unlist(list_type1_snr2$SCA_P$RMSE$Group2),
                                         unlist(list_type1_snr2$GICA$RMSE$Group2),
                                         unlist(list_type1_snr2$DSCA$RMSE$Group2),
                                         unlist(list_type1_snr2$DGICA$RMSE$Group2),
                                         unlist(list_type1_snr2$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_type1_snr2$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr2_type1$SNR <- 0.5
RMSE_G2_snr3_type1 <- data.frame(value=c(unlist(list_type1_snr3$GRIDY$RMSE$Group2),
                                         unlist(list_type1_snr3$SCA_P$RMSE$Group2),
                                         unlist(list_type1_snr3$GICA$RMSE$Group2),
                                         unlist(list_type1_snr3$DSCA$RMSE$Group2),
                                         unlist(list_type1_snr3$DGICA$RMSE$Group2),
                                         unlist(list_type1_snr3$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_type1_snr3$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr3_type1$SNR <- 0.75
RMSE_G2_snr4_type1 <- data.frame(value=c(unlist(list_type1_snr4$GRIDY$RMSE$Group2),
                                         unlist(list_type1_snr4$SCA_P$RMSE$Group2),
                                         unlist(list_type1_snr4$GICA$RMSE$Group2),
                                         unlist(list_type1_snr4$DSCA$RMSE$Group2),
                                         unlist(list_type1_snr4$DGICA$RMSE$Group2),
                                         unlist(list_type1_snr4$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_type1_snr4$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr4_type1$SNR <- 1.25
RMSE_G2_snr5_type1 <- data.frame(value=c(unlist(list_type1_snr5$GRIDY$RMSE$Group2),
                                         unlist(list_type1_snr5$SCA_P$RMSE$Group2),
                                         unlist(list_type1_snr5$GICA$RMSE$Group2),
                                         unlist(list_type1_snr5$DSCA$RMSE$Group2),
                                         unlist(list_type1_snr5$DGICA$RMSE$Group2),
                                         unlist(list_type1_snr5$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_type1_snr5$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr5_type1$SNR <- 2
RMSE_G2_snr6_type1 <- data.frame(value=c(unlist(list_type1_snr6$GRIDY$RMSE$Group2),
                                         unlist(list_type1_snr6$SCA_P$RMSE$Group2),
                                         unlist(list_type1_snr6$GICA$RMSE$Group2),
                                         unlist(list_type1_snr6$DSCA$RMSE$Group2),
                                         unlist(list_type1_snr6$DGICA$RMSE$Group2),
                                         unlist(list_type1_snr6$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_type1_snr6$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr6_type1$SNR <- 4
RMSE_G2_type1 <- rbind(RMSE_G2_snr1_type1,RMSE_G2_snr2_type1,RMSE_G2_snr3_type1,
                       RMSE_G2_snr4_type1,RMSE_G2_snr5_type1,RMSE_G2_snr6_type1)
RMSE_G2_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                  "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

SNR_list <- c("0.25","0.5","0.75","1.25","2","4")
RMSE_J_type1$SNR <- factor(RMSE_J_type1$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
RMSE_G1_type1$SNR <- factor(RMSE_G1_type1$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
RMSE_G2_type1$SNR <- factor(RMSE_G2_type1$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)

method_list <- c("GRIDY","SCA-P","GICA","DSCA","DGICA","Unfitted_SCA-PF2","Unfitted_GICA")
RMSE_J_type1$method <- factor(RMSE_J_type1$method,levels=method_list,labels=method_list)
RMSE_G1_type1$method <- factor(RMSE_G1_type1$method,levels=method_list,labels=method_list)
RMSE_G2_type1$method <- factor(RMSE_G2_type1$method,levels=method_list,labels=method_list)

pl_RMSE_J_type1 <- ggplot(RMSE_J_type1,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name="RMSE (Joint)") +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_RMSE_G1_type1 <- ggplot(RMSE_G1_type1,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name="RMSE (Group1)") +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_RMSE_G2_type1 <- ggplot(RMSE_G2_type1,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name="RMSE (Group2)") +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))


# CC_B:
CC_B_J_snr1_type1 <- data.frame(value=c(unlist(list_type1_snr1$GRIDY$CC_B$Joint),
                                        unlist(list_type1_snr1$SCA_P$CC_B$Joint),
                                        unlist(list_type1_snr1$GICA$CC_B$Joint),
                                        unlist(list_type1_snr1$DSCA$CC_B$Joint),
                                        unlist(list_type1_snr1$DGICA$CC_B$Joint),
                                        unlist(list_type1_snr1$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_type1_snr1$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr1_type1$SNR <- 0.25
CC_B_J_snr2_type1 <- data.frame(value=c(unlist(list_type1_snr2$GRIDY$CC_B$Joint),
                                        unlist(list_type1_snr2$SCA_P$CC_B$Joint),
                                        unlist(list_type1_snr2$GICA$CC_B$Joint),
                                        unlist(list_type1_snr2$DSCA$CC_B$Joint),
                                        unlist(list_type1_snr2$DGICA$CC_B$Joint),
                                        unlist(list_type1_snr2$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_type1_snr2$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr2_type1$SNR <- 0.5
CC_B_J_snr3_type1 <- data.frame(value=c(unlist(list_type1_snr3$GRIDY$CC_B$Joint),
                                        unlist(list_type1_snr3$SCA_P$CC_B$Joint),
                                        unlist(list_type1_snr3$GICA$CC_B$Joint),
                                        unlist(list_type1_snr3$DSCA$CC_B$Joint),
                                        unlist(list_type1_snr3$DGICA$CC_B$Joint),
                                        unlist(list_type1_snr3$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_type1_snr3$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr3_type1$SNR <- 0.75
CC_B_J_snr4_type1 <- data.frame(value=c(unlist(list_type1_snr4$GRIDY$CC_B$Joint),
                                        unlist(list_type1_snr4$SCA_P$CC_B$Joint),
                                        unlist(list_type1_snr4$GICA$CC_B$Joint),
                                        unlist(list_type1_snr4$DSCA$CC_B$Joint),
                                        unlist(list_type1_snr4$DGICA$CC_B$Joint),
                                        unlist(list_type1_snr4$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_type1_snr4$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr4_type1$SNR <- 1.25
CC_B_J_snr5_type1 <- data.frame(value=c(unlist(list_type1_snr5$GRIDY$CC_B$Joint),
                                        unlist(list_type1_snr5$SCA_P$CC_B$Joint),
                                        unlist(list_type1_snr5$GICA$CC_B$Joint),
                                        unlist(list_type1_snr5$DSCA$CC_B$Joint),
                                        unlist(list_type1_snr5$DGICA$CC_B$Joint),
                                        unlist(list_type1_snr5$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_type1_snr5$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr5_type1$SNR <- 2
CC_B_J_snr6_type1 <- data.frame(value=c(unlist(list_type1_snr6$GRIDY$CC_B$Joint),
                                        unlist(list_type1_snr6$SCA_P$CC_B$Joint),
                                        unlist(list_type1_snr6$GICA$CC_B$Joint),
                                        unlist(list_type1_snr6$DSCA$CC_B$Joint),
                                        unlist(list_type1_snr6$DGICA$CC_B$Joint),
                                        unlist(list_type1_snr6$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_type1_snr6$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr6_type1$SNR <- 4
CC_B_J_type1 <- rbind(CC_B_J_snr1_type1,CC_B_J_snr2_type1,CC_B_J_snr3_type1,
                      CC_B_J_snr4_type1,CC_B_J_snr5_type1,CC_B_J_snr6_type1)
CC_B_J_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                 "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

CC_B_G1_snr1_type1 <- data.frame(value=c(unlist(list_type1_snr1$GRIDY$CC_B$Group1),
                                         unlist(list_type1_snr1$SCA_P$CC_B$Group1),
                                         unlist(list_type1_snr1$GICA$CC_B$Group1),
                                         unlist(list_type1_snr1$DSCA$CC_B$Group1),
                                         unlist(list_type1_snr1$DGICA$CC_B$Group1),
                                         unlist(list_type1_snr1$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_type1_snr1$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr1_type1$SNR <- 0.25
CC_B_G1_snr2_type1 <- data.frame(value=c(unlist(list_type1_snr2$GRIDY$CC_B$Group1),
                                         unlist(list_type1_snr2$SCA_P$CC_B$Group1),
                                         unlist(list_type1_snr2$GICA$CC_B$Group1),
                                         unlist(list_type1_snr2$DSCA$CC_B$Group1),
                                         unlist(list_type1_snr2$DGICA$CC_B$Group1),
                                         unlist(list_type1_snr2$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_type1_snr2$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr2_type1$SNR <- 0.5
CC_B_G1_snr3_type1 <- data.frame(value=c(unlist(list_type1_snr3$GRIDY$CC_B$Group1),
                                         unlist(list_type1_snr3$SCA_P$CC_B$Group1),
                                         unlist(list_type1_snr3$GICA$CC_B$Group1),
                                         unlist(list_type1_snr3$DSCA$CC_B$Group1),
                                         unlist(list_type1_snr3$DGICA$CC_B$Group1),
                                         unlist(list_type1_snr3$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_type1_snr3$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr3_type1$SNR <- 0.75
CC_B_G1_snr4_type1 <- data.frame(value=c(unlist(list_type1_snr4$GRIDY$CC_B$Group1),
                                         unlist(list_type1_snr4$SCA_P$CC_B$Group1),
                                         unlist(list_type1_snr4$GICA$CC_B$Group1),
                                         unlist(list_type1_snr4$DSCA$CC_B$Group1),
                                         unlist(list_type1_snr4$DGICA$CC_B$Group1),
                                         unlist(list_type1_snr4$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_type1_snr4$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr4_type1$SNR <- 1.25
CC_B_G1_snr5_type1 <- data.frame(value=c(unlist(list_type1_snr5$GRIDY$CC_B$Group1),
                                         unlist(list_type1_snr5$SCA_P$CC_B$Group1),
                                         unlist(list_type1_snr5$GICA$CC_B$Group1),
                                         unlist(list_type1_snr5$DSCA$CC_B$Group1),
                                         unlist(list_type1_snr5$DGICA$CC_B$Group1),
                                         unlist(list_type1_snr5$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_type1_snr5$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr5_type1$SNR <- 2
CC_B_G1_snr6_type1 <- data.frame(value=c(unlist(list_type1_snr6$GRIDY$CC_B$Group1),
                                         unlist(list_type1_snr6$SCA_P$CC_B$Group1),
                                         unlist(list_type1_snr6$GICA$CC_B$Group1),
                                         unlist(list_type1_snr6$DSCA$CC_B$Group1),
                                         unlist(list_type1_snr6$DGICA$CC_B$Group1),
                                         unlist(list_type1_snr6$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_type1_snr6$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr6_type1$SNR <- 4
CC_B_G1_type1 <- rbind(CC_B_G1_snr1_type1,CC_B_G1_snr2_type1,CC_B_G1_snr3_type1,
                       CC_B_G1_snr4_type1,CC_B_G1_snr5_type1,CC_B_G1_snr6_type1)
CC_B_G1_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                  "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


CC_B_G2_snr1_type1 <- data.frame(value=c(unlist(list_type1_snr1$GRIDY$CC_B$Group2),
                                         unlist(list_type1_snr1$SCA_P$CC_B$Group2),
                                         unlist(list_type1_snr1$GICA$CC_B$Group2),
                                         unlist(list_type1_snr1$DSCA$CC_B$Group2),
                                         unlist(list_type1_snr1$DGICA$CC_B$Group2),
                                         unlist(list_type1_snr1$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_type1_snr1$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr1_type1$SNR <- 0.25
CC_B_G2_snr2_type1 <- data.frame(value=c(unlist(list_type1_snr2$GRIDY$CC_B$Group2),
                                         unlist(list_type1_snr2$SCA_P$CC_B$Group2),
                                         unlist(list_type1_snr2$GICA$CC_B$Group2),
                                         unlist(list_type1_snr2$DSCA$CC_B$Group2),
                                         unlist(list_type1_snr2$DGICA$CC_B$Group2),
                                         unlist(list_type1_snr2$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_type1_snr2$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr2_type1$SNR <- 0.5
CC_B_G2_snr3_type1 <- data.frame(value=c(unlist(list_type1_snr3$GRIDY$CC_B$Group2),
                                         unlist(list_type1_snr3$SCA_P$CC_B$Group2),
                                         unlist(list_type1_snr3$GICA$CC_B$Group2),
                                         unlist(list_type1_snr3$DSCA$CC_B$Group2),
                                         unlist(list_type1_snr3$DGICA$CC_B$Group2),
                                         unlist(list_type1_snr3$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_type1_snr3$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr3_type1$SNR <- 0.75
CC_B_G2_snr4_type1 <- data.frame(value=c(unlist(list_type1_snr4$GRIDY$CC_B$Group2),
                                         unlist(list_type1_snr4$SCA_P$CC_B$Group2),
                                         unlist(list_type1_snr4$GICA$CC_B$Group2),
                                         unlist(list_type1_snr4$DSCA$CC_B$Group2),
                                         unlist(list_type1_snr4$DGICA$CC_B$Group2),
                                         unlist(list_type1_snr4$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_type1_snr4$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr4_type1$SNR <- 1.25
CC_B_G2_snr5_type1 <- data.frame(value=c(unlist(list_type1_snr5$GRIDY$CC_B$Group2),
                                         unlist(list_type1_snr5$SCA_P$CC_B$Group2),
                                         unlist(list_type1_snr5$GICA$CC_B$Group2),
                                         unlist(list_type1_snr5$DSCA$CC_B$Group2),
                                         unlist(list_type1_snr5$DGICA$CC_B$Group2),
                                         unlist(list_type1_snr5$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_type1_snr5$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr5_type1$SNR <- 2
CC_B_G2_snr6_type1 <- data.frame(value=c(unlist(list_type1_snr6$GRIDY$CC_B$Group2),
                                         unlist(list_type1_snr6$SCA_P$CC_B$Group2),
                                         unlist(list_type1_snr6$GICA$CC_B$Group2),
                                         unlist(list_type1_snr6$DSCA$CC_B$Group2),
                                         unlist(list_type1_snr6$DGICA$CC_B$Group2),
                                         unlist(list_type1_snr6$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_type1_snr6$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr6_type1$SNR <- 4
CC_B_G2_type1 <- rbind(CC_B_G2_snr1_type1,CC_B_G2_snr2_type1,CC_B_G2_snr3_type1,
                       CC_B_G2_snr4_type1,CC_B_G2_snr5_type1,CC_B_G2_snr6_type1)
CC_B_G2_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                  "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

SNR_list <- c("0.25","0.5","0.75","1.25","2","4")
CC_B_J_type1$SNR <- factor(CC_B_J_type1$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
CC_B_G1_type1$SNR <- factor(CC_B_G1_type1$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
CC_B_G2_type1$SNR <- factor(CC_B_G2_type1$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)

method_list <- c("GRIDY","SCA-P","GICA","DSCA","DGICA","Unfitted_SCA-PF2","Unfitted_GICA")
CC_B_J_type1$method <- factor(CC_B_J_type1$method,levels=method_list,labels=method_list)
CC_B_G1_type1$method <- factor(CC_B_G1_type1$method,levels=method_list,labels=method_list)
CC_B_G2_type1$method <- factor(CC_B_G2_type1$method,levels=method_list,labels=method_list)

pl_CC_B_J_type1 <- ggplot(CC_B_J_type1,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_B$ (Joint)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_B_G1_type1 <- ggplot(CC_B_G1_type1,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_B$ (Group1)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_B_G2_type1 <- ggplot(CC_B_G2_type1,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_B$ (Group2)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))


# CC_F:
CC_F_J_snr1_type1 <- data.frame(value=c(unlist(list_type1_snr1$GRIDY$CC_F$Joint),
                                        unlist(list_type1_snr1$SCA_P$CC_F$Joint),
                                        unlist(list_type1_snr1$GICA$CC_F$Joint),
                                        unlist(list_type1_snr1$DSCA$CC_F$Joint),
                                        unlist(list_type1_snr1$DGICA$CC_F$Joint),
                                        unlist(list_type1_snr1$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_type1_snr1$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr1_type1$SNR <- 0.25
CC_F_J_snr2_type1 <- data.frame(value=c(unlist(list_type1_snr2$GRIDY$CC_F$Joint),
                                        unlist(list_type1_snr2$SCA_P$CC_F$Joint),
                                        unlist(list_type1_snr2$GICA$CC_F$Joint),
                                        unlist(list_type1_snr2$DSCA$CC_F$Joint),
                                        unlist(list_type1_snr2$DGICA$CC_F$Joint),
                                        unlist(list_type1_snr2$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_type1_snr2$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr2_type1$SNR <- 0.5
CC_F_J_snr3_type1 <- data.frame(value=c(unlist(list_type1_snr3$GRIDY$CC_F$Joint),
                                        unlist(list_type1_snr3$SCA_P$CC_F$Joint),
                                        unlist(list_type1_snr3$GICA$CC_F$Joint),
                                        unlist(list_type1_snr3$DSCA$CC_F$Joint),
                                        unlist(list_type1_snr3$DGICA$CC_F$Joint),
                                        unlist(list_type1_snr3$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_type1_snr3$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr3_type1$SNR <- 0.75
CC_F_J_snr4_type1 <- data.frame(value=c(unlist(list_type1_snr4$GRIDY$CC_F$Joint),
                                        unlist(list_type1_snr4$SCA_P$CC_F$Joint),
                                        unlist(list_type1_snr4$GICA$CC_F$Joint),
                                        unlist(list_type1_snr4$DSCA$CC_F$Joint),
                                        unlist(list_type1_snr4$DGICA$CC_F$Joint),
                                        unlist(list_type1_snr4$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_type1_snr4$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr4_type1$SNR <- 1.25
CC_F_J_snr5_type1 <- data.frame(value=c(unlist(list_type1_snr5$GRIDY$CC_F$Joint),
                                        unlist(list_type1_snr5$SCA_P$CC_F$Joint),
                                        unlist(list_type1_snr5$GICA$CC_F$Joint),
                                        unlist(list_type1_snr5$DSCA$CC_F$Joint),
                                        unlist(list_type1_snr5$DGICA$CC_F$Joint),
                                        unlist(list_type1_snr5$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_type1_snr5$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr5_type1$SNR <- 2
CC_F_J_snr6_type1 <- data.frame(value=c(unlist(list_type1_snr6$GRIDY$CC_F$Joint),
                                        unlist(list_type1_snr6$SCA_P$CC_F$Joint),
                                        unlist(list_type1_snr6$GICA$CC_F$Joint),
                                        unlist(list_type1_snr6$DSCA$CC_F$Joint),
                                        unlist(list_type1_snr6$DGICA$CC_F$Joint),
                                        unlist(list_type1_snr6$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_type1_snr6$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr6_type1$SNR <- 4
CC_F_J_type1 <- rbind(CC_F_J_snr1_type1,CC_F_J_snr2_type1,CC_F_J_snr3_type1,
                      CC_F_J_snr4_type1,CC_F_J_snr5_type1,CC_F_J_snr6_type1)
CC_F_J_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                 "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

CC_F_G1_snr1_type1 <- data.frame(value=c(unlist(list_type1_snr1$GRIDY$CC_F$Group1),
                                         unlist(list_type1_snr1$SCA_P$CC_F$Group1),
                                         unlist(list_type1_snr1$GICA$CC_F$Group1),
                                         unlist(list_type1_snr1$DSCA$CC_F$Group1),
                                         unlist(list_type1_snr1$DGICA$CC_F$Group1),
                                         unlist(list_type1_snr1$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_type1_snr1$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr1_type1$SNR <- 0.25
CC_F_G1_snr2_type1 <- data.frame(value=c(unlist(list_type1_snr2$GRIDY$CC_F$Group1),
                                         unlist(list_type1_snr2$SCA_P$CC_F$Group1),
                                         unlist(list_type1_snr2$GICA$CC_F$Group1),
                                         unlist(list_type1_snr2$DSCA$CC_F$Group1),
                                         unlist(list_type1_snr2$DGICA$CC_F$Group1),
                                         unlist(list_type1_snr2$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_type1_snr2$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr2_type1$SNR <- 0.5
CC_F_G1_snr3_type1 <- data.frame(value=c(unlist(list_type1_snr3$GRIDY$CC_F$Group1),
                                         unlist(list_type1_snr3$SCA_P$CC_F$Group1),
                                         unlist(list_type1_snr3$GICA$CC_F$Group1),
                                         unlist(list_type1_snr3$DSCA$CC_F$Group1),
                                         unlist(list_type1_snr3$DGICA$CC_F$Group1),
                                         unlist(list_type1_snr3$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_type1_snr3$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr3_type1$SNR <- 0.75
CC_F_G1_snr4_type1 <- data.frame(value=c(unlist(list_type1_snr4$GRIDY$CC_F$Group1),
                                         unlist(list_type1_snr4$SCA_P$CC_F$Group1),
                                         unlist(list_type1_snr4$GICA$CC_F$Group1),
                                         unlist(list_type1_snr4$DSCA$CC_F$Group1),
                                         unlist(list_type1_snr4$DGICA$CC_F$Group1),
                                         unlist(list_type1_snr4$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_type1_snr4$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr4_type1$SNR <- 1.25
CC_F_G1_snr5_type1 <- data.frame(value=c(unlist(list_type1_snr5$GRIDY$CC_F$Group1),
                                         unlist(list_type1_snr5$SCA_P$CC_F$Group1),
                                         unlist(list_type1_snr5$GICA$CC_F$Group1),
                                         unlist(list_type1_snr5$DSCA$CC_F$Group1),
                                         unlist(list_type1_snr5$DGICA$CC_F$Group1),
                                         unlist(list_type1_snr5$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_type1_snr5$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr5_type1$SNR <- 2
CC_F_G1_snr6_type1 <- data.frame(value=c(unlist(list_type1_snr6$GRIDY$CC_F$Group1),
                                         unlist(list_type1_snr6$SCA_P$CC_F$Group1),
                                         unlist(list_type1_snr6$GICA$CC_F$Group1),
                                         unlist(list_type1_snr6$DSCA$CC_F$Group1),
                                         unlist(list_type1_snr6$DGICA$CC_F$Group1),
                                         unlist(list_type1_snr6$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_type1_snr6$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr6_type1$SNR <- 4
CC_F_G1_type1 <- rbind(CC_F_G1_snr1_type1,CC_F_G1_snr2_type1,CC_F_G1_snr3_type1,
                       CC_F_G1_snr4_type1,CC_F_G1_snr5_type1,CC_F_G1_snr6_type1)
CC_F_G1_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                  "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


CC_F_G2_snr1_type1 <- data.frame(value=c(unlist(list_type1_snr1$GRIDY$CC_F$Group2),
                                         unlist(list_type1_snr1$SCA_P$CC_F$Group2),
                                         unlist(list_type1_snr1$GICA$CC_F$Group2),
                                         unlist(list_type1_snr1$DSCA$CC_F$Group2),
                                         unlist(list_type1_snr1$DGICA$CC_F$Group2),
                                         unlist(list_type1_snr1$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_type1_snr1$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr1_type1$SNR <- 0.25
CC_F_G2_snr2_type1 <- data.frame(value=c(unlist(list_type1_snr2$GRIDY$CC_F$Group2),
                                         unlist(list_type1_snr2$SCA_P$CC_F$Group2),
                                         unlist(list_type1_snr2$GICA$CC_F$Group2),
                                         unlist(list_type1_snr2$DSCA$CC_F$Group2),
                                         unlist(list_type1_snr2$DGICA$CC_F$Group2),
                                         unlist(list_type1_snr2$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_type1_snr2$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr2_type1$SNR <- 0.5
CC_F_G2_snr3_type1 <- data.frame(value=c(unlist(list_type1_snr3$GRIDY$CC_F$Group2),
                                         unlist(list_type1_snr3$SCA_P$CC_F$Group2),
                                         unlist(list_type1_snr3$GICA$CC_F$Group2),
                                         unlist(list_type1_snr3$DSCA$CC_F$Group2),
                                         unlist(list_type1_snr3$DGICA$CC_F$Group2),
                                         unlist(list_type1_snr3$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_type1_snr3$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr3_type1$SNR <- 0.75
CC_F_G2_snr4_type1 <- data.frame(value=c(unlist(list_type1_snr4$GRIDY$CC_F$Group2),
                                         unlist(list_type1_snr4$SCA_P$CC_F$Group2),
                                         unlist(list_type1_snr4$GICA$CC_F$Group2),
                                         unlist(list_type1_snr4$DSCA$CC_F$Group2),
                                         unlist(list_type1_snr4$DGICA$CC_F$Group2),
                                         unlist(list_type1_snr4$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_type1_snr4$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr4_type1$SNR <- 1.25
CC_F_G2_snr5_type1 <- data.frame(value=c(unlist(list_type1_snr5$GRIDY$CC_F$Group2),
                                         unlist(list_type1_snr5$SCA_P$CC_F$Group2),
                                         unlist(list_type1_snr5$GICA$CC_F$Group2),
                                         unlist(list_type1_snr5$DSCA$CC_F$Group2),
                                         unlist(list_type1_snr5$DGICA$CC_F$Group2),
                                         unlist(list_type1_snr5$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_type1_snr5$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr5_type1$SNR <- 2
CC_F_G2_snr6_type1 <- data.frame(value=c(unlist(list_type1_snr6$GRIDY$CC_F$Group2),
                                         unlist(list_type1_snr6$SCA_P$CC_F$Group2),
                                         unlist(list_type1_snr6$GICA$CC_F$Group2),
                                         unlist(list_type1_snr6$DSCA$CC_F$Group2),
                                         unlist(list_type1_snr6$DGICA$CC_F$Group2),
                                         unlist(list_type1_snr6$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_type1_snr6$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr6_type1$SNR <- 4
CC_F_G2_type1 <- rbind(CC_F_G2_snr1_type1,CC_F_G2_snr2_type1,CC_F_G2_snr3_type1,
                       CC_F_G2_snr4_type1,CC_F_G2_snr5_type1,CC_F_G2_snr6_type1)
CC_F_G2_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                  "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

SNR_list <- c("0.25","0.5","0.75","1.25","2","4")
CC_F_J_type1$SNR <- factor(CC_F_J_type1$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),SNR_list)
CC_F_G1_type1$SNR <- factor(CC_F_G1_type1$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),SNR_list)
CC_F_G2_type1$SNR <- factor(CC_F_G2_type1$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),SNR_list)

method_list <- c("GRIDY","SCA-P","GICA","DSCA","DGICA","Unfitted_SCA-PF2","Unfitted_GICA")
CC_F_J_type1$method <- factor(CC_F_J_type1$method,levels=method_list,labels=method_list)
CC_F_G1_type1$method <- factor(CC_F_G1_type1$method,levels=method_list,labels=method_list)
CC_F_G2_type1$method <- factor(CC_F_G2_type1$method,levels=method_list,labels=method_list)

pl_CC_F_J_type1 <- ggplot(CC_F_J_type1,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_F$ (Joint)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_F_G1_type1 <- ggplot(CC_F_G1_type1,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_F$ (Group1)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_F_G2_type1 <- ggplot(CC_F_G2_type1,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_F$ (Group2)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))



# CC_PSI:
CC_PSI_J_snr1_type1 <- data.frame(value=c(unlist(list_type1_snr1$GRIDY$CC_PSI$Joint),
                                          unlist(list_type1_snr1$SCA_P$CC_PSI$Joint),
                                          unlist(list_type1_snr1$GICA$CC_PSI$Joint),
                                          unlist(list_type1_snr1$DSCA$CC_PSI$Joint),
                                          unlist(list_type1_snr1$DGICA$CC_PSI$Joint),
                                          unlist(list_type1_snr1$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_type1_snr1$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr1_type1$SNR <- 0.25
CC_PSI_J_snr2_type1 <- data.frame(value=c(unlist(list_type1_snr2$GRIDY$CC_PSI$Joint),
                                          unlist(list_type1_snr2$SCA_P$CC_PSI$Joint),
                                          unlist(list_type1_snr2$GICA$CC_PSI$Joint),
                                          unlist(list_type1_snr2$DSCA$CC_PSI$Joint),
                                          unlist(list_type1_snr2$DGICA$CC_PSI$Joint),
                                          unlist(list_type1_snr2$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_type1_snr2$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr2_type1$SNR <- 0.5
CC_PSI_J_snr3_type1 <- data.frame(value=c(unlist(list_type1_snr3$GRIDY$CC_PSI$Joint),
                                          unlist(list_type1_snr3$SCA_P$CC_PSI$Joint),
                                          unlist(list_type1_snr3$GICA$CC_PSI$Joint),
                                          unlist(list_type1_snr3$DSCA$CC_PSI$Joint),
                                          unlist(list_type1_snr3$DGICA$CC_PSI$Joint),
                                          unlist(list_type1_snr3$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_type1_snr3$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr3_type1$SNR <- 0.75
CC_PSI_J_snr4_type1 <- data.frame(value=c(unlist(list_type1_snr4$GRIDY$CC_PSI$Joint),
                                          unlist(list_type1_snr4$SCA_P$CC_PSI$Joint),
                                          unlist(list_type1_snr4$GICA$CC_PSI$Joint),
                                          unlist(list_type1_snr4$DSCA$CC_PSI$Joint),
                                          unlist(list_type1_snr4$DGICA$CC_PSI$Joint),
                                          unlist(list_type1_snr4$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_type1_snr4$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr4_type1$SNR <- 1.25
CC_PSI_J_snr5_type1 <- data.frame(value=c(unlist(list_type1_snr5$GRIDY$CC_PSI$Joint),
                                          unlist(list_type1_snr5$SCA_P$CC_PSI$Joint),
                                          unlist(list_type1_snr5$GICA$CC_PSI$Joint),
                                          unlist(list_type1_snr5$DSCA$CC_PSI$Joint),
                                          unlist(list_type1_snr5$DGICA$CC_PSI$Joint),
                                          unlist(list_type1_snr5$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_type1_snr5$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr5_type1$SNR <- 2
CC_PSI_J_snr6_type1 <- data.frame(value=c(unlist(list_type1_snr6$GRIDY$CC_PSI$Joint),
                                          unlist(list_type1_snr6$SCA_P$CC_PSI$Joint),
                                          unlist(list_type1_snr6$GICA$CC_PSI$Joint),
                                          unlist(list_type1_snr6$DSCA$CC_PSI$Joint),
                                          unlist(list_type1_snr6$DGICA$CC_PSI$Joint),
                                          unlist(list_type1_snr6$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_type1_snr6$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr6_type1$SNR <- 4
CC_PSI_J_type1 <- rbind(CC_PSI_J_snr1_type1,CC_PSI_J_snr2_type1,CC_PSI_J_snr3_type1,
                        CC_PSI_J_snr4_type1,CC_PSI_J_snr5_type1,CC_PSI_J_snr6_type1)
CC_PSI_J_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                   "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

CC_PSI_G1_snr1_type1 <- data.frame(value=c(unlist(list_type1_snr1$GRIDY$CC_PSI$Group1),
                                           unlist(list_type1_snr1$SCA_P$CC_PSI$Group1),
                                           unlist(list_type1_snr1$GICA$CC_PSI$Group1),
                                           unlist(list_type1_snr1$DSCA$CC_PSI$Group1),
                                           unlist(list_type1_snr1$DGICA$CC_PSI$Group1),
                                           unlist(list_type1_snr1$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_type1_snr1$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr1_type1$SNR <- 0.25
CC_PSI_G1_snr2_type1 <- data.frame(value=c(unlist(list_type1_snr2$GRIDY$CC_PSI$Group1),
                                           unlist(list_type1_snr2$SCA_P$CC_PSI$Group1),
                                           unlist(list_type1_snr2$GICA$CC_PSI$Group1),
                                           unlist(list_type1_snr2$DSCA$CC_PSI$Group1),
                                           unlist(list_type1_snr2$DGICA$CC_PSI$Group1),
                                           unlist(list_type1_snr2$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_type1_snr2$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr2_type1$SNR <- 0.5
CC_PSI_G1_snr3_type1 <- data.frame(value=c(unlist(list_type1_snr3$GRIDY$CC_PSI$Group1),
                                           unlist(list_type1_snr3$SCA_P$CC_PSI$Group1),
                                           unlist(list_type1_snr3$GICA$CC_PSI$Group1),
                                           unlist(list_type1_snr3$DSCA$CC_PSI$Group1),
                                           unlist(list_type1_snr3$DGICA$CC_PSI$Group1),
                                           unlist(list_type1_snr3$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_type1_snr3$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr3_type1$SNR <- 0.75
CC_PSI_G1_snr4_type1 <- data.frame(value=c(unlist(list_type1_snr4$GRIDY$CC_PSI$Group1),
                                           unlist(list_type1_snr4$SCA_P$CC_PSI$Group1),
                                           unlist(list_type1_snr4$GICA$CC_PSI$Group1),
                                           unlist(list_type1_snr4$DSCA$CC_PSI$Group1),
                                           unlist(list_type1_snr4$DGICA$CC_PSI$Group1),
                                           unlist(list_type1_snr4$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_type1_snr4$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr4_type1$SNR <- 1.25
CC_PSI_G1_snr5_type1 <- data.frame(value=c(unlist(list_type1_snr5$GRIDY$CC_PSI$Group1),
                                           unlist(list_type1_snr5$SCA_P$CC_PSI$Group1),
                                           unlist(list_type1_snr5$GICA$CC_PSI$Group1),
                                           unlist(list_type1_snr5$DSCA$CC_PSI$Group1),
                                           unlist(list_type1_snr5$DGICA$CC_PSI$Group1),
                                           unlist(list_type1_snr5$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_type1_snr5$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr5_type1$SNR <- 2
CC_PSI_G1_snr6_type1 <- data.frame(value=c(unlist(list_type1_snr6$GRIDY$CC_PSI$Group1),
                                           unlist(list_type1_snr6$SCA_P$CC_PSI$Group1),
                                           unlist(list_type1_snr6$GICA$CC_PSI$Group1),
                                           unlist(list_type1_snr6$DSCA$CC_PSI$Group1),
                                           unlist(list_type1_snr6$DGICA$CC_PSI$Group1),
                                           unlist(list_type1_snr6$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_type1_snr6$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr6_type1$SNR <- 4
CC_PSI_G1_type1 <- rbind(CC_PSI_G1_snr1_type1,CC_PSI_G1_snr2_type1,CC_PSI_G1_snr3_type1,
                         CC_PSI_G1_snr4_type1,CC_PSI_G1_snr5_type1,CC_PSI_G1_snr6_type1)
CC_PSI_G1_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                    "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


CC_PSI_G2_snr1_type1 <- data.frame(value=c(unlist(list_type1_snr1$GRIDY$CC_PSI$Group2),
                                           unlist(list_type1_snr1$SCA_P$CC_PSI$Group2),
                                           unlist(list_type1_snr1$GICA$CC_PSI$Group2),
                                           unlist(list_type1_snr1$DSCA$CC_PSI$Group2),
                                           unlist(list_type1_snr1$DGICA$CC_PSI$Group2),
                                           unlist(list_type1_snr1$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_type1_snr1$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr1_type1$SNR <- 0.25
CC_PSI_G2_snr2_type1 <- data.frame(value=c(unlist(list_type1_snr2$GRIDY$CC_PSI$Group2),
                                           unlist(list_type1_snr2$SCA_P$CC_PSI$Group2),
                                           unlist(list_type1_snr2$GICA$CC_PSI$Group2),
                                           unlist(list_type1_snr2$DSCA$CC_PSI$Group2),
                                           unlist(list_type1_snr2$DGICA$CC_PSI$Group2),
                                           unlist(list_type1_snr2$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_type1_snr2$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr2_type1$SNR <- 0.5
CC_PSI_G2_snr3_type1 <- data.frame(value=c(unlist(list_type1_snr3$GRIDY$CC_PSI$Group2),
                                           unlist(list_type1_snr3$SCA_P$CC_PSI$Group2),
                                           unlist(list_type1_snr3$GICA$CC_PSI$Group2),
                                           unlist(list_type1_snr3$DSCA$CC_PSI$Group2),
                                           unlist(list_type1_snr3$DGICA$CC_PSI$Group2),
                                           unlist(list_type1_snr3$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_type1_snr3$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr3_type1$SNR <- 0.75
CC_PSI_G2_snr4_type1 <- data.frame(value=c(unlist(list_type1_snr4$GRIDY$CC_PSI$Group2),
                                           unlist(list_type1_snr4$SCA_P$CC_PSI$Group2),
                                           unlist(list_type1_snr4$GICA$CC_PSI$Group2),
                                           unlist(list_type1_snr4$DSCA$CC_PSI$Group2),
                                           unlist(list_type1_snr4$DGICA$CC_PSI$Group2),
                                           unlist(list_type1_snr4$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_type1_snr4$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr4_type1$SNR <- 1.25
CC_PSI_G2_snr5_type1 <- data.frame(value=c(unlist(list_type1_snr5$GRIDY$CC_PSI$Group2),
                                           unlist(list_type1_snr5$SCA_P$CC_PSI$Group2),
                                           unlist(list_type1_snr5$GICA$CC_PSI$Group2),
                                           unlist(list_type1_snr5$DSCA$CC_PSI$Group2),
                                           unlist(list_type1_snr5$DGICA$CC_PSI$Group2),
                                           unlist(list_type1_snr5$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_type1_snr5$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr5_type1$SNR <- 2
CC_PSI_G2_snr6_type1 <- data.frame(value=c(unlist(list_type1_snr6$GRIDY$CC_PSI$Group2),
                                           unlist(list_type1_snr6$SCA_P$CC_PSI$Group2),
                                           unlist(list_type1_snr6$GICA$CC_PSI$Group2),
                                           unlist(list_type1_snr6$DSCA$CC_PSI$Group2),
                                           unlist(list_type1_snr6$DGICA$CC_PSI$Group2),
                                           unlist(list_type1_snr6$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_type1_snr6$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr6_type1$SNR <- 4
CC_PSI_G2_type1 <- rbind(CC_PSI_G2_snr1_type1,CC_PSI_G2_snr2_type1,CC_PSI_G2_snr3_type1,
                         CC_PSI_G2_snr4_type1,CC_PSI_G2_snr5_type1,CC_PSI_G2_snr6_type1)
CC_PSI_G2_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                    "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

SNR_list <- c("0.25","0.5","0.75","1.25","2","4")
CC_PSI_J_type1$SNR <- factor(CC_PSI_J_type1$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
CC_PSI_G1_type1$SNR <- factor(CC_PSI_G1_type1$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
CC_PSI_G2_type1$SNR <- factor(CC_PSI_G2_type1$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)

method_list <- c("GRIDY","SCA-P","GICA","DSCA","DGICA","Unfitted_SCA-PF2","Unfitted_GICA")
CC_PSI_J_type1$method <- factor(CC_PSI_J_type1$method,levels=method_list,labels=method_list)
CC_PSI_G1_type1$method <- factor(CC_PSI_G1_type1$method,levels=method_list,labels=method_list)
CC_PSI_G2_type1$method <- factor(CC_PSI_G2_type1$method,levels=method_list,labels=method_list)

pl_CC_PSI_J_type1 <- ggplot(CC_PSI_J_type1,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_{\\Psi}$ (Joint)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_PSI_G1_type1 <- ggplot(CC_PSI_G1_type1,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_{\\Psi}$ (Group1)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_PSI_G2_type1 <- ggplot(CC_PSI_G2_type1,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_{\\Psi}$ (Group2)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))


# Plotting:
mylegend_type1 <- g_legend(pl_R2_J_type1)
grid.arrange(arrangeGrob(pl_R2_J_type1 + theme(legend.position="none"),
                         pl_R2_G1_type1 + theme(legend.position="none"), 
                         pl_R2_G2_type1 + theme(legend.position="none"),
                         pl_RMSE_J_type1 + theme(legend.position="none"),
                         pl_RMSE_G1_type1 + theme(legend.position="none"), 
                         pl_RMSE_G2_type1 + theme(legend.position="none"),
                         pl_CC_B_J_type1 + theme(legend.position="none"),
                         pl_CC_B_G1_type1 + theme(legend.position="none"), 
                         pl_CC_B_G2_type1 + theme(legend.position="none"),
                         pl_CC_F_J_type1 + theme(legend.position="none"),
                         pl_CC_F_G1_type1 + theme(legend.position="none"), 
                         pl_CC_F_G2_type1 + theme(legend.position="none"),
                         pl_CC_PSI_J_type1 + theme(legend.position="none"),
                         pl_CC_PSI_G1_type1 + theme(legend.position="none"), 
                         pl_CC_PSI_G2_type1 + theme(legend.position="none"),
                         nrow=5,ncol=3), mylegend_type1,
             nrow = 2, heights = c(5,0.2))


#-----------------------------------------------------------------------------#
# type 2:
#-----------------------------------------------------------------------------#
# R2:
R2_J_snr1_type2 <- data.frame(value=c(unlist(list_type2_snr1$GRIDY$R2$Joint),
                                      unlist(list_type2_snr1$SCA_P$R2$Joint),
                                      unlist(list_type2_snr1$GICA$R2$Joint),
                                      unlist(list_type2_snr1$DSCA$R2$Joint),
                                      unlist(list_type2_snr1$DGICA$R2$Joint),
                                      unlist(list_type2_snr1$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_type2_snr1$Unfitted_GICA$R2$Joint)))
R2_J_snr1_type2$SNR <- 0.25
R2_J_snr2_type2 <- data.frame(value=c(unlist(list_type2_snr2$GRIDY$R2$Joint),
                                      unlist(list_type2_snr2$SCA_P$R2$Joint),
                                      unlist(list_type2_snr2$GICA$R2$Joint),
                                      unlist(list_type2_snr2$DSCA$R2$Joint),
                                      unlist(list_type2_snr2$DGICA$R2$Joint),
                                      unlist(list_type2_snr2$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_type2_snr2$Unfitted_GICA$R2$Joint)))
R2_J_snr2_type2$SNR <- 0.5
R2_J_snr3_type2 <- data.frame(value=c(unlist(list_type2_snr3$GRIDY$R2$Joint),
                                      unlist(list_type2_snr3$SCA_P$R2$Joint),
                                      unlist(list_type2_snr3$GICA$R2$Joint),
                                      unlist(list_type2_snr3$DSCA$R2$Joint),
                                      unlist(list_type2_snr3$DGICA$R2$Joint),
                                      unlist(list_type2_snr3$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_type2_snr3$Unfitted_GICA$R2$Joint)))
R2_J_snr3_type2$SNR <- 0.75
R2_J_snr4_type2 <- data.frame(value=c(unlist(list_type2_snr4$GRIDY$R2$Joint),
                                      unlist(list_type2_snr4$SCA_P$R2$Joint),
                                      unlist(list_type2_snr4$GICA$R2$Joint),
                                      unlist(list_type2_snr4$DSCA$R2$Joint),
                                      unlist(list_type2_snr4$DGICA$R2$Joint),
                                      unlist(list_type2_snr4$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_type2_snr4$Unfitted_GICA$R2$Joint)))
R2_J_snr4_type2$SNR <- 1.25
R2_J_snr5_type2 <- data.frame(value=c(unlist(list_type2_snr5$GRIDY$R2$Joint),
                                      unlist(list_type2_snr5$SCA_P$R2$Joint),
                                      unlist(list_type2_snr5$GICA$R2$Joint),
                                      unlist(list_type2_snr5$DSCA$R2$Joint),
                                      unlist(list_type2_snr5$DGICA$R2$Joint),
                                      unlist(list_type2_snr5$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_type2_snr5$Unfitted_GICA$R2$Joint)))
R2_J_snr5_type2$SNR <- 2
R2_J_snr6_type2 <- data.frame(value=c(unlist(list_type2_snr6$GRIDY$R2$Joint),
                                      unlist(list_type2_snr6$SCA_P$R2$Joint),
                                      unlist(list_type2_snr6$GICA$R2$Joint),
                                      unlist(list_type2_snr6$DSCA$R2$Joint),
                                      unlist(list_type2_snr6$DGICA$R2$Joint),
                                      unlist(list_type2_snr6$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_type2_snr6$Unfitted_GICA$R2$Joint)))
R2_J_snr6_type2$SNR <- 4
R2_J_type2 <- rbind(R2_J_snr1_type2,R2_J_snr2_type2,R2_J_snr3_type2,
                    R2_J_snr4_type2,R2_J_snr5_type2,R2_J_snr6_type2)
R2_J_type2$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                               "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

R2_G1_snr1_type2 <- data.frame(value=c(unlist(list_type2_snr1$GRIDY$R2$Group1),
                                       unlist(list_type2_snr1$SCA_P$R2$Group1),
                                       unlist(list_type2_snr1$GICA$R2$Group1),
                                       unlist(list_type2_snr1$DSCA$R2$Group1),
                                       unlist(list_type2_snr1$DGICA$R2$Group1),
                                       unlist(list_type2_snr1$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_type2_snr1$Unfitted_GICA$R2$Group1)))
R2_G1_snr1_type2$SNR <- 0.25
R2_G1_snr2_type2 <- data.frame(value=c(unlist(list_type2_snr2$GRIDY$R2$Group1),
                                       unlist(list_type2_snr2$SCA_P$R2$Group1),
                                       unlist(list_type2_snr2$GICA$R2$Group1),
                                       unlist(list_type2_snr2$DSCA$R2$Group1),
                                       unlist(list_type2_snr2$DGICA$R2$Group1),
                                       unlist(list_type2_snr2$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_type2_snr2$Unfitted_GICA$R2$Group1)))
R2_G1_snr2_type2$SNR <- 0.5
R2_G1_snr3_type2 <- data.frame(value=c(unlist(list_type2_snr3$GRIDY$R2$Group1),
                                       unlist(list_type2_snr3$SCA_P$R2$Group1),
                                       unlist(list_type2_snr3$GICA$R2$Group1),
                                       unlist(list_type2_snr3$DSCA$R2$Group1),
                                       unlist(list_type2_snr3$DGICA$R2$Group1),
                                       unlist(list_type2_snr3$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_type2_snr3$Unfitted_GICA$R2$Group1)))
R2_G1_snr3_type2$SNR <- 0.75
R2_G1_snr4_type2 <- data.frame(value=c(unlist(list_type2_snr4$GRIDY$R2$Group1),
                                       unlist(list_type2_snr4$SCA_P$R2$Group1),
                                       unlist(list_type2_snr4$GICA$R2$Group1),
                                       unlist(list_type2_snr4$DSCA$R2$Group1),
                                       unlist(list_type2_snr4$DGICA$R2$Group1),
                                       unlist(list_type2_snr4$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_type2_snr4$Unfitted_GICA$R2$Group1)))
R2_G1_snr4_type2$SNR <- 1.25
R2_G1_snr5_type2 <- data.frame(value=c(unlist(list_type2_snr5$GRIDY$R2$Group1),
                                       unlist(list_type2_snr5$SCA_P$R2$Group1),
                                       unlist(list_type2_snr5$GICA$R2$Group1),
                                       unlist(list_type2_snr5$DSCA$R2$Group1),
                                       unlist(list_type2_snr5$DGICA$R2$Group1),
                                       unlist(list_type2_snr5$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_type2_snr5$Unfitted_GICA$R2$Group1)))
R2_G1_snr5_type2$SNR <- 2
R2_G1_snr6_type2 <- data.frame(value=c(unlist(list_type2_snr6$GRIDY$R2$Group1),
                                       unlist(list_type2_snr6$SCA_P$R2$Group1),
                                       unlist(list_type2_snr6$GICA$R2$Group1),
                                       unlist(list_type2_snr6$DSCA$R2$Group1),
                                       unlist(list_type2_snr6$DGICA$R2$Group1),
                                       unlist(list_type2_snr6$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_type2_snr6$Unfitted_GICA$R2$Group1)))
R2_G1_snr6_type2$SNR <- 4
R2_G1_type2 <- rbind(R2_G1_snr1_type2,R2_G1_snr2_type2,R2_G1_snr3_type2,
                     R2_G1_snr4_type2,R2_G1_snr5_type2,R2_G1_snr6_type2)
R2_G1_type2$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


R2_G2_snr1_type2 <- data.frame(value=c(unlist(list_type2_snr1$GRIDY$R2$Group2),
                                       unlist(list_type2_snr1$SCA_P$R2$Group2),
                                       unlist(list_type2_snr1$GICA$R2$Group2),
                                       unlist(list_type2_snr1$DSCA$R2$Group2),
                                       unlist(list_type2_snr1$DGICA$R2$Group2),
                                       unlist(list_type2_snr1$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_type2_snr1$Unfitted_GICA$R2$Group2)))
R2_G2_snr1_type2$SNR <- 0.25
R2_G2_snr2_type2 <- data.frame(value=c(unlist(list_type2_snr2$GRIDY$R2$Group2),
                                       unlist(list_type2_snr2$SCA_P$R2$Group2),
                                       unlist(list_type2_snr2$GICA$R2$Group2),
                                       unlist(list_type2_snr2$DSCA$R2$Group2),
                                       unlist(list_type2_snr2$DGICA$R2$Group2),
                                       unlist(list_type2_snr2$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_type2_snr2$Unfitted_GICA$R2$Group2)))
R2_G2_snr2_type2$SNR <- 0.5
R2_G2_snr3_type2 <- data.frame(value=c(unlist(list_type2_snr3$GRIDY$R2$Group2),
                                       unlist(list_type2_snr3$SCA_P$R2$Group2),
                                       unlist(list_type2_snr3$GICA$R2$Group2),
                                       unlist(list_type2_snr3$DSCA$R2$Group2),
                                       unlist(list_type2_snr3$DGICA$R2$Group2),
                                       unlist(list_type2_snr3$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_type2_snr3$Unfitted_GICA$R2$Group2)))
R2_G2_snr3_type2$SNR <- 0.75
R2_G2_snr4_type2 <- data.frame(value=c(unlist(list_type2_snr4$GRIDY$R2$Group2),
                                       unlist(list_type2_snr4$SCA_P$R2$Group2),
                                       unlist(list_type2_snr4$GICA$R2$Group2),
                                       unlist(list_type2_snr4$DSCA$R2$Group2),
                                       unlist(list_type2_snr4$DGICA$R2$Group2),
                                       unlist(list_type2_snr4$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_type2_snr4$Unfitted_GICA$R2$Group2)))
R2_G2_snr4_type2$SNR <- 1.25
R2_G2_snr5_type2 <- data.frame(value=c(unlist(list_type2_snr5$GRIDY$R2$Group2),
                                       unlist(list_type2_snr5$SCA_P$R2$Group2),
                                       unlist(list_type2_snr5$GICA$R2$Group2),
                                       unlist(list_type2_snr5$DSCA$R2$Group2),
                                       unlist(list_type2_snr5$DGICA$R2$Group2),
                                       unlist(list_type2_snr5$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_type2_snr5$Unfitted_GICA$R2$Group2)))
R2_G2_snr5_type2$SNR <- 2
R2_G2_snr6_type2 <- data.frame(value=c(unlist(list_type2_snr6$GRIDY$R2$Group2),
                                       unlist(list_type2_snr6$SCA_P$R2$Group2),
                                       unlist(list_type2_snr6$GICA$R2$Group2),
                                       unlist(list_type2_snr6$DSCA$R2$Group2),
                                       unlist(list_type2_snr6$DGICA$R2$Group2),
                                       unlist(list_type2_snr6$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_type2_snr6$Unfitted_GICA$R2$Group2)))
R2_G2_snr6_type2$SNR <- 4
R2_G2_type2 <- rbind(R2_G2_snr1_type2,R2_G2_snr2_type2,R2_G2_snr3_type2,
                     R2_G2_snr4_type2,R2_G2_snr5_type2,R2_G2_snr6_type2)
R2_G2_type2$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


SNR_list <- c("0.25","0.5","0.75","1.25","2","4")
R2_J_type2$SNR <- factor(R2_J_type2$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
R2_G1_type2$SNR <- factor(R2_G1_type2$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
R2_G2_type2$SNR <- factor(R2_G2_type2$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)

method_list <- c("GRIDY","SCA-P","GICA","DSCA","DGICA","Unfitted_SCA-PF2","Unfitted_GICA")
R2_J_type2$method <- factor(R2_J_type2$method,levels=method_list,labels=method_list)
R2_G1_type2$method <- factor(R2_G1_type2$method,levels=method_list,labels=method_list)
R2_G2_type2$method <- factor(R2_G2_type2$method,levels=method_list,labels=method_list)

pl_R2_J_type2 <- ggplot(R2_J_type2,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$R^2$ (Joint)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_R2_G1_type2 <- ggplot(R2_G1_type2,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$R^2$ (Group1)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_R2_G2_type2 <- ggplot(R2_G2_type2,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$R^2$ (Group2)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))


# RMSE:
RMSE_J_snr1_type2 <- data.frame(value=c(unlist(list_type2_snr1$GRIDY$RMSE$Joint),
                                        unlist(list_type2_snr1$SCA_P$RMSE$Joint),
                                        unlist(list_type2_snr1$GICA$RMSE$Joint),
                                        unlist(list_type2_snr1$DSCA$RMSE$Joint),
                                        unlist(list_type2_snr1$DGICA$RMSE$Joint),
                                        unlist(list_type2_snr1$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_type2_snr1$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr1_type2$SNR <- 0.25
RMSE_J_snr2_type2 <- data.frame(value=c(unlist(list_type2_snr2$GRIDY$RMSE$Joint),
                                        unlist(list_type2_snr2$SCA_P$RMSE$Joint),
                                        unlist(list_type2_snr2$GICA$RMSE$Joint),
                                        unlist(list_type2_snr2$DSCA$RMSE$Joint),
                                        unlist(list_type2_snr2$DGICA$RMSE$Joint),
                                        unlist(list_type2_snr2$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_type2_snr2$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr2_type2$SNR <- 0.5
RMSE_J_snr3_type2 <- data.frame(value=c(unlist(list_type2_snr3$GRIDY$RMSE$Joint),
                                        unlist(list_type2_snr3$SCA_P$RMSE$Joint),
                                        unlist(list_type2_snr3$GICA$RMSE$Joint),
                                        unlist(list_type2_snr3$DSCA$RMSE$Joint),
                                        unlist(list_type2_snr3$DGICA$RMSE$Joint),
                                        unlist(list_type2_snr3$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_type2_snr3$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr3_type2$SNR <- 0.75
RMSE_J_snr4_type2 <- data.frame(value=c(unlist(list_type2_snr4$GRIDY$RMSE$Joint),
                                        unlist(list_type2_snr4$SCA_P$RMSE$Joint),
                                        unlist(list_type2_snr4$GICA$RMSE$Joint),
                                        unlist(list_type2_snr4$DSCA$RMSE$Joint),
                                        unlist(list_type2_snr4$DGICA$RMSE$Joint),
                                        unlist(list_type2_snr4$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_type2_snr4$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr4_type2$SNR <- 1.25
RMSE_J_snr5_type2 <- data.frame(value=c(unlist(list_type2_snr5$GRIDY$RMSE$Joint),
                                        unlist(list_type2_snr5$SCA_P$RMSE$Joint),
                                        unlist(list_type2_snr5$GICA$RMSE$Joint),
                                        unlist(list_type2_snr5$DSCA$RMSE$Joint),
                                        unlist(list_type2_snr5$DGICA$RMSE$Joint),
                                        unlist(list_type2_snr5$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_type2_snr5$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr5_type2$SNR <- 2
RMSE_J_snr6_type2 <- data.frame(value=c(unlist(list_type2_snr6$GRIDY$RMSE$Joint),
                                        unlist(list_type2_snr6$SCA_P$RMSE$Joint),
                                        unlist(list_type2_snr6$GICA$RMSE$Joint),
                                        unlist(list_type2_snr6$DSCA$RMSE$Joint),
                                        unlist(list_type2_snr6$DGICA$RMSE$Joint),
                                        unlist(list_type2_snr6$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_type2_snr6$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr6_type2$SNR <- 4
RMSE_J_type2 <- rbind(RMSE_J_snr1_type2,RMSE_J_snr2_type2,RMSE_J_snr3_type2,
                      RMSE_J_snr4_type2,RMSE_J_snr5_type2,RMSE_J_snr6_type2)
RMSE_J_type2$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                 "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

RMSE_G1_snr1_type2 <- data.frame(value=c(unlist(list_type2_snr1$GRIDY$RMSE$Group1),
                                         unlist(list_type2_snr1$SCA_P$RMSE$Group1),
                                         unlist(list_type2_snr1$GICA$RMSE$Group1),
                                         unlist(list_type2_snr1$DSCA$RMSE$Group1),
                                         unlist(list_type2_snr1$DGICA$RMSE$Group1),
                                         unlist(list_type2_snr1$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_type2_snr1$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr1_type2$SNR <- 0.25
RMSE_G1_snr2_type2 <- data.frame(value=c(unlist(list_type2_snr2$GRIDY$RMSE$Group1),
                                         unlist(list_type2_snr2$SCA_P$RMSE$Group1),
                                         unlist(list_type2_snr2$GICA$RMSE$Group1),
                                         unlist(list_type2_snr2$DSCA$RMSE$Group1),
                                         unlist(list_type2_snr2$DGICA$RMSE$Group1),
                                         unlist(list_type2_snr2$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_type2_snr2$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr2_type2$SNR <- 0.5
RMSE_G1_snr3_type2 <- data.frame(value=c(unlist(list_type2_snr3$GRIDY$RMSE$Group1),
                                         unlist(list_type2_snr3$SCA_P$RMSE$Group1),
                                         unlist(list_type2_snr3$GICA$RMSE$Group1),
                                         unlist(list_type2_snr3$DSCA$RMSE$Group1),
                                         unlist(list_type2_snr3$DGICA$RMSE$Group1),
                                         unlist(list_type2_snr3$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_type2_snr3$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr3_type2$SNR <- 0.75
RMSE_G1_snr4_type2 <- data.frame(value=c(unlist(list_type2_snr4$GRIDY$RMSE$Group1),
                                         unlist(list_type2_snr4$SCA_P$RMSE$Group1),
                                         unlist(list_type2_snr4$GICA$RMSE$Group1),
                                         unlist(list_type2_snr4$DSCA$RMSE$Group1),
                                         unlist(list_type2_snr4$DGICA$RMSE$Group1),
                                         unlist(list_type2_snr4$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_type2_snr4$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr4_type2$SNR <- 1.25
RMSE_G1_snr5_type2 <- data.frame(value=c(unlist(list_type2_snr5$GRIDY$RMSE$Group1),
                                         unlist(list_type2_snr5$SCA_P$RMSE$Group1),
                                         unlist(list_type2_snr5$GICA$RMSE$Group1),
                                         unlist(list_type2_snr5$DSCA$RMSE$Group1),
                                         unlist(list_type2_snr5$DGICA$RMSE$Group1),
                                         unlist(list_type2_snr5$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_type2_snr5$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr5_type2$SNR <- 2
RMSE_G1_snr6_type2 <- data.frame(value=c(unlist(list_type2_snr6$GRIDY$RMSE$Group1),
                                         unlist(list_type2_snr6$SCA_P$RMSE$Group1),
                                         unlist(list_type2_snr6$GICA$RMSE$Group1),
                                         unlist(list_type2_snr6$DSCA$RMSE$Group1),
                                         unlist(list_type2_snr6$DGICA$RMSE$Group1),
                                         unlist(list_type2_snr6$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_type2_snr6$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr6_type2$SNR <- 4
RMSE_G1_type2 <- rbind(RMSE_G1_snr1_type2,RMSE_G1_snr2_type2,RMSE_G1_snr3_type2,
                       RMSE_G1_snr4_type2,RMSE_G1_snr5_type2,RMSE_G1_snr6_type2)
RMSE_G1_type2$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                  "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


RMSE_G2_snr1_type2 <- data.frame(value=c(unlist(list_type2_snr1$GRIDY$RMSE$Group2),
                                         unlist(list_type2_snr1$SCA_P$RMSE$Group2),
                                         unlist(list_type2_snr1$GICA$RMSE$Group2),
                                         unlist(list_type2_snr1$DSCA$RMSE$Group2),
                                         unlist(list_type2_snr1$DGICA$RMSE$Group2),
                                         unlist(list_type2_snr1$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_type2_snr1$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr1_type2$SNR <- 0.25
RMSE_G2_snr2_type2 <- data.frame(value=c(unlist(list_type2_snr2$GRIDY$RMSE$Group2),
                                         unlist(list_type2_snr2$SCA_P$RMSE$Group2),
                                         unlist(list_type2_snr2$GICA$RMSE$Group2),
                                         unlist(list_type2_snr2$DSCA$RMSE$Group2),
                                         unlist(list_type2_snr2$DGICA$RMSE$Group2),
                                         unlist(list_type2_snr2$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_type2_snr2$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr2_type2$SNR <- 0.5
RMSE_G2_snr3_type2 <- data.frame(value=c(unlist(list_type2_snr3$GRIDY$RMSE$Group2),
                                         unlist(list_type2_snr3$SCA_P$RMSE$Group2),
                                         unlist(list_type2_snr3$GICA$RMSE$Group2),
                                         unlist(list_type2_snr3$DSCA$RMSE$Group2),
                                         unlist(list_type2_snr3$DGICA$RMSE$Group2),
                                         unlist(list_type2_snr3$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_type2_snr3$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr3_type2$SNR <- 0.75
RMSE_G2_snr4_type2 <- data.frame(value=c(unlist(list_type2_snr4$GRIDY$RMSE$Group2),
                                         unlist(list_type2_snr4$SCA_P$RMSE$Group2),
                                         unlist(list_type2_snr4$GICA$RMSE$Group2),
                                         unlist(list_type2_snr4$DSCA$RMSE$Group2),
                                         unlist(list_type2_snr4$DGICA$RMSE$Group2),
                                         unlist(list_type2_snr4$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_type2_snr4$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr4_type2$SNR <- 1.25
RMSE_G2_snr5_type2 <- data.frame(value=c(unlist(list_type2_snr5$GRIDY$RMSE$Group2),
                                         unlist(list_type2_snr5$SCA_P$RMSE$Group2),
                                         unlist(list_type2_snr5$GICA$RMSE$Group2),
                                         unlist(list_type2_snr5$DSCA$RMSE$Group2),
                                         unlist(list_type2_snr5$DGICA$RMSE$Group2),
                                         unlist(list_type2_snr5$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_type2_snr5$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr5_type2$SNR <- 2
RMSE_G2_snr6_type2 <- data.frame(value=c(unlist(list_type2_snr6$GRIDY$RMSE$Group2),
                                         unlist(list_type2_snr6$SCA_P$RMSE$Group2),
                                         unlist(list_type2_snr6$GICA$RMSE$Group2),
                                         unlist(list_type2_snr6$DSCA$RMSE$Group2),
                                         unlist(list_type2_snr6$DGICA$RMSE$Group2),
                                         unlist(list_type2_snr6$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_type2_snr6$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr6_type2$SNR <- 4
RMSE_G2_type2 <- rbind(RMSE_G2_snr1_type2,RMSE_G2_snr2_type2,RMSE_G2_snr3_type2,
                       RMSE_G2_snr4_type2,RMSE_G2_snr5_type2,RMSE_G2_snr6_type2)
RMSE_G2_type2$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                  "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

SNR_list <- c("0.25","0.5","0.75","1.25","2","4")
RMSE_J_type2$SNR <- factor(RMSE_J_type2$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
RMSE_G1_type2$SNR <- factor(RMSE_G1_type2$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
RMSE_G2_type2$SNR <- factor(RMSE_G2_type2$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)

method_list <- c("GRIDY","SCA-P","GICA","DSCA","DGICA","Unfitted_SCA-PF2","Unfitted_GICA")
RMSE_J_type2$method <- factor(RMSE_J_type2$method,levels=method_list,labels=method_list)
RMSE_G1_type2$method <- factor(RMSE_G1_type2$method,levels=method_list,labels=method_list)
RMSE_G2_type2$method <- factor(RMSE_G2_type2$method,levels=method_list,labels=method_list)

pl_RMSE_J_type2 <- ggplot(RMSE_J_type2,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name="RMSE (Joint)") +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_RMSE_G1_type2 <- ggplot(RMSE_G1_type2,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name="RMSE (Group1)") +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_RMSE_G2_type2 <- ggplot(RMSE_G2_type2,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name="RMSE (Group2)") +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))


# CC_B:
CC_B_J_snr1_type2 <- data.frame(value=c(unlist(list_type2_snr1$GRIDY$CC_B$Joint),
                                        unlist(list_type2_snr1$SCA_P$CC_B$Joint),
                                        unlist(list_type2_snr1$GICA$CC_B$Joint),
                                        unlist(list_type2_snr1$DSCA$CC_B$Joint),
                                        unlist(list_type2_snr1$DGICA$CC_B$Joint),
                                        unlist(list_type2_snr1$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_type2_snr1$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr1_type2$SNR <- 0.25
CC_B_J_snr2_type2 <- data.frame(value=c(unlist(list_type2_snr2$GRIDY$CC_B$Joint),
                                        unlist(list_type2_snr2$SCA_P$CC_B$Joint),
                                        unlist(list_type2_snr2$GICA$CC_B$Joint),
                                        unlist(list_type2_snr2$DSCA$CC_B$Joint),
                                        unlist(list_type2_snr2$DGICA$CC_B$Joint),
                                        unlist(list_type2_snr2$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_type2_snr2$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr2_type2$SNR <- 0.5
CC_B_J_snr3_type2 <- data.frame(value=c(unlist(list_type2_snr3$GRIDY$CC_B$Joint),
                                        unlist(list_type2_snr3$SCA_P$CC_B$Joint),
                                        unlist(list_type2_snr3$GICA$CC_B$Joint),
                                        unlist(list_type2_snr3$DSCA$CC_B$Joint),
                                        unlist(list_type2_snr3$DGICA$CC_B$Joint),
                                        unlist(list_type2_snr3$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_type2_snr3$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr3_type2$SNR <- 0.75
CC_B_J_snr4_type2 <- data.frame(value=c(unlist(list_type2_snr4$GRIDY$CC_B$Joint),
                                        unlist(list_type2_snr4$SCA_P$CC_B$Joint),
                                        unlist(list_type2_snr4$GICA$CC_B$Joint),
                                        unlist(list_type2_snr4$DSCA$CC_B$Joint),
                                        unlist(list_type2_snr4$DGICA$CC_B$Joint),
                                        unlist(list_type2_snr4$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_type2_snr4$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr4_type2$SNR <- 1.25
CC_B_J_snr5_type2 <- data.frame(value=c(unlist(list_type2_snr5$GRIDY$CC_B$Joint),
                                        unlist(list_type2_snr5$SCA_P$CC_B$Joint),
                                        unlist(list_type2_snr5$GICA$CC_B$Joint),
                                        unlist(list_type2_snr5$DSCA$CC_B$Joint),
                                        unlist(list_type2_snr5$DGICA$CC_B$Joint),
                                        unlist(list_type2_snr5$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_type2_snr5$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr5_type2$SNR <- 2
CC_B_J_snr6_type2 <- data.frame(value=c(unlist(list_type2_snr6$GRIDY$CC_B$Joint),
                                        unlist(list_type2_snr6$SCA_P$CC_B$Joint),
                                        unlist(list_type2_snr6$GICA$CC_B$Joint),
                                        unlist(list_type2_snr6$DSCA$CC_B$Joint),
                                        unlist(list_type2_snr6$DGICA$CC_B$Joint),
                                        unlist(list_type2_snr6$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_type2_snr6$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr6_type2$SNR <- 4
CC_B_J_type2 <- rbind(CC_B_J_snr1_type2,CC_B_J_snr2_type2,CC_B_J_snr3_type2,
                      CC_B_J_snr4_type2,CC_B_J_snr5_type2,CC_B_J_snr6_type2)
CC_B_J_type2$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                 "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

CC_B_G1_snr1_type2 <- data.frame(value=c(unlist(list_type2_snr1$GRIDY$CC_B$Group1),
                                         unlist(list_type2_snr1$SCA_P$CC_B$Group1),
                                         unlist(list_type2_snr1$GICA$CC_B$Group1),
                                         unlist(list_type2_snr1$DSCA$CC_B$Group1),
                                         unlist(list_type2_snr1$DGICA$CC_B$Group1),
                                         unlist(list_type2_snr1$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_type2_snr1$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr1_type2$SNR <- 0.25
CC_B_G1_snr2_type2 <- data.frame(value=c(unlist(list_type2_snr2$GRIDY$CC_B$Group1),
                                         unlist(list_type2_snr2$SCA_P$CC_B$Group1),
                                         unlist(list_type2_snr2$GICA$CC_B$Group1),
                                         unlist(list_type2_snr2$DSCA$CC_B$Group1),
                                         unlist(list_type2_snr2$DGICA$CC_B$Group1),
                                         unlist(list_type2_snr2$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_type2_snr2$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr2_type2$SNR <- 0.5
CC_B_G1_snr3_type2 <- data.frame(value=c(unlist(list_type2_snr3$GRIDY$CC_B$Group1),
                                         unlist(list_type2_snr3$SCA_P$CC_B$Group1),
                                         unlist(list_type2_snr3$GICA$CC_B$Group1),
                                         unlist(list_type2_snr3$DSCA$CC_B$Group1),
                                         unlist(list_type2_snr3$DGICA$CC_B$Group1),
                                         unlist(list_type2_snr3$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_type2_snr3$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr3_type2$SNR <- 0.75
CC_B_G1_snr4_type2 <- data.frame(value=c(unlist(list_type2_snr4$GRIDY$CC_B$Group1),
                                         unlist(list_type2_snr4$SCA_P$CC_B$Group1),
                                         unlist(list_type2_snr4$GICA$CC_B$Group1),
                                         unlist(list_type2_snr4$DSCA$CC_B$Group1),
                                         unlist(list_type2_snr4$DGICA$CC_B$Group1),
                                         unlist(list_type2_snr4$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_type2_snr4$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr4_type2$SNR <- 1.25
CC_B_G1_snr5_type2 <- data.frame(value=c(unlist(list_type2_snr5$GRIDY$CC_B$Group1),
                                         unlist(list_type2_snr5$SCA_P$CC_B$Group1),
                                         unlist(list_type2_snr5$GICA$CC_B$Group1),
                                         unlist(list_type2_snr5$DSCA$CC_B$Group1),
                                         unlist(list_type2_snr5$DGICA$CC_B$Group1),
                                         unlist(list_type2_snr5$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_type2_snr5$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr5_type2$SNR <- 2
CC_B_G1_snr6_type2 <- data.frame(value=c(unlist(list_type2_snr6$GRIDY$CC_B$Group1),
                                         unlist(list_type2_snr6$SCA_P$CC_B$Group1),
                                         unlist(list_type2_snr6$GICA$CC_B$Group1),
                                         unlist(list_type2_snr6$DSCA$CC_B$Group1),
                                         unlist(list_type2_snr6$DGICA$CC_B$Group1),
                                         unlist(list_type2_snr6$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_type2_snr6$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr6_type2$SNR <- 4
CC_B_G1_type2 <- rbind(CC_B_G1_snr1_type2,CC_B_G1_snr2_type2,CC_B_G1_snr3_type2,
                       CC_B_G1_snr4_type2,CC_B_G1_snr5_type2,CC_B_G1_snr6_type2)
CC_B_G1_type2$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                  "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


CC_B_G2_snr1_type2 <- data.frame(value=c(unlist(list_type2_snr1$GRIDY$CC_B$Group2),
                                         unlist(list_type2_snr1$SCA_P$CC_B$Group2),
                                         unlist(list_type2_snr1$GICA$CC_B$Group2),
                                         unlist(list_type2_snr1$DSCA$CC_B$Group2),
                                         unlist(list_type2_snr1$DGICA$CC_B$Group2),
                                         unlist(list_type2_snr1$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_type2_snr1$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr1_type2$SNR <- 0.25
CC_B_G2_snr2_type2 <- data.frame(value=c(unlist(list_type2_snr2$GRIDY$CC_B$Group2),
                                         unlist(list_type2_snr2$SCA_P$CC_B$Group2),
                                         unlist(list_type2_snr2$GICA$CC_B$Group2),
                                         unlist(list_type2_snr2$DSCA$CC_B$Group2),
                                         unlist(list_type2_snr2$DGICA$CC_B$Group2),
                                         unlist(list_type2_snr2$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_type2_snr2$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr2_type2$SNR <- 0.5
CC_B_G2_snr3_type2 <- data.frame(value=c(unlist(list_type2_snr3$GRIDY$CC_B$Group2),
                                         unlist(list_type2_snr3$SCA_P$CC_B$Group2),
                                         unlist(list_type2_snr3$GICA$CC_B$Group2),
                                         unlist(list_type2_snr3$DSCA$CC_B$Group2),
                                         unlist(list_type2_snr3$DGICA$CC_B$Group2),
                                         unlist(list_type2_snr3$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_type2_snr3$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr3_type2$SNR <- 0.75
CC_B_G2_snr4_type2 <- data.frame(value=c(unlist(list_type2_snr4$GRIDY$CC_B$Group2),
                                         unlist(list_type2_snr4$SCA_P$CC_B$Group2),
                                         unlist(list_type2_snr4$GICA$CC_B$Group2),
                                         unlist(list_type2_snr4$DSCA$CC_B$Group2),
                                         unlist(list_type2_snr4$DGICA$CC_B$Group2),
                                         unlist(list_type2_snr4$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_type2_snr4$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr4_type2$SNR <- 1.25
CC_B_G2_snr5_type2 <- data.frame(value=c(unlist(list_type2_snr5$GRIDY$CC_B$Group2),
                                         unlist(list_type2_snr5$SCA_P$CC_B$Group2),
                                         unlist(list_type2_snr5$GICA$CC_B$Group2),
                                         unlist(list_type2_snr5$DSCA$CC_B$Group2),
                                         unlist(list_type2_snr5$DGICA$CC_B$Group2),
                                         unlist(list_type2_snr5$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_type2_snr5$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr5_type2$SNR <- 2
CC_B_G2_snr6_type2 <- data.frame(value=c(unlist(list_type2_snr6$GRIDY$CC_B$Group2),
                                         unlist(list_type2_snr6$SCA_P$CC_B$Group2),
                                         unlist(list_type2_snr6$GICA$CC_B$Group2),
                                         unlist(list_type2_snr6$DSCA$CC_B$Group2),
                                         unlist(list_type2_snr6$DGICA$CC_B$Group2),
                                         unlist(list_type2_snr6$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_type2_snr6$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr6_type2$SNR <- 4
CC_B_G2_type2 <- rbind(CC_B_G2_snr1_type2,CC_B_G2_snr2_type2,CC_B_G2_snr3_type2,
                       CC_B_G2_snr4_type2,CC_B_G2_snr5_type2,CC_B_G2_snr6_type2)
CC_B_G2_type2$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                  "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

SNR_list <- c("0.25","0.5","0.75","1.25","2","4")
CC_B_J_type2$SNR <- factor(CC_B_J_type2$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
CC_B_G1_type2$SNR <- factor(CC_B_G1_type2$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
CC_B_G2_type2$SNR <- factor(CC_B_G2_type2$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)

method_list <- c("GRIDY","SCA-P","GICA","DSCA","DGICA","Unfitted_SCA-PF2","Unfitted_GICA")
CC_B_J_type2$method <- factor(CC_B_J_type2$method,levels=method_list,labels=method_list)
CC_B_G1_type2$method <- factor(CC_B_G1_type2$method,levels=method_list,labels=method_list)
CC_B_G2_type2$method <- factor(CC_B_G2_type2$method,levels=method_list,labels=method_list)

pl_CC_B_J_type2 <- ggplot(CC_B_J_type2,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_B$ (Joint)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_B_G1_type2 <- ggplot(CC_B_G1_type2,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_B$ (Group1)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_B_G2_type2 <- ggplot(CC_B_G2_type2,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_B$ (Group2)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))


# CC_F:
CC_F_J_snr1_type2 <- data.frame(value=c(unlist(list_type2_snr1$GRIDY$CC_F$Joint),
                                        unlist(list_type2_snr1$SCA_P$CC_F$Joint),
                                        unlist(list_type2_snr1$GICA$CC_F$Joint),
                                        unlist(list_type2_snr1$DSCA$CC_F$Joint),
                                        unlist(list_type2_snr1$DGICA$CC_F$Joint),
                                        unlist(list_type2_snr1$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_type2_snr1$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr1_type2$SNR <- 0.25
CC_F_J_snr2_type2 <- data.frame(value=c(unlist(list_type2_snr2$GRIDY$CC_F$Joint),
                                        unlist(list_type2_snr2$SCA_P$CC_F$Joint),
                                        unlist(list_type2_snr2$GICA$CC_F$Joint),
                                        unlist(list_type2_snr2$DSCA$CC_F$Joint),
                                        unlist(list_type2_snr2$DGICA$CC_F$Joint),
                                        unlist(list_type2_snr2$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_type2_snr2$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr2_type2$SNR <- 0.5
CC_F_J_snr3_type2 <- data.frame(value=c(unlist(list_type2_snr3$GRIDY$CC_F$Joint),
                                        unlist(list_type2_snr3$SCA_P$CC_F$Joint),
                                        unlist(list_type2_snr3$GICA$CC_F$Joint),
                                        unlist(list_type2_snr3$DSCA$CC_F$Joint),
                                        unlist(list_type2_snr3$DGICA$CC_F$Joint),
                                        unlist(list_type2_snr3$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_type2_snr3$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr3_type2$SNR <- 0.75
CC_F_J_snr4_type2 <- data.frame(value=c(unlist(list_type2_snr4$GRIDY$CC_F$Joint),
                                        unlist(list_type2_snr4$SCA_P$CC_F$Joint),
                                        unlist(list_type2_snr4$GICA$CC_F$Joint),
                                        unlist(list_type2_snr4$DSCA$CC_F$Joint),
                                        unlist(list_type2_snr4$DGICA$CC_F$Joint),
                                        unlist(list_type2_snr4$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_type2_snr4$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr4_type2$SNR <- 1.25
CC_F_J_snr5_type2 <- data.frame(value=c(unlist(list_type2_snr5$GRIDY$CC_F$Joint),
                                        unlist(list_type2_snr5$SCA_P$CC_F$Joint),
                                        unlist(list_type2_snr5$GICA$CC_F$Joint),
                                        unlist(list_type2_snr5$DSCA$CC_F$Joint),
                                        unlist(list_type2_snr5$DGICA$CC_F$Joint),
                                        unlist(list_type2_snr5$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_type2_snr5$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr5_type2$SNR <- 2
CC_F_J_snr6_type2 <- data.frame(value=c(unlist(list_type2_snr6$GRIDY$CC_F$Joint),
                                        unlist(list_type2_snr6$SCA_P$CC_F$Joint),
                                        unlist(list_type2_snr6$GICA$CC_F$Joint),
                                        unlist(list_type2_snr6$DSCA$CC_F$Joint),
                                        unlist(list_type2_snr6$DGICA$CC_F$Joint),
                                        unlist(list_type2_snr6$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_type2_snr6$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr6_type2$SNR <- 4
CC_F_J_type2 <- rbind(CC_F_J_snr1_type2,CC_F_J_snr2_type2,CC_F_J_snr3_type2,
                      CC_F_J_snr4_type2,CC_F_J_snr5_type2,CC_F_J_snr6_type2)
CC_F_J_type2$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                 "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

CC_F_G1_snr1_type2 <- data.frame(value=c(unlist(list_type2_snr1$GRIDY$CC_F$Group1),
                                         unlist(list_type2_snr1$SCA_P$CC_F$Group1),
                                         unlist(list_type2_snr1$GICA$CC_F$Group1),
                                         unlist(list_type2_snr1$DSCA$CC_F$Group1),
                                         unlist(list_type2_snr1$DGICA$CC_F$Group1),
                                         unlist(list_type2_snr1$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_type2_snr1$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr1_type2$SNR <- 0.25
CC_F_G1_snr2_type2 <- data.frame(value=c(unlist(list_type2_snr2$GRIDY$CC_F$Group1),
                                         unlist(list_type2_snr2$SCA_P$CC_F$Group1),
                                         unlist(list_type2_snr2$GICA$CC_F$Group1),
                                         unlist(list_type2_snr2$DSCA$CC_F$Group1),
                                         unlist(list_type2_snr2$DGICA$CC_F$Group1),
                                         unlist(list_type2_snr2$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_type2_snr2$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr2_type2$SNR <- 0.5
CC_F_G1_snr3_type2 <- data.frame(value=c(unlist(list_type2_snr3$GRIDY$CC_F$Group1),
                                         unlist(list_type2_snr3$SCA_P$CC_F$Group1),
                                         unlist(list_type2_snr3$GICA$CC_F$Group1),
                                         unlist(list_type2_snr3$DSCA$CC_F$Group1),
                                         unlist(list_type2_snr3$DGICA$CC_F$Group1),
                                         unlist(list_type2_snr3$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_type2_snr3$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr3_type2$SNR <- 0.75
CC_F_G1_snr4_type2 <- data.frame(value=c(unlist(list_type2_snr4$GRIDY$CC_F$Group1),
                                         unlist(list_type2_snr4$SCA_P$CC_F$Group1),
                                         unlist(list_type2_snr4$GICA$CC_F$Group1),
                                         unlist(list_type2_snr4$DSCA$CC_F$Group1),
                                         unlist(list_type2_snr4$DGICA$CC_F$Group1),
                                         unlist(list_type2_snr4$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_type2_snr4$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr4_type2$SNR <- 1.25
CC_F_G1_snr5_type2 <- data.frame(value=c(unlist(list_type2_snr5$GRIDY$CC_F$Group1),
                                         unlist(list_type2_snr5$SCA_P$CC_F$Group1),
                                         unlist(list_type2_snr5$GICA$CC_F$Group1),
                                         unlist(list_type2_snr5$DSCA$CC_F$Group1),
                                         unlist(list_type2_snr5$DGICA$CC_F$Group1),
                                         unlist(list_type2_snr5$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_type2_snr5$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr5_type2$SNR <- 2
CC_F_G1_snr6_type2 <- data.frame(value=c(unlist(list_type2_snr6$GRIDY$CC_F$Group1),
                                         unlist(list_type2_snr6$SCA_P$CC_F$Group1),
                                         unlist(list_type2_snr6$GICA$CC_F$Group1),
                                         unlist(list_type2_snr6$DSCA$CC_F$Group1),
                                         unlist(list_type2_snr6$DGICA$CC_F$Group1),
                                         unlist(list_type2_snr6$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_type2_snr6$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr6_type2$SNR <- 4
CC_F_G1_type2 <- rbind(CC_F_G1_snr1_type2,CC_F_G1_snr2_type2,CC_F_G1_snr3_type2,
                       CC_F_G1_snr4_type2,CC_F_G1_snr5_type2,CC_F_G1_snr6_type2)
CC_F_G1_type2$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                  "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


CC_F_G2_snr1_type2 <- data.frame(value=c(unlist(list_type2_snr1$GRIDY$CC_F$Group2),
                                         unlist(list_type2_snr1$SCA_P$CC_F$Group2),
                                         unlist(list_type2_snr1$GICA$CC_F$Group2),
                                         unlist(list_type2_snr1$DSCA$CC_F$Group2),
                                         unlist(list_type2_snr1$DGICA$CC_F$Group2),
                                         unlist(list_type2_snr1$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_type2_snr1$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr1_type2$SNR <- 0.25
CC_F_G2_snr2_type2 <- data.frame(value=c(unlist(list_type2_snr2$GRIDY$CC_F$Group2),
                                         unlist(list_type2_snr2$SCA_P$CC_F$Group2),
                                         unlist(list_type2_snr2$GICA$CC_F$Group2),
                                         unlist(list_type2_snr2$DSCA$CC_F$Group2),
                                         unlist(list_type2_snr2$DGICA$CC_F$Group2),
                                         unlist(list_type2_snr2$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_type2_snr2$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr2_type2$SNR <- 0.5
CC_F_G2_snr3_type2 <- data.frame(value=c(unlist(list_type2_snr3$GRIDY$CC_F$Group2),
                                         unlist(list_type2_snr3$SCA_P$CC_F$Group2),
                                         unlist(list_type2_snr3$GICA$CC_F$Group2),
                                         unlist(list_type2_snr3$DSCA$CC_F$Group2),
                                         unlist(list_type2_snr3$DGICA$CC_F$Group2),
                                         unlist(list_type2_snr3$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_type2_snr3$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr3_type2$SNR <- 0.75
CC_F_G2_snr4_type2 <- data.frame(value=c(unlist(list_type2_snr4$GRIDY$CC_F$Group2),
                                         unlist(list_type2_snr4$SCA_P$CC_F$Group2),
                                         unlist(list_type2_snr4$GICA$CC_F$Group2),
                                         unlist(list_type2_snr4$DSCA$CC_F$Group2),
                                         unlist(list_type2_snr4$DGICA$CC_F$Group2),
                                         unlist(list_type2_snr4$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_type2_snr4$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr4_type2$SNR <- 1.25
CC_F_G2_snr5_type2 <- data.frame(value=c(unlist(list_type2_snr5$GRIDY$CC_F$Group2),
                                         unlist(list_type2_snr5$SCA_P$CC_F$Group2),
                                         unlist(list_type2_snr5$GICA$CC_F$Group2),
                                         unlist(list_type2_snr5$DSCA$CC_F$Group2),
                                         unlist(list_type2_snr5$DGICA$CC_F$Group2),
                                         unlist(list_type2_snr5$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_type2_snr5$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr5_type2$SNR <- 2
CC_F_G2_snr6_type2 <- data.frame(value=c(unlist(list_type2_snr6$GRIDY$CC_F$Group2),
                                         unlist(list_type2_snr6$SCA_P$CC_F$Group2),
                                         unlist(list_type2_snr6$GICA$CC_F$Group2),
                                         unlist(list_type2_snr6$DSCA$CC_F$Group2),
                                         unlist(list_type2_snr6$DGICA$CC_F$Group2),
                                         unlist(list_type2_snr6$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_type2_snr6$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr6_type2$SNR <- 4
CC_F_G2_type2 <- rbind(CC_F_G2_snr1_type2,CC_F_G2_snr2_type2,CC_F_G2_snr3_type2,
                       CC_F_G2_snr4_type2,CC_F_G2_snr5_type2,CC_F_G2_snr6_type2)
CC_F_G2_type2$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                  "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

SNR_list <- c("0.25","0.5","0.75","1.25","2","4")
CC_F_J_type2$SNR <- factor(CC_F_J_type2$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),SNR_list)
CC_F_G1_type2$SNR <- factor(CC_F_G1_type2$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),SNR_list)
CC_F_G2_type2$SNR <- factor(CC_F_G2_type2$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),SNR_list)

method_list <- c("GRIDY","SCA-P","GICA","DSCA","DGICA","Unfitted_SCA-PF2","Unfitted_GICA")
CC_F_J_type2$method <- factor(CC_F_J_type2$method,levels=method_list,labels=method_list)
CC_F_G1_type2$method <- factor(CC_F_G1_type2$method,levels=method_list,labels=method_list)
CC_F_G2_type2$method <- factor(CC_F_G2_type2$method,levels=method_list,labels=method_list)

pl_CC_F_J_type2 <- ggplot(CC_F_J_type2,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_F$ (Joint)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_F_G1_type2 <- ggplot(CC_F_G1_type2,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_F$ (Group1)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_F_G2_type2 <- ggplot(CC_F_G2_type2,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_F$ (Group2)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))



# CC_PSI:
CC_PSI_J_snr1_type2 <- data.frame(value=c(unlist(list_type2_snr1$GRIDY$CC_PSI$Joint),
                                          unlist(list_type2_snr1$SCA_P$CC_PSI$Joint),
                                          unlist(list_type2_snr1$GICA$CC_PSI$Joint),
                                          unlist(list_type2_snr1$DSCA$CC_PSI$Joint),
                                          unlist(list_type2_snr1$DGICA$CC_PSI$Joint),
                                          unlist(list_type2_snr1$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_type2_snr1$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr1_type2$SNR <- 0.25
CC_PSI_J_snr2_type2 <- data.frame(value=c(unlist(list_type2_snr2$GRIDY$CC_PSI$Joint),
                                          unlist(list_type2_snr2$SCA_P$CC_PSI$Joint),
                                          unlist(list_type2_snr2$GICA$CC_PSI$Joint),
                                          unlist(list_type2_snr2$DSCA$CC_PSI$Joint),
                                          unlist(list_type2_snr2$DGICA$CC_PSI$Joint),
                                          unlist(list_type2_snr2$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_type2_snr2$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr2_type2$SNR <- 0.5
CC_PSI_J_snr3_type2 <- data.frame(value=c(unlist(list_type2_snr3$GRIDY$CC_PSI$Joint),
                                          unlist(list_type2_snr3$SCA_P$CC_PSI$Joint),
                                          unlist(list_type2_snr3$GICA$CC_PSI$Joint),
                                          unlist(list_type2_snr3$DSCA$CC_PSI$Joint),
                                          unlist(list_type2_snr3$DGICA$CC_PSI$Joint),
                                          unlist(list_type2_snr3$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_type2_snr3$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr3_type2$SNR <- 0.75
CC_PSI_J_snr4_type2 <- data.frame(value=c(unlist(list_type2_snr4$GRIDY$CC_PSI$Joint),
                                          unlist(list_type2_snr4$SCA_P$CC_PSI$Joint),
                                          unlist(list_type2_snr4$GICA$CC_PSI$Joint),
                                          unlist(list_type2_snr4$DSCA$CC_PSI$Joint),
                                          unlist(list_type2_snr4$DGICA$CC_PSI$Joint),
                                          unlist(list_type2_snr4$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_type2_snr4$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr4_type2$SNR <- 1.25
CC_PSI_J_snr5_type2 <- data.frame(value=c(unlist(list_type2_snr5$GRIDY$CC_PSI$Joint),
                                          unlist(list_type2_snr5$SCA_P$CC_PSI$Joint),
                                          unlist(list_type2_snr5$GICA$CC_PSI$Joint),
                                          unlist(list_type2_snr5$DSCA$CC_PSI$Joint),
                                          unlist(list_type2_snr5$DGICA$CC_PSI$Joint),
                                          unlist(list_type2_snr5$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_type2_snr5$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr5_type2$SNR <- 2
CC_PSI_J_snr6_type2 <- data.frame(value=c(unlist(list_type2_snr6$GRIDY$CC_PSI$Joint),
                                          unlist(list_type2_snr6$SCA_P$CC_PSI$Joint),
                                          unlist(list_type2_snr6$GICA$CC_PSI$Joint),
                                          unlist(list_type2_snr6$DSCA$CC_PSI$Joint),
                                          unlist(list_type2_snr6$DGICA$CC_PSI$Joint),
                                          unlist(list_type2_snr6$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_type2_snr6$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr6_type2$SNR <- 4
CC_PSI_J_type2 <- rbind(CC_PSI_J_snr1_type2,CC_PSI_J_snr2_type2,CC_PSI_J_snr3_type2,
                        CC_PSI_J_snr4_type2,CC_PSI_J_snr5_type2,CC_PSI_J_snr6_type2)
CC_PSI_J_type2$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                   "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

CC_PSI_G1_snr1_type2 <- data.frame(value=c(unlist(list_type2_snr1$GRIDY$CC_PSI$Group1),
                                           unlist(list_type2_snr1$SCA_P$CC_PSI$Group1),
                                           unlist(list_type2_snr1$GICA$CC_PSI$Group1),
                                           unlist(list_type2_snr1$DSCA$CC_PSI$Group1),
                                           unlist(list_type2_snr1$DGICA$CC_PSI$Group1),
                                           unlist(list_type2_snr1$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_type2_snr1$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr1_type2$SNR <- 0.25
CC_PSI_G1_snr2_type2 <- data.frame(value=c(unlist(list_type2_snr2$GRIDY$CC_PSI$Group1),
                                           unlist(list_type2_snr2$SCA_P$CC_PSI$Group1),
                                           unlist(list_type2_snr2$GICA$CC_PSI$Group1),
                                           unlist(list_type2_snr2$DSCA$CC_PSI$Group1),
                                           unlist(list_type2_snr2$DGICA$CC_PSI$Group1),
                                           unlist(list_type2_snr2$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_type2_snr2$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr2_type2$SNR <- 0.5
CC_PSI_G1_snr3_type2 <- data.frame(value=c(unlist(list_type2_snr3$GRIDY$CC_PSI$Group1),
                                           unlist(list_type2_snr3$SCA_P$CC_PSI$Group1),
                                           unlist(list_type2_snr3$GICA$CC_PSI$Group1),
                                           unlist(list_type2_snr3$DSCA$CC_PSI$Group1),
                                           unlist(list_type2_snr3$DGICA$CC_PSI$Group1),
                                           unlist(list_type2_snr3$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_type2_snr3$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr3_type2$SNR <- 0.75
CC_PSI_G1_snr4_type2 <- data.frame(value=c(unlist(list_type2_snr4$GRIDY$CC_PSI$Group1),
                                           unlist(list_type2_snr4$SCA_P$CC_PSI$Group1),
                                           unlist(list_type2_snr4$GICA$CC_PSI$Group1),
                                           unlist(list_type2_snr4$DSCA$CC_PSI$Group1),
                                           unlist(list_type2_snr4$DGICA$CC_PSI$Group1),
                                           unlist(list_type2_snr4$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_type2_snr4$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr4_type2$SNR <- 1.25
CC_PSI_G1_snr5_type2 <- data.frame(value=c(unlist(list_type2_snr5$GRIDY$CC_PSI$Group1),
                                           unlist(list_type2_snr5$SCA_P$CC_PSI$Group1),
                                           unlist(list_type2_snr5$GICA$CC_PSI$Group1),
                                           unlist(list_type2_snr5$DSCA$CC_PSI$Group1),
                                           unlist(list_type2_snr5$DGICA$CC_PSI$Group1),
                                           unlist(list_type2_snr5$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_type2_snr5$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr5_type2$SNR <- 2
CC_PSI_G1_snr6_type2 <- data.frame(value=c(unlist(list_type2_snr6$GRIDY$CC_PSI$Group1),
                                           unlist(list_type2_snr6$SCA_P$CC_PSI$Group1),
                                           unlist(list_type2_snr6$GICA$CC_PSI$Group1),
                                           unlist(list_type2_snr6$DSCA$CC_PSI$Group1),
                                           unlist(list_type2_snr6$DGICA$CC_PSI$Group1),
                                           unlist(list_type2_snr6$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_type2_snr6$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr6_type2$SNR <- 4
CC_PSI_G1_type2 <- rbind(CC_PSI_G1_snr1_type2,CC_PSI_G1_snr2_type2,CC_PSI_G1_snr3_type2,
                         CC_PSI_G1_snr4_type2,CC_PSI_G1_snr5_type2,CC_PSI_G1_snr6_type2)
CC_PSI_G1_type2$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                    "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


CC_PSI_G2_snr1_type2 <- data.frame(value=c(unlist(list_type2_snr1$GRIDY$CC_PSI$Group2),
                                           unlist(list_type2_snr1$SCA_P$CC_PSI$Group2),
                                           unlist(list_type2_snr1$GICA$CC_PSI$Group2),
                                           unlist(list_type2_snr1$DSCA$CC_PSI$Group2),
                                           unlist(list_type2_snr1$DGICA$CC_PSI$Group2),
                                           unlist(list_type2_snr1$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_type2_snr1$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr1_type2$SNR <- 0.25
CC_PSI_G2_snr2_type2 <- data.frame(value=c(unlist(list_type2_snr2$GRIDY$CC_PSI$Group2),
                                           unlist(list_type2_snr2$SCA_P$CC_PSI$Group2),
                                           unlist(list_type2_snr2$GICA$CC_PSI$Group2),
                                           unlist(list_type2_snr2$DSCA$CC_PSI$Group2),
                                           unlist(list_type2_snr2$DGICA$CC_PSI$Group2),
                                           unlist(list_type2_snr2$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_type2_snr2$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr2_type2$SNR <- 0.5
CC_PSI_G2_snr3_type2 <- data.frame(value=c(unlist(list_type2_snr3$GRIDY$CC_PSI$Group2),
                                           unlist(list_type2_snr3$SCA_P$CC_PSI$Group2),
                                           unlist(list_type2_snr3$GICA$CC_PSI$Group2),
                                           unlist(list_type2_snr3$DSCA$CC_PSI$Group2),
                                           unlist(list_type2_snr3$DGICA$CC_PSI$Group2),
                                           unlist(list_type2_snr3$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_type2_snr3$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr3_type2$SNR <- 0.75
CC_PSI_G2_snr4_type2 <- data.frame(value=c(unlist(list_type2_snr4$GRIDY$CC_PSI$Group2),
                                           unlist(list_type2_snr4$SCA_P$CC_PSI$Group2),
                                           unlist(list_type2_snr4$GICA$CC_PSI$Group2),
                                           unlist(list_type2_snr4$DSCA$CC_PSI$Group2),
                                           unlist(list_type2_snr4$DGICA$CC_PSI$Group2),
                                           unlist(list_type2_snr4$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_type2_snr4$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr4_type2$SNR <- 1.25
CC_PSI_G2_snr5_type2 <- data.frame(value=c(unlist(list_type2_snr5$GRIDY$CC_PSI$Group2),
                                           unlist(list_type2_snr5$SCA_P$CC_PSI$Group2),
                                           unlist(list_type2_snr5$GICA$CC_PSI$Group2),
                                           unlist(list_type2_snr5$DSCA$CC_PSI$Group2),
                                           unlist(list_type2_snr5$DGICA$CC_PSI$Group2),
                                           unlist(list_type2_snr5$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_type2_snr5$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr5_type2$SNR <- 2
CC_PSI_G2_snr6_type2 <- data.frame(value=c(unlist(list_type2_snr6$GRIDY$CC_PSI$Group2),
                                           unlist(list_type2_snr6$SCA_P$CC_PSI$Group2),
                                           unlist(list_type2_snr6$GICA$CC_PSI$Group2),
                                           unlist(list_type2_snr6$DSCA$CC_PSI$Group2),
                                           unlist(list_type2_snr6$DGICA$CC_PSI$Group2),
                                           unlist(list_type2_snr6$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_type2_snr6$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr6_type2$SNR <- 4
CC_PSI_G2_type2 <- rbind(CC_PSI_G2_snr1_type2,CC_PSI_G2_snr2_type2,CC_PSI_G2_snr3_type2,
                         CC_PSI_G2_snr4_type2,CC_PSI_G2_snr5_type2,CC_PSI_G2_snr6_type2)
CC_PSI_G2_type2$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                    "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

SNR_list <- c("0.25","0.5","0.75","1.25","2","4")
CC_PSI_J_type2$SNR <- factor(CC_PSI_J_type2$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
CC_PSI_G1_type2$SNR <- factor(CC_PSI_G1_type2$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)
CC_PSI_G2_type2$SNR <- factor(CC_PSI_G2_type2$SNR,levels=c(0.25,0.5,0.75,1.25,2,4),labels=SNR_list)

method_list <- c("GRIDY","SCA-P","GICA","DSCA","DGICA","Unfitted_SCA-PF2","Unfitted_GICA")
CC_PSI_J_type2$method <- factor(CC_PSI_J_type2$method,levels=method_list,labels=method_list)
CC_PSI_G1_type2$method <- factor(CC_PSI_G1_type2$method,levels=method_list,labels=method_list)
CC_PSI_G2_type2$method <- factor(CC_PSI_G2_type2$method,levels=method_list,labels=method_list)

pl_CC_PSI_J_type2 <- ggplot(CC_PSI_J_type2,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_{\\Psi}$ (Joint)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_PSI_G1_type2 <- ggplot(CC_PSI_G1_type2,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_{\\Psi}$ (Group1)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_PSI_G2_type2 <- ggplot(CC_PSI_G2_type2,aes(x=SNR,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="SNR",labels=c("0.25","0.5","0.75","1.25","2","4")) +
  scale_y_continuous(name=TeX("$CC_{\\Psi}$ (Group2)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))


# Plotting:
mylegend_type2 <- g_legend(pl_R2_J_type2)
grid.arrange(arrangeGrob(pl_R2_J_type2 + theme(legend.position="none"),
                         pl_R2_G1_type2 + theme(legend.position="none"), 
                         pl_R2_G2_type2 + theme(legend.position="none"),
                         pl_RMSE_J_type2 + theme(legend.position="none"),
                         pl_RMSE_G1_type2 + theme(legend.position="none"), 
                         pl_RMSE_G2_type2 + theme(legend.position="none"),
                         pl_CC_B_J_type2 + theme(legend.position="none"),
                         pl_CC_B_G1_type2 + theme(legend.position="none"), 
                         pl_CC_B_G2_type2 + theme(legend.position="none"),
                         pl_CC_F_J_type2 + theme(legend.position="none"),
                         pl_CC_F_G1_type2 + theme(legend.position="none"), 
                         pl_CC_F_G2_type2 + theme(legend.position="none"),
                         pl_CC_PSI_J_type2 + theme(legend.position="none"),
                         pl_CC_PSI_G1_type2 + theme(legend.position="none"), 
                         pl_CC_PSI_G2_type2 + theme(legend.position="none"),
                         nrow=5,ncol=3), mylegend_type2,
             nrow = 2, heights = c(5,0.2))
