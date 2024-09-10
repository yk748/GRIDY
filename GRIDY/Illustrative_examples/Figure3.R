#-----------------------------------------------------------------------------#
#   File name : Figure3.R    												  
#
#   Project : "Group integrative dynamic factor models 
#             with application to multiple subject brain connectivity"
#
#   Maintainer : Younghoon Kim                     
#
#   Date : Sep. 1st, 2024
#
#   Purpose : Code for producing results for plotting Figure 3 in Illustrative example 
#             section
#
#   R version 4.0.5 (2021-03-31)                                     
#
#   Input data file : /Illustrative_examples/result/sim2_snr-_type1.RData
# 
#   Output data file : /Illustrative_examples/figures/Figure3.pdf
#
#   Required R packages : ggplot2_3.4.0, latex2exp_0.9.6, ggpubr_0.5.0, RColorBrewer_1.1-3,
#                         reshape2_1.4.4, and gridExtra_2.3
#-----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#
# Loading simulated data
# -----------------------------------------------------------------------------#
source(paste0(dirname(getwd()),"/","library_simulation.R"))

# Defining function for a legend 
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

load("./result/sim2_snr1_type1.RData"); 
list_snr_1_type1 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./result/sim2_snr2_type1.RData"); 
list_snr_2_type1 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./result/sim2_snr3_type1.RData"); 
list_snr_3_type1 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./result/sim2_snr4_type1.RData"); 
list_snr_4_type1 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./result/sim2_snr5_type1.RData"); 
list_snr_5_type1 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./result/sim2_snr6_type1.RData"); 
list_snr_6_type1 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                         DSCA=list_DSCA,DGICA=list_DGICA,
                         Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)


#-----------------------------------------------------------------------------#
# R2
#-----------------------------------------------------------------------------#
R2_J_snr1_type1 <- data.frame(value=c(unlist(list_snr_1_type1$GRIDY$R2$Joint),
                                      unlist(list_snr_1_type1$SCA_P$R2$Joint),
                                      unlist(list_snr_1_type1$GICA$R2$Joint),
                                      unlist(list_snr_1_type1$DSCA$R2$Joint),
                                      unlist(list_snr_1_type1$DGICA$R2$Joint),
                                      unlist(list_snr_1_type1$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_snr_1_type1$Unfitted_GICA$R2$Joint)))
R2_J_snr1_type1$SNR <- 0.25
R2_J_snr2_type1 <- data.frame(value=c(unlist(list_snr_2_type1$GRIDY$R2$Joint),
                                      unlist(list_snr_2_type1$SCA_P$R2$Joint),
                                      unlist(list_snr_2_type1$GICA$R2$Joint),
                                      unlist(list_snr_2_type1$DSCA$R2$Joint),
                                      unlist(list_snr_2_type1$DGICA$R2$Joint),
                                      unlist(list_snr_2_type1$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_snr_2_type1$Unfitted_GICA$R2$Joint)))
R2_J_snr2_type1$SNR <- 0.5
R2_J_snr3_type1 <- data.frame(value=c(unlist(list_snr_3_type1$GRIDY$R2$Joint),
                                      unlist(list_snr_3_type1$SCA_P$R2$Joint),
                                      unlist(list_snr_3_type1$GICA$R2$Joint),
                                      unlist(list_snr_3_type1$DSCA$R2$Joint),
                                      unlist(list_snr_3_type1$DGICA$R2$Joint),
                                      unlist(list_snr_3_type1$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_snr_3_type1$Unfitted_GICA$R2$Joint)))
R2_J_snr3_type1$SNR <- 0.75
R2_J_snr4_type1 <- data.frame(value=c(unlist(list_snr_4_type1$GRIDY$R2$Joint),
                                      unlist(list_snr_4_type1$SCA_P$R2$Joint),
                                      unlist(list_snr_4_type1$GICA$R2$Joint),
                                      unlist(list_snr_4_type1$DSCA$R2$Joint),
                                      unlist(list_snr_4_type1$DGICA$R2$Joint),
                                      unlist(list_snr_4_type1$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_snr_4_type1$Unfitted_GICA$R2$Joint)))
R2_J_snr4_type1$SNR <- 1.25
R2_J_snr5_type1 <- data.frame(value=c(unlist(list_snr_5_type1$GRIDY$R2$Joint),
                                      unlist(list_snr_5_type1$SCA_P$R2$Joint),
                                      unlist(list_snr_5_type1$GICA$R2$Joint),
                                      unlist(list_snr_5_type1$DSCA$R2$Joint),
                                      unlist(list_snr_5_type1$DGICA$R2$Joint),
                                      unlist(list_snr_5_type1$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_snr_5_type1$Unfitted_GICA$R2$Joint)))
R2_J_snr5_type1$SNR <- 2
R2_J_snr6_type1 <- data.frame(value=c(unlist(list_snr_6_type1$GRIDY$R2$Joint),
                                      unlist(list_snr_6_type1$SCA_P$R2$Joint),
                                      unlist(list_snr_6_type1$GICA$R2$Joint),
                                      unlist(list_snr_6_type1$DSCA$R2$Joint),
                                      unlist(list_snr_6_type1$DGICA$R2$Joint),
                                      unlist(list_snr_6_type1$Unfitted_SCA_PF2$R2$Joint),
                                      unlist(list_snr_6_type1$Unfitted_GICA$R2$Joint)))
R2_J_snr6_type1$SNR <- 4
R2_J_type1 <- rbind(R2_J_snr1_type1,R2_J_snr2_type1,R2_J_snr3_type1,
                    R2_J_snr4_type1,R2_J_snr5_type1,R2_J_snr6_type1)
R2_J_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                               "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

R2_G1_snr1_type1 <- data.frame(value=c(unlist(list_snr_1_type1$GRIDY$R2$Group1),
                                       unlist(list_snr_1_type1$SCA_P$R2$Group1),
                                       unlist(list_snr_1_type1$GICA$R2$Group1),
                                       unlist(list_snr_1_type1$DSCA$R2$Group1),
                                       unlist(list_snr_1_type1$DGICA$R2$Group1),
                                       unlist(list_snr_1_type1$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_snr_1_type1$Unfitted_GICA$R2$Group1)))
R2_G1_snr1_type1$SNR <- 0.25
R2_G1_snr2_type1 <- data.frame(value=c(unlist(list_snr_2_type1$GRIDY$R2$Group1),
                                       unlist(list_snr_2_type1$SCA_P$R2$Group1),
                                       unlist(list_snr_2_type1$GICA$R2$Group1),
                                       unlist(list_snr_2_type1$DSCA$R2$Group1),
                                       unlist(list_snr_2_type1$DGICA$R2$Group1),
                                       unlist(list_snr_2_type1$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_snr_2_type1$Unfitted_GICA$R2$Group1)))
R2_G1_snr2_type1$SNR <- 0.5
R2_G1_snr3_type1 <- data.frame(value=c(unlist(list_snr_3_type1$GRIDY$R2$Group1),
                                       unlist(list_snr_3_type1$SCA_P$R2$Group1),
                                       unlist(list_snr_3_type1$GICA$R2$Group1),
                                       unlist(list_snr_3_type1$DSCA$R2$Group1),
                                       unlist(list_snr_3_type1$DGICA$R2$Group1),
                                       unlist(list_snr_3_type1$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_snr_3_type1$Unfitted_GICA$R2$Group1)))
R2_G1_snr3_type1$SNR <- 0.75
R2_G1_snr4_type1 <- data.frame(value=c(unlist(list_snr_4_type1$GRIDY$R2$Group1),
                                       unlist(list_snr_4_type1$SCA_P$R2$Group1),
                                       unlist(list_snr_4_type1$GICA$R2$Group1),
                                       unlist(list_snr_4_type1$DSCA$R2$Group1),
                                       unlist(list_snr_4_type1$DGICA$R2$Group1),
                                       unlist(list_snr_4_type1$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_snr_4_type1$Unfitted_GICA$R2$Group1)))
R2_G1_snr4_type1$SNR <- 1.25
R2_G1_snr5_type1 <- data.frame(value=c(unlist(list_snr_5_type1$GRIDY$R2$Group1),
                                       unlist(list_snr_5_type1$SCA_P$R2$Group1),
                                       unlist(list_snr_5_type1$GICA$R2$Group1),
                                       unlist(list_snr_5_type1$DSCA$R2$Group1),
                                       unlist(list_snr_5_type1$DGICA$R2$Group1),
                                       unlist(list_snr_5_type1$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_snr_5_type1$Unfitted_GICA$R2$Group1)))
R2_G1_snr5_type1$SNR <- 2
R2_G1_snr6_type1 <- data.frame(value=c(unlist(list_snr_6_type1$GRIDY$R2$Group1),
                                       unlist(list_snr_6_type1$SCA_P$R2$Group1),
                                       unlist(list_snr_6_type1$GICA$R2$Group1),
                                       unlist(list_snr_6_type1$DSCA$R2$Group1),
                                       unlist(list_snr_6_type1$DGICA$R2$Group1),
                                       unlist(list_snr_6_type1$Unfitted_SCA_PF2$R2$Group1),
                                       unlist(list_snr_6_type1$Unfitted_GICA$R2$Group1)))
R2_G1_snr6_type1$SNR <- 4
R2_G1_type1 <- rbind(R2_G1_snr1_type1,R2_G1_snr2_type1,R2_G1_snr3_type1,
                     R2_G1_snr4_type1,R2_G1_snr5_type1,R2_G1_snr6_type1)
R2_G1_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


R2_G2_snr1_type1 <- data.frame(value=c(unlist(list_snr_1_type1$GRIDY$R2$Group2),
                                       unlist(list_snr_1_type1$SCA_P$R2$Group2),
                                       unlist(list_snr_1_type1$GICA$R2$Group2),
                                       unlist(list_snr_1_type1$DSCA$R2$Group2),
                                       unlist(list_snr_1_type1$DGICA$R2$Group2),
                                       unlist(list_snr_1_type1$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_snr_1_type1$Unfitted_GICA$R2$Group2)))
R2_G2_snr1_type1$SNR <- 0.25
R2_G2_snr2_type1 <- data.frame(value=c(unlist(list_snr_2_type1$GRIDY$R2$Group2),
                                       unlist(list_snr_2_type1$SCA_P$R2$Group2),
                                       unlist(list_snr_2_type1$GICA$R2$Group2),
                                       unlist(list_snr_2_type1$DSCA$R2$Group2),
                                       unlist(list_snr_2_type1$DGICA$R2$Group2),
                                       unlist(list_snr_2_type1$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_snr_2_type1$Unfitted_GICA$R2$Group2)))
R2_G2_snr2_type1$SNR <- 0.5
R2_G2_snr3_type1 <- data.frame(value=c(unlist(list_snr_3_type1$GRIDY$R2$Group2),
                                       unlist(list_snr_3_type1$SCA_P$R2$Group2),
                                       unlist(list_snr_3_type1$GICA$R2$Group2),
                                       unlist(list_snr_3_type1$DSCA$R2$Group2),
                                       unlist(list_snr_3_type1$DGICA$R2$Group2),
                                       unlist(list_snr_3_type1$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_snr_3_type1$Unfitted_GICA$R2$Group2)))
R2_G2_snr3_type1$SNR <- 0.75
R2_G2_snr4_type1 <- data.frame(value=c(unlist(list_snr_4_type1$GRIDY$R2$Group2),
                                       unlist(list_snr_4_type1$SCA_P$R2$Group2),
                                       unlist(list_snr_4_type1$GICA$R2$Group2),
                                       unlist(list_snr_4_type1$DSCA$R2$Group2),
                                       unlist(list_snr_4_type1$DGICA$R2$Group2),
                                       unlist(list_snr_4_type1$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_snr_4_type1$Unfitted_GICA$R2$Group2)))
R2_G2_snr4_type1$SNR <- 1.25
R2_G2_snr5_type1 <- data.frame(value=c(unlist(list_snr_5_type1$GRIDY$R2$Group2),
                                       unlist(list_snr_5_type1$SCA_P$R2$Group2),
                                       unlist(list_snr_5_type1$GICA$R2$Group2),
                                       unlist(list_snr_5_type1$DSCA$R2$Group2),
                                       unlist(list_snr_5_type1$DGICA$R2$Group2),
                                       unlist(list_snr_5_type1$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_snr_5_type1$Unfitted_GICA$R2$Group2)))
R2_G2_snr5_type1$SNR <- 2
R2_G2_snr6_type1 <- data.frame(value=c(unlist(list_snr_6_type1$GRIDY$R2$Group2),
                                       unlist(list_snr_6_type1$SCA_P$R2$Group2),
                                       unlist(list_snr_6_type1$GICA$R2$Group2),
                                       unlist(list_snr_6_type1$DSCA$R2$Group2),
                                       unlist(list_snr_6_type1$DGICA$R2$Group2),
                                       unlist(list_snr_6_type1$Unfitted_SCA_PF2$R2$Group2),
                                       unlist(list_snr_6_type1$Unfitted_GICA$R2$Group2)))
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


#-----------------------------------------------------------------------------#
# RMSE
#-----------------------------------------------------------------------------#
RMSE_J_snr1_type1 <- data.frame(value=c(unlist(list_snr_1_type1$GRIDY$RMSE$Joint),
                                        unlist(list_snr_1_type1$SCA_P$RMSE$Joint),
                                        unlist(list_snr_1_type1$GICA$RMSE$Joint),
                                        unlist(list_snr_1_type1$DSCA$RMSE$Joint),
                                        unlist(list_snr_1_type1$DGICA$RMSE$Joint),
                                        unlist(list_snr_1_type1$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_snr_1_type1$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr1_type1$SNR <- 0.25
RMSE_J_snr2_type1 <- data.frame(value=c(unlist(list_snr_2_type1$GRIDY$RMSE$Joint),
                                        unlist(list_snr_2_type1$SCA_P$RMSE$Joint),
                                        unlist(list_snr_2_type1$GICA$RMSE$Joint),
                                        unlist(list_snr_2_type1$DSCA$RMSE$Joint),
                                        unlist(list_snr_2_type1$DGICA$RMSE$Joint),
                                        unlist(list_snr_2_type1$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_snr_2_type1$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr2_type1$SNR <- 0.5
RMSE_J_snr3_type1 <- data.frame(value=c(unlist(list_snr_3_type1$GRIDY$RMSE$Joint),
                                        unlist(list_snr_3_type1$SCA_P$RMSE$Joint),
                                        unlist(list_snr_3_type1$GICA$RMSE$Joint),
                                        unlist(list_snr_3_type1$DSCA$RMSE$Joint),
                                        unlist(list_snr_3_type1$DGICA$RMSE$Joint),
                                        unlist(list_snr_3_type1$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_snr_3_type1$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr3_type1$SNR <- 0.75
RMSE_J_snr4_type1 <- data.frame(value=c(unlist(list_snr_4_type1$GRIDY$RMSE$Joint),
                                        unlist(list_snr_4_type1$SCA_P$RMSE$Joint),
                                        unlist(list_snr_4_type1$GICA$RMSE$Joint),
                                        unlist(list_snr_4_type1$DSCA$RMSE$Joint),
                                        unlist(list_snr_4_type1$DGICA$RMSE$Joint),
                                        unlist(list_snr_4_type1$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_snr_4_type1$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr4_type1$SNR <- 1.25
RMSE_J_snr5_type1 <- data.frame(value=c(unlist(list_snr_5_type1$GRIDY$RMSE$Joint),
                                        unlist(list_snr_5_type1$SCA_P$RMSE$Joint),
                                        unlist(list_snr_5_type1$GICA$RMSE$Joint),
                                        unlist(list_snr_5_type1$DSCA$RMSE$Joint),
                                        unlist(list_snr_5_type1$DGICA$RMSE$Joint),
                                        unlist(list_snr_5_type1$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_snr_5_type1$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr5_type1$SNR <- 2
RMSE_J_snr6_type1 <- data.frame(value=c(unlist(list_snr_6_type1$GRIDY$RMSE$Joint),
                                        unlist(list_snr_6_type1$SCA_P$RMSE$Joint),
                                        unlist(list_snr_6_type1$GICA$RMSE$Joint),
                                        unlist(list_snr_6_type1$DSCA$RMSE$Joint),
                                        unlist(list_snr_6_type1$DGICA$RMSE$Joint),
                                        unlist(list_snr_6_type1$Unfitted_SCA_PF2$RMSE$Joint),
                                        unlist(list_snr_6_type1$Unfitted_GICA$RMSE$Joint)))
RMSE_J_snr6_type1$SNR <- 4
RMSE_J_type1 <- rbind(RMSE_J_snr1_type1,RMSE_J_snr2_type1,RMSE_J_snr3_type1,
                      RMSE_J_snr4_type1,RMSE_J_snr5_type1,RMSE_J_snr6_type1)
RMSE_J_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                 "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

RMSE_G1_snr1_type1 <- data.frame(value=c(unlist(list_snr_1_type1$GRIDY$RMSE$Group1),
                                         unlist(list_snr_1_type1$SCA_P$RMSE$Group1),
                                         unlist(list_snr_1_type1$GICA$RMSE$Group1),
                                         unlist(list_snr_1_type1$DSCA$RMSE$Group1),
                                         unlist(list_snr_1_type1$DGICA$RMSE$Group1),
                                         unlist(list_snr_1_type1$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_snr_1_type1$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr1_type1$SNR <- 0.25
RMSE_G1_snr2_type1 <- data.frame(value=c(unlist(list_snr_2_type1$GRIDY$RMSE$Group1),
                                         unlist(list_snr_2_type1$SCA_P$RMSE$Group1),
                                         unlist(list_snr_2_type1$GICA$RMSE$Group1),
                                         unlist(list_snr_2_type1$DSCA$RMSE$Group1),
                                         unlist(list_snr_2_type1$DGICA$RMSE$Group1),
                                         unlist(list_snr_2_type1$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_snr_2_type1$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr2_type1$SNR <- 0.5
RMSE_G1_snr3_type1 <- data.frame(value=c(unlist(list_snr_3_type1$GRIDY$RMSE$Group1),
                                         unlist(list_snr_3_type1$SCA_P$RMSE$Group1),
                                         unlist(list_snr_3_type1$GICA$RMSE$Group1),
                                         unlist(list_snr_3_type1$DSCA$RMSE$Group1),
                                         unlist(list_snr_3_type1$DGICA$RMSE$Group1),
                                         unlist(list_snr_3_type1$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_snr_3_type1$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr3_type1$SNR <- 0.75
RMSE_G1_snr4_type1 <- data.frame(value=c(unlist(list_snr_4_type1$GRIDY$RMSE$Group1),
                                         unlist(list_snr_4_type1$SCA_P$RMSE$Group1),
                                         unlist(list_snr_4_type1$GICA$RMSE$Group1),
                                         unlist(list_snr_4_type1$DSCA$RMSE$Group1),
                                         unlist(list_snr_4_type1$DGICA$RMSE$Group1),
                                         unlist(list_snr_4_type1$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_snr_4_type1$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr4_type1$SNR <- 1.25
RMSE_G1_snr5_type1 <- data.frame(value=c(unlist(list_snr_5_type1$GRIDY$RMSE$Group1),
                                         unlist(list_snr_5_type1$SCA_P$RMSE$Group1),
                                         unlist(list_snr_5_type1$GICA$RMSE$Group1),
                                         unlist(list_snr_5_type1$DSCA$RMSE$Group1),
                                         unlist(list_snr_5_type1$DGICA$RMSE$Group1),
                                         unlist(list_snr_5_type1$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_snr_5_type1$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr5_type1$SNR <- 2
RMSE_G1_snr6_type1 <- data.frame(value=c(unlist(list_snr_6_type1$GRIDY$RMSE$Group1),
                                         unlist(list_snr_6_type1$SCA_P$RMSE$Group1),
                                         unlist(list_snr_6_type1$GICA$RMSE$Group1),
                                         unlist(list_snr_6_type1$DSCA$RMSE$Group1),
                                         unlist(list_snr_6_type1$DGICA$RMSE$Group1),
                                         unlist(list_snr_6_type1$Unfitted_SCA_PF2$RMSE$Group1),
                                         unlist(list_snr_6_type1$Unfitted_GICA$RMSE$Group1)))
RMSE_G1_snr6_type1$SNR <- 4
RMSE_G1_type1 <- rbind(RMSE_G1_snr1_type1,RMSE_G1_snr2_type1,RMSE_G1_snr3_type1,
                       RMSE_G1_snr4_type1,RMSE_G1_snr5_type1,RMSE_G1_snr6_type1)
RMSE_G1_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                  "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


RMSE_G2_snr1_type1 <- data.frame(value=c(unlist(list_snr_1_type1$GRIDY$RMSE$Group2),
                                         unlist(list_snr_1_type1$SCA_P$RMSE$Group2),
                                         unlist(list_snr_1_type1$GICA$RMSE$Group2),
                                         unlist(list_snr_1_type1$DSCA$RMSE$Group2),
                                         unlist(list_snr_1_type1$DGICA$RMSE$Group2),
                                         unlist(list_snr_1_type1$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_snr_1_type1$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr1_type1$SNR <- 0.25
RMSE_G2_snr2_type1 <- data.frame(value=c(unlist(list_snr_2_type1$GRIDY$RMSE$Group2),
                                         unlist(list_snr_2_type1$SCA_P$RMSE$Group2),
                                         unlist(list_snr_2_type1$GICA$RMSE$Group2),
                                         unlist(list_snr_2_type1$DSCA$RMSE$Group2),
                                         unlist(list_snr_2_type1$DGICA$RMSE$Group2),
                                         unlist(list_snr_2_type1$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_snr_2_type1$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr2_type1$SNR <- 0.5
RMSE_G2_snr3_type1 <- data.frame(value=c(unlist(list_snr_3_type1$GRIDY$RMSE$Group2),
                                         unlist(list_snr_3_type1$SCA_P$RMSE$Group2),
                                         unlist(list_snr_3_type1$GICA$RMSE$Group2),
                                         unlist(list_snr_3_type1$DSCA$RMSE$Group2),
                                         unlist(list_snr_3_type1$DGICA$RMSE$Group2),
                                         unlist(list_snr_3_type1$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_snr_3_type1$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr3_type1$SNR <- 0.75
RMSE_G2_snr4_type1 <- data.frame(value=c(unlist(list_snr_4_type1$GRIDY$RMSE$Group2),
                                         unlist(list_snr_4_type1$SCA_P$RMSE$Group2),
                                         unlist(list_snr_4_type1$GICA$RMSE$Group2),
                                         unlist(list_snr_4_type1$DSCA$RMSE$Group2),
                                         unlist(list_snr_4_type1$DGICA$RMSE$Group2),
                                         unlist(list_snr_4_type1$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_snr_4_type1$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr4_type1$SNR <- 1.25
RMSE_G2_snr5_type1 <- data.frame(value=c(unlist(list_snr_5_type1$GRIDY$RMSE$Group2),
                                         unlist(list_snr_5_type1$SCA_P$RMSE$Group2),
                                         unlist(list_snr_5_type1$GICA$RMSE$Group2),
                                         unlist(list_snr_5_type1$DSCA$RMSE$Group2),
                                         unlist(list_snr_5_type1$DGICA$RMSE$Group2),
                                         unlist(list_snr_5_type1$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_snr_5_type1$Unfitted_GICA$RMSE$Group2)))
RMSE_G2_snr5_type1$SNR <- 2
RMSE_G2_snr6_type1 <- data.frame(value=c(unlist(list_snr_6_type1$GRIDY$RMSE$Group2),
                                         unlist(list_snr_6_type1$SCA_P$RMSE$Group2),
                                         unlist(list_snr_6_type1$GICA$RMSE$Group2),
                                         unlist(list_snr_6_type1$DSCA$RMSE$Group2),
                                         unlist(list_snr_6_type1$DGICA$RMSE$Group2),
                                         unlist(list_snr_6_type1$Unfitted_SCA_PF2$RMSE$Group2),
                                         unlist(list_snr_6_type1$Unfitted_GICA$RMSE$Group2)))
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


#-----------------------------------------------------------------------------#
# CC_B
#-----------------------------------------------------------------------------#
CC_B_J_snr1_type1 <- data.frame(value=c(unlist(list_snr_1_type1$GRIDY$CC_B$Joint),
                                        unlist(list_snr_1_type1$SCA_P$CC_B$Joint),
                                        unlist(list_snr_1_type1$GICA$CC_B$Joint),
                                        unlist(list_snr_1_type1$DSCA$CC_B$Joint),
                                        unlist(list_snr_1_type1$DGICA$CC_B$Joint),
                                        unlist(list_snr_1_type1$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_snr_1_type1$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr1_type1$SNR <- 0.25
CC_B_J_snr2_type1 <- data.frame(value=c(unlist(list_snr_2_type1$GRIDY$CC_B$Joint),
                                        unlist(list_snr_2_type1$SCA_P$CC_B$Joint),
                                        unlist(list_snr_2_type1$GICA$CC_B$Joint),
                                        unlist(list_snr_2_type1$DSCA$CC_B$Joint),
                                        unlist(list_snr_2_type1$DGICA$CC_B$Joint),
                                        unlist(list_snr_2_type1$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_snr_2_type1$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr2_type1$SNR <- 0.5
CC_B_J_snr3_type1 <- data.frame(value=c(unlist(list_snr_3_type1$GRIDY$CC_B$Joint),
                                        unlist(list_snr_3_type1$SCA_P$CC_B$Joint),
                                        unlist(list_snr_3_type1$GICA$CC_B$Joint),
                                        unlist(list_snr_3_type1$DSCA$CC_B$Joint),
                                        unlist(list_snr_3_type1$DGICA$CC_B$Joint),
                                        unlist(list_snr_3_type1$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_snr_3_type1$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr3_type1$SNR <- 0.75
CC_B_J_snr4_type1 <- data.frame(value=c(unlist(list_snr_4_type1$GRIDY$CC_B$Joint),
                                        unlist(list_snr_4_type1$SCA_P$CC_B$Joint),
                                        unlist(list_snr_4_type1$GICA$CC_B$Joint),
                                        unlist(list_snr_4_type1$DSCA$CC_B$Joint),
                                        unlist(list_snr_4_type1$DGICA$CC_B$Joint),
                                        unlist(list_snr_4_type1$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_snr_4_type1$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr4_type1$SNR <- 1.25
CC_B_J_snr5_type1 <- data.frame(value=c(unlist(list_snr_5_type1$GRIDY$CC_B$Joint),
                                        unlist(list_snr_5_type1$SCA_P$CC_B$Joint),
                                        unlist(list_snr_5_type1$GICA$CC_B$Joint),
                                        unlist(list_snr_5_type1$DSCA$CC_B$Joint),
                                        unlist(list_snr_5_type1$DGICA$CC_B$Joint),
                                        unlist(list_snr_5_type1$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_snr_5_type1$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr5_type1$SNR <- 2
CC_B_J_snr6_type1 <- data.frame(value=c(unlist(list_snr_6_type1$GRIDY$CC_B$Joint),
                                        unlist(list_snr_6_type1$SCA_P$CC_B$Joint),
                                        unlist(list_snr_6_type1$GICA$CC_B$Joint),
                                        unlist(list_snr_6_type1$DSCA$CC_B$Joint),
                                        unlist(list_snr_6_type1$DGICA$CC_B$Joint),
                                        unlist(list_snr_6_type1$Unfitted_SCA_PF2$CC_B$Joint),
                                        unlist(list_snr_6_type1$Unfitted_GICA$CC_B$Joint)))
CC_B_J_snr6_type1$SNR <- 4
CC_B_J_type1 <- rbind(CC_B_J_snr1_type1,CC_B_J_snr2_type1,CC_B_J_snr3_type1,
                      CC_B_J_snr4_type1,CC_B_J_snr5_type1,CC_B_J_snr6_type1)
CC_B_J_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                 "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

CC_B_G1_snr1_type1 <- data.frame(value=c(unlist(list_snr_1_type1$GRIDY$CC_B$Group1),
                                         unlist(list_snr_1_type1$SCA_P$CC_B$Group1),
                                         unlist(list_snr_1_type1$GICA$CC_B$Group1),
                                         unlist(list_snr_1_type1$DSCA$CC_B$Group1),
                                         unlist(list_snr_1_type1$DGICA$CC_B$Group1),
                                         unlist(list_snr_1_type1$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_snr_1_type1$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr1_type1$SNR <- 0.25
CC_B_G1_snr2_type1 <- data.frame(value=c(unlist(list_snr_2_type1$GRIDY$CC_B$Group1),
                                         unlist(list_snr_2_type1$SCA_P$CC_B$Group1),
                                         unlist(list_snr_2_type1$GICA$CC_B$Group1),
                                         unlist(list_snr_2_type1$DSCA$CC_B$Group1),
                                         unlist(list_snr_2_type1$DGICA$CC_B$Group1),
                                         unlist(list_snr_2_type1$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_snr_2_type1$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr2_type1$SNR <- 0.5
CC_B_G1_snr3_type1 <- data.frame(value=c(unlist(list_snr_3_type1$GRIDY$CC_B$Group1),
                                         unlist(list_snr_3_type1$SCA_P$CC_B$Group1),
                                         unlist(list_snr_3_type1$GICA$CC_B$Group1),
                                         unlist(list_snr_3_type1$DSCA$CC_B$Group1),
                                         unlist(list_snr_3_type1$DGICA$CC_B$Group1),
                                         unlist(list_snr_3_type1$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_snr_3_type1$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr3_type1$SNR <- 0.75
CC_B_G1_snr4_type1 <- data.frame(value=c(unlist(list_snr_4_type1$GRIDY$CC_B$Group1),
                                         unlist(list_snr_4_type1$SCA_P$CC_B$Group1),
                                         unlist(list_snr_4_type1$GICA$CC_B$Group1),
                                         unlist(list_snr_4_type1$DSCA$CC_B$Group1),
                                         unlist(list_snr_4_type1$DGICA$CC_B$Group1),
                                         unlist(list_snr_4_type1$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_snr_4_type1$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr4_type1$SNR <- 1.25
CC_B_G1_snr5_type1 <- data.frame(value=c(unlist(list_snr_5_type1$GRIDY$CC_B$Group1),
                                         unlist(list_snr_5_type1$SCA_P$CC_B$Group1),
                                         unlist(list_snr_5_type1$GICA$CC_B$Group1),
                                         unlist(list_snr_5_type1$DSCA$CC_B$Group1),
                                         unlist(list_snr_5_type1$DGICA$CC_B$Group1),
                                         unlist(list_snr_5_type1$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_snr_5_type1$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr5_type1$SNR <- 2
CC_B_G1_snr6_type1 <- data.frame(value=c(unlist(list_snr_6_type1$GRIDY$CC_B$Group1),
                                         unlist(list_snr_6_type1$SCA_P$CC_B$Group1),
                                         unlist(list_snr_6_type1$GICA$CC_B$Group1),
                                         unlist(list_snr_6_type1$DSCA$CC_B$Group1),
                                         unlist(list_snr_6_type1$DGICA$CC_B$Group1),
                                         unlist(list_snr_6_type1$Unfitted_SCA_PF2$CC_B$Group1),
                                         unlist(list_snr_6_type1$Unfitted_GICA$CC_B$Group1)))
CC_B_G1_snr6_type1$SNR <- 4
CC_B_G1_type1 <- rbind(CC_B_G1_snr1_type1,CC_B_G1_snr2_type1,CC_B_G1_snr3_type1,
                       CC_B_G1_snr4_type1,CC_B_G1_snr5_type1,CC_B_G1_snr6_type1)
CC_B_G1_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                  "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


CC_B_G2_snr1_type1 <- data.frame(value=c(unlist(list_snr_1_type1$GRIDY$CC_B$Group2),
                                         unlist(list_snr_1_type1$SCA_P$CC_B$Group2),
                                         unlist(list_snr_1_type1$GICA$CC_B$Group2),
                                         unlist(list_snr_1_type1$DSCA$CC_B$Group2),
                                         unlist(list_snr_1_type1$DGICA$CC_B$Group2),
                                         unlist(list_snr_1_type1$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_snr_1_type1$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr1_type1$SNR <- 0.25
CC_B_G2_snr2_type1 <- data.frame(value=c(unlist(list_snr_2_type1$GRIDY$CC_B$Group2),
                                         unlist(list_snr_2_type1$SCA_P$CC_B$Group2),
                                         unlist(list_snr_2_type1$GICA$CC_B$Group2),
                                         unlist(list_snr_2_type1$DSCA$CC_B$Group2),
                                         unlist(list_snr_2_type1$DGICA$CC_B$Group2),
                                         unlist(list_snr_2_type1$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_snr_2_type1$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr2_type1$SNR <- 0.5
CC_B_G2_snr3_type1 <- data.frame(value=c(unlist(list_snr_3_type1$GRIDY$CC_B$Group2),
                                         unlist(list_snr_3_type1$SCA_P$CC_B$Group2),
                                         unlist(list_snr_3_type1$GICA$CC_B$Group2),
                                         unlist(list_snr_3_type1$DSCA$CC_B$Group2),
                                         unlist(list_snr_3_type1$DGICA$CC_B$Group2),
                                         unlist(list_snr_3_type1$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_snr_3_type1$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr3_type1$SNR <- 0.75
CC_B_G2_snr4_type1 <- data.frame(value=c(unlist(list_snr_4_type1$GRIDY$CC_B$Group2),
                                         unlist(list_snr_4_type1$SCA_P$CC_B$Group2),
                                         unlist(list_snr_4_type1$GICA$CC_B$Group2),
                                         unlist(list_snr_4_type1$DSCA$CC_B$Group2),
                                         unlist(list_snr_4_type1$DGICA$CC_B$Group2),
                                         unlist(list_snr_4_type1$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_snr_4_type1$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr4_type1$SNR <- 1.25
CC_B_G2_snr5_type1 <- data.frame(value=c(unlist(list_snr_5_type1$GRIDY$CC_B$Group2),
                                         unlist(list_snr_5_type1$SCA_P$CC_B$Group2),
                                         unlist(list_snr_5_type1$GICA$CC_B$Group2),
                                         unlist(list_snr_5_type1$DSCA$CC_B$Group2),
                                         unlist(list_snr_5_type1$DGICA$CC_B$Group2),
                                         unlist(list_snr_5_type1$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_snr_5_type1$Unfitted_GICA$CC_B$Group2)))
CC_B_G2_snr5_type1$SNR <- 2
CC_B_G2_snr6_type1 <- data.frame(value=c(unlist(list_snr_6_type1$GRIDY$CC_B$Group2),
                                         unlist(list_snr_6_type1$SCA_P$CC_B$Group2),
                                         unlist(list_snr_6_type1$GICA$CC_B$Group2),
                                         unlist(list_snr_6_type1$DSCA$CC_B$Group2),
                                         unlist(list_snr_6_type1$DGICA$CC_B$Group2),
                                         unlist(list_snr_6_type1$Unfitted_SCA_PF2$CC_B$Group2),
                                         unlist(list_snr_6_type1$Unfitted_GICA$CC_B$Group2)))
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


# -----------------------------------------------------------------------------#
# CC_F
# -----------------------------------------------------------------------------#
CC_F_J_snr1_type1 <- data.frame(value=c(unlist(list_snr_1_type1$GRIDY$CC_F$Joint),
                                        unlist(list_snr_1_type1$SCA_P$CC_F$Joint),
                                        unlist(list_snr_1_type1$GICA$CC_F$Joint),
                                        unlist(list_snr_1_type1$DSCA$CC_F$Joint),
                                        unlist(list_snr_1_type1$DGICA$CC_F$Joint),
                                        unlist(list_snr_1_type1$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_snr_1_type1$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr1_type1$SNR <- 0.25
CC_F_J_snr2_type1 <- data.frame(value=c(unlist(list_snr_2_type1$GRIDY$CC_F$Joint),
                                        unlist(list_snr_2_type1$SCA_P$CC_F$Joint),
                                        unlist(list_snr_2_type1$GICA$CC_F$Joint),
                                        unlist(list_snr_2_type1$DSCA$CC_F$Joint),
                                        unlist(list_snr_2_type1$DGICA$CC_F$Joint),
                                        unlist(list_snr_2_type1$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_snr_2_type1$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr2_type1$SNR <- 0.5
CC_F_J_snr3_type1 <- data.frame(value=c(unlist(list_snr_3_type1$GRIDY$CC_F$Joint),
                                        unlist(list_snr_3_type1$SCA_P$CC_F$Joint),
                                        unlist(list_snr_3_type1$GICA$CC_F$Joint),
                                        unlist(list_snr_3_type1$DSCA$CC_F$Joint),
                                        unlist(list_snr_3_type1$DGICA$CC_F$Joint),
                                        unlist(list_snr_3_type1$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_snr_3_type1$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr3_type1$SNR <- 0.75
CC_F_J_snr4_type1 <- data.frame(value=c(unlist(list_snr_4_type1$GRIDY$CC_F$Joint),
                                        unlist(list_snr_4_type1$SCA_P$CC_F$Joint),
                                        unlist(list_snr_4_type1$GICA$CC_F$Joint),
                                        unlist(list_snr_4_type1$DSCA$CC_F$Joint),
                                        unlist(list_snr_4_type1$DGICA$CC_F$Joint),
                                        unlist(list_snr_4_type1$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_snr_4_type1$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr4_type1$SNR <- 1.25
CC_F_J_snr5_type1 <- data.frame(value=c(unlist(list_snr_5_type1$GRIDY$CC_F$Joint),
                                        unlist(list_snr_5_type1$SCA_P$CC_F$Joint),
                                        unlist(list_snr_5_type1$GICA$CC_F$Joint),
                                        unlist(list_snr_5_type1$DSCA$CC_F$Joint),
                                        unlist(list_snr_5_type1$DGICA$CC_F$Joint),
                                        unlist(list_snr_5_type1$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_snr_5_type1$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr5_type1$SNR <- 2
CC_F_J_snr6_type1 <- data.frame(value=c(unlist(list_snr_6_type1$GRIDY$CC_F$Joint),
                                        unlist(list_snr_6_type1$SCA_P$CC_F$Joint),
                                        unlist(list_snr_6_type1$GICA$CC_F$Joint),
                                        unlist(list_snr_6_type1$DSCA$CC_F$Joint),
                                        unlist(list_snr_6_type1$DGICA$CC_F$Joint),
                                        unlist(list_snr_6_type1$Unfitted_SCA_PF2$CC_F$Joint),
                                        unlist(list_snr_6_type1$Unfitted_GICA$CC_F$Joint)))
CC_F_J_snr6_type1$SNR <- 4
CC_F_J_type1 <- rbind(CC_F_J_snr1_type1,CC_F_J_snr2_type1,CC_F_J_snr3_type1,
                      CC_F_J_snr4_type1,CC_F_J_snr5_type1,CC_F_J_snr6_type1)
CC_F_J_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                 "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

CC_F_G1_snr1_type1 <- data.frame(value=c(unlist(list_snr_1_type1$GRIDY$CC_F$Group1),
                                         unlist(list_snr_1_type1$SCA_P$CC_F$Group1),
                                         unlist(list_snr_1_type1$GICA$CC_F$Group1),
                                         unlist(list_snr_1_type1$DSCA$CC_F$Group1),
                                         unlist(list_snr_1_type1$DGICA$CC_F$Group1),
                                         unlist(list_snr_1_type1$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_snr_1_type1$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr1_type1$SNR <- 0.25
CC_F_G1_snr2_type1 <- data.frame(value=c(unlist(list_snr_2_type1$GRIDY$CC_F$Group1),
                                         unlist(list_snr_2_type1$SCA_P$CC_F$Group1),
                                         unlist(list_snr_2_type1$GICA$CC_F$Group1),
                                         unlist(list_snr_2_type1$DSCA$CC_F$Group1),
                                         unlist(list_snr_2_type1$DGICA$CC_F$Group1),
                                         unlist(list_snr_2_type1$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_snr_2_type1$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr2_type1$SNR <- 0.5
CC_F_G1_snr3_type1 <- data.frame(value=c(unlist(list_snr_3_type1$GRIDY$CC_F$Group1),
                                         unlist(list_snr_3_type1$SCA_P$CC_F$Group1),
                                         unlist(list_snr_3_type1$GICA$CC_F$Group1),
                                         unlist(list_snr_3_type1$DSCA$CC_F$Group1),
                                         unlist(list_snr_3_type1$DGICA$CC_F$Group1),
                                         unlist(list_snr_3_type1$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_snr_3_type1$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr3_type1$SNR <- 0.75
CC_F_G1_snr4_type1 <- data.frame(value=c(unlist(list_snr_4_type1$GRIDY$CC_F$Group1),
                                         unlist(list_snr_4_type1$SCA_P$CC_F$Group1),
                                         unlist(list_snr_4_type1$GICA$CC_F$Group1),
                                         unlist(list_snr_4_type1$DSCA$CC_F$Group1),
                                         unlist(list_snr_4_type1$DGICA$CC_F$Group1),
                                         unlist(list_snr_4_type1$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_snr_4_type1$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr4_type1$SNR <- 1.25
CC_F_G1_snr5_type1 <- data.frame(value=c(unlist(list_snr_5_type1$GRIDY$CC_F$Group1),
                                         unlist(list_snr_5_type1$SCA_P$CC_F$Group1),
                                         unlist(list_snr_5_type1$GICA$CC_F$Group1),
                                         unlist(list_snr_5_type1$DSCA$CC_F$Group1),
                                         unlist(list_snr_5_type1$DGICA$CC_F$Group1),
                                         unlist(list_snr_5_type1$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_snr_5_type1$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr5_type1$SNR <- 2
CC_F_G1_snr6_type1 <- data.frame(value=c(unlist(list_snr_6_type1$GRIDY$CC_F$Group1),
                                         unlist(list_snr_6_type1$SCA_P$CC_F$Group1),
                                         unlist(list_snr_6_type1$GICA$CC_F$Group1),
                                         unlist(list_snr_6_type1$DSCA$CC_F$Group1),
                                         unlist(list_snr_6_type1$DGICA$CC_F$Group1),
                                         unlist(list_snr_6_type1$Unfitted_SCA_PF2$CC_F$Group1),
                                         unlist(list_snr_6_type1$Unfitted_GICA$CC_F$Group1)))
CC_F_G1_snr6_type1$SNR <- 4
CC_F_G1_type1 <- rbind(CC_F_G1_snr1_type1,CC_F_G1_snr2_type1,CC_F_G1_snr3_type1,
                       CC_F_G1_snr4_type1,CC_F_G1_snr5_type1,CC_F_G1_snr6_type1)
CC_F_G1_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                  "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


CC_F_G2_snr1_type1 <- data.frame(value=c(unlist(list_snr_1_type1$GRIDY$CC_F$Group2),
                                         unlist(list_snr_1_type1$SCA_P$CC_F$Group2),
                                         unlist(list_snr_1_type1$GICA$CC_F$Group2),
                                         unlist(list_snr_1_type1$DSCA$CC_F$Group2),
                                         unlist(list_snr_1_type1$DGICA$CC_F$Group2),
                                         unlist(list_snr_1_type1$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_snr_1_type1$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr1_type1$SNR <- 0.25
CC_F_G2_snr2_type1 <- data.frame(value=c(unlist(list_snr_2_type1$GRIDY$CC_F$Group2),
                                         unlist(list_snr_2_type1$SCA_P$CC_F$Group2),
                                         unlist(list_snr_2_type1$GICA$CC_F$Group2),
                                         unlist(list_snr_2_type1$DSCA$CC_F$Group2),
                                         unlist(list_snr_2_type1$DGICA$CC_F$Group2),
                                         unlist(list_snr_2_type1$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_snr_2_type1$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr2_type1$SNR <- 0.5
CC_F_G2_snr3_type1 <- data.frame(value=c(unlist(list_snr_3_type1$GRIDY$CC_F$Group2),
                                         unlist(list_snr_3_type1$SCA_P$CC_F$Group2),
                                         unlist(list_snr_3_type1$GICA$CC_F$Group2),
                                         unlist(list_snr_3_type1$DSCA$CC_F$Group2),
                                         unlist(list_snr_3_type1$DGICA$CC_F$Group2),
                                         unlist(list_snr_3_type1$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_snr_3_type1$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr3_type1$SNR <- 0.75
CC_F_G2_snr4_type1 <- data.frame(value=c(unlist(list_snr_4_type1$GRIDY$CC_F$Group2),
                                         unlist(list_snr_4_type1$SCA_P$CC_F$Group2),
                                         unlist(list_snr_4_type1$GICA$CC_F$Group2),
                                         unlist(list_snr_4_type1$DSCA$CC_F$Group2),
                                         unlist(list_snr_4_type1$DGICA$CC_F$Group2),
                                         unlist(list_snr_4_type1$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_snr_4_type1$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr4_type1$SNR <- 1.25
CC_F_G2_snr5_type1 <- data.frame(value=c(unlist(list_snr_5_type1$GRIDY$CC_F$Group2),
                                         unlist(list_snr_5_type1$SCA_P$CC_F$Group2),
                                         unlist(list_snr_5_type1$GICA$CC_F$Group2),
                                         unlist(list_snr_5_type1$DSCA$CC_F$Group2),
                                         unlist(list_snr_5_type1$DGICA$CC_F$Group2),
                                         unlist(list_snr_5_type1$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_snr_5_type1$Unfitted_GICA$CC_F$Group2)))
CC_F_G2_snr5_type1$SNR <- 2
CC_F_G2_snr6_type1 <- data.frame(value=c(unlist(list_snr_6_type1$GRIDY$CC_F$Group2),
                                         unlist(list_snr_6_type1$SCA_P$CC_F$Group2),
                                         unlist(list_snr_6_type1$GICA$CC_F$Group2),
                                         unlist(list_snr_6_type1$DSCA$CC_F$Group2),
                                         unlist(list_snr_6_type1$DGICA$CC_F$Group2),
                                         unlist(list_snr_6_type1$Unfitted_SCA_PF2$CC_F$Group2),
                                         unlist(list_snr_6_type1$Unfitted_GICA$CC_F$Group2)))
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



# -----------------------------------------------------------------------------#
# CC_PSI
# -----------------------------------------------------------------------------#
CC_PSI_J_snr1_type1 <- data.frame(value=c(unlist(list_snr_1_type1$GRIDY$CC_PSI$Joint),
                                          unlist(list_snr_1_type1$SCA_P$CC_PSI$Joint),
                                          unlist(list_snr_1_type1$GICA$CC_PSI$Joint),
                                          unlist(list_snr_1_type1$DSCA$CC_PSI$Joint),
                                          unlist(list_snr_1_type1$DGICA$CC_PSI$Joint),
                                          unlist(list_snr_1_type1$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_snr_1_type1$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr1_type1$SNR <- 0.25
CC_PSI_J_snr2_type1 <- data.frame(value=c(unlist(list_snr_2_type1$GRIDY$CC_PSI$Joint),
                                          unlist(list_snr_2_type1$SCA_P$CC_PSI$Joint),
                                          unlist(list_snr_2_type1$GICA$CC_PSI$Joint),
                                          unlist(list_snr_2_type1$DSCA$CC_PSI$Joint),
                                          unlist(list_snr_2_type1$DGICA$CC_PSI$Joint),
                                          unlist(list_snr_2_type1$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_snr_2_type1$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr2_type1$SNR <- 0.5
CC_PSI_J_snr3_type1 <- data.frame(value=c(unlist(list_snr_3_type1$GRIDY$CC_PSI$Joint),
                                          unlist(list_snr_3_type1$SCA_P$CC_PSI$Joint),
                                          unlist(list_snr_3_type1$GICA$CC_PSI$Joint),
                                          unlist(list_snr_3_type1$DSCA$CC_PSI$Joint),
                                          unlist(list_snr_3_type1$DGICA$CC_PSI$Joint),
                                          unlist(list_snr_3_type1$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_snr_3_type1$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr3_type1$SNR <- 0.75
CC_PSI_J_snr4_type1 <- data.frame(value=c(unlist(list_snr_4_type1$GRIDY$CC_PSI$Joint),
                                          unlist(list_snr_4_type1$SCA_P$CC_PSI$Joint),
                                          unlist(list_snr_4_type1$GICA$CC_PSI$Joint),
                                          unlist(list_snr_4_type1$DSCA$CC_PSI$Joint),
                                          unlist(list_snr_4_type1$DGICA$CC_PSI$Joint),
                                          unlist(list_snr_4_type1$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_snr_4_type1$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr4_type1$SNR <- 1.25
CC_PSI_J_snr5_type1 <- data.frame(value=c(unlist(list_snr_5_type1$GRIDY$CC_PSI$Joint),
                                          unlist(list_snr_5_type1$SCA_P$CC_PSI$Joint),
                                          unlist(list_snr_5_type1$GICA$CC_PSI$Joint),
                                          unlist(list_snr_5_type1$DSCA$CC_PSI$Joint),
                                          unlist(list_snr_5_type1$DGICA$CC_PSI$Joint),
                                          unlist(list_snr_5_type1$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_snr_5_type1$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr5_type1$SNR <- 2
CC_PSI_J_snr6_type1 <- data.frame(value=c(unlist(list_snr_6_type1$GRIDY$CC_PSI$Joint),
                                          unlist(list_snr_6_type1$SCA_P$CC_PSI$Joint),
                                          unlist(list_snr_6_type1$GICA$CC_PSI$Joint),
                                          unlist(list_snr_6_type1$DSCA$CC_PSI$Joint),
                                          unlist(list_snr_6_type1$DGICA$CC_PSI$Joint),
                                          unlist(list_snr_6_type1$Unfitted_SCA_PF2$CC_PSI$Joint),
                                          unlist(list_snr_6_type1$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_J_snr6_type1$SNR <- 4
CC_PSI_J_type1 <- rbind(CC_PSI_J_snr1_type1,CC_PSI_J_snr2_type1,CC_PSI_J_snr3_type1,
                        CC_PSI_J_snr4_type1,CC_PSI_J_snr5_type1,CC_PSI_J_snr6_type1)
CC_PSI_J_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                   "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)

CC_PSI_G1_snr1_type1 <- data.frame(value=c(unlist(list_snr_1_type1$GRIDY$CC_PSI$Group1),
                                           unlist(list_snr_1_type1$SCA_P$CC_PSI$Group1),
                                           unlist(list_snr_1_type1$GICA$CC_PSI$Group1),
                                           unlist(list_snr_1_type1$DSCA$CC_PSI$Group1),
                                           unlist(list_snr_1_type1$DGICA$CC_PSI$Group1),
                                           unlist(list_snr_1_type1$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_snr_1_type1$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr1_type1$SNR <- 0.25
CC_PSI_G1_snr2_type1 <- data.frame(value=c(unlist(list_snr_2_type1$GRIDY$CC_PSI$Group1),
                                           unlist(list_snr_2_type1$SCA_P$CC_PSI$Group1),
                                           unlist(list_snr_2_type1$GICA$CC_PSI$Group1),
                                           unlist(list_snr_2_type1$DSCA$CC_PSI$Group1),
                                           unlist(list_snr_2_type1$DGICA$CC_PSI$Group1),
                                           unlist(list_snr_2_type1$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_snr_2_type1$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr2_type1$SNR <- 0.5
CC_PSI_G1_snr3_type1 <- data.frame(value=c(unlist(list_snr_3_type1$GRIDY$CC_PSI$Group1),
                                           unlist(list_snr_3_type1$SCA_P$CC_PSI$Group1),
                                           unlist(list_snr_3_type1$GICA$CC_PSI$Group1),
                                           unlist(list_snr_3_type1$DSCA$CC_PSI$Group1),
                                           unlist(list_snr_3_type1$DGICA$CC_PSI$Group1),
                                           unlist(list_snr_3_type1$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_snr_3_type1$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr3_type1$SNR <- 0.75
CC_PSI_G1_snr4_type1 <- data.frame(value=c(unlist(list_snr_4_type1$GRIDY$CC_PSI$Group1),
                                           unlist(list_snr_4_type1$SCA_P$CC_PSI$Group1),
                                           unlist(list_snr_4_type1$GICA$CC_PSI$Group1),
                                           unlist(list_snr_4_type1$DSCA$CC_PSI$Group1),
                                           unlist(list_snr_4_type1$DGICA$CC_PSI$Group1),
                                           unlist(list_snr_4_type1$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_snr_4_type1$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr4_type1$SNR <- 1.25
CC_PSI_G1_snr5_type1 <- data.frame(value=c(unlist(list_snr_5_type1$GRIDY$CC_PSI$Group1),
                                           unlist(list_snr_5_type1$SCA_P$CC_PSI$Group1),
                                           unlist(list_snr_5_type1$GICA$CC_PSI$Group1),
                                           unlist(list_snr_5_type1$DSCA$CC_PSI$Group1),
                                           unlist(list_snr_5_type1$DGICA$CC_PSI$Group1),
                                           unlist(list_snr_5_type1$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_snr_5_type1$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr5_type1$SNR <- 2
CC_PSI_G1_snr6_type1 <- data.frame(value=c(unlist(list_snr_6_type1$GRIDY$CC_PSI$Group1),
                                           unlist(list_snr_6_type1$SCA_P$CC_PSI$Group1),
                                           unlist(list_snr_6_type1$GICA$CC_PSI$Group1),
                                           unlist(list_snr_6_type1$DSCA$CC_PSI$Group1),
                                           unlist(list_snr_6_type1$DGICA$CC_PSI$Group1),
                                           unlist(list_snr_6_type1$Unfitted_SCA_PF2$CC_PSI$Group1),
                                           unlist(list_snr_6_type1$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G1_snr6_type1$SNR <- 4
CC_PSI_G1_type1 <- rbind(CC_PSI_G1_snr1_type1,CC_PSI_G1_snr2_type1,CC_PSI_G1_snr3_type1,
                         CC_PSI_G1_snr4_type1,CC_PSI_G1_snr5_type1,CC_PSI_G1_snr6_type1)
CC_PSI_G1_type1$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                    "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),6)


CC_PSI_G2_snr1_type1 <- data.frame(value=c(unlist(list_snr_1_type1$GRIDY$CC_PSI$Group2),
                                           unlist(list_snr_1_type1$SCA_P$CC_PSI$Group2),
                                           unlist(list_snr_1_type1$GICA$CC_PSI$Group2),
                                           unlist(list_snr_1_type1$DSCA$CC_PSI$Group2),
                                           unlist(list_snr_1_type1$DGICA$CC_PSI$Group2),
                                           unlist(list_snr_1_type1$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_snr_1_type1$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr1_type1$SNR <- 0.25
CC_PSI_G2_snr2_type1 <- data.frame(value=c(unlist(list_snr_2_type1$GRIDY$CC_PSI$Group2),
                                           unlist(list_snr_2_type1$SCA_P$CC_PSI$Group2),
                                           unlist(list_snr_2_type1$GICA$CC_PSI$Group2),
                                           unlist(list_snr_2_type1$DSCA$CC_PSI$Group2),
                                           unlist(list_snr_2_type1$DGICA$CC_PSI$Group2),
                                           unlist(list_snr_2_type1$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_snr_2_type1$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr2_type1$SNR <- 0.5
CC_PSI_G2_snr3_type1 <- data.frame(value=c(unlist(list_snr_3_type1$GRIDY$CC_PSI$Group2),
                                           unlist(list_snr_3_type1$SCA_P$CC_PSI$Group2),
                                           unlist(list_snr_3_type1$GICA$CC_PSI$Group2),
                                           unlist(list_snr_3_type1$DSCA$CC_PSI$Group2),
                                           unlist(list_snr_3_type1$DGICA$CC_PSI$Group2),
                                           unlist(list_snr_3_type1$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_snr_3_type1$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr3_type1$SNR <- 0.75
CC_PSI_G2_snr4_type1 <- data.frame(value=c(unlist(list_snr_4_type1$GRIDY$CC_PSI$Group2),
                                           unlist(list_snr_4_type1$SCA_P$CC_PSI$Group2),
                                           unlist(list_snr_4_type1$GICA$CC_PSI$Group2),
                                           unlist(list_snr_4_type1$DSCA$CC_PSI$Group2),
                                           unlist(list_snr_4_type1$DGICA$CC_PSI$Group2),
                                           unlist(list_snr_4_type1$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_snr_4_type1$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr4_type1$SNR <- 1.25
CC_PSI_G2_snr5_type1 <- data.frame(value=c(unlist(list_snr_5_type1$GRIDY$CC_PSI$Group2),
                                           unlist(list_snr_5_type1$SCA_P$CC_PSI$Group2),
                                           unlist(list_snr_5_type1$GICA$CC_PSI$Group2),
                                           unlist(list_snr_5_type1$DSCA$CC_PSI$Group2),
                                           unlist(list_snr_5_type1$DGICA$CC_PSI$Group2),
                                           unlist(list_snr_5_type1$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_snr_5_type1$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_G2_snr5_type1$SNR <- 2
CC_PSI_G2_snr6_type1 <- data.frame(value=c(unlist(list_snr_6_type1$GRIDY$CC_PSI$Group2),
                                           unlist(list_snr_6_type1$SCA_P$CC_PSI$Group2),
                                           unlist(list_snr_6_type1$GICA$CC_PSI$Group2),
                                           unlist(list_snr_6_type1$DSCA$CC_PSI$Group2),
                                           unlist(list_snr_6_type1$DGICA$CC_PSI$Group2),
                                           unlist(list_snr_6_type1$Unfitted_SCA_PF2$CC_PSI$Group2),
                                           unlist(list_snr_6_type1$Unfitted_GICA$CC_PSI$Group2)))
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



# -----------------------------------------------------------------------------#
# Saving Figure 3
# -----------------------------------------------------------------------------#
mylegend_type1 <- g_legend(pl_R2_J_type1)

pdf(file = "./figures/Figure3.pdf", width = 11.5, height = 10.5)
par(mar=c(0,0,0,0))
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
dev.off()
