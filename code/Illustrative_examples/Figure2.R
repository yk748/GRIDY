#-----------------------------------------------------------------------------#
# Code for producing results for plotting Figure 2 in Illustrative example section
#-----------------------------------------------------------------------------#


# -----------------------------------------------------------------------------#
# Loading simulated data
# -----------------------------------------------------------------------------#
# Packages required
library(ggplot2)
library(latex2exp)
library(ggpubr)
library(RColorBrewer)
library(reshape2)
library(gridExtra)

# Defining function for a legend 
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

load("./result/sim1_d200T100K50_type2.RData")
list_type2_d200T100K50 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                               DSCA=list_DSCA,DGICA=list_DGICA,
                               Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./result/sim1_d400T100K50_type2.RData")
list_type2_d400T100K50 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                               DSCA=list_DSCA,DGICA=list_DGICA,
                               Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)

load("./result/sim1_d100T200K50_type2.RData")
list_type2_d100T200K50 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                               DSCA=list_DSCA,DGICA=list_DGICA,
                               Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./result/sim1_d100T400K50_type2.RData")
list_type2_d100T400K50 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                               DSCA=list_DSCA,DGICA=list_DGICA,
                               Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)

load("./result/sim1_d100T200K10_type2.RData")
list_type2_d100T200K10 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                               DSCA=list_DSCA,DGICA=list_DGICA,
                               Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)
load("./result/sim1_d100T200K100_type2.RData")
list_type2_d100T200K100 <- list(GRIDY=list_GRIDY,SCA_P=list_SCA_P,GICA=list_GICA,
                                DSCA=list_DSCA,DGICA=list_DGICA,
                                Unfitted_SCA_PF2=list_Unfitted_SCA_PF2,Unfitted_GICA=list_Unfitted_GICA)

#-----------------------------------------------------------------------------#
# R2
#-----------------------------------------------------------------------------#
R2_J_type2_d200T100K50 <- data.frame(value=c(unlist(list_type2_d200T100K50$GRIDY$R2$Joint),
                                             unlist(list_type2_d200T100K50$SCA_P$R2$Joint),
                                             unlist(list_type2_d200T100K50$GICA$R2$Joint),
                                             unlist(list_type2_d200T100K50$DSCA$R2$Joint),
                                             unlist(list_type2_d200T100K50$DGICA$R2$Joint),
                                             unlist(list_type2_d200T100K50$Unfitted_SCA_PF2$R2$Joint),
                                             unlist(list_type2_d200T100K50$Unfitted_GICA$R2$Joint)))
R2_G1_type2_d200T100K50 <- data.frame(value=c(unlist(list_type2_d200T100K50$GRIDY$R2$Group1),
                                              unlist(list_type2_d200T100K50$SCA_P$R2$Group1),
                                              unlist(list_type2_d200T100K50$GICA$R2$Group1),
                                              unlist(list_type2_d200T100K50$DSCA$R2$Group1),
                                              unlist(list_type2_d200T100K50$DGICA$R2$Group1),
                                              unlist(list_type2_d200T100K50$Unfitted_SCA_PF2$R2$Group1),
                                              unlist(list_type2_d200T100K50$Unfitted_GICA$R2$Group1)))
R2_G2_type2_d200T100K50 <- data.frame(value=c(unlist(list_type2_d200T100K50$GRIDY$R2$Group2),
                                              unlist(list_type2_d200T100K50$SCA_P$R2$Group2),
                                              unlist(list_type2_d200T100K50$GICA$R2$Group2),
                                              unlist(list_type2_d200T100K50$DSCA$R2$Group2),
                                              unlist(list_type2_d200T100K50$DGICA$R2$Group2),
                                              unlist(list_type2_d200T100K50$Unfitted_SCA_PF2$R2$Group2),
                                              unlist(list_type2_d200T100K50$Unfitted_GICA$R2$Group2)))
R2_type2_d200T100K50 <- rbind(R2_J_type2_d200T100K50,R2_G1_type2_d200T100K50,R2_G2_type2_d200T100K50)
R2_type2_d200T100K50$size <- c("d=200")
R2_type2_d200T100K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                         "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
R2_type2_d200T100K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))


R2_J_type2_d400T100K50 <- data.frame(value=c(unlist(list_type2_d400T100K50$GRIDY$R2$Joint),
                                             unlist(list_type2_d400T100K50$SCA_P$R2$Joint),
                                             unlist(list_type2_d400T100K50$GICA$R2$Joint),
                                             unlist(list_type2_d400T100K50$DSCA$R2$Joint),
                                             unlist(list_type2_d400T100K50$DGICA$R2$Joint),
                                             unlist(list_type2_d400T100K50$Unfitted_SCA_PF2$R2$Joint),
                                             unlist(list_type2_d400T100K50$Unfitted_GICA$R2$Joint)))
R2_G1_type2_d400T100K50 <- data.frame(value=c(unlist(list_type2_d400T100K50$GRIDY$R2$Group1),
                                              unlist(list_type2_d400T100K50$SCA_P$R2$Group1),
                                              unlist(list_type2_d400T100K50$GICA$R2$Group1),
                                              unlist(list_type2_d400T100K50$DSCA$R2$Group1),
                                              unlist(list_type2_d400T100K50$DGICA$R2$Group1),
                                              unlist(list_type2_d400T100K50$Unfitted_SCA_PF2$R2$Group1),
                                              unlist(list_type2_d400T100K50$Unfitted_GICA$R2$Group1)))
R2_G2_type2_d400T100K50 <- data.frame(value=c(unlist(list_type2_d400T100K50$GRIDY$R2$Group2),
                                              unlist(list_type2_d400T100K50$SCA_P$R2$Group2),
                                              unlist(list_type2_d400T100K50$GICA$R2$Group2),
                                              unlist(list_type2_d400T100K50$DSCA$R2$Group2),
                                              unlist(list_type2_d400T100K50$DGICA$R2$Group2),
                                              unlist(list_type2_d400T100K50$Unfitted_SCA_PF2$R2$Group2),
                                              unlist(list_type2_d400T100K50$Unfitted_GICA$R2$Group2)))
R2_type2_d400T100K50 <- rbind(R2_J_type2_d400T100K50,R2_G1_type2_d400T100K50,R2_G2_type2_d400T100K50)
R2_type2_d400T100K50$size <- c("d=400")
R2_type2_d400T100K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                         "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
R2_type2_d400T100K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



R2_J_type2_d100T200K50 <- data.frame(value=c(unlist(list_type2_d100T200K50$GRIDY$R2$Joint),
                                             unlist(list_type2_d100T200K50$SCA_P$R2$Joint),
                                             unlist(list_type2_d100T200K50$GICA$R2$Joint),
                                             unlist(list_type2_d100T200K50$DSCA$R2$Joint),
                                             unlist(list_type2_d100T200K50$DGICA$R2$Joint),
                                             unlist(list_type2_d100T200K50$Unfitted_SCA_PF2$R2$Joint),
                                             unlist(list_type2_d100T200K50$Unfitted_GICA$R2$Joint)))
R2_G1_type2_d100T200K50 <- data.frame(value=c(unlist(list_type2_d100T200K50$GRIDY$R2$Group1),
                                              unlist(list_type2_d100T200K50$SCA_P$R2$Group1),
                                              unlist(list_type2_d100T200K50$GICA$R2$Group1),
                                              unlist(list_type2_d100T200K50$DSCA$R2$Group1),
                                              unlist(list_type2_d100T200K50$DGICA$R2$Group1),
                                              unlist(list_type2_d100T200K50$Unfitted_SCA_PF2$R2$Group1),
                                              unlist(list_type2_d100T200K50$Unfitted_GICA$R2$Group1)))
R2_G2_type2_d100T200K50 <- data.frame(value=c(unlist(list_type2_d100T200K50$GRIDY$R2$Group2),
                                              unlist(list_type2_d100T200K50$SCA_P$R2$Group2),
                                              unlist(list_type2_d100T200K50$GICA$R2$Group2),
                                              unlist(list_type2_d100T200K50$DSCA$R2$Group2),
                                              unlist(list_type2_d100T200K50$DGICA$R2$Group2),
                                              unlist(list_type2_d100T200K50$Unfitted_SCA_PF2$R2$Group2),
                                              unlist(list_type2_d100T200K50$Unfitted_GICA$R2$Group2)))
R2_type2_d100T200K50 <- rbind(R2_J_type2_d100T200K50,R2_G1_type2_d100T200K50,R2_G2_type2_d100T200K50)
R2_type2_d100T200K50$size <- c("T=200")
R2_type2_d100T200K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                         "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
R2_type2_d100T200K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



R2_J_type2_d100T400K50 <- data.frame(value=c(unlist(list_type2_d100T400K50$GRIDY$R2$Joint),
                                             unlist(list_type2_d100T400K50$SCA_P$R2$Joint),
                                             unlist(list_type2_d100T400K50$GICA$R2$Joint),
                                             unlist(list_type2_d100T400K50$DSCA$R2$Joint),
                                             unlist(list_type2_d100T400K50$DGICA$R2$Joint),
                                             unlist(list_type2_d100T400K50$Unfitted_SCA_PF2$R2$Joint),
                                             unlist(list_type2_d100T400K50$Unfitted_GICA$R2$Joint)))
R2_G1_type2_d100T400K50 <- data.frame(value=c(unlist(list_type2_d100T400K50$GRIDY$R2$Group1),
                                              unlist(list_type2_d100T400K50$SCA_P$R2$Group1),
                                              unlist(list_type2_d100T400K50$GICA$R2$Group1),
                                              unlist(list_type2_d100T400K50$DSCA$R2$Group1),
                                              unlist(list_type2_d100T400K50$DGICA$R2$Group1),
                                              unlist(list_type2_d100T400K50$Unfitted_SCA_PF2$R2$Group1),
                                              unlist(list_type2_d100T400K50$Unfitted_GICA$R2$Group1)))
R2_G2_type2_d100T400K50 <- data.frame(value=c(unlist(list_type2_d100T400K50$GRIDY$R2$Group2),
                                              unlist(list_type2_d100T400K50$SCA_P$R2$Group2),
                                              unlist(list_type2_d100T400K50$GICA$R2$Group2),
                                              unlist(list_type2_d100T400K50$DSCA$R2$Group2),
                                              unlist(list_type2_d100T400K50$DGICA$R2$Group2),
                                              unlist(list_type2_d100T400K50$Unfitted_SCA_PF2$R2$Group2),
                                              unlist(list_type2_d100T400K50$Unfitted_GICA$R2$Group2)))
R2_type2_d100T400K50 <- rbind(R2_J_type2_d100T400K50,R2_G1_type2_d100T400K50,R2_G2_type2_d100T400K50)
R2_type2_d100T400K50$size <- c("T=400")
R2_type2_d100T400K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                         "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
R2_type2_d100T400K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



R2_J_type2_d100T200K10 <- data.frame(value=c(unlist(list_type2_d100T200K10$GRIDY$R2$Joint),
                                             unlist(list_type2_d100T200K10$SCA_P$R2$Joint),
                                             unlist(list_type2_d100T200K10$GICA$R2$Joint),
                                             unlist(list_type2_d100T200K10$DSCA$R2$Joint),
                                             unlist(list_type2_d100T200K10$DGICA$R2$Joint),
                                             unlist(list_type2_d100T200K10$Unfitted_SCA_PF2$R2$Joint),
                                             unlist(list_type2_d100T200K10$Unfitted_GICA$R2$Joint)))
R2_G1_type2_d100T200K10 <- data.frame(value=c(unlist(list_type2_d100T200K10$GRIDY$R2$Group1),
                                              unlist(list_type2_d100T200K10$SCA_P$R2$Group1),
                                              unlist(list_type2_d100T200K10$GICA$R2$Group1),
                                              unlist(list_type2_d100T200K10$DSCA$R2$Group1),
                                              unlist(list_type2_d100T200K10$DGICA$R2$Group1),
                                              unlist(list_type2_d100T200K10$Unfitted_SCA_PF2$R2$Group1),
                                              unlist(list_type2_d100T200K10$Unfitted_GICA$R2$Group1)))
R2_G2_type2_d100T200K10 <- data.frame(value=c(unlist(list_type2_d100T200K10$GRIDY$R2$Group2),
                                              unlist(list_type2_d100T200K10$SCA_P$R2$Group2),
                                              unlist(list_type2_d100T200K10$GICA$R2$Group2),
                                              unlist(list_type2_d100T200K10$DSCA$R2$Group2),
                                              unlist(list_type2_d100T200K10$DGICA$R2$Group2),
                                              unlist(list_type2_d100T200K10$Unfitted_SCA_PF2$R2$Group2),
                                              unlist(list_type2_d100T200K10$Unfitted_GICA$R2$Group2)))
R2_type2_d100T200K10 <- rbind(R2_J_type2_d100T200K10,R2_G1_type2_d100T200K10,R2_G2_type2_d100T200K10)
R2_type2_d100T200K10$size <- c("2K=20")
R2_type2_d100T200K10$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                         "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
R2_type2_d100T200K10$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))


R2_J_type2_d100T200K100 <- data.frame(value=c(unlist(list_type2_d100T200K100$GRIDY$R2$Joint),
                                              unlist(list_type2_d100T200K100$SCA_P$R2$Joint),
                                              unlist(list_type2_d100T200K100$GICA$R2$Joint),
                                              unlist(list_type2_d100T200K100$DSCA$R2$Joint),
                                              unlist(list_type2_d100T200K100$DGICA$R2$Joint),
                                              unlist(list_type2_d100T200K100$Unfitted_SCA_PF2$R2$Joint),
                                              unlist(list_type2_d100T200K100$Unfitted_GICA$R2$Joint)))
R2_G1_type2_d100T200K100 <- data.frame(value=c(unlist(list_type2_d100T200K100$GRIDY$R2$Group1),
                                               unlist(list_type2_d100T200K100$SCA_P$R2$Group1),
                                               unlist(list_type2_d100T200K100$GICA$R2$Group1),
                                               unlist(list_type2_d100T200K100$DSCA$R2$Group1),
                                               unlist(list_type2_d100T200K100$DGICA$R2$Group1),
                                               unlist(list_type2_d100T200K100$Unfitted_SCA_PF2$R2$Group1),
                                               unlist(list_type2_d100T200K100$Unfitted_GICA$R2$Group1)))
R2_G2_type2_d100T200K100 <- data.frame(value=c(unlist(list_type2_d100T200K100$GRIDY$R2$Group2),
                                               unlist(list_type2_d100T200K100$SCA_P$R2$Group2),
                                               unlist(list_type2_d100T200K100$GICA$R2$Group2),
                                               unlist(list_type2_d100T200K100$DSCA$R2$Group2),
                                               unlist(list_type2_d100T200K100$DGICA$R2$Group2),
                                               unlist(list_type2_d100T200K100$Unfitted_SCA_PF2$R2$Group2),
                                               unlist(list_type2_d100T200K100$Unfitted_GICA$R2$Group2)))
R2_type2_d100T200K100 <- rbind(R2_J_type2_d100T200K100,R2_G1_type2_d100T200K100,R2_G2_type2_d100T200K100)
R2_type2_d100T200K100$size <- c("2K=200")
R2_type2_d100T200K100$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                          "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
R2_type2_d100T200K100$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



R2_J_type2 <- rbind(R2_type2_d200T100K50[which(R2_type2_d200T100K50$structure == "Joint"),],
                    R2_type2_d400T100K50[which(R2_type2_d400T100K50$structure == "Joint"),],
                    R2_type2_d100T200K50[which(R2_type2_d100T200K50$structure == "Joint"),],
                    R2_type2_d100T400K50[which(R2_type2_d100T200K50$structure == "Joint"),],
                    R2_type2_d100T200K10[which(R2_type2_d100T200K10$structure == "Joint"),],
                    R2_type2_d100T200K100[which(R2_type2_d100T200K100$structure == "Joint"),])
R2_G1_type2 <- rbind(R2_type2_d200T100K50[which(R2_type2_d200T100K50$structure == "Group Individual1"),],
                     R2_type2_d400T100K50[which(R2_type2_d400T100K50$structure == "Group Individual1"),],
                     R2_type2_d100T200K50[which(R2_type2_d100T200K50$structure == "Group Individual1"),],
                     R2_type2_d100T400K50[which(R2_type2_d100T200K50$structure == "Group Individual1"),],
                     R2_type2_d100T200K10[which(R2_type2_d100T200K10$structure == "Group Individual1"),],
                     R2_type2_d100T200K100[which(R2_type2_d100T200K100$structure == "Group Individual1"),])
R2_G2_type2 <- rbind(R2_type2_d200T100K50[which(R2_type2_d200T100K50$structure == "Group Individual2"),],
                     R2_type2_d400T100K50[which(R2_type2_d400T100K50$structure == "Group Individual2"),],
                     R2_type2_d100T200K50[which(R2_type2_d100T200K50$structure == "Group Individual2"),],
                     R2_type2_d100T400K50[which(R2_type2_d100T200K50$structure == "Group Individual2"),],
                     R2_type2_d100T200K10[which(R2_type2_d100T200K10$structure == "Group Individual2"),],
                     R2_type2_d100T200K100[which(R2_type2_d100T200K100$structure == "Group Individual2"),])


method_list <- c("GRIDY","SCA-P","GICA","DSCA","DGICA","Unfitted_SCA-PF2","Unfitted_GICA")
R2_J_type2$method <- factor(R2_J_type2$method,levels=method_list,labels=method_list)
R2_G1_type2$method <- factor(R2_G1_type2$method,levels=method_list,labels=method_list)
R2_G2_type2$method <- factor(R2_G2_type2$method,levels=method_list,labels=method_list)

size_list <- c("d=200","d=400","T=200","T=400","2K=20","2K=200")
R2_J_type2$size<- factor(R2_J_type2$size,levels=size_list,labels=size_list)
R2_G1_type2$size <- factor(R2_G1_type2$size,levels=size_list,labels=size_list)
R2_G2_type2$size <- factor(R2_G2_type2$size,levels=size_list,labels=size_list)


pl_R2_J_type2 <- ggplot(R2_J_type2,aes(x=size,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="Size") +
  scale_y_continuous(name=TeX("$R^2$ (Joint)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_R2_G1_type2 <- ggplot(R2_G1_type2,aes(x=size,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="Size") +
  scale_y_continuous(name=TeX("$R^2$ (Group 1)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_R2_G2_type2 <- ggplot(R2_G2_type2,aes(x=size,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="Size") +
  scale_y_continuous(name=TeX("$R^2$ (Group 2)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))


#-----------------------------------------------------------------------------#
# RMSE
#-----------------------------------------------------------------------------#
RMSE_J_type2_d200T100K50 <- data.frame(value=c(unlist(list_type2_d200T100K50$GRIDY$RMSE$Joint),
                                               unlist(list_type2_d200T100K50$SCA_P$RMSE$Joint),
                                               unlist(list_type2_d200T100K50$GICA$RMSE$Joint),
                                               unlist(list_type2_d200T100K50$DSCA$RMSE$Joint),
                                               unlist(list_type2_d200T100K50$DGICA$RMSE$Joint),
                                               unlist(list_type2_d200T100K50$Unfitted_SCA_PF2$RMSE$Joint),
                                               unlist(list_type2_d200T100K50$Unfitted_GICA$RMSE$Joint)))
RMSE_G1_type2_d200T100K50 <- data.frame(value=c(unlist(list_type2_d200T100K50$GRIDY$RMSE$Group1),
                                                unlist(list_type2_d200T100K50$SCA_P$RMSE$Group1),
                                                unlist(list_type2_d200T100K50$GICA$RMSE$Group1),
                                                unlist(list_type2_d200T100K50$DSCA$RMSE$Group1),
                                                unlist(list_type2_d200T100K50$DGICA$RMSE$Group1),
                                                unlist(list_type2_d200T100K50$Unfitted_SCA_PF2$RMSE$Group1),
                                                unlist(list_type2_d200T100K50$Unfitted_GICA$RMSE$Group1)))
RMSE_G2_type2_d200T100K50 <- data.frame(value=c(unlist(list_type2_d200T100K50$GRIDY$RMSE$Group2),
                                                unlist(list_type2_d200T100K50$SCA_P$RMSE$Group2),
                                                unlist(list_type2_d200T100K50$GICA$RMSE$Group2),
                                                unlist(list_type2_d200T100K50$DSCA$RMSE$Group2),
                                                unlist(list_type2_d200T100K50$DGICA$RMSE$Group2),
                                                unlist(list_type2_d200T100K50$Unfitted_SCA_PF2$RMSE$Group2),
                                                unlist(list_type2_d200T100K50$Unfitted_GICA$RMSE$Group2)))
RMSE_type2_d200T100K50 <- rbind(RMSE_J_type2_d200T100K50,RMSE_G1_type2_d200T100K50,RMSE_G2_type2_d200T100K50)
RMSE_type2_d200T100K50$size <- c("d=200")
RMSE_type2_d200T100K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                           "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
RMSE_type2_d200T100K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))


RMSE_J_type2_d400T100K50 <- data.frame(value=c(unlist(list_type2_d400T100K50$GRIDY$RMSE$Joint),
                                               unlist(list_type2_d400T100K50$SCA_P$RMSE$Joint),
                                               unlist(list_type2_d400T100K50$GICA$RMSE$Joint),
                                               unlist(list_type2_d400T100K50$DSCA$RMSE$Joint),
                                               unlist(list_type2_d400T100K50$DGICA$RMSE$Joint),
                                               unlist(list_type2_d400T100K50$Unfitted_SCA_PF2$RMSE$Joint),
                                               unlist(list_type2_d400T100K50$Unfitted_GICA$RMSE$Joint)))
RMSE_G1_type2_d400T100K50 <- data.frame(value=c(unlist(list_type2_d400T100K50$GRIDY$RMSE$Group1),
                                                unlist(list_type2_d400T100K50$SCA_P$RMSE$Group1),
                                                unlist(list_type2_d400T100K50$GICA$RMSE$Group1),
                                                unlist(list_type2_d400T100K50$DSCA$RMSE$Group1),
                                                unlist(list_type2_d400T100K50$DGICA$RMSE$Group1),
                                                unlist(list_type2_d400T100K50$Unfitted_SCA_PF2$RMSE$Group1),
                                                unlist(list_type2_d400T100K50$Unfitted_GICA$RMSE$Group1)))
RMSE_G2_type2_d400T100K50 <- data.frame(value=c(unlist(list_type2_d400T100K50$GRIDY$RMSE$Group2),
                                                unlist(list_type2_d400T100K50$SCA_P$RMSE$Group2),
                                                unlist(list_type2_d400T100K50$GICA$RMSE$Group2),
                                                unlist(list_type2_d400T100K50$DSCA$RMSE$Group2),
                                                unlist(list_type2_d400T100K50$DGICA$RMSE$Group2),
                                                unlist(list_type2_d400T100K50$Unfitted_SCA_PF2$RMSE$Group2),
                                                unlist(list_type2_d400T100K50$Unfitted_GICA$RMSE$Group2)))
RMSE_type2_d400T100K50 <- rbind(RMSE_J_type2_d400T100K50,RMSE_G1_type2_d400T100K50,RMSE_G2_type2_d400T100K50)
RMSE_type2_d400T100K50$size <- c("d=400")
RMSE_type2_d400T100K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                           "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
RMSE_type2_d400T100K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



RMSE_J_type2_d100T200K50 <- data.frame(value=c(unlist(list_type2_d100T200K50$GRIDY$RMSE$Joint),
                                               unlist(list_type2_d100T200K50$SCA_P$RMSE$Joint),
                                               unlist(list_type2_d100T200K50$GICA$RMSE$Joint),
                                               unlist(list_type2_d100T200K50$DSCA$RMSE$Joint),
                                               unlist(list_type2_d100T200K50$DGICA$RMSE$Joint),
                                               unlist(list_type2_d100T200K50$Unfitted_SCA_PF2$RMSE$Joint),
                                               unlist(list_type2_d100T200K50$Unfitted_GICA$RMSE$Joint)))
RMSE_G1_type2_d100T200K50 <- data.frame(value=c(unlist(list_type2_d100T200K50$GRIDY$RMSE$Group1),
                                                unlist(list_type2_d100T200K50$SCA_P$RMSE$Group1),
                                                unlist(list_type2_d100T200K50$GICA$RMSE$Group1),
                                                unlist(list_type2_d100T200K50$DSCA$RMSE$Group1),
                                                unlist(list_type2_d100T200K50$DGICA$RMSE$Group1),
                                                unlist(list_type2_d100T200K50$Unfitted_SCA_PF2$RMSE$Group1),
                                                unlist(list_type2_d100T200K50$Unfitted_GICA$RMSE$Group1)))
RMSE_G2_type2_d100T200K50 <- data.frame(value=c(unlist(list_type2_d100T200K50$GRIDY$RMSE$Group2),
                                                unlist(list_type2_d100T200K50$SCA_P$RMSE$Group2),
                                                unlist(list_type2_d100T200K50$GICA$RMSE$Group2),
                                                unlist(list_type2_d100T200K50$DSCA$RMSE$Group2),
                                                unlist(list_type2_d100T200K50$DGICA$RMSE$Group2),
                                                unlist(list_type2_d100T200K50$Unfitted_SCA_PF2$RMSE$Group2),
                                                unlist(list_type2_d100T200K50$Unfitted_GICA$RMSE$Group2)))
RMSE_type2_d100T200K50 <- rbind(RMSE_J_type2_d100T200K50,RMSE_G1_type2_d100T200K50,RMSE_G2_type2_d100T200K50)
RMSE_type2_d100T200K50$size <- c("T=200")
RMSE_type2_d100T200K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                           "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
RMSE_type2_d100T200K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



RMSE_J_type2_d100T400K50 <- data.frame(value=c(unlist(list_type2_d100T400K50$GRIDY$RMSE$Joint),
                                               unlist(list_type2_d100T400K50$SCA_P$RMSE$Joint),
                                               unlist(list_type2_d100T400K50$GICA$RMSE$Joint),
                                               unlist(list_type2_d100T400K50$DSCA$RMSE$Joint),
                                               unlist(list_type2_d100T400K50$DGICA$RMSE$Joint),
                                               unlist(list_type2_d100T400K50$Unfitted_SCA_PF2$RMSE$Joint),
                                               unlist(list_type2_d100T400K50$Unfitted_GICA$RMSE$Joint)))
RMSE_G1_type2_d100T400K50 <- data.frame(value=c(unlist(list_type2_d100T400K50$GRIDY$RMSE$Group1),
                                                unlist(list_type2_d100T400K50$SCA_P$RMSE$Group1),
                                                unlist(list_type2_d100T400K50$GICA$RMSE$Group1),
                                                unlist(list_type2_d100T400K50$DSCA$RMSE$Group1),
                                                unlist(list_type2_d100T400K50$DGICA$RMSE$Group1),
                                                unlist(list_type2_d100T400K50$Unfitted_SCA_PF2$RMSE$Group1),
                                                unlist(list_type2_d100T400K50$Unfitted_GICA$RMSE$Group1)))
RMSE_G2_type2_d100T400K50 <- data.frame(value=c(unlist(list_type2_d100T400K50$GRIDY$RMSE$Group2),
                                                unlist(list_type2_d100T400K50$SCA_P$RMSE$Group2),
                                                unlist(list_type2_d100T400K50$GICA$RMSE$Group2),
                                                unlist(list_type2_d100T400K50$DSCA$RMSE$Group2),
                                                unlist(list_type2_d100T400K50$DGICA$RMSE$Group2),
                                                unlist(list_type2_d100T400K50$Unfitted_SCA_PF2$RMSE$Group2),
                                                unlist(list_type2_d100T400K50$Unfitted_GICA$RMSE$Group2)))
RMSE_type2_d100T400K50 <- rbind(RMSE_J_type2_d100T400K50,RMSE_G1_type2_d100T400K50,RMSE_G2_type2_d100T400K50)
RMSE_type2_d100T400K50$size <- c("T=400")
RMSE_type2_d100T400K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                           "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
RMSE_type2_d100T400K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



RMSE_J_type2_d100T200K10 <- data.frame(value=c(unlist(list_type2_d100T200K10$GRIDY$RMSE$Joint),
                                               unlist(list_type2_d100T200K10$SCA_P$RMSE$Joint),
                                               unlist(list_type2_d100T200K10$GICA$RMSE$Joint),
                                               unlist(list_type2_d100T200K10$DSCA$RMSE$Joint),
                                               unlist(list_type2_d100T200K10$DGICA$RMSE$Joint),
                                               unlist(list_type2_d100T200K10$Unfitted_SCA_PF2$RMSE$Joint),
                                               unlist(list_type2_d100T200K10$Unfitted_GICA$RMSE$Joint)))
RMSE_G1_type2_d100T200K10 <- data.frame(value=c(unlist(list_type2_d100T200K10$GRIDY$RMSE$Group1),
                                                unlist(list_type2_d100T200K10$SCA_P$RMSE$Group1),
                                                unlist(list_type2_d100T200K10$GICA$RMSE$Group1),
                                                unlist(list_type2_d100T200K10$DSCA$RMSE$Group1),
                                                unlist(list_type2_d100T200K10$DGICA$RMSE$Group1),
                                                unlist(list_type2_d100T200K10$Unfitted_SCA_PF2$RMSE$Group1),
                                                unlist(list_type2_d100T200K10$Unfitted_GICA$RMSE$Group1)))
RMSE_G2_type2_d100T200K10 <- data.frame(value=c(unlist(list_type2_d100T200K10$GRIDY$RMSE$Group2),
                                                unlist(list_type2_d100T200K10$SCA_P$RMSE$Group2),
                                                unlist(list_type2_d100T200K10$GICA$RMSE$Group2),
                                                unlist(list_type2_d100T200K10$DSCA$RMSE$Group2),
                                                unlist(list_type2_d100T200K10$DGICA$RMSE$Group2),
                                                unlist(list_type2_d100T200K10$Unfitted_SCA_PF2$RMSE$Group2),
                                                unlist(list_type2_d100T200K10$Unfitted_GICA$RMSE$Group2)))
RMSE_type2_d100T200K10 <- rbind(RMSE_J_type2_d100T200K10,RMSE_G1_type2_d100T200K10,RMSE_G2_type2_d100T200K10)
RMSE_type2_d100T200K10$size <- c("2K=20")
RMSE_type2_d100T200K10$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                           "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
RMSE_type2_d100T200K10$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))


RMSE_J_type2_d100T200K100 <- data.frame(value=c(unlist(list_type2_d100T200K100$GRIDY$RMSE$Joint),
                                                unlist(list_type2_d100T200K100$SCA_P$RMSE$Joint),
                                                unlist(list_type2_d100T200K100$GICA$RMSE$Joint),
                                                unlist(list_type2_d100T200K100$DSCA$RMSE$Joint),
                                                unlist(list_type2_d100T200K100$DGICA$RMSE$Joint),
                                                unlist(list_type2_d100T200K100$Unfitted_SCA_PF2$RMSE$Joint),
                                                unlist(list_type2_d100T200K100$Unfitted_GICA$RMSE$Joint)))
RMSE_G1_type2_d100T200K100 <- data.frame(value=c(unlist(list_type2_d100T200K100$GRIDY$RMSE$Group1),
                                                 unlist(list_type2_d100T200K100$SCA_P$RMSE$Group1),
                                                 unlist(list_type2_d100T200K100$GICA$RMSE$Group1),
                                                 unlist(list_type2_d100T200K100$DSCA$RMSE$Group1),
                                                 unlist(list_type2_d100T200K100$DGICA$RMSE$Group1),
                                                 unlist(list_type2_d100T200K100$Unfitted_SCA_PF2$RMSE$Group1),
                                                 unlist(list_type2_d100T200K100$Unfitted_GICA$RMSE$Group1)))
RMSE_G2_type2_d100T200K100 <- data.frame(value=c(unlist(list_type2_d100T200K100$GRIDY$RMSE$Group2),
                                                 unlist(list_type2_d100T200K100$SCA_P$RMSE$Group2),
                                                 unlist(list_type2_d100T200K100$GICA$RMSE$Group2),
                                                 unlist(list_type2_d100T200K100$DSCA$RMSE$Group2),
                                                 unlist(list_type2_d100T200K100$DGICA$RMSE$Group2),
                                                 unlist(list_type2_d100T200K100$Unfitted_SCA_PF2$RMSE$Group2),
                                                 unlist(list_type2_d100T200K100$Unfitted_GICA$RMSE$Group2)))
RMSE_type2_d100T200K100 <- rbind(RMSE_J_type2_d100T200K100,RMSE_G1_type2_d100T200K100,RMSE_G2_type2_d100T200K100)
RMSE_type2_d100T200K100$size <- c("2K=200")
RMSE_type2_d100T200K100$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                            "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
RMSE_type2_d100T200K100$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



RMSE_J_type2 <- rbind(RMSE_type2_d200T100K50[which(RMSE_type2_d200T100K50$structure == "Joint"),],
                      RMSE_type2_d400T100K50[which(RMSE_type2_d400T100K50$structure == "Joint"),],
                      RMSE_type2_d100T200K50[which(RMSE_type2_d100T200K50$structure == "Joint"),],
                      RMSE_type2_d100T400K50[which(RMSE_type2_d100T200K50$structure == "Joint"),],
                      RMSE_type2_d100T200K10[which(RMSE_type2_d100T200K10$structure == "Joint"),],
                      RMSE_type2_d100T200K100[which(RMSE_type2_d100T200K100$structure == "Joint"),])
RMSE_G1_type2 <- rbind(RMSE_type2_d200T100K50[which(RMSE_type2_d200T100K50$structure == "Group Individual1"),],
                       RMSE_type2_d400T100K50[which(RMSE_type2_d400T100K50$structure == "Group Individual1"),],
                       RMSE_type2_d100T200K50[which(RMSE_type2_d100T200K50$structure == "Group Individual1"),],
                       RMSE_type2_d100T400K50[which(RMSE_type2_d100T200K50$structure == "Group Individual1"),],
                       RMSE_type2_d100T200K10[which(RMSE_type2_d100T200K10$structure == "Group Individual1"),],
                       RMSE_type2_d100T200K100[which(RMSE_type2_d100T200K100$structure == "Group Individual1"),])
RMSE_G2_type2 <- rbind(RMSE_type2_d200T100K50[which(RMSE_type2_d200T100K50$structure == "Group Individual2"),],
                       RMSE_type2_d400T100K50[which(RMSE_type2_d400T100K50$structure == "Group Individual2"),],
                       RMSE_type2_d100T200K50[which(RMSE_type2_d100T200K50$structure == "Group Individual2"),],
                       RMSE_type2_d100T400K50[which(RMSE_type2_d100T200K50$structure == "Group Individual2"),],
                       RMSE_type2_d100T200K10[which(RMSE_type2_d100T200K10$structure == "Group Individual2"),],
                       RMSE_type2_d100T200K100[which(RMSE_type2_d100T200K100$structure == "Group Individual2"),])


method_list <- c("GRIDY","SCA-P","GICA","DSCA","DGICA","Unfitted_SCA-PF2","Unfitted_GICA")
RMSE_J_type2$method <- factor(RMSE_J_type2$method,levels=method_list,labels=method_list)
RMSE_G1_type2$method <- factor(RMSE_G1_type2$method,levels=method_list,labels=method_list)
RMSE_G2_type2$method <- factor(RMSE_G2_type2$method,levels=method_list,labels=method_list)

size_list <- c("d=200","d=400","T=200","T=400","2K=20","2K=200")
RMSE_J_type2$size<- factor(RMSE_J_type2$size,levels=size_list,labels=size_list)
RMSE_G1_type2$size <- factor(RMSE_G1_type2$size,levels=size_list,labels=size_list)
RMSE_G2_type2$size <- factor(RMSE_G2_type2$size,levels=size_list,labels=size_list)


pl_RMSE_J_type2 <- ggplot(RMSE_J_type2,aes(x=size,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="Size") +
  scale_y_continuous(name="RMSE (Joint)") +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_RMSE_G1_type2 <- ggplot(RMSE_G1_type2,aes(x=size,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="Size") +
  scale_y_continuous(name="RMSE (Group1)") +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_RMSE_G2_type2 <- ggplot(RMSE_G2_type2,aes(x=size,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="Size") +
  scale_y_continuous(name="RMSE (Group2)") +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))



#-----------------------------------------------------------------------------#
# CC_B
#-----------------------------------------------------------------------------#
CC_B_J_type2_d200T100K50 <- data.frame(value=c(unlist(list_type2_d200T100K50$GRIDY$CC_B$Joint),
                                               unlist(list_type2_d200T100K50$SCA_P$CC_B$Joint),
                                               unlist(list_type2_d200T100K50$GICA$CC_B$Joint),
                                               unlist(list_type2_d200T100K50$DSCA$CC_B$Joint),
                                               unlist(list_type2_d200T100K50$DGICA$CC_B$Joint),
                                               unlist(list_type2_d200T100K50$Unfitted_SCA_PF2$CC_B$Joint),
                                               unlist(list_type2_d200T100K50$Unfitted_GICA$CC_B$Joint)))
CC_B_G1_type2_d200T100K50 <- data.frame(value=c(unlist(list_type2_d200T100K50$GRIDY$CC_B$Group1),
                                                unlist(list_type2_d200T100K50$SCA_P$CC_B$Group1),
                                                unlist(list_type2_d200T100K50$GICA$CC_B$Group1),
                                                unlist(list_type2_d200T100K50$DSCA$CC_B$Group1),
                                                unlist(list_type2_d200T100K50$DGICA$CC_B$Group1),
                                                unlist(list_type2_d200T100K50$Unfitted_SCA_PF2$CC_B$Group1),
                                                unlist(list_type2_d200T100K50$Unfitted_GICA$CC_B$Group1)))
CC_B_G2_type2_d200T100K50 <- data.frame(value=c(unlist(list_type2_d200T100K50$GRIDY$CC_B$Group2),
                                                unlist(list_type2_d200T100K50$SCA_P$CC_B$Group2),
                                                unlist(list_type2_d200T100K50$GICA$CC_B$Group2),
                                                unlist(list_type2_d200T100K50$DSCA$CC_B$Group2),
                                                unlist(list_type2_d200T100K50$DGICA$CC_B$Group2),
                                                unlist(list_type2_d200T100K50$Unfitted_SCA_PF2$CC_B$Group2),
                                                unlist(list_type2_d200T100K50$Unfitted_GICA$CC_B$Group2)))
CC_B_type2_d200T100K50 <- rbind(CC_B_J_type2_d200T100K50,CC_B_G1_type2_d200T100K50,CC_B_G2_type2_d200T100K50)
CC_B_type2_d200T100K50$size <- c("d=200")
CC_B_type2_d200T100K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                           "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_B_type2_d200T100K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))


CC_B_J_type2_d400T100K50 <- data.frame(value=c(unlist(list_type2_d400T100K50$GRIDY$CC_B$Joint),
                                               unlist(list_type2_d400T100K50$SCA_P$CC_B$Joint),
                                               unlist(list_type2_d400T100K50$GICA$CC_B$Joint),
                                               unlist(list_type2_d400T100K50$DSCA$CC_B$Joint),
                                               unlist(list_type2_d400T100K50$DGICA$CC_B$Joint),
                                               unlist(list_type2_d400T100K50$Unfitted_SCA_PF2$CC_B$Joint),
                                               unlist(list_type2_d400T100K50$Unfitted_GICA$CC_B$Joint)))
CC_B_G1_type2_d400T100K50 <- data.frame(value=c(unlist(list_type2_d400T100K50$GRIDY$CC_B$Group1),
                                                unlist(list_type2_d400T100K50$SCA_P$CC_B$Group1),
                                                unlist(list_type2_d400T100K50$GICA$CC_B$Group1),
                                                unlist(list_type2_d400T100K50$DSCA$CC_B$Group1),
                                                unlist(list_type2_d400T100K50$DGICA$CC_B$Group1),
                                                unlist(list_type2_d400T100K50$Unfitted_SCA_PF2$CC_B$Group1),
                                                unlist(list_type2_d400T100K50$Unfitted_GICA$CC_B$Group1)))
CC_B_G2_type2_d400T100K50 <- data.frame(value=c(unlist(list_type2_d400T100K50$GRIDY$CC_B$Group2),
                                                unlist(list_type2_d400T100K50$SCA_P$CC_B$Group2),
                                                unlist(list_type2_d400T100K50$GICA$CC_B$Group2),
                                                unlist(list_type2_d400T100K50$DSCA$CC_B$Group2),
                                                unlist(list_type2_d400T100K50$DGICA$CC_B$Group2),
                                                unlist(list_type2_d400T100K50$Unfitted_SCA_PF2$CC_B$Group2),
                                                unlist(list_type2_d400T100K50$Unfitted_GICA$CC_B$Group2)))
CC_B_type2_d400T100K50 <- rbind(CC_B_J_type2_d400T100K50,CC_B_G1_type2_d400T100K50,CC_B_G2_type2_d400T100K50)
CC_B_type2_d400T100K50$size <- c("d=400")
CC_B_type2_d400T100K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                           "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_B_type2_d400T100K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



CC_B_J_type2_d100T200K50 <- data.frame(value=c(unlist(list_type2_d100T200K50$GRIDY$CC_B$Joint),
                                               unlist(list_type2_d100T200K50$SCA_P$CC_B$Joint),
                                               unlist(list_type2_d100T200K50$GICA$CC_B$Joint),
                                               unlist(list_type2_d100T200K50$DSCA$CC_B$Joint),
                                               unlist(list_type2_d100T200K50$DGICA$CC_B$Joint),
                                               unlist(list_type2_d100T200K50$Unfitted_SCA_PF2$CC_B$Joint),
                                               unlist(list_type2_d100T200K50$Unfitted_GICA$CC_B$Joint)))
CC_B_G1_type2_d100T200K50 <- data.frame(value=c(unlist(list_type2_d100T200K50$GRIDY$CC_B$Group1),
                                                unlist(list_type2_d100T200K50$SCA_P$CC_B$Group1),
                                                unlist(list_type2_d100T200K50$GICA$CC_B$Group1),
                                                unlist(list_type2_d100T200K50$DSCA$CC_B$Group1),
                                                unlist(list_type2_d100T200K50$DGICA$CC_B$Group1),
                                                unlist(list_type2_d100T200K50$Unfitted_SCA_PF2$CC_B$Group1),
                                                unlist(list_type2_d100T200K50$Unfitted_GICA$CC_B$Group1)))
CC_B_G2_type2_d100T200K50 <- data.frame(value=c(unlist(list_type2_d100T200K50$GRIDY$CC_B$Group2),
                                                unlist(list_type2_d100T200K50$SCA_P$CC_B$Group2),
                                                unlist(list_type2_d100T200K50$GICA$CC_B$Group2),
                                                unlist(list_type2_d100T200K50$DSCA$CC_B$Group2),
                                                unlist(list_type2_d100T200K50$DGICA$CC_B$Group2),
                                                unlist(list_type2_d100T200K50$Unfitted_SCA_PF2$CC_B$Group2),
                                                unlist(list_type2_d100T200K50$Unfitted_GICA$CC_B$Group2)))
CC_B_type2_d100T200K50 <- rbind(CC_B_J_type2_d100T200K50,CC_B_G1_type2_d100T200K50,CC_B_G2_type2_d100T200K50)
CC_B_type2_d100T200K50$size <- c("T=200")
CC_B_type2_d100T200K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                           "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_B_type2_d100T200K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



CC_B_J_type2_d100T400K50 <- data.frame(value=c(unlist(list_type2_d100T400K50$GRIDY$CC_B$Joint),
                                               unlist(list_type2_d100T400K50$SCA_P$CC_B$Joint),
                                               unlist(list_type2_d100T400K50$GICA$CC_B$Joint),
                                               unlist(list_type2_d100T400K50$DSCA$CC_B$Joint),
                                               unlist(list_type2_d100T400K50$DGICA$CC_B$Joint),
                                               unlist(list_type2_d100T400K50$Unfitted_SCA_PF2$CC_B$Joint),
                                               unlist(list_type2_d100T400K50$Unfitted_GICA$CC_B$Joint)))
CC_B_G1_type2_d100T400K50 <- data.frame(value=c(unlist(list_type2_d100T400K50$GRIDY$CC_B$Group1),
                                                unlist(list_type2_d100T400K50$SCA_P$CC_B$Group1),
                                                unlist(list_type2_d100T400K50$GICA$CC_B$Group1),
                                                unlist(list_type2_d100T400K50$DSCA$CC_B$Group1),
                                                unlist(list_type2_d100T400K50$DGICA$CC_B$Group1),
                                                unlist(list_type2_d100T400K50$Unfitted_SCA_PF2$CC_B$Group1),
                                                unlist(list_type2_d100T400K50$Unfitted_GICA$CC_B$Group1)))
CC_B_G2_type2_d100T400K50 <- data.frame(value=c(unlist(list_type2_d100T400K50$GRIDY$CC_B$Group2),
                                                unlist(list_type2_d100T400K50$SCA_P$CC_B$Group2),
                                                unlist(list_type2_d100T400K50$GICA$CC_B$Group2),
                                                unlist(list_type2_d100T400K50$DSCA$CC_B$Group2),
                                                unlist(list_type2_d100T400K50$DGICA$CC_B$Group2),
                                                unlist(list_type2_d100T400K50$Unfitted_SCA_PF2$CC_B$Group2),
                                                unlist(list_type2_d100T400K50$Unfitted_GICA$CC_B$Group2)))
CC_B_type2_d100T400K50 <- rbind(CC_B_J_type2_d100T400K50,CC_B_G1_type2_d100T400K50,CC_B_G2_type2_d100T400K50)
CC_B_type2_d100T400K50$size <- c("T=400")
CC_B_type2_d100T400K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                           "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_B_type2_d100T400K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



CC_B_J_type2_d100T200K10 <- data.frame(value=c(unlist(list_type2_d100T200K10$GRIDY$CC_B$Joint),
                                               unlist(list_type2_d100T200K10$SCA_P$CC_B$Joint),
                                               unlist(list_type2_d100T200K10$GICA$CC_B$Joint),
                                               unlist(list_type2_d100T200K10$DSCA$CC_B$Joint),
                                               unlist(list_type2_d100T200K10$DGICA$CC_B$Joint),
                                               unlist(list_type2_d100T200K10$Unfitted_SCA_PF2$CC_B$Joint),
                                               unlist(list_type2_d100T200K10$Unfitted_GICA$CC_B$Joint)))
CC_B_G1_type2_d100T200K10 <- data.frame(value=c(unlist(list_type2_d100T200K10$GRIDY$CC_B$Group1),
                                                unlist(list_type2_d100T200K10$SCA_P$CC_B$Group1),
                                                unlist(list_type2_d100T200K10$GICA$CC_B$Group1),
                                                unlist(list_type2_d100T200K10$DSCA$CC_B$Group1),
                                                unlist(list_type2_d100T200K10$DGICA$CC_B$Group1),
                                                unlist(list_type2_d100T200K10$Unfitted_SCA_PF2$CC_B$Group1),
                                                unlist(list_type2_d100T200K10$Unfitted_GICA$CC_B$Group1)))
CC_B_G2_type2_d100T200K10 <- data.frame(value=c(unlist(list_type2_d100T200K10$GRIDY$CC_B$Group2),
                                                unlist(list_type2_d100T200K10$SCA_P$CC_B$Group2),
                                                unlist(list_type2_d100T200K10$GICA$CC_B$Group2),
                                                unlist(list_type2_d100T200K10$DSCA$CC_B$Group2),
                                                unlist(list_type2_d100T200K10$DGICA$CC_B$Group2),
                                                unlist(list_type2_d100T200K10$Unfitted_SCA_PF2$CC_B$Group2),
                                                unlist(list_type2_d100T200K10$Unfitted_GICA$CC_B$Group2)))
CC_B_type2_d100T200K10 <- rbind(CC_B_J_type2_d100T200K10,CC_B_G1_type2_d100T200K10,CC_B_G2_type2_d100T200K10)
CC_B_type2_d100T200K10$size <- c("2K=20")
CC_B_type2_d100T200K10$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                           "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_B_type2_d100T200K10$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))


CC_B_J_type2_d100T200K100 <- data.frame(value=c(unlist(list_type2_d100T200K100$GRIDY$CC_B$Joint),
                                                unlist(list_type2_d100T200K100$SCA_P$CC_B$Joint),
                                                unlist(list_type2_d100T200K100$GICA$CC_B$Joint),
                                                unlist(list_type2_d100T200K100$DSCA$CC_B$Joint),
                                                unlist(list_type2_d100T200K100$DGICA$CC_B$Joint),
                                                unlist(list_type2_d100T200K100$Unfitted_SCA_PF2$CC_B$Joint),
                                                unlist(list_type2_d100T200K100$Unfitted_GICA$CC_B$Joint)))
CC_B_G1_type2_d100T200K100 <- data.frame(value=c(unlist(list_type2_d100T200K100$GRIDY$CC_B$Group1),
                                                 unlist(list_type2_d100T200K100$SCA_P$CC_B$Group1),
                                                 unlist(list_type2_d100T200K100$GICA$CC_B$Group1),
                                                 unlist(list_type2_d100T200K100$DSCA$CC_B$Group1),
                                                 unlist(list_type2_d100T200K100$DGICA$CC_B$Group1),
                                                 unlist(list_type2_d100T200K100$Unfitted_SCA_PF2$CC_B$Group1),
                                                 unlist(list_type2_d100T200K100$Unfitted_GICA$CC_B$Group1)))
CC_B_G2_type2_d100T200K100 <- data.frame(value=c(unlist(list_type2_d100T200K100$GRIDY$CC_B$Group2),
                                                 unlist(list_type2_d100T200K100$SCA_P$CC_B$Group2),
                                                 unlist(list_type2_d100T200K100$GICA$CC_B$Group2),
                                                 unlist(list_type2_d100T200K100$DSCA$CC_B$Group2),
                                                 unlist(list_type2_d100T200K100$DGICA$CC_B$Group2),
                                                 unlist(list_type2_d100T200K100$Unfitted_SCA_PF2$CC_B$Group2),
                                                 unlist(list_type2_d100T200K100$Unfitted_GICA$CC_B$Group2)))
CC_B_type2_d100T200K100 <- rbind(CC_B_J_type2_d100T200K100,CC_B_G1_type2_d100T200K100,CC_B_G2_type2_d100T200K100)
CC_B_type2_d100T200K100$size <- c("2K=200")
CC_B_type2_d100T200K100$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                            "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_B_type2_d100T200K100$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



CC_B_J_type2 <- rbind(CC_B_type2_d200T100K50[which(CC_B_type2_d200T100K50$structure == "Joint"),],
                      CC_B_type2_d400T100K50[which(CC_B_type2_d400T100K50$structure == "Joint"),],
                      CC_B_type2_d100T200K50[which(CC_B_type2_d100T200K50$structure == "Joint"),],
                      CC_B_type2_d100T400K50[which(CC_B_type2_d100T200K50$structure == "Joint"),],
                      CC_B_type2_d100T200K10[which(CC_B_type2_d100T200K10$structure == "Joint"),],
                      CC_B_type2_d100T200K100[which(CC_B_type2_d100T200K100$structure == "Joint"),])
CC_B_G1_type2 <- rbind(CC_B_type2_d200T100K50[which(CC_B_type2_d200T100K50$structure == "Group Individual1"),],
                       CC_B_type2_d400T100K50[which(CC_B_type2_d400T100K50$structure == "Group Individual1"),],
                       CC_B_type2_d100T200K50[which(CC_B_type2_d100T200K50$structure == "Group Individual1"),],
                       CC_B_type2_d100T400K50[which(CC_B_type2_d100T200K50$structure == "Group Individual1"),],
                       CC_B_type2_d100T200K10[which(CC_B_type2_d100T200K10$structure == "Group Individual1"),],
                       CC_B_type2_d100T200K100[which(CC_B_type2_d100T200K100$structure == "Group Individual1"),])
CC_B_G2_type2 <- rbind(CC_B_type2_d200T100K50[which(CC_B_type2_d200T100K50$structure == "Group Individual2"),],
                       CC_B_type2_d400T100K50[which(CC_B_type2_d400T100K50$structure == "Group Individual2"),],
                       CC_B_type2_d100T200K50[which(CC_B_type2_d100T200K50$structure == "Group Individual2"),],
                       CC_B_type2_d100T400K50[which(CC_B_type2_d100T200K50$structure == "Group Individual2"),],
                       CC_B_type2_d100T200K10[which(CC_B_type2_d100T200K10$structure == "Group Individual2"),],
                       CC_B_type2_d100T200K100[which(CC_B_type2_d100T200K100$structure == "Group Individual2"),])


method_list <- c("GRIDY","SCA-P","GICA","DSCA","DGICA","Unfitted_SCA-PF2","Unfitted_GICA")
CC_B_J_type2$method <- factor(CC_B_J_type2$method,levels=method_list,labels=method_list)
CC_B_G1_type2$method <- factor(CC_B_G1_type2$method,levels=method_list,labels=method_list)
CC_B_G2_type2$method <- factor(CC_B_G2_type2$method,levels=method_list,labels=method_list)

size_list <- c("d=200","d=400","T=200","T=400","2K=20","2K=200")
CC_B_J_type2$size<- factor(CC_B_J_type2$size,levels=size_list,labels=size_list)
CC_B_G1_type2$size <- factor(CC_B_G1_type2$size,levels=size_list,labels=size_list)
CC_B_G2_type2$size <- factor(CC_B_G2_type2$size,levels=size_list,labels=size_list)


pl_CC_B_J_type2 <- ggplot(CC_B_J_type2,aes(x=size,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="Size") +
  scale_y_continuous(name=TeX("$CC_B$ (Joint)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_B_G1_type2 <- ggplot(CC_B_G1_type2,aes(x=size,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="Size") +
  scale_y_continuous(name=TeX("$CC_B$ (Group 1)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_B_G2_type2 <- ggplot(CC_B_G2_type2,aes(x=size,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="Size") +
  scale_y_continuous(name=TeX("$CC_B$ (Group 2)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))



# -----------------------------------------------------------------------------#
# CC_F
# -----------------------------------------------------------------------------#
CC_F_J_type2_d200T100K50 <- data.frame(value=c(unlist(list_type2_d200T100K50$GRIDY$CC_F$Joint),
                                               unlist(list_type2_d200T100K50$SCA_P$CC_F$Joint),
                                               unlist(list_type2_d200T100K50$GICA$CC_F$Joint),
                                               unlist(list_type2_d200T100K50$DSCA$CC_F$Joint),
                                               unlist(list_type2_d200T100K50$DGICA$CC_F$Joint),
                                               unlist(list_type2_d200T100K50$Unfitted_SCA_PF2$CC_F$Joint),
                                               unlist(list_type2_d200T100K50$Unfitted_GICA$CC_F$Joint)))
CC_F_G1_type2_d200T100K50 <- data.frame(value=c(unlist(list_type2_d200T100K50$GRIDY$CC_F$Group1),
                                                unlist(list_type2_d200T100K50$SCA_P$CC_F$Group1),
                                                unlist(list_type2_d200T100K50$GICA$CC_F$Group1),
                                                unlist(list_type2_d200T100K50$DSCA$CC_F$Group1),
                                                unlist(list_type2_d200T100K50$DGICA$CC_F$Group1),
                                                unlist(list_type2_d200T100K50$Unfitted_SCA_PF2$CC_F$Group1),
                                                unlist(list_type2_d200T100K50$Unfitted_GICA$CC_F$Group1)))
CC_F_G2_type2_d200T100K50 <- data.frame(value=c(unlist(list_type2_d200T100K50$GRIDY$CC_F$Group2),
                                                unlist(list_type2_d200T100K50$SCA_P$CC_F$Group2),
                                                unlist(list_type2_d200T100K50$GICA$CC_F$Group2),
                                                unlist(list_type2_d200T100K50$DSCA$CC_F$Group2),
                                                unlist(list_type2_d200T100K50$DGICA$CC_F$Group2),
                                                unlist(list_type2_d200T100K50$Unfitted_SCA_PF2$CC_F$Group2),
                                                unlist(list_type2_d200T100K50$Unfitted_GICA$CC_F$Group2)))
CC_F_type2_d200T100K50 <- rbind(CC_F_J_type2_d200T100K50,CC_F_G1_type2_d200T100K50,CC_F_G2_type2_d200T100K50)
CC_F_type2_d200T100K50$size <- c("d=200")
CC_F_type2_d200T100K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                           "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_F_type2_d200T100K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))


CC_F_J_type2_d400T100K50 <- data.frame(value=c(unlist(list_type2_d400T100K50$GRIDY$CC_F$Joint),
                                               unlist(list_type2_d400T100K50$SCA_P$CC_F$Joint),
                                               unlist(list_type2_d400T100K50$GICA$CC_F$Joint),
                                               unlist(list_type2_d400T100K50$DSCA$CC_F$Joint),
                                               unlist(list_type2_d400T100K50$DGICA$CC_F$Joint),
                                               unlist(list_type2_d400T100K50$Unfitted_SCA_PF2$CC_F$Joint),
                                               unlist(list_type2_d400T100K50$Unfitted_GICA$CC_F$Joint)))
CC_F_G1_type2_d400T100K50 <- data.frame(value=c(unlist(list_type2_d400T100K50$GRIDY$CC_F$Group1),
                                                unlist(list_type2_d400T100K50$SCA_P$CC_F$Group1),
                                                unlist(list_type2_d400T100K50$GICA$CC_F$Group1),
                                                unlist(list_type2_d400T100K50$DSCA$CC_F$Group1),
                                                unlist(list_type2_d400T100K50$DGICA$CC_F$Group1),
                                                unlist(list_type2_d400T100K50$Unfitted_SCA_PF2$CC_F$Group1),
                                                unlist(list_type2_d400T100K50$Unfitted_GICA$CC_F$Group1)))
CC_F_G2_type2_d400T100K50 <- data.frame(value=c(unlist(list_type2_d400T100K50$GRIDY$CC_F$Group2),
                                                unlist(list_type2_d400T100K50$SCA_P$CC_F$Group2),
                                                unlist(list_type2_d400T100K50$GICA$CC_F$Group2),
                                                unlist(list_type2_d400T100K50$DSCA$CC_F$Group2),
                                                unlist(list_type2_d400T100K50$DGICA$CC_F$Group2),
                                                unlist(list_type2_d400T100K50$Unfitted_SCA_PF2$CC_F$Group2),
                                                unlist(list_type2_d400T100K50$Unfitted_GICA$CC_F$Group2)))
CC_F_type2_d400T100K50 <- rbind(CC_F_J_type2_d400T100K50,CC_F_G1_type2_d400T100K50,CC_F_G2_type2_d400T100K50)
CC_F_type2_d400T100K50$size <- c("d=400")
CC_F_type2_d400T100K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                           "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_F_type2_d400T100K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



CC_F_J_type2_d100T200K50 <- data.frame(value=c(unlist(list_type2_d100T200K50$GRIDY$CC_F$Joint),
                                               unlist(list_type2_d100T200K50$SCA_P$CC_F$Joint),
                                               unlist(list_type2_d100T200K50$GICA$CC_F$Joint),
                                               unlist(list_type2_d100T200K50$DSCA$CC_F$Joint),
                                               unlist(list_type2_d100T200K50$DGICA$CC_F$Joint),
                                               unlist(list_type2_d100T200K50$Unfitted_SCA_PF2$CC_F$Joint),
                                               unlist(list_type2_d100T200K50$Unfitted_GICA$CC_F$Joint)))
CC_F_G1_type2_d100T200K50 <- data.frame(value=c(unlist(list_type2_d100T200K50$GRIDY$CC_F$Group1),
                                                unlist(list_type2_d100T200K50$SCA_P$CC_F$Group1),
                                                unlist(list_type2_d100T200K50$GICA$CC_F$Group1),
                                                unlist(list_type2_d100T200K50$DSCA$CC_F$Group1),
                                                unlist(list_type2_d100T200K50$DGICA$CC_F$Group1),
                                                unlist(list_type2_d100T200K50$Unfitted_SCA_PF2$CC_F$Group1),
                                                unlist(list_type2_d100T200K50$Unfitted_GICA$CC_F$Group1)))
CC_F_G2_type2_d100T200K50 <- data.frame(value=c(unlist(list_type2_d100T200K50$GRIDY$CC_F$Group2),
                                                unlist(list_type2_d100T200K50$SCA_P$CC_F$Group2),
                                                unlist(list_type2_d100T200K50$GICA$CC_F$Group2),
                                                unlist(list_type2_d100T200K50$DSCA$CC_F$Group2),
                                                unlist(list_type2_d100T200K50$DGICA$CC_F$Group2),
                                                unlist(list_type2_d100T200K50$Unfitted_SCA_PF2$CC_F$Group2),
                                                unlist(list_type2_d100T200K50$Unfitted_GICA$CC_F$Group2)))
CC_F_type2_d100T200K50 <- rbind(CC_F_J_type2_d100T200K50,CC_F_G1_type2_d100T200K50,CC_F_G2_type2_d100T200K50)
CC_F_type2_d100T200K50$size <- c("T=200")
CC_F_type2_d100T200K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                           "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_F_type2_d100T200K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



CC_F_J_type2_d100T400K50 <- data.frame(value=c(unlist(list_type2_d100T400K50$GRIDY$CC_F$Joint),
                                               unlist(list_type2_d100T400K50$SCA_P$CC_F$Joint),
                                               unlist(list_type2_d100T400K50$GICA$CC_F$Joint),
                                               unlist(list_type2_d100T400K50$DSCA$CC_F$Joint),
                                               unlist(list_type2_d100T400K50$DGICA$CC_F$Joint),
                                               unlist(list_type2_d100T400K50$Unfitted_SCA_PF2$CC_F$Joint),
                                               unlist(list_type2_d100T400K50$Unfitted_GICA$CC_F$Joint)))
CC_F_G1_type2_d100T400K50 <- data.frame(value=c(unlist(list_type2_d100T400K50$GRIDY$CC_F$Group1),
                                                unlist(list_type2_d100T400K50$SCA_P$CC_F$Group1),
                                                unlist(list_type2_d100T400K50$GICA$CC_F$Group1),
                                                unlist(list_type2_d100T400K50$DSCA$CC_F$Group1),
                                                unlist(list_type2_d100T400K50$DGICA$CC_F$Group1),
                                                unlist(list_type2_d100T400K50$Unfitted_SCA_PF2$CC_F$Group1),
                                                unlist(list_type2_d100T400K50$Unfitted_GICA$CC_F$Group1)))
CC_F_G2_type2_d100T400K50 <- data.frame(value=c(unlist(list_type2_d100T400K50$GRIDY$CC_F$Group2),
                                                unlist(list_type2_d100T400K50$SCA_P$CC_F$Group2),
                                                unlist(list_type2_d100T400K50$GICA$CC_F$Group2),
                                                unlist(list_type2_d100T400K50$DSCA$CC_F$Group2),
                                                unlist(list_type2_d100T400K50$DGICA$CC_F$Group2),
                                                unlist(list_type2_d100T400K50$Unfitted_SCA_PF2$CC_F$Group2),
                                                unlist(list_type2_d100T400K50$Unfitted_GICA$CC_F$Group2)))
CC_F_type2_d100T400K50 <- rbind(CC_F_J_type2_d100T400K50,CC_F_G1_type2_d100T400K50,CC_F_G2_type2_d100T400K50)
CC_F_type2_d100T400K50$size <- c("T=400")
CC_F_type2_d100T400K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                           "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_F_type2_d100T400K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



CC_F_J_type2_d100T200K10 <- data.frame(value=c(unlist(list_type2_d100T200K10$GRIDY$CC_F$Joint),
                                               unlist(list_type2_d100T200K10$SCA_P$CC_F$Joint),
                                               unlist(list_type2_d100T200K10$GICA$CC_F$Joint),
                                               unlist(list_type2_d100T200K10$DSCA$CC_F$Joint),
                                               unlist(list_type2_d100T200K10$DGICA$CC_F$Joint),
                                               unlist(list_type2_d100T200K10$Unfitted_SCA_PF2$CC_F$Joint),
                                               unlist(list_type2_d100T200K10$Unfitted_GICA$CC_F$Joint)))
CC_F_G1_type2_d100T200K10 <- data.frame(value=c(unlist(list_type2_d100T200K10$GRIDY$CC_F$Group1),
                                                unlist(list_type2_d100T200K10$SCA_P$CC_F$Group1),
                                                unlist(list_type2_d100T200K10$GICA$CC_F$Group1),
                                                unlist(list_type2_d100T200K10$DSCA$CC_F$Group1),
                                                unlist(list_type2_d100T200K10$DGICA$CC_F$Group1),
                                                unlist(list_type2_d100T200K10$Unfitted_SCA_PF2$CC_F$Group1),
                                                unlist(list_type2_d100T200K10$Unfitted_GICA$CC_F$Group1)))
CC_F_G2_type2_d100T200K10 <- data.frame(value=c(unlist(list_type2_d100T200K10$GRIDY$CC_F$Group2),
                                                unlist(list_type2_d100T200K10$SCA_P$CC_F$Group2),
                                                unlist(list_type2_d100T200K10$GICA$CC_F$Group2),
                                                unlist(list_type2_d100T200K10$DSCA$CC_F$Group2),
                                                unlist(list_type2_d100T200K10$DGICA$CC_F$Group2),
                                                unlist(list_type2_d100T200K10$Unfitted_SCA_PF2$CC_F$Group2),
                                                unlist(list_type2_d100T200K10$Unfitted_GICA$CC_F$Group2)))
CC_F_type2_d100T200K10 <- rbind(CC_F_J_type2_d100T200K10,CC_F_G1_type2_d100T200K10,CC_F_G2_type2_d100T200K10)
CC_F_type2_d100T200K10$size <- c("2K=20")
CC_F_type2_d100T200K10$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                           "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_F_type2_d100T200K10$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))


CC_F_J_type2_d100T200K100 <- data.frame(value=c(unlist(list_type2_d100T200K100$GRIDY$CC_F$Joint),
                                                unlist(list_type2_d100T200K100$SCA_P$CC_F$Joint),
                                                unlist(list_type2_d100T200K100$GICA$CC_F$Joint),
                                                unlist(list_type2_d100T200K100$DSCA$CC_F$Joint),
                                                unlist(list_type2_d100T200K100$DGICA$CC_F$Joint),
                                                unlist(list_type2_d100T200K100$Unfitted_SCA_PF2$CC_F$Joint),
                                                unlist(list_type2_d100T200K100$Unfitted_GICA$CC_F$Joint)))
CC_F_G1_type2_d100T200K100 <- data.frame(value=c(unlist(list_type2_d100T200K100$GRIDY$CC_F$Group1),
                                                 unlist(list_type2_d100T200K100$SCA_P$CC_F$Group1),
                                                 unlist(list_type2_d100T200K100$GICA$CC_F$Group1),
                                                 unlist(list_type2_d100T200K100$DSCA$CC_F$Group1),
                                                 unlist(list_type2_d100T200K100$DGICA$CC_F$Group1),
                                                 unlist(list_type2_d100T200K100$Unfitted_SCA_PF2$CC_F$Group1),
                                                 unlist(list_type2_d100T200K100$Unfitted_GICA$CC_F$Group1)))
CC_F_G2_type2_d100T200K100 <- data.frame(value=c(unlist(list_type2_d100T200K100$GRIDY$CC_F$Group2),
                                                 unlist(list_type2_d100T200K100$SCA_P$CC_F$Group2),
                                                 unlist(list_type2_d100T200K100$GICA$CC_F$Group2),
                                                 unlist(list_type2_d100T200K100$DSCA$CC_F$Group2),
                                                 unlist(list_type2_d100T200K100$DGICA$CC_F$Group2),
                                                 unlist(list_type2_d100T200K100$Unfitted_SCA_PF2$CC_F$Group2),
                                                 unlist(list_type2_d100T200K100$Unfitted_GICA$CC_F$Group2)))
CC_F_type2_d100T200K100 <- rbind(CC_F_J_type2_d100T200K100,CC_F_G1_type2_d100T200K100,CC_F_G2_type2_d100T200K100)
CC_F_type2_d100T200K100$size <- c("2K=200")
CC_F_type2_d100T200K100$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                            "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_F_type2_d100T200K100$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



CC_F_J_type2 <- rbind(CC_F_type2_d200T100K50[which(CC_F_type2_d200T100K50$structure == "Joint"),],
                      CC_F_type2_d400T100K50[which(CC_F_type2_d400T100K50$structure == "Joint"),],
                      CC_F_type2_d100T200K50[which(CC_F_type2_d100T200K50$structure == "Joint"),],
                      CC_F_type2_d100T400K50[which(CC_F_type2_d100T200K50$structure == "Joint"),],
                      CC_F_type2_d100T200K10[which(CC_F_type2_d100T200K10$structure == "Joint"),],
                      CC_F_type2_d100T200K100[which(CC_F_type2_d100T200K100$structure == "Joint"),])
CC_F_G1_type2 <- rbind(CC_F_type2_d200T100K50[which(CC_F_type2_d200T100K50$structure == "Group Individual1"),],
                       CC_F_type2_d400T100K50[which(CC_F_type2_d400T100K50$structure == "Group Individual1"),],
                       CC_F_type2_d100T200K50[which(CC_F_type2_d100T200K50$structure == "Group Individual1"),],
                       CC_F_type2_d100T400K50[which(CC_F_type2_d100T200K50$structure == "Group Individual1"),],
                       CC_F_type2_d100T200K10[which(CC_F_type2_d100T200K10$structure == "Group Individual1"),],
                       CC_F_type2_d100T200K100[which(CC_F_type2_d100T200K100$structure == "Group Individual1"),])
CC_F_G2_type2 <- rbind(CC_F_type2_d200T100K50[which(CC_F_type2_d200T100K50$structure == "Group Individual2"),],
                       CC_F_type2_d400T100K50[which(CC_F_type2_d400T100K50$structure == "Group Individual2"),],
                       CC_F_type2_d100T200K50[which(CC_F_type2_d100T200K50$structure == "Group Individual2"),],
                       CC_F_type2_d100T400K50[which(CC_F_type2_d100T200K50$structure == "Group Individual2"),],
                       CC_F_type2_d100T200K10[which(CC_F_type2_d100T200K10$structure == "Group Individual2"),],
                       CC_F_type2_d100T200K100[which(CC_F_type2_d100T200K100$structure == "Group Individual2"),])


method_list <- c("GRIDY","SCA-P","GICA","DSCA","DGICA","Unfitted_SCA-PF2","Unfitted_GICA")
CC_F_J_type2$method <- factor(CC_F_J_type2$method,levels=method_list,labels=method_list)
CC_F_G1_type2$method <- factor(CC_F_G1_type2$method,levels=method_list,labels=method_list)
CC_F_G2_type2$method <- factor(CC_F_G2_type2$method,levels=method_list,labels=method_list)

size_list <- c("d=200","d=400","T=200","T=400","2K=20","2K=200")
CC_F_J_type2$size<- factor(CC_F_J_type2$size,levels=size_list,labels=size_list)
CC_F_G1_type2$size <- factor(CC_F_G1_type2$size,levels=size_list,labels=size_list)
CC_F_G2_type2$size <- factor(CC_F_G2_type2$size,levels=size_list,labels=size_list)


pl_CC_F_J_type2 <- ggplot(CC_F_J_type2,aes(x=size,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="Size") +
  scale_y_continuous(name=TeX("$CC_F$ (Joint)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_F_G1_type2 <- ggplot(CC_F_G1_type2,aes(x=size,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="Size") +
  scale_y_continuous(name=TeX("$CC_F$ (Group 1)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_F_G2_type2 <- ggplot(CC_F_G2_type2,aes(x=size,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="Size") +
  scale_y_continuous(name=TeX("$CC_F$ (Group 2)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))



# -----------------------------------------------------------------------------#
# CC_PSI
# -----------------------------------------------------------------------------#
CC_PSI_J_type2_d200T100K50 <- data.frame(value=c(unlist(list_type2_d200T100K50$GRIDY$CC_PSI$Joint),
                                                 unlist(list_type2_d200T100K50$SCA_P$CC_PSI$Joint),
                                                 unlist(list_type2_d200T100K50$GICA$CC_PSI$Joint),
                                                 unlist(list_type2_d200T100K50$DSCA$CC_PSI$Joint),
                                                 unlist(list_type2_d200T100K50$DGICA$CC_PSI$Joint),
                                                 unlist(list_type2_d200T100K50$Unfitted_SCA_PF2$CC_PSI$Joint),
                                                 unlist(list_type2_d200T100K50$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_G1_type2_d200T100K50 <- data.frame(value=c(unlist(list_type2_d200T100K50$GRIDY$CC_PSI$Group1),
                                                  unlist(list_type2_d200T100K50$SCA_P$CC_PSI$Group1),
                                                  unlist(list_type2_d200T100K50$GICA$CC_PSI$Group1),
                                                  unlist(list_type2_d200T100K50$DSCA$CC_PSI$Group1),
                                                  unlist(list_type2_d200T100K50$DGICA$CC_PSI$Group1),
                                                  unlist(list_type2_d200T100K50$Unfitted_SCA_PF2$CC_PSI$Group1),
                                                  unlist(list_type2_d200T100K50$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G2_type2_d200T100K50 <- data.frame(value=c(unlist(list_type2_d200T100K50$GRIDY$CC_PSI$Group2),
                                                  unlist(list_type2_d200T100K50$SCA_P$CC_PSI$Group2),
                                                  unlist(list_type2_d200T100K50$GICA$CC_PSI$Group2),
                                                  unlist(list_type2_d200T100K50$DSCA$CC_PSI$Group2),
                                                  unlist(list_type2_d200T100K50$DGICA$CC_PSI$Group2),
                                                  unlist(list_type2_d200T100K50$Unfitted_SCA_PF2$CC_PSI$Group2),
                                                  unlist(list_type2_d200T100K50$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_type2_d200T100K50 <- rbind(CC_PSI_J_type2_d200T100K50,CC_PSI_G1_type2_d200T100K50,CC_PSI_G2_type2_d200T100K50)
CC_PSI_type2_d200T100K50$size <- c("d=200")
CC_PSI_type2_d200T100K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                             "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_PSI_type2_d200T100K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))


CC_PSI_J_type2_d400T100K50 <- data.frame(value=c(unlist(list_type2_d400T100K50$GRIDY$CC_PSI$Joint),
                                                 unlist(list_type2_d400T100K50$SCA_P$CC_PSI$Joint),
                                                 unlist(list_type2_d400T100K50$GICA$CC_PSI$Joint),
                                                 unlist(list_type2_d400T100K50$DSCA$CC_PSI$Joint),
                                                 unlist(list_type2_d400T100K50$DGICA$CC_PSI$Joint),
                                                 unlist(list_type2_d400T100K50$Unfitted_SCA_PF2$CC_PSI$Joint),
                                                 unlist(list_type2_d400T100K50$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_G1_type2_d400T100K50 <- data.frame(value=c(unlist(list_type2_d400T100K50$GRIDY$CC_PSI$Group1),
                                                  unlist(list_type2_d400T100K50$SCA_P$CC_PSI$Group1),
                                                  unlist(list_type2_d400T100K50$GICA$CC_PSI$Group1),
                                                  unlist(list_type2_d400T100K50$DSCA$CC_PSI$Group1),
                                                  unlist(list_type2_d400T100K50$DGICA$CC_PSI$Group1),
                                                  unlist(list_type2_d400T100K50$Unfitted_SCA_PF2$CC_PSI$Group1),
                                                  unlist(list_type2_d400T100K50$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G2_type2_d400T100K50 <- data.frame(value=c(unlist(list_type2_d400T100K50$GRIDY$CC_PSI$Group2),
                                                  unlist(list_type2_d400T100K50$SCA_P$CC_PSI$Group2),
                                                  unlist(list_type2_d400T100K50$GICA$CC_PSI$Group2),
                                                  unlist(list_type2_d400T100K50$DSCA$CC_PSI$Group2),
                                                  unlist(list_type2_d400T100K50$DGICA$CC_PSI$Group2),
                                                  unlist(list_type2_d400T100K50$Unfitted_SCA_PF2$CC_PSI$Group2),
                                                  unlist(list_type2_d400T100K50$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_type2_d400T100K50 <- rbind(CC_PSI_J_type2_d400T100K50,CC_PSI_G1_type2_d400T100K50,CC_PSI_G2_type2_d400T100K50)
CC_PSI_type2_d400T100K50$size <- c("d=400")
CC_PSI_type2_d400T100K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                             "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_PSI_type2_d400T100K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



CC_PSI_J_type2_d100T200K50 <- data.frame(value=c(unlist(list_type2_d100T200K50$GRIDY$CC_PSI$Joint),
                                                 unlist(list_type2_d100T200K50$SCA_P$CC_PSI$Joint),
                                                 unlist(list_type2_d100T200K50$GICA$CC_PSI$Joint),
                                                 unlist(list_type2_d100T200K50$DSCA$CC_PSI$Joint),
                                                 unlist(list_type2_d100T200K50$DGICA$CC_PSI$Joint),
                                                 unlist(list_type2_d100T200K50$Unfitted_SCA_PF2$CC_PSI$Joint),
                                                 unlist(list_type2_d100T200K50$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_G1_type2_d100T200K50 <- data.frame(value=c(unlist(list_type2_d100T200K50$GRIDY$CC_PSI$Group1),
                                                  unlist(list_type2_d100T200K50$SCA_P$CC_PSI$Group1),
                                                  unlist(list_type2_d100T200K50$GICA$CC_PSI$Group1),
                                                  unlist(list_type2_d100T200K50$DSCA$CC_PSI$Group1),
                                                  unlist(list_type2_d100T200K50$DGICA$CC_PSI$Group1),
                                                  unlist(list_type2_d100T200K50$Unfitted_SCA_PF2$CC_PSI$Group1),
                                                  unlist(list_type2_d100T200K50$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G2_type2_d100T200K50 <- data.frame(value=c(unlist(list_type2_d100T200K50$GRIDY$CC_PSI$Group2),
                                                  unlist(list_type2_d100T200K50$SCA_P$CC_PSI$Group2),
                                                  unlist(list_type2_d100T200K50$GICA$CC_PSI$Group2),
                                                  unlist(list_type2_d100T200K50$DSCA$CC_PSI$Group2),
                                                  unlist(list_type2_d100T200K50$DGICA$CC_PSI$Group2),
                                                  unlist(list_type2_d100T200K50$Unfitted_SCA_PF2$CC_PSI$Group2),
                                                  unlist(list_type2_d100T200K50$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_type2_d100T200K50 <- rbind(CC_PSI_J_type2_d100T200K50,CC_PSI_G1_type2_d100T200K50,CC_PSI_G2_type2_d100T200K50)
CC_PSI_type2_d100T200K50$size <- c("T=200")
CC_PSI_type2_d100T200K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                             "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_PSI_type2_d100T200K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



CC_PSI_J_type2_d100T400K50 <- data.frame(value=c(unlist(list_type2_d100T400K50$GRIDY$CC_PSI$Joint),
                                                 unlist(list_type2_d100T400K50$SCA_P$CC_PSI$Joint),
                                                 unlist(list_type2_d100T400K50$GICA$CC_PSI$Joint),
                                                 unlist(list_type2_d100T400K50$DSCA$CC_PSI$Joint),
                                                 unlist(list_type2_d100T400K50$DGICA$CC_PSI$Joint),
                                                 unlist(list_type2_d100T400K50$Unfitted_SCA_PF2$CC_PSI$Joint),
                                                 unlist(list_type2_d100T400K50$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_G1_type2_d100T400K50 <- data.frame(value=c(unlist(list_type2_d100T400K50$GRIDY$CC_PSI$Group1),
                                                  unlist(list_type2_d100T400K50$SCA_P$CC_PSI$Group1),
                                                  unlist(list_type2_d100T400K50$GICA$CC_PSI$Group1),
                                                  unlist(list_type2_d100T400K50$DSCA$CC_PSI$Group1),
                                                  unlist(list_type2_d100T400K50$DGICA$CC_PSI$Group1),
                                                  unlist(list_type2_d100T400K50$Unfitted_SCA_PF2$CC_PSI$Group1),
                                                  unlist(list_type2_d100T400K50$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G2_type2_d100T400K50 <- data.frame(value=c(unlist(list_type2_d100T400K50$GRIDY$CC_PSI$Group2),
                                                  unlist(list_type2_d100T400K50$SCA_P$CC_PSI$Group2),
                                                  unlist(list_type2_d100T400K50$GICA$CC_PSI$Group2),
                                                  unlist(list_type2_d100T400K50$DSCA$CC_PSI$Group2),
                                                  unlist(list_type2_d100T400K50$DGICA$CC_PSI$Group2),
                                                  unlist(list_type2_d100T400K50$Unfitted_SCA_PF2$CC_PSI$Group2),
                                                  unlist(list_type2_d100T400K50$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_type2_d100T400K50 <- rbind(CC_PSI_J_type2_d100T400K50,CC_PSI_G1_type2_d100T400K50,CC_PSI_G2_type2_d100T400K50)
CC_PSI_type2_d100T400K50$size <- c("T=400")
CC_PSI_type2_d100T400K50$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                             "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_PSI_type2_d100T400K50$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



CC_PSI_J_type2_d100T200K10 <- data.frame(value=c(unlist(list_type2_d100T200K10$GRIDY$CC_PSI$Joint),
                                                 unlist(list_type2_d100T200K10$SCA_P$CC_PSI$Joint),
                                                 unlist(list_type2_d100T200K10$GICA$CC_PSI$Joint),
                                                 unlist(list_type2_d100T200K10$DSCA$CC_PSI$Joint),
                                                 unlist(list_type2_d100T200K10$DGICA$CC_PSI$Joint),
                                                 unlist(list_type2_d100T200K10$Unfitted_SCA_PF2$CC_PSI$Joint),
                                                 unlist(list_type2_d100T200K10$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_G1_type2_d100T200K10 <- data.frame(value=c(unlist(list_type2_d100T200K10$GRIDY$CC_PSI$Group1),
                                                  unlist(list_type2_d100T200K10$SCA_P$CC_PSI$Group1),
                                                  unlist(list_type2_d100T200K10$GICA$CC_PSI$Group1),
                                                  unlist(list_type2_d100T200K10$DSCA$CC_PSI$Group1),
                                                  unlist(list_type2_d100T200K10$DGICA$CC_PSI$Group1),
                                                  unlist(list_type2_d100T200K10$Unfitted_SCA_PF2$CC_PSI$Group1),
                                                  unlist(list_type2_d100T200K10$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G2_type2_d100T200K10 <- data.frame(value=c(unlist(list_type2_d100T200K10$GRIDY$CC_PSI$Group2),
                                                  unlist(list_type2_d100T200K10$SCA_P$CC_PSI$Group2),
                                                  unlist(list_type2_d100T200K10$GICA$CC_PSI$Group2),
                                                  unlist(list_type2_d100T200K10$DSCA$CC_PSI$Group2),
                                                  unlist(list_type2_d100T200K10$DGICA$CC_PSI$Group2),
                                                  unlist(list_type2_d100T200K10$Unfitted_SCA_PF2$CC_PSI$Group2),
                                                  unlist(list_type2_d100T200K10$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_type2_d100T200K10 <- rbind(CC_PSI_J_type2_d100T200K10,CC_PSI_G1_type2_d100T200K10,CC_PSI_G2_type2_d100T200K10)
CC_PSI_type2_d100T200K10$size <- c("2K=20")
CC_PSI_type2_d100T200K10$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                             "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_PSI_type2_d100T200K10$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))


CC_PSI_J_type2_d100T200K100 <- data.frame(value=c(unlist(list_type2_d100T200K100$GRIDY$CC_PSI$Joint),
                                                  unlist(list_type2_d100T200K100$SCA_P$CC_PSI$Joint),
                                                  unlist(list_type2_d100T200K100$GICA$CC_PSI$Joint),
                                                  unlist(list_type2_d100T200K100$DSCA$CC_PSI$Joint),
                                                  unlist(list_type2_d100T200K100$DGICA$CC_PSI$Joint),
                                                  unlist(list_type2_d100T200K100$Unfitted_SCA_PF2$CC_PSI$Joint),
                                                  unlist(list_type2_d100T200K100$Unfitted_GICA$CC_PSI$Joint)))
CC_PSI_G1_type2_d100T200K100 <- data.frame(value=c(unlist(list_type2_d100T200K100$GRIDY$CC_PSI$Group1),
                                                   unlist(list_type2_d100T200K100$SCA_P$CC_PSI$Group1),
                                                   unlist(list_type2_d100T200K100$GICA$CC_PSI$Group1),
                                                   unlist(list_type2_d100T200K100$DSCA$CC_PSI$Group1),
                                                   unlist(list_type2_d100T200K100$DGICA$CC_PSI$Group1),
                                                   unlist(list_type2_d100T200K100$Unfitted_SCA_PF2$CC_PSI$Group1),
                                                   unlist(list_type2_d100T200K100$Unfitted_GICA$CC_PSI$Group1)))
CC_PSI_G2_type2_d100T200K100 <- data.frame(value=c(unlist(list_type2_d100T200K100$GRIDY$CC_PSI$Group2),
                                                   unlist(list_type2_d100T200K100$SCA_P$CC_PSI$Group2),
                                                   unlist(list_type2_d100T200K100$GICA$CC_PSI$Group2),
                                                   unlist(list_type2_d100T200K100$DSCA$CC_PSI$Group2),
                                                   unlist(list_type2_d100T200K100$DGICA$CC_PSI$Group2),
                                                   unlist(list_type2_d100T200K100$Unfitted_SCA_PF2$CC_PSI$Group2),
                                                   unlist(list_type2_d100T200K100$Unfitted_GICA$CC_PSI$Group2)))
CC_PSI_type2_d100T200K100 <- rbind(CC_PSI_J_type2_d100T200K100,CC_PSI_G1_type2_d100T200K100,CC_PSI_G2_type2_d100T200K100)
CC_PSI_type2_d100T200K100$size <- c("2K=200")
CC_PSI_type2_d100T200K100$method <- rep(rep(c("GRIDY","SCA-P","GICA","DSCA","DGICA",
                                              "Unfitted_SCA-PF2","Unfitted_GICA"),each=100),3)
CC_PSI_type2_d100T200K100$structure <- c(rep("Joint",700),rep("Group Individual1",700),rep("Group Individual2",700))



CC_PSI_J_type2 <- rbind(CC_PSI_type2_d200T100K50[which(CC_PSI_type2_d200T100K50$structure == "Joint"),],
                        CC_PSI_type2_d400T100K50[which(CC_PSI_type2_d400T100K50$structure == "Joint"),],
                        CC_PSI_type2_d100T200K50[which(CC_PSI_type2_d100T200K50$structure == "Joint"),],
                        CC_PSI_type2_d100T400K50[which(CC_PSI_type2_d100T200K50$structure == "Joint"),],
                        CC_PSI_type2_d100T200K10[which(CC_PSI_type2_d100T200K10$structure == "Joint"),],
                        CC_PSI_type2_d100T200K100[which(CC_PSI_type2_d100T200K100$structure == "Joint"),])
CC_PSI_G1_type2 <- rbind(CC_PSI_type2_d200T100K50[which(CC_PSI_type2_d200T100K50$structure == "Group Individual1"),],
                         CC_PSI_type2_d400T100K50[which(CC_PSI_type2_d400T100K50$structure == "Group Individual1"),],
                         CC_PSI_type2_d100T200K50[which(CC_PSI_type2_d100T200K50$structure == "Group Individual1"),],
                         CC_PSI_type2_d100T400K50[which(CC_PSI_type2_d100T200K50$structure == "Group Individual1"),],
                         CC_PSI_type2_d100T200K10[which(CC_PSI_type2_d100T200K10$structure == "Group Individual1"),],
                         CC_PSI_type2_d100T200K100[which(CC_PSI_type2_d100T200K100$structure == "Group Individual1"),])
CC_PSI_G2_type2 <- rbind(CC_PSI_type2_d200T100K50[which(CC_PSI_type2_d200T100K50$structure == "Group Individual2"),],
                         CC_PSI_type2_d400T100K50[which(CC_PSI_type2_d400T100K50$structure == "Group Individual2"),],
                         CC_PSI_type2_d100T200K50[which(CC_PSI_type2_d100T200K50$structure == "Group Individual2"),],
                         CC_PSI_type2_d100T400K50[which(CC_PSI_type2_d100T200K50$structure == "Group Individual2"),],
                         CC_PSI_type2_d100T200K10[which(CC_PSI_type2_d100T200K10$structure == "Group Individual2"),],
                         CC_PSI_type2_d100T200K100[which(CC_PSI_type2_d100T200K100$structure == "Group Individual2"),])


method_list <- c("GRIDY","SCA-P","GICA","DSCA","DGICA","Unfitted_SCA-PF2","Unfitted_GICA")
CC_PSI_J_type2$method <- factor(CC_PSI_J_type2$method,levels=method_list,labels=method_list)
CC_PSI_G1_type2$method <- factor(CC_PSI_G1_type2$method,levels=method_list,labels=method_list)
CC_PSI_G2_type2$method <- factor(CC_PSI_G2_type2$method,levels=method_list,labels=method_list)

size_list <- c("d=200","d=400","T=200","T=400","2K=20","2K=200")
CC_PSI_J_type2$size<- factor(CC_PSI_J_type2$size,levels=size_list,labels=size_list)
CC_PSI_G1_type2$size <- factor(CC_PSI_G1_type2$size,levels=size_list,labels=size_list)
CC_PSI_G2_type2$size <- factor(CC_PSI_G2_type2$size,levels=size_list,labels=size_list)


pl_CC_PSI_J_type2 <- ggplot(CC_PSI_J_type2,aes(x=size,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="Size") +
  scale_y_continuous(name=TeX("$CC_{\\Psi}$ (Joint)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_PSI_G1_type2 <- ggplot(CC_PSI_G1_type2,aes(x=size,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="Size") +
  scale_y_continuous(name=TeX("$CC_{\\Psi}$ (Group 1)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
pl_CC_PSI_G2_type2 <- ggplot(CC_PSI_G2_type2,aes(x=size,y=value,color=method,fill=method)) +
  geom_boxplot(size=0.1,alpha=1)  +
  scale_fill_brewer(name="Method",palette="Set1") +
  scale_color_brewer(name="Method",palette="Set1") +
  scale_x_discrete(name="Size") +
  scale_y_continuous(name=TeX("$CC_{\\Psi}$ (Group 2)")) +
  theme(axis.text.y=element_text(size=7),legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))




# -----------------------------------------------------------------------------#
# Saving Figure 2
# -----------------------------------------------------------------------------#
mylegend_type2 <- g_legend(pl_R2_J_type2)

pdf(file = "./Figure2.pdf", width = 11.5, height = 10.5)
par(mar=c(0,0,0,0))
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
dev.off()


