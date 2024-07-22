#-----------------------------------------------------------------------------#
# Code for producing Figure S2 in Supplementary material
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
# Creating Figure S2
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


# Plotting: This produces Figure S2
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


# -----------------------------------------------------------------------------#
# Saving Figure S2
# -----------------------------------------------------------------------------#
pdf(file = "./FigureS2.pdf", width = 11.5, height = 10.5)
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
dev.off()

