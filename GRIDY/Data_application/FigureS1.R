#-----------------------------------------------------------------------------#
#   File name : FigureS1.R    												  
#
#   Project : "Group integrative dynamic factor models 
#             with application to multiple subject brain connectivity"
#
#   Maintainer : Younghoon Kim                     
#
#   Date : Sep. 1st, 2024
#
#   Purpose : code for producing results for plotting Figure S1 in Supplemental material
#
#   R version 4.0.5 (2021-03-31)                                      
#
#   Input data file : /Data_application/result/result_intermediate.RData
# 
#   Output data file : /Data_application/figures/FigureS1.pdf
#
#   Required R packages : ggplot2_3.4.0, latex2exp_0.9.6, ggpubr_0.5.0, RColorBrewer_1.1-3,
#                         reshape2_1.4.4, and gridExtra_2.3
#-----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
# Preparation
# -----------------------------------------------------------------------------#
# Load source code from the parent directory
source(paste0(dirname(getwd()),"/","library_application.R"))

# Load data from the previous step
load("./result/result_intermediate.RData")

# -----------------------------------------------------------------------------#
# Creating Figure S1
# -----------------------------------------------------------------------------#
# rank table:
rank_table <- data.frame(rank = unique(df_index$rank)[order(unique(df_index$rank))],
                         frequency = tabulate(match(r_hat,unique(df_index$rank)[order(unique(df_index$rank))])))

rank_table_G1 <- data.frame(rank = unique(df_index$rank[which(df_index$group == 1)])[order(unique(df_index$rank[which(df_index$group == 1)]))],
                            frequency = tabulate(match(r_hat[which(df_index$group == 1)],
                                                       unique(df_index$rank[which(df_index$group == 1)])[order(unique(df_index$rank[which(df_index$group == 1)]))])))

rank_table_G2 <- data.frame(rank = unique(df_index$rank[which(df_index$group == 2)])[order(unique(df_index$rank[which(df_index$group == 2)]))],
                            frequency = tabulate(match(r_hat[which(df_index$group == 2)],
                                                       unique(df_index$rank[which(df_index$group == 2)])[order(unique(df_index$rank[which(df_index$group == 2)]))])))


# Plotting: This produces Figure S1
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
  geom_step(data=rank_table,aes(x=rank, y=cumsum(frequency/dim(df_hist)[1]*250)),alpha=0.5,linetype="dashed")

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


# -----------------------------------------------------------------------------#
# Saving Figure S1
# -----------------------------------------------------------------------------#
plot_tmp <- ggarrange(G_hist,G1_hist,G2_hist,nrow=3,common.legend=TRUE,legend="bottom")

pdf(file = "./figures/FigureS1.pdf", width = 11.5, height = 10.5)
par(mar=c(0,0,0,0))
plot_tmp
dev.off()

