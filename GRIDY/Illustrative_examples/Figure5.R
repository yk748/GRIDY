#-----------------------------------------------------------------------------#
#   File name : Figure5.R    												  
#
#   Project : "Group integrative dynamic factor models 
#             with application to multiple subject brain connectivity"
#
#   Maintainer : Younghoon Kim                     
#
#   Date : Sep. 1st, 2024
#
#   Purpose : code for producing results for plotting Figure 5 in Illustrative example 
#             section
#
#   R version 4.0.5 (2021-03-31)                                    
#
#   Input data file : /Illustrative_examples/result/sim3_set-_d-T-K-.RData
# 
#   Output data file : /Illustrative_examples/figures/Figure5.pdf
#
#   Required R packages : ggplot2_3.4.0, latex2exp_0.9.6, ggpubr_0.5.0, RColorBrewer_1.1-3,
#                         reshape2_1.4.4, and gridExtra_2.3
#-----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#
# Loading simulated data
# -----------------------------------------------------------------------------#
source(paste0(dirname(getwd()),"/","library_simulation.R"))

initial_rank <- c("(1,1)","(1,2)","(1,3)","(2,2)","(1,4)",
                  "(2,3)","(2,4)","(3,3)","(3,4)","(4,4)")

load("./result/sim3_set1_d100T100K50.RData"); r_hat_c1_d100T100K50 <- na.omit(r_hat)
load("./result/sim3_set2_d100T100K50.RData"); r_hat_c2_d100T100K50 <- na.omit(r_hat)
load("./result/sim3_set3_d100T100K50.RData"); r_hat_c3_d100T100K50 <- na.omit(r_hat)
load("./result/sim3_set4_d100T100K50.RData"); r_hat_c4_d100T100K50 <- na.omit(r_hat)
load("./result/sim3_set5_d100T100K50.RData"); r_hat_c5_d100T100K50 <- na.omit(r_hat)
load("./result/sim3_set6_d100T100K50.RData"); r_hat_c6_d100T100K50 <- na.omit(r_hat)
colnames(r_hat_c1_d100T100K50) <- colnames(r_hat_c2_d100T100K50) <- colnames(r_hat_c3_d100T100K50) <- colnames(r_hat_c4_d100T100K50) <- colnames(r_hat_c5_d100T100K50) <- colnames(r_hat_c6_d100T100K50) <- initial_rank

load("./result/sim3_set1_d200T100K50.RData"); r_hat_c1_d200T100K50 <- na.omit(r_hat)
load("./result/sim3_set2_d200T100K50.RData"); r_hat_c2_d200T100K50 <- na.omit(r_hat)
load("./result/sim3_set3_d200T100K50.RData"); r_hat_c3_d200T100K50 <- na.omit(r_hat)
load("./result/sim3_set4_d200T100K50.RData"); r_hat_c4_d200T100K50 <- na.omit(r_hat)
load("./result/sim3_set5_d200T100K50.RData"); r_hat_c5_d200T100K50 <- na.omit(r_hat)
load("./result/sim3_set6_d200T100K50.RData"); r_hat_c6_d200T100K50 <- na.omit(r_hat)
colnames(r_hat_c1_d200T100K50) <- colnames(r_hat_c2_d200T100K50) <- colnames(r_hat_c3_d200T100K50) <- colnames(r_hat_c4_d200T100K50) <- colnames(r_hat_c5_d200T100K50) <- colnames(r_hat_c6_d200T100K50) <- initial_rank

load("./result/sim3_set1_d300T100K50.RData"); r_hat_c1_d300T100K50 <- na.omit(r_hat)
load("./result/sim3_set2_d300T100K50.RData"); r_hat_c2_d300T100K50 <- na.omit(r_hat)
load("./result/sim3_set3_d300T100K50.RData"); r_hat_c3_d300T100K50 <- na.omit(r_hat)
load("./result/sim3_set4_d300T100K50.RData"); r_hat_c4_d300T100K50 <- na.omit(r_hat)
load("./result/sim3_set5_d300T100K50.RData"); r_hat_c5_d300T100K50 <- na.omit(r_hat)
load("./result/sim3_set6_d300T100K50.RData"); r_hat_c6_d300T100K50 <- na.omit(r_hat)
colnames(r_hat_c1_d300T100K50) <- colnames(r_hat_c2_d300T100K50) <- colnames(r_hat_c3_d300T100K50) <- colnames(r_hat_c4_d300T100K50) <- colnames(r_hat_c5_d300T100K50) <- colnames(r_hat_c6_d300T100K50) <- initial_rank

load("./result/sim3_set1_d100T200K50.RData"); r_hat_c1_d100T200K50 <- na.omit(r_hat)
load("./result/sim3_set2_d100T200K50.RData"); r_hat_c2_d100T200K50 <- na.omit(r_hat)
load("./result/sim3_set3_d100T200K50.RData"); r_hat_c3_d100T200K50 <- na.omit(r_hat)
load("./result/sim3_set4_d100T200K50.RData"); r_hat_c4_d100T200K50 <- na.omit(r_hat)
load("./result/sim3_set5_d100T200K50.RData"); r_hat_c5_d100T200K50 <- na.omit(r_hat)
load("./result/sim3_set6_d100T200K50.RData"); r_hat_c6_d100T200K50 <- na.omit(r_hat)
colnames(r_hat_c1_d100T200K50) <- colnames(r_hat_c2_d100T200K50) <- colnames(r_hat_c3_d100T200K50) <- colnames(r_hat_c4_d100T200K50) <- colnames(r_hat_c5_d100T200K50) <- colnames(r_hat_c6_d100T200K50) <- initial_rank

load("./result/sim3_set1_d100T300K50.RData"); r_hat_c1_d100T300K50 <- na.omit(r_hat)
load("./result/sim3_set2_d100T300K50.RData"); r_hat_c2_d100T300K50 <- na.omit(r_hat)
load("./result/sim3_set3_d100T300K50.RData"); r_hat_c3_d100T300K50 <- na.omit(r_hat)
load("./result/sim3_set4_d100T300K50.RData"); r_hat_c4_d100T300K50 <- na.omit(r_hat)
load("./result/sim3_set5_d100T300K50.RData"); r_hat_c5_d100T300K50 <- na.omit(r_hat)
load("./result/sim3_set6_d100T300K50.RData"); r_hat_c6_d100T300K50 <- na.omit(r_hat)
colnames(r_hat_c1_d100T300K50) <- colnames(r_hat_c2_d100T300K50) <- colnames(r_hat_c3_d100T300K50) <- colnames(r_hat_c4_d100T300K50) <- colnames(r_hat_c5_d100T300K50) <- colnames(r_hat_c6_d100T300K50) <- initial_rank


df_c1_d100T100K50 <- melt(r_hat_c1_d100T100K50)[,-1]; df_c1_d100T100K50$snr <- 0.25
df_c2_d100T100K50 <- melt(r_hat_c2_d100T100K50)[,-1]; df_c2_d100T100K50$snr <- 0.5
df_c3_d100T100K50 <- melt(r_hat_c3_d100T100K50)[,-1]; df_c3_d100T100K50$snr <- 0.75
df_c4_d100T100K50 <- melt(r_hat_c4_d100T100K50)[,-1]; df_c4_d100T100K50$snr <- 1.25
df_c5_d100T100K50 <- melt(r_hat_c5_d100T100K50)[,-1]; df_c5_d100T100K50$snr <- 2
df_c6_d100T100K50 <- melt(r_hat_c6_d100T100K50)[,-1]; df_c6_d100T100K50$snr <- 4
df_d100T100K50 <- rbind(df_c1_d100T100K50,df_c2_d100T100K50,df_c3_d100T100K50,
                        df_c4_d100T100K50,df_c5_d100T100K50,df_c6_d100T100K50)
df_d100T100K50$size <- 1

df_c1_d200T100K50 <- melt(r_hat_c1_d200T100K50)[,-1]; df_c1_d200T100K50$snr <- 0.25
df_c2_d200T100K50 <- melt(r_hat_c2_d200T100K50)[,-1]; df_c2_d200T100K50$snr <- 0.5
df_c3_d200T100K50 <- melt(r_hat_c3_d200T100K50)[,-1]; df_c3_d200T100K50$snr <- 0.75
df_c4_d200T100K50 <- melt(r_hat_c4_d200T100K50)[,-1]; df_c4_d200T100K50$snr <- 1.25
df_c5_d200T100K50 <- melt(r_hat_c5_d200T100K50)[,-1]; df_c5_d200T100K50$snr <- 2
df_c6_d200T100K50 <- melt(r_hat_c6_d200T100K50)[,-1]; df_c6_d200T100K50$snr <- 4
df_d200T100K50 <- rbind(df_c1_d200T100K50,df_c2_d200T100K50,df_c3_d200T100K50,
                        df_c4_d200T100K50,df_c5_d200T100K50,df_c6_d200T100K50)
df_d200T100K50$size <- 2

df_c1_d300T100K50 <- melt(r_hat_c1_d300T100K50)[,-1]; df_c1_d300T100K50$snr <- 0.25
df_c2_d300T100K50 <- melt(r_hat_c2_d300T100K50)[,-1]; df_c2_d300T100K50$snr <- 0.5
df_c3_d300T100K50 <- melt(r_hat_c3_d300T100K50)[,-1]; df_c3_d300T100K50$snr <- 0.75
df_c4_d300T100K50 <- melt(r_hat_c4_d300T100K50)[,-1]; df_c4_d300T100K50$snr <- 1.25
df_c5_d300T100K50 <- melt(r_hat_c5_d300T100K50)[,-1]; df_c5_d300T100K50$snr <- 2
df_c6_d300T100K50 <- melt(r_hat_c6_d300T100K50)[,-1]; df_c6_d300T100K50$snr <- 4
df_d300T100K50 <- rbind(df_c1_d300T100K50,df_c2_d300T100K50,df_c3_d300T100K50,
                        df_c4_d300T100K50,df_c5_d300T100K50,df_c6_d300T100K50)
df_d300T100K50$size <- 3

df_c1_d100T200K50 <- melt(r_hat_c1_d100T200K50)[,-1]; df_c1_d100T200K50$snr <- 0.25
df_c2_d100T200K50 <- melt(r_hat_c2_d100T200K50)[,-1]; df_c2_d100T200K50$snr <- 0.5
df_c3_d100T200K50 <- melt(r_hat_c3_d100T200K50)[,-1]; df_c3_d100T200K50$snr <- 0.75
df_c4_d100T200K50 <- melt(r_hat_c4_d100T200K50)[,-1]; df_c4_d100T200K50$snr <- 1.25
df_c5_d100T200K50 <- melt(r_hat_c5_d100T200K50)[,-1]; df_c5_d100T200K50$snr <- 2
df_c6_d100T200K50 <- melt(r_hat_c6_d100T200K50)[,-1]; df_c6_d100T200K50$snr <- 4
df_d100T200K50 <- rbind(df_c1_d100T200K50,df_c2_d100T200K50,df_c3_d100T200K50,
                        df_c4_d100T200K50,df_c5_d100T200K50,df_c6_d100T200K50)
df_d100T200K50$size <- 4

df_c1_d100T300K50 <- melt(r_hat_c1_d100T300K50)[,-1]; df_c1_d100T300K50$snr <- 0.25
df_c2_d100T300K50 <- melt(r_hat_c2_d100T300K50)[,-1]; df_c2_d100T300K50$snr <- 0.5
df_c3_d100T300K50 <- melt(r_hat_c3_d100T300K50)[,-1]; df_c3_d100T300K50$snr <- 0.75
df_c4_d100T300K50 <- melt(r_hat_c4_d100T300K50)[,-1]; df_c4_d100T300K50$snr <- 1.25
df_c5_d100T300K50 <- melt(r_hat_c5_d100T300K50)[,-1]; df_c5_d100T300K50$snr <- 2
df_c6_d100T300K50 <- melt(r_hat_c6_d100T300K50)[,-1]; df_c6_d100T300K50$snr <- 4
df_d100T300K50 <- rbind(df_c1_d100T300K50,df_c2_d100T300K50,df_c3_d100T300K50,
                        df_c4_d100T300K50,df_c5_d100T300K50,df_c6_d100T300K50)
df_d100T300K50$size <- 5

df_rank <- rbind(df_d100T100K50,df_d200T100K50,df_d300T100K50,df_d100T200K50,df_d100T300K50)
df_rank$size <- factor(df_rank$size,
                       levels=c(1,2,4,3,5),
                       labels=c("(100,100)","(200,100)","(100,200)","(300,100)","(100,300)"))
df_rank$snr <- factor(df_rank$snr,
                      levels=c(0.25,0.50,0.75,1.25,2.00,4.00),
                      labels=c("SNR=0.25","SNR=0.5","SNR=0.75","SNR=1.25","SNR=2","SNR=4"))

mycolors <- colorRampPalette(brewer.pal(9,"Set1"))(10)
size.labs <- c("(100,100)","(200,100)","(100,200)","(300,100)","(100,300)")
names(size.labs) <- c("(100,100)","(200,100)","(100,200)","(300,100)","(100,300)")
snr.labs <- c("SNR=0.25","SNR=0.5","SNR=0.75","SNR=1.25","SNR=2","SNR=4")
names(snr.labs) <- c("SNR=0.25","SNR=0.5","SNR=0.75","SNR=1.25","SNR=2","SNR=4")

# -----------------------------------------------------------------------------#
# Saving Figure 5
# -----------------------------------------------------------------------------#
pdf(file = "./figures/Figure5.pdf", width = 11.5, height = 10.5)
par(mar=c(0,0,0,0))
ggplot(df_rank,
       aes(x=as.factor(value),fill=Var2)) +
  geom_bar(alpha=1,color="white",aes(y=after_stat(count)),
           position='stack',width = 1,color = "NA") +
  scale_fill_manual(name=TeX("$(r_{J},r_{G})$"),
                    values=mycolors,labels=initial_rank) +
  scale_x_discrete(name=TeX("Estimated initial rank $\\hat{r}$")) +
  scale_y_continuous(name="Frequency") + 
  facet_grid(snr~size, 
             labeller = labeller(size = size.labs, snr = snr.labs)) + 
  theme(axis.text.y=element_text(size=7),legend.position = "bottom") + 
  guides(fill = guide_legend(nrow = 1)) 
dev.off()


