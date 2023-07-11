setwd("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/")

library(data.table)
library(tidyverse)
all <- fread("ALL_sample.tsv")
pc_unrelated <- fread("grm_unrelated.pc")
pc_related <- fread("grm_related.projections")

pcs <- bind_rows(pc_unrelated, pc_related) %>% select(-FID)

data <- all %>% inner_join(pcs, by="IID")

data <- data %>% mutate(PCAoutlier = case_when(PC1 > 0.04 & PC2 > 0.025 ~ TRUE,
                                               PC3 < -0.08 ~ TRUE,
                                               PC5 < -0.08 ~ TRUE,
                                               TRUE ~ FALSE))


data <- data %>% mutate(Status = case_when(group == "Control" ~ "Control",
                                           group == "sIPF" ~ "sporadic IPF",
                                           group == "fIPF" ~ "FPF (IPF)",
                                           TRUE ~ "FPF (non-IPF)"))

p1 <- ggplot(data, aes(x=PC1, y=PC2, col=Status)) + geom_point() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + theme(axis.text.x = element_text(size=20,face="bold"),
                                                              axis.text.y=element_text(size=20,face="bold"),
                                                              axis.title=element_text(size=20,face="bold"), 
                                                              strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=20,face="bold"),
                                                              legend.title = element_text(size=20,face="bold"), 
                                                              legend.text = element_text(size=20,face="bold"), 
                                                              plot.title = element_text(size = 20),
                                                              strip.text = element_text(size=25)) 
p2 <- ggplot(data, aes(x=PC1, y=PC3, col=Status)) + geom_point() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + theme(axis.text.x = element_text(size=20,face="bold"),
                                                           axis.text.y=element_text(size=20,face="bold"),
                                                           axis.title=element_text(size=20,face="bold"), 
                                                           strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=20,face="bold"),
                                                           legend.title = element_text(size=20,face="bold"), 
                                                           legend.text = element_text(size=20,face="bold"), 
                                                           plot.title = element_text(size = 20),
                                                           strip.text = element_text(size=25)) 
p3 <- ggplot(data, aes(x=PC4, y=PC5, col=Status)) + geom_point() + 
  theme_bw()+ scale_colour_brewer(palette = "Dark2") + theme(axis.text.x = element_text(size=20,face="bold"),
                                                           axis.text.y=element_text(size=20,face="bold"),
                                                           axis.title=element_text(size=20,face="bold"), 
                                                           strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=20,face="bold"),
                                                           legend.title = element_text(size=20,face="bold"), 
                                                           legend.text = element_text(size=20,face="bold"), 
                                                           plot.title = element_text(size = 20),
                                                           strip.text = element_text(size=25)) 


p4 <- ggplot(data, aes(x=PC1, y=PC2, col=PCAoutlier)) + geom_point() + 
  theme_bw() + scale_colour_brewer(palette = "Set1") + theme(axis.text.x = element_text(size=20,face="bold"),
                                                              axis.text.y=element_text(size=20,face="bold"),
                                                              axis.title=element_text(size=20,face="bold"), 
                                                              strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=20,face="bold"),
                                                              legend.title = element_text(size=20,face="bold"), 
                                                              legend.text = element_text(size=20,face="bold"), 
                                                              plot.title = element_text(size = 20),
                                                              strip.text = element_text(size=25)) 
p5 <- ggplot(data, aes(x=PC1, y=PC3, col=PCAoutlier)) + geom_point() + 
  theme_bw() + scale_colour_brewer(palette = "Set1") + theme(axis.text.x = element_text(size=20,face="bold"),
                                                              axis.text.y=element_text(size=20,face="bold"),
                                                              axis.title=element_text(size=20,face="bold"), 
                                                              strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=20,face="bold"),
                                                              legend.title = element_text(size=20,face="bold"), 
                                                              legend.text = element_text(size=20,face="bold"), 
                                                              plot.title = element_text(size = 20),
                                                              strip.text = element_text(size=25)) 
p6 <- ggplot(data, aes(x=PC4, y=PC5, col=PCAoutlier)) + geom_point() + 
  theme_bw()+ scale_colour_brewer(palette = "Set1") + theme(axis.text.x = element_text(size=20,face="bold"),
                                                             axis.text.y=element_text(size=20,face="bold"),
                                                             axis.title=element_text(size=20,face="bold"), 
                                                             strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=20,face="bold"),
                                                             legend.title = element_text(size=20,face="bold"), 
                                                             legend.text = element_text(size=20,face="bold"), 
                                                             plot.title = element_text(size = 20),
                                                             strip.text = element_text(size=25)) 

library(gridExtra)
png(file="../repo/ILD_GWAS/05.PCA/PCAplot.png", width = 1000, height = 1000)
grid.arrange(arrangeGrob(p1, p2, p3, ncol=1), arrangeGrob(p4, p5, p6, ncol=1),
             widths = c(10, 9))
#grid.arrange(arrangeGrob(p1, p2, p3, ncol=3))
dev.off()

write.table(data, "ALL_sample_PCs.tsv", quote=F, col.names = T, row.names = F, sep="\t")
