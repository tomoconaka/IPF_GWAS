library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggtext)
d <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/ALL_meta_IPF.b38_filtered.txt.gz")
d$Pvalue <- as.numeric(d$Pvalue)
d <- d %>% filter(Pvalue < 0.05)
d <- d %>% mutate(Pvalue = ifelse(Pvalue < 1e-300, 1e-300, Pvalue))
data_cum <- d %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(POS)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(CHR, bp_add)

gwas_data <- d %>% 
  inner_join(data_cum, by = "CHR") %>% 
  mutate(bp_cum = POS + bp_add)

axis_set <- gwas_data %>% 
  group_by(CHR) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  filter(Pvalue == min(Pvalue)) %>% 
  mutate(ylim = abs(floor(log10(Pvalue))) + 2) %>% 
  pull(ylim)

sig <- 5e-8

leadvariant <- c("1:9109660", "1:150579566", "1:155199564", "1:163372193","1:214484297", 
                 "2:64652534",
                 "3:44804157", "3:169774313", 
                 "4:88963935",
                 "5:1282299", "5:169588475", 
                 "6:7562999", "6:27698141","6:32436600","6:43385242","6:122431405", 
                 "7:1936821", "7:100032719",
                 "8:119931204", 
                 "9:106717987", 
                 "10:941334","10:103920828", 
                 "11:1219991",
                 "12:101892202",
                 "13:82382479", "13:112880670",
                 "15:40426335", "15:85742593", 
                 "16:112241", "16:67960834", 
                 "17:43370681", "17:46253848", "17:76029575", 
                 "18:666625", 
                 "19:4717660", "19:5847989",
                 "20:63585923"
)


gwas_data <- gwas_data %>% mutate(shape = ifelse(MarkerName %in% leadvariant, 18, 16),
                                  color = case_when(CHR == 1	& POS >= 150579566 - 500000 & POS <= 150579566 + 500000 ~ "#d7191c",#EUR specific
                                                    CHR == 1	& POS >= 155199564 - 500000 & POS <= 155199564 + 500000 ~ "#d7191c",#EUR specific
                                                    CHR == 1	& POS >= 163372193 - 500000 & POS <= 163372193 + 500000 ~ "#d7191c",#EUR specific
                                                    CHR == 2	& POS >= 64652534 - 500000 & POS <= 64652534 + 500000 ~ "#d7191c",#EUR and EAS
                                                    #CHR == 6	& POS >= 27698141 - 500000 & POS <= 27698141 + 500000 ~ "#ca0020",#EUR specific
                                                    CHR == 6	& POS >= 32436600 - 500000 & POS <= 32436600 + 500000 ~ "#d7191c",#EUR and EAS
                                                    #CHR == 6	& POS >= 122431405 - 500000 & POS <= 122431405 + 500000 ~ "#ca0020",#EUR specific
                                                    #CHR == 10	& POS >= 941334 - 500000 & POS <= 941334 + 500000 ~ "#ca0020",#EUR specific
                                                    CHR == 12	& POS >= 101892202 - 500000 & POS <= 101892202 + 500000 ~ "#d7191c",#EUR specific
                                                    CHR == 13	& POS >= 82382479 - 500000 & POS <= 82382479 + 500000 ~ "#1b7837",#EAS specific
                                                    CHR == 17	& POS >= 43370681 - 500000 & POS <= 43370681 + 500000 ~ "#d7191c",#EUR and EAS
                                                    CHR == 17	& POS >= 76029575 - 500000 & POS <= 76029575 + 500000 ~ "#d7191c",#EUR and EAS
                                                    CHR == 18	& POS >= 666625 - 500000 & POS <= 666625 + 500000 ~ "#d7191c",#EUR specific
                                                    CHR %in% c(1,3,5,7,9,11,13,15,17,19,21,23) ~ "#276FBF",
                                                    TRUE ~ "#183059"))

gwas_data <- gwas_data %>% mutate(size = case_when(MarkerName %in% leadvariant ~ 3,
                                                   color %in% c("#1b7837", "#d7191c") ~ 2, 
                                                   TRUE ~ 1))


##e4007f
manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(Pvalue), 
                                  color = color, size = size, shape = shape)) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
  scale_y_log10() + scale_color_identity() + scale_size_identity() +
  #scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$CHR)))) +
  #scale_size_continuous(range = c(0.5,3)) +
  scale_shape_identity() +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

print(manhplot)
ggsave(paste0("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/Fig1a_Manhattan.png"), 
       width=40, height=15, units = "cm", dpi=300)



d <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/EUR_meta_IPF.b38_filtered.txt.gz")
d$Pvalue <- as.numeric(d$Pvalue)
d <- d %>% filter(Pvalue < 0.05)
d <- d %>% mutate(Pvalue = ifelse(Pvalue < 1e-300, 1e-300, Pvalue))
data_cum <- d %>% 
  group_by(`#CHR`) %>% 
  summarise(max_bp = max(POS)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(`#CHR`, bp_add)

gwas_data <- d %>% 
  inner_join(data_cum, by = "#CHR") %>% 
  mutate(bp_cum = POS + bp_add)

axis_set <- gwas_data %>% 
  group_by(`#CHR`) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  filter(Pvalue == min(Pvalue)) %>% 
  mutate(ylim = abs(floor(log10(Pvalue))) + 2) %>% 
  pull(ylim)

sig <- 5e-8

leadvariant <- c("1:9109660", "1:150579566", "1:155199564", "1:163372193","1:214484297", 
                 "2:64652534",
                 "3:44804157", "3:169774313", 
                 "4:88963935",
                 "5:1282299", "5:169588475", 
                 "6:7562999", "6:27698141","6:32436600","6:43385242","6:122431405", 
                 "7:1936821", "7:100032719",
                 "8:119931204", 
                 "9:106717987", 
                 "10:941334","10:103920828", 
                 "11:1219991",
                 "12:101892202",
                 "13:82382479", "13:112880670",
                 "15:40426335", "15:85742593", 
                 "16:112241", "16:67960834", 
                 "17:43370681", "17:46253848", "17:76029575", 
                 "18:666625", 
                 "19:4717660", "19:5847989",
                 "20:63585923"
)


gwas_data <- gwas_data %>% mutate(shape = ifelse(MarkerName %in% leadvariant, 18, 16),
                                  size = ifelse(MarkerName %in% leadvariant, 3, 0.5),
                                  color = case_when(`#CHR` == 6	& POS >= 27698141 - 500000 & POS <= 27698141 + 500000 ~ "#d7191c",#EUR specific
                                                    `#CHR` == 6	& POS >= 122431405 - 500000 & POS <= 122431405 + 500000 ~ "#d7191c",#EUR specific
                                                    `#CHR` == 10	& POS >= 941334 - 500000 & POS <= 941334 + 500000 ~ "#d7191c",#EUR specific
                                                    `#CHR` %in% c(1,3,5,7,9,11,13,15,17,19,21,23) ~ "#276FBF",
                                                    TRUE ~ "#183059"))

gwas_data <- gwas_data %>% mutate(size = case_when(MarkerName %in% leadvariant ~ 3,
                                                   color %in% c("#1b7837", "#d7191c") ~ 2, 
                                                   TRUE ~ 1))

##e4007f
manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(Pvalue), 
                                  color = color, size = size, shape = shape)) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$`#CHR`, breaks = axis_set$center) +
  scale_y_log10() + scale_color_identity() + scale_size_identity() +
  #scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$`#CHR`)))) +
  #scale_size_continuous(range = c(0.5,3)) +
  scale_shape_identity() + #geom_text_repel() + 
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

print(manhplot)
ggsave(paste0("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/Fig1b_Manhattan.png"), 
       width=40, height=15, units = "cm", dpi=300)







