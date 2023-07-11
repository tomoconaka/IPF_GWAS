
ukb <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/UKB_scratch/category_UKB_assoc.summary.tsv")
ukb <- ukb %>% filter(outcome %in% c("IPF", "RA_ILD", "SjS_ILD",
                                     "SSc_ILD", "ANCA_ILD")) %>% mutate(Beta = case_when(variable == "AGE" ~ 10*Beta,
                                                                                         TRUE ~ Beta),
                                                                        se = ifelse(variable == "AGE", 10*se, se))

ukb <- ukb %>% mutate(OR = exp(Beta),
                      LL = exp(Beta + qnorm(0.025)*se),
                      UL = exp(Beta + qnorm(0.975)*se),
                      p = pval) %>% filter(!(grepl("PC", variable))) %>% 
  select(variable, OR, LL, UL, p, Beta, se, Ncases, Ncontrols, outcome)

ukb <- ukb %>% mutate(cohort = "UKB")

ALL <- ukb %>% mutate(variable = case_when(variable %in% c("HighPRS") ~ "High IPF-PRS",
                                                          variable %in% c("telomere_gene", "tel_gene") ~ "Rare telomere variants carrier",
                                                          variable %in% c("ShortTS")~ "Short LTL",
                                                          variable %in% c("age.x", "AGE") ~ "10 year increase in age at recruitment",
                                                          variable %in% c("geneticSex.xMale", "SEX.y") ~ "Male",
                                                          variable %in% c("smokingever", "smokingEver") ~ "Ever smoker"))


ALL$OR[14] <- NA
ALL$LL[14] <- NA
ALL$UL[14] <- NA

for(i in unique(ALL$outcome)[-1]){
  ra <- ALL  %>% filter(outcome %in% c("IPF", i))
  cov <- unique(ra$variable)
  for(j in c(1:length(cov))){
    tmp <- ra %>% filter(variable == cov[j]) 
    m1 <- meta::metagen(Beta,
                        se,
                        data=tmp,
                        studlab=paste(cohort),
                        comb.fixed = TRUE,
                        comb.random = TRUE,
                        prediction = FALSE,
                        sm="OR")
    ALL$hetero_p[ALL$variable == cov[j] & ALL$outcome == i] <- m1$pval.Q
  }
}


library(forcats)
ALL %>% mutate(name = fct_relevel(variable, "Rare telomere variants carrier", "High IPF-PRS", "Short LTL",
                                  "10 year increase in age at recruitment", "Male", "Ever smoker")) %>%
  mutate(ILD = fct_relevel(outcome, "ANCA_ILD","SjS_ILD","SSc_ILD","RA_ILD", "IPF")) %>%
  ggplot(aes(x=ILD, y=OR, ymin=LL, ymax=UL, color=ILD)) +
  geom_pointrange(aes(col=ILD), lwd=0.8) + geom_hline(aes(fill=ILD), yintercept =1, linetype=2) +
  xlab("") + ylab("Odds Ratio (95% Confidence Interval)") +  geom_text(aes(label=round(OR,1), y=OR, col=outcome, hjust = 0.5, vjust = 1.5), size=5) +
  geom_errorbar(aes(ymin=LL, ymax=UL, col=ILD), width=0.1, cex=1) + #ylim(0.8,3.5)+
  facet_wrap(~name,  strip.position = 'left', nrow = 10) + theme_minimal() +scale_y_log10() +
  scale_color_manual(labels = c("ANCA_ILD","SjS_ILD","SSc_ILD","RA_ILD", "IPF"), values=c("#d73027","#fc8d59", "#fee090", "#91bfdb", "#4575b4")) + labs(color="Type of ILDs") +
  theme(plot.title=element_text(size=15),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"), 
        legend.text = element_text(size=15,face="bold"), 
        legend.title = element_text(size=15,face="bold"), 
        strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=15,face="bold"))+
  coord_flip() +guides(col = guide_legend(reverse = TRUE))

ggsave(paste0("../repo/IPF_GWAS/Figures/Fig4a.png"), width=25, height=30, units = "cm", dpi=300)

ALL <- ALL %>% select(variable, OR, LL, UL, p, hetero_p, Ncases, Ncontrols, outcome)
ALL %>% openxlsx::write.xlsx("../repo/IPF_GWAS/Tables/NonIPF_susceptibility.xlsx")



ukb <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/UKB_scratch/category_RAILD_survival.tsv")
ukb <- ukb %>% mutate(coef = case_when(variable == "AGE" ~ 10*coef,
                                       TRUE ~ coef),
                      se.coef. = ifelse(variable == "AGE", 10*se.coef., se.coef.))

ukb <- ukb %>% mutate(HR = exp(coef),
                      LL = exp(coef + qnorm(0.025)*se.coef.),
                      UL = exp(coef + qnorm(0.975)*se.coef.),
                      p = Pr...z..) %>% filter(!(grepl("PC", variable))) %>% 
  select(variable, HR, LL, UL, p, Ncases, Ntotal, hetero_p)

ALL <- ukb %>% mutate(variable = case_when(variable %in% c("HighPRS") ~ "High IPF-PRS",
                                           variable %in% c("telomere_gene", "tel_gene") ~ "Rare telomere variants carrier",
                                           variable %in% c("ShortTS")~ "Short LTL",
                                           variable %in% c("age.x", "AGE") ~ "10 year increase in age at recruitment",
                                           variable %in% c("geneticSex.xMale", "SEX.x") ~ "Male",
                                           variable %in% c("smokingever", "smokingEver") ~ "Ever smoker"))

library(forcats)
ALL$outcome <- "RA_ILD"
ALL %>% mutate(name = fct_relevel(variable, "Rare telomere variants carrier", "High IPF-PRS", "Short LTL",
                                  "10 year increase in age at recruitment", "Male", "Ever smoker")) %>%
  ggplot(aes(x=outcome, y=HR, ymin=LL, ymax=UL, color="#91bfdb")) +
  geom_pointrange(aes(col="#91bfdb"), lwd=0.8) + geom_hline(aes(fill="#91bfdb"), yintercept =1, linetype=2) +
  xlab("") + ylab("Hazard Ratio (95% Confidence Interval)") +  geom_text(aes(label=round(HR,1), y=HR, col="#91bfdb", hjust = 0.5, vjust = 1.5), size=5) +
  geom_errorbar(aes(ymin=LL, ymax=UL, col="#91bfdb"), width=0.1, cex=1) + #ylim(0.8,3.5)+
  facet_wrap(~name,  strip.position = 'left', nrow = 10) + theme_minimal() +scale_y_log10() +
  scale_color_identity() + 
  theme(plot.title=element_text(size=15),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"), 
        legend.text = element_text(size=15,face="bold"), 
        legend.title = element_text(size=15,face="bold"), 
        strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=15,face="bold"))+
  coord_flip() +guides(col = guide_legend(reverse = TRUE))

ggsave(paste0("../repo/IPF_GWAS/Figures/Fig4b.png"), width=25, height=10, units = "cm", dpi=300)

ALL <- ALL %>% select(variable, HR, LL, UL, p, hetero_p, Ncases, Ntotal)
ALL %>% openxlsx::write.xlsx("../repo/IPF_GWAS/Tables/RA_survival.xlsx")




