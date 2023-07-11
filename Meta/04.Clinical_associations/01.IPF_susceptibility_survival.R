

kyoto <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/category_IPF.susceptibility.assoc.tsv")
ukb <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/UKB_scratch/category_UKB_assoc.summary.tsv")

kyoto <- kyoto %>% mutate(Estimate = case_when(variable == "age.x" ~ 10*Estimate,
                                               TRUE ~ Estimate),
                          Std..Error = ifelse(variable == "age.x", 10*Std..Error, Std..Error))
kyoto <- kyoto %>% mutate(OR = exp(Estimate),
                          LL = exp(Estimate + qnorm(0.025)*Std..Error),
                          UL = exp(Estimate + qnorm(0.975)*Std..Error),
                          p = Pr...z..) %>% filter(!(grepl("PC", variable))) %>% 
  mutate(Beta = Estimate, se = Std..Error) %>% 
  select(variable, OR, LL, UL, p, Beta, se)

kyoto <- kyoto %>% mutate(cohort = "Kyoto-ILD")
kyoto <- kyoto %>% mutate(Ncases = 106, Ncontrols = 2664)
ukb <- ukb %>% filter(outcome == "IPF") %>% mutate(Beta = case_when(variable == "AGE" ~ 10*Beta,
                                                                    TRUE ~ Beta),
                                                   se = ifelse(variable == "AGE", 10*se, se))

ukb <- ukb %>% mutate(OR = exp(Beta),
                      LL = exp(Beta + qnorm(0.025)*se),
                      UL = exp(Beta + qnorm(0.975)*se),
                      p = pval) %>% filter(!(grepl("PC", variable))) %>% 
  select(variable, OR, LL, UL, p, Beta, se, Ncases, Ncontrols)

ukb <- ukb %>% mutate(cohort = "UKB")

ALL <- bind_rows(kyoto, ukb)
ALL <- ALL %>% mutate(variable = case_when(variable %in% c("HighPRS") ~ "High IPF-PRS",
                                           variable %in% c("telomere_gene", "tel_gene") ~ "Rare telomere variants carrier",
                                           variable %in% c("ShortTS")~ "Short LTL",
                                           variable %in% c("age.x", "AGE") ~ "10 year increase in age at recruitment",
                                           variable %in% c("geneticSex.xMale", "SEX.y") ~ "Male",
                                           variable %in% c("smokingever", "smokingEver") ~ "Ever smoker"))

cov <- unique(ALL$variable)
for(j in c(1:length(cov))){
  tmp <- ALL %>% filter(variable == cov[j]) 
  m1 <- meta::metagen(Beta,
                      se,
                      data=tmp,
                      studlab=paste(cohort),
                      comb.fixed = TRUE,
                      comb.random = TRUE,
                      prediction = FALSE,
                      sm="OR")
  ALL$hetero_p[ALL$variable == cov[j]] <- m1$pval.Q
}


library(forcats)
ALL %>% mutate(name = fct_relevel(variable, "Rare telomere variants carrier", "High IPF-PRS", "Short LTL",
                                  "10 year increase in age at recruitment", "Male", "Ever smoker")) %>%
  ggplot(aes(x=cohort, y=OR, ymin=LL, ymax=UL, color=cohort)) +
  geom_pointrange(aes(col=cohort), lwd=0.8) + geom_hline(aes(fill=cohort), yintercept =1, linetype=2) +
  xlab("") + ylab("Odds Ratio (95% Confidence Interval)") + geom_text(aes(label=round(OR,1), y=OR, col=cohort, hjust = 0.5, vjust = 1.5), size=4) +
  geom_errorbar(aes(ymin=LL, ymax=UL, col=cohort), width=0.1, cex=1) + #ylim(0.8,3.5)+
  facet_wrap(~name,  strip.position = 'left', nrow = 10) + theme_minimal() +scale_y_log10() +
  scale_color_manual(labels = c("Kyoto-ILD", "UKB"), values=c("#e4007f", "#0068b7")) + labs(color="Cohort") +
  theme(plot.title=element_text(size=15),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"), 
        legend.text = element_text(size=15,face="bold"), 
        legend.title = element_text(size=15,face="bold"), 
        strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=15,face="bold"))+
  coord_flip() +guides(col = guide_legend(reverse = TRUE))

ggsave(paste0("../repo/IPF_GWAS/Figures/Fig3a.png"), width=25, height=15, units = "cm", dpi=300)

ALL <- ALL %>% select(variable, OR, LL, UL, p, cohort, hetero_p, Ncases, Ncontrols)
ALL %>% openxlsx::write.xlsx("../repo/IPF_GWAS/Tables/IPF_susceptibility.xlsx")

####################SURVIVAL
kyoto <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/category_IPF.survival.assoc.tsv")
ukb <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/UKB_scratch/category_IPF_survival.tsv")

kyoto <- kyoto %>% mutate(coef = case_when(variable == "age.x" ~ 10*coef,
                                           TRUE ~ coef),
                          se.coef. = ifelse(variable == "age.x", 10*se.coef., se.coef.))

kyoto <- kyoto %>% mutate(HR = exp(coef),
                          LL = exp(coef + qnorm(0.025)*se.coef.),
                          UL = exp(coef + qnorm(0.975)*se.coef.),
                          p = Pr...z..) %>% filter(!(grepl("PC", variable))) %>% 
  mutate(Beta = coef, se = se.coef.) %>% 
  select(variable, HR, LL, UL, p, Beta, se, Ncases, Ncontrols) %>% rename(Ntotal = Ncontrols)

kyoto <- kyoto %>% mutate(cohort = "Kyoto-ILD")

ukb <- ukb %>% mutate(Beta = case_when(variable == "AGE" ~ 10*coef,
                                       TRUE ~ coef),
                      se = ifelse(variable == "AGE", 10*se.coef., se.coef.))

ukb <- ukb %>% mutate(HR = exp(Beta),
                      LL = exp(Beta + qnorm(0.025)*se),
                      UL = exp(Beta + qnorm(0.975)*se),
                      p = Pr...z..) %>% filter(!(grepl("PC", variable))) %>% 
  select(variable, HR, LL, UL, p, Beta, se,  Ncases, Ntotal)

ukb <- ukb %>% mutate(cohort = "UKB")

ALL <- bind_rows(kyoto, ukb)
ALL <- ALL %>% mutate(variable = case_when(variable %in% c("HighPRS") ~ "High IPF-PRS",
                                           variable %in% c("telomere_gene", "tel_gene") ~ "Rare telomere variants carrier",
                                           variable %in% c("ShortTS")~ "Short LTL",
                                           variable %in% c("age.x", "AGE") ~ "10 year increase in age at recruitment",
                                           variable %in% c("geneticSex.xMale", "SEX.x") ~ "Male",
                                           variable %in% c("smokingever", "smokingEver") ~ "Ever smoker"))

cov <- unique(ALL$variable)
for(j in c(1:length(cov))){
  tmp <- ALL %>% filter(variable == cov[j]) 
  m1 <- meta::metagen(Beta,
                      se,
                      data=tmp,
                      studlab=paste(cohort),
                      comb.fixed = TRUE,
                      comb.random = TRUE,
                      prediction = FALSE,
                      sm="HR")
  ALL$hetero_p[ALL$variable == cov[j]] <- m1$pval.Q
}


library(forcats)
ALL %>% mutate(name = fct_relevel(variable, "Rare telomere variants carrier", "High IPF-PRS", "Short LTL",
                                  "10 year increase in age at recruitment", "Male", "Ever smoker")) %>%
  ggplot(aes(x=cohort, y=HR, ymin=LL, ymax=UL, color=cohort)) +
  geom_pointrange(aes(col=cohort), lwd=0.8) + geom_hline(aes(fill=cohort), yintercept =1, linetype=2) +
  xlab("") + ylab("Hazard Ratio (95% Confidence Interval)") + geom_text(aes(label=round(HR,1), y=HR, col=cohort, hjust = 0.5, vjust = 1.5), size=4) +
  geom_errorbar(aes(ymin=LL, ymax=UL, col=cohort), width=0.1, cex=1) + scale_y_log10() +
  facet_wrap(~name,  strip.position = 'left', nrow = 10) + theme_minimal() +
  scale_color_manual(labels = c("Kyoto-ILD", "UKB"), values=c("#e4007f", "#0068b7")) + labs(color="Cohort") +
  theme(plot.title=element_text(size=15),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"), 
        legend.text = element_text(size=15,face="bold"), 
        legend.title = element_text(size=15,face="bold"), 
        strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=15,face="bold"))+
  coord_flip() +guides(col = guide_legend(reverse = TRUE))

ggsave(paste0("../repo/IPF_GWAS/Figures/Fig3b.png"), width=25, height=15, units = "cm", dpi=300)

ALL <- ALL %>% select(variable, HR, LL, UL, p, cohort, hetero_p, Ncases, Ntotal)
ALL %>% openxlsx::write.xlsx("../repo/IPF_GWAS/Tables/IPF_survival.xlsx")
