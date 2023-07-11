ukb <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/UKB_scratch/UKB_assoc.summary.tsv")

ukb <- ukb %>% filter(outcome %in% c("IPF", "RA_ILD", "SjS_ILD",
                                     "SSc_ILD", "ANCA_ILD")) %>% mutate(Beta = case_when(variable == "AGE" ~ 10*Beta,
                                                                    variable == "TS_sd" ~ -1*Beta,
                                                                    TRUE ~ Beta),
                                                   se = ifelse(variable == "AGE", 10*se, se))

ukb <- ukb %>% mutate(OR = exp(Beta),
                      LL = exp(Beta + qnorm(0.025)*se),
                      UL = exp(Beta + qnorm(0.975)*se),
                      p = pval) %>% filter(!(grepl("PC", variable))) %>% 
  select(variable, OR, LL, UL, p, Beta, se, Ncases, Ncontrols, outcome)

ukb <- ukb %>% mutate(cohort = "UKB")

ALL <- ukb %>% mutate(variable = case_when(variable %in% c("PRS", "prs") ~ "1SD increase in PRS",
                                           variable %in% c("telomere_gene", "tel_gene") ~ "Rare telomere variants carrier",
                                           variable %in% c("tel_sd", "TS_sd")~ "1SD decrease in age/sex adjusted LTL",
                                           variable %in% c("age.x", "AGE") ~ "10 year increase in age at recruitment",
                                           variable %in% c("geneticSex.xMale", "SEX.y") ~ "Male",
                                           variable %in% c("smokingever", "smokingEver") ~ "Ever smoker"))

ALL$OR[14] <- NA
ALL$LL[14] <- NA
ALL$UL[14] <- NA


library(forcats)
ALL %>% mutate(name = fct_relevel(variable, "Rare telomere variants carrier", "1SD increase in PRS", "1SD decrease in age/sex adjusted LTL",
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

ggsave(paste0("../repo/IPF_GWAS/Figures/Fig4.png"), width=25, height=30, units = "cm", dpi=300)

ALL <- ukb %>% mutate(variable = case_when(variable %in% c("PRS", "prs") ~ "1SD increase in PRS",
                                           variable %in% c("telomere_gene", "tel_gene") ~ "Rare telomere variants carrier",
                                           variable %in% c("tel_sd", "TS_sd")~ "1SD decrease in age/sex adjusted LTL",
                                           variable %in% c("age.x", "AGE") ~ "10 year increase in age",
                                           variable %in% c("geneticSex.xMale", "SEX.y") ~ "Male",
                                           variable %in% c("smokingever", "smokingEver") ~ "Ever smoker"))

ALL <- ALL %>% select(variable, OR, LL, UL, p, Ncases, Ncontrols, outcome)
ALL %>% openxlsx::write.xlsx("../repo/IPF_GWAS/Tables/NonIPF_susceptibility.xlsx")

