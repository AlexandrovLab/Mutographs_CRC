# Author: Marcos Diaz-Gay
# Date: Dec 24, 2024
# RStudio

library(tidyverse)
library(scales)
library(ggbeeswarm)
library(rstatix)
library(ggpubr)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')
met = met %>%
    filter(Status == 'MSS')

ancestry = read.csv('../../../../../../../../../Restricted/Mutographs/CRC_Ancestry_2024-11-20/Ancestry/admixture_pc.csv')
anc_prop = ancestry %>%
    select(Germline_ID, EAS, AFR, EUR)
ancestry = ancestry %>%
    select(Germline_ID, PC1, PC2, PC3, PC4, PC5) %>%
    mutate(donor_id = paste0(str_sub(Germline_ID,1,7),'a')) %>%
    select(-Germline_ID)

met = met %>%
    left_join(ancestry)



anc_prop$ancestry = NA
for (i in 1:nrow(anc_prop)){
    row = anc_prop[i,2:4]
    max_anc = names(which.max(row))
    anc_selected = ifelse(sum(row>0.2) >1,
                          # paste0(max_anc,'_admix'),
                          'ADMIX',
                          max_anc)
    anc_prop$ancestry[i] = anc_selected
}
anc_prop = anc_prop %>%
    mutate(donor_id = paste0(str_sub(Germline_ID,1,7),'a')) %>%
    select(-Germline_ID)

met = met %>%
    left_join(anc_prop)


met$sex = factor(met$sex,
                    levels = c('Male','Female'))
met$tumorsite_group = factor(met$tumorsite_group,
                    levels = c('Proximal colon', 'Distal colon', 'Rectum'))
met$country = factor(met$country)
met$country = relevel(met$country, ref = 'Brazil')
met$age_eo = ifelse(met$age_eo == '0-49', '0-49 (n=97)',
                    '50+ (n=705)')

## TMB SBS
model = lm(log10(TMB_SBS) ~ age_eo + sex + tumorsite_group + PC1 + PC2 + PC3 + PC4 + PC5 + purity, data = met)
summary(model)
p_value_adj = summary(model)$coefficients['age_eo50+ (n=705)', 'Pr(>|t|)']

met$age_eo = factor(met$age_eo,
                    levels = c('0-49 (n=97)','50+ (n=705)'))
met %>%
    mutate(analysis = 'SBS') %>%
ggplot() +
    aes(x=age_eo, y = TMB_SBS,
        fill = age_eo) +
    facet_wrap(.~analysis) +
    geom_quasirandom(aes(col = age_eo), size=1) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    stat_pwc(
             aes(label = ifelse(
                 after_stat(p) > 0,
                 sprintf("p = %5.3f", p_value_adj),
                 ))) +
    scale_y_continuous(transform = 'log10',
        labels = scales::label_comma(),
        expand = c(0.05,0,0.1,0)
        ) +
    scale_fill_manual(values = scales::seq_gradient_pal("#8C3FC0", "#73C03F", "Lab")(seq(0,1,length.out=2))) +
    scale_colour_manual(values = scales::seq_gradient_pal("#8C3FC0", "#73C03F", "Lab")(seq(0,1,length.out=2))) +
    theme_bw() +
    labs(title = 'MSS molecular subgroup',
         subtitle = 'Adjusted by sex, tumor subsite,\npurity, and genetic ancestry',
         x = 'Age of diagnosis',
         y = 'Number of mutations') +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          strip.text = element_text(size = 14,face = 'bold',hjust = 0.5),
          strip.background = element_blank(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    guides(fill = 'none',
           col = 'none')
# 380 x 450

fold_change_SBS = 10^summary(model)$coefficients['age_eo50+ (n=705)', 'Estimate']


## TMB ID
model = lm(log10(TMB_ID) ~ age_eo + sex + tumorsite_group + PC1 + PC2 + PC3 + PC4 + PC5 + purity, data = met)
summary(model)
p_value_adj = summary(model)$coefficients['age_eo50+ (n=705)', 'Pr(>|t|)']

met %>%
    mutate(analysis = 'ID') %>%
    ggplot() +
    aes(x=age_eo, y = TMB_ID,
        fill = age_eo) +
    facet_wrap(.~analysis) +
    geom_quasirandom(aes(col = age_eo), size=1) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    stat_pwc(
        aes(label = ifelse(
            after_stat(p) > 0,
            sprintf("p = %5.3f", p_value_adj),
        ))) +
    scale_y_continuous(transform = 'log10',
                       labels = scales::label_comma(),
                       expand = c(0.05,0,0.1,0)
    ) +
    scale_fill_manual(values = scales::seq_gradient_pal("#8C3FC0", "#73C03F", "Lab")(seq(0,1,length.out=2))) +
    scale_colour_manual(values = scales::seq_gradient_pal("#8C3FC0", "#73C03F", "Lab")(seq(0,1,length.out=2))) +
    theme_bw() +
    labs(title = '',
         subtitle = '\n',
         x = 'Age of diagnosis',
         y = '') +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14,face = 'bold',hjust = 0.5),
          strip.text = element_text(size = 14,face = 'bold',hjust = 0.5),
          strip.background = element_blank(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    guides(fill = 'none',
           col = 'none')
# 380 x 450

fold_change_ID = 10^summary(model)$coefficients['age_eo50+ (n=705)', 'Estimate']

