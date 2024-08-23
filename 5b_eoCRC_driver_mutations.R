# Author: Marcos Diaz-Gay
# Date: Aug 5, 2024
# RStudio

library(tidyverse)
library(scales)
library(ggbeeswarm)
library(rstatix)
library(ggpubr)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')
met = met %>%
    filter(Status == 'MSS')
met$sex = factor(met$sex,
                    levels = c('Male','Female'))
met$tumorsite_group = factor(met$tumorsite_group,
                    levels = c('Proximal colon', 'Distal colon', 'Rectum'))
met$country = factor(met$country)
met$country = relevel(met$country, ref = 'Brazil')
met$age_eo = ifelse(met$age_eo == '0-49', '0-49 (n=97)',
                    '50+ (n=705)')

## TMB SBS
model = lm(total_driver_mutations ~ age_eo + sex + tumorsite_group + country + purity, data = met)
summary(model)
p_value_adj = summary(model)$coefficients['age_eo50+ (n=705)', 'Pr(>|t|)']

met$age_eo = factor(met$age_eo,
                    levels = c('0-49 (n=97)','50+ (n=705)'))
met %>%
    mutate(analysis = 'Total driver mutations') %>%
ggplot() +
    aes(x=age_eo, y = total_driver_mutations,
        fill = age_eo) +
    facet_wrap(.~analysis) +
    geom_quasirandom(aes(col = age_eo), size=1) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    stat_pwc(
             aes(label = ifelse(
                 after_stat(p) > 0,
                 sprintf("p = %5.1e", p_value_adj),
                 ))) +
    scale_y_continuous(#transform = 'log10',
        labels = scales::label_comma(),
        expand = c(0.05,0,0.1,0)
        ) +
    scale_fill_manual(values = scales::seq_gradient_pal("#8C3FC0", "#73C03F", "Lab")(seq(0,1,length.out=2))) +
    scale_colour_manual(values = scales::seq_gradient_pal("#8C3FC0", "#73C03F", "Lab")(seq(0,1,length.out=2))) +
    theme_bw() +
    labs(title = 'Age of onset distribution',
         subtitle = 'Adjusted by sex, country, tumor subsite,\nand purity',
         x = 'Age of onset',
         y = 'Number of driver mutations') +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          strip.text = element_text(size = 14,face = 'bold',hjust = 0.5),
          strip.background = element_blank(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    guides(fill = 'none',
           col = 'none')
# 400 x 500

(fold_change = (summary(model)$coefficients['age_eo50+ (n=705)', 'Estimate'] +
        mean(met$total_driver_mutations[met$age_eo=='0-49 (n=97)'])) /
        mean(met$total_driver_mutations[met$age_eo=='0-49 (n=97)']) )
