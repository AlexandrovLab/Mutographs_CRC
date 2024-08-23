# Author: Marcos Diaz-Gay
# Date: Aug 4, 2024
# RStudio

library(tidyverse)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')
met = met %>%
    filter(Status %in% c('MSS', 'MSI')) %>%
    mutate(Status = factor(Status, levels = c('MSS', 'MSI')),
           sex = factor(sex),
           sex = relevel(sex, ref = 'Male'),
           country = factor(country),
           country = relevel(country, ref = 'Brazil'),
           tumorsite_group = factor(tumorsite_group),
           tumorsite_group = relevel(tumorsite_group,
                                     ref = 'Proximal colon'))

signature = 'SBS18_c'
signature = 'SBS5_c'
signature = 'ID1_c'
signature = 'ID2_c'

(FC = mean(met[met$Status=='MSI',signature]) /
    mean(met[met$Status=='MSS',signature]))

(sum = summary(lm(get(signature) ~ Status + age_diag + sex + tumorsite_group + country + purity,
           data = met)))
sum$coefficients['StatusMSI','Pr(>|t|)']

wilcox.test(get(signature) ~ Status, met)
