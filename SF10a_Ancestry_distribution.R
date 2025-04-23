# Author: Marcos Diaz-Gay
# Date: Dec 24, 2024
# RStudio

library(tidyverse)
library(scales)
library(ggrepel)
library(broom)
library(logistf)

###################
status = 'MSS'
# status = 'MSI'
variable = 'ancestry'
###################

if (variable == 'country') {
    title = 'Country enrichment'
    subtitle = 'Adjusted by age, sex, tumor subsite, and purity'
}
if (variable == 'ancestry') {
    title = 'Genetic ancestry enrichment'
    subtitle = 'Adjusted by age, sex, tumor subsite, and purity'
}

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')

df = met #%>% filter(Status == status)
df$country[df$country=='Czech Republic'] = 'Czechia'

ancestry = read.csv('../../../../../../../../../Restricted/Mutographs/CRC_Ancestry_2024-11-20/Ancestry/admixture_pc.csv')
anc_prop = ancestry %>%
    select(Germline_ID, EAS, AFR, EUR)

ancestry = ancestry %>%
    select(Germline_ID, PC1, PC2, PC3, PC4, PC5) %>%
    mutate(donor_id = paste0(str_sub(Germline_ID,1,7),'a')) %>%
    select(-Germline_ID)

df = df %>%
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

df = df %>%
    left_join(anc_prop)

df = df %>%
    filter(!is.na(ancestry))


# Heatmap

df_country = df %>%
    group_by(country) %>%
    summarise(all_country = n())

df %>%
    group_by(ancestry,country) %>%
    summarise(count = n()) %>%
    left_join(df_country) %>%
    mutate(rel_count = count / all_country,
           col_label = ifelse(rel_count>0.5,
                              'white','black')) %>%
ggplot() +
    aes(y = reorder(ancestry, desc(ancestry)), x = country, fill = rel_count) +
    geom_tile() +
    geom_text(aes(label = paste0(count, ' (',
                                 round(rel_count * 100,1),'%)'),
                  col = col_label)) +
    scale_fill_viridis_c(direction = -1, option = 'G') +
    scale_color_manual(values = c('black','white')) +
    theme_bw() +
    labs(x = '', y = '',
         title = 'Genetic ancestry distribution by country') +
    guides(col = 'none',
           fill = 'none') +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(size = 13),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 14, face = 'bold'),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 13))
# Exported 1300 x 400
