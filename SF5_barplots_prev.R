# Author: Marcos Diaz-Gay
# Date: Dec 30, 2024
# RStudio

library(tidyverse)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')
met$country[met$country=='Czech Republic'] = 'Czechia'

met = met %>%
    filter(Status == 'MSS')

all_sigsint = grep('_c_rel$', colnames(met), value = T)
sigsint = all_sigsint[colSums(apply(met[,all_sigsint],2,is.na))!=nrow(met)]

################################################################################
##### SBS
vartype = 'SBS'

sigs = sigsint[grep('SBS', sigsint)]

all_country = met %>%
    group_by(country) %>%
    summarise(total_country = n())
prev_sigs_by_country = met %>%
    group_by(country) %>%
    summarise(across(all_of(sigs), ~ sum(. > 0))) %>%
    pivot_longer(cols = -country, names_to = "signature", values_to = "prev") %>%
    left_join(all_country) %>%
    mutate(rel_prev = prev / total_country)


prev_sigs_by_country %>%
    filter(signature != 'SBS_Unknown_c_rel') %>%
    mutate(signature_label = sapply(strsplit(signature,'_c'),'[',1)) %>%
ggplot() +
    aes(x = reorder(signature_label, desc(rel_prev)), y = rel_prev, fill = country) +
    geom_col(position = 'dodge') +
    labs(y = 'Proportion of cases with the signaturee', x= '',
         fill = '',
         title = 'Prevalence of SBS mutational signatures by country') +
    theme_bw() +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(size = 13),#, angle = 45, hjust = 1),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 14, face = 'bold'),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 13),
          legend.position = 'top') +
    scale_fill_manual(values = c(
        "blue","#B9E769", "#59A14F","#76B7B2",
        "#FF9DA7", "#9C755F", "#A5D8F3",
         "#EDC948","#E15759",
        "#B07AA1",  'grey90', 'grey30'
    )) +
    guides(fill = guide_legend(nrow = 1))
# 1600 x 400


################################################################################
##### ID / DBS / CN / SV
# vartype = 'ID'
# vartype = 'DBS'
# vartype = 'CN'
vartype = 'SV'

sigs = sigsint[grep(vartype, sigsint)]

if (vartype %in% c('ID', 'DBS')) {
    all_country = met %>%
        group_by(country) %>%
        summarise(total_country = n())
    prev_sigs_by_country = met %>%
        group_by(country) %>%
        summarise(across(all_of(sigs), ~ sum(. > 0))) %>%
        pivot_longer(cols = -country, names_to = "signature", values_to = "prev") %>%
        left_join(all_country) %>%
        mutate(rel_prev = prev / total_country)
} else {
    
    all_country = met %>%
        filter(!is.na(CN1_c)) %>%
        group_by(country) %>%
        summarise(total_country = n())
    prev_sigs_by_country = met %>%
        group_by(country) %>%
        # summarise(across(all_of(sigs), ~ sum(. > 0))) %>%
        summarise(across(all_of(sigs), ~ sum(. > 0, na.rm = T))) %>%
        pivot_longer(cols = -country, names_to = "signature", values_to = "prev") %>%
        left_join(all_country) %>%
        mutate(rel_prev = prev / total_country)
    
}

prev_sigs_by_country %>%
    filter(signature != paste0(vartype, '_Unknown_c_rel')) %>%
    mutate(signature_label = sapply(strsplit(signature,'_c'),'[',1)) %>%
    ggplot() +
    aes(x = reorder(signature_label, desc(rel_prev)), y = rel_prev, fill = country) +
    geom_col(position = 'dodge') +
    labs(y = 'Proportion of cases\nwith the signature', x= '',
         fill = '',
         title = paste0('Prevalence of ',
                        vartype,' mutational signatures by country')) +
    theme_bw() +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(size = 13),#, angle = 45, hjust = 1),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 14, face = 'bold'),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 13),
          legend.position = 'top') +
    scale_fill_manual(values = c(
        "blue","#B9E769", "#59A14F","#76B7B2",
        "#FF9DA7", "#9C755F", "#A5D8F3",
        "#EDC948","#E15759",
        "#B07AA1",  'grey90', 'grey30'
    )) +
    guides(fill = guide_legend(nrow = 2))
# 800 x 400
