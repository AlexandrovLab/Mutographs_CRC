# Author: Marcos Diaz-Gay
# Date: Aug 5, 2024
# RStudio

library(tidyverse)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv', 
                 check.names = F)
met$country[met$country == 'Czech Republic'] = 'Czechia'

median_no_zeros = function(vector){
    vector[vector==0] = NA
    return(median(vector, na.rm = T))
}


df = met %>% filter(Status == 'MSS')

# Drivers
d_mss = read.delim('../../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/genes_with_driver_muts_MSS.tsv')
mss_driver_genes = colnames(d_mss)[-1]

d_rec_mss = read.delim('../../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/genes_with_recurrent_driver_muts_MSS.tsv',
                       check.names = F)
threshold_rec_muts = 10
rec_muts_to_consider = names(colSums(d_rec_mss[,-1])[colSums(d_rec_mss[,-1]) > threshold_rec_muts])
mss_rec_driver_muts = rec_muts_to_consider

sigsint = c(mss_driver_genes, mss_rec_driver_muts)

dfplot = df %>%
    group_by(age_group) %>%
    summarise(across(all_of(sigsint),median_no_zeros)) %>%
    pivot_longer(cols = sigsint, names_to = 'Signature', values_to = 'Median_activity')

dfplot2 = df %>%
    group_by(age_group) %>%
    summarise(across(all_of(sigsint),~sum(. > 0, na.rm=T))) %>%
    pivot_longer(cols = sigsint, names_to = 'Signature', values_to = 'Prevalence')

dfplot = merge(dfplot, dfplot2)

dfplot3 = df %>%
    group_by(age_group) %>%
    summarise(Total_cases = n()) %>%
    mutate(vartype_group ='DR_RDM')

dfplot = dfplot %>%
    mutate(vartype = ifelse(str_detect(Signature,'_c'),'RDM','DR'),
           vartype_group = 'DR_RDM')

dfplot = merge(dfplot, dfplot3, all.x = T)

dfplot = dfplot %>%
    mutate(Prevalence_perc = Prevalence / Total_cases,
           Signature_label = Signature
    )

dfplot$Signature_label = factor(dfplot$Signature_label,
                                levels = gtools::mixedsort(unique(dfplot$Signature_label),
                                                           decreasing = T))

age_group_order = unique(dfplot$age_group)
dfplot$age_group = factor(dfplot$age_group, levels = age_group_order)

dfplot$vartype_plot = dfplot$vartype
dfplot$vartype_plot[dfplot$vartype=='DR'] = 'Cancer driver genes'
dfplot$vartype_plot[dfplot$vartype=='RDM'] = 'Recurrent driver mutations'


dfplot$vartype_plot = factor(dfplot$vartype_plot,
                             levels = c('Cancer driver genes',
                                        'Recurrent driver mutations'))

dfplot = dfplot %>%
    mutate(tendency = ifelse(Signature_label == 'APC',
                             'Enriched in older patients',
                             'No significant enrichment'))


dfplot %>%
    filter(vartype == 'DR') %>%
ggplot() +
    aes(x=age_group, y = Prevalence_perc,
        group = Signature_label,
        col = tendency) +
    facet_grid(.~vartype_plot) +
    # geom_point(size = 3) +
    geom_line(aes(size = tendency,
                  col = tendency),
              # key_glyph = "point"
              ) +
    geom_text_repel(aes(label = ifelse(age_group=='70+' & tendency=='No significant enrichment',
                                        as.character(Signature_label),NA),
                        ),
                    show.legend = F) +
    geom_label_repel(aes(label = ifelse(age_group=='70+' & tendency!='No significant enrichment',
                                        as.character(Signature_label),NA)),
                     nudge_x = 0.35, min.segment.length = 1, fontface='bold',
                     show.legend = F) +
    theme_bw() +
    labs(x = 'Age of onset',
         y='Prevalence',
         color = 'Age of onset enrichment'
    ) +
    guides(
        color = 'none',
        size = 'none'
        ) +
    scale_size_manual(values = c(1,0.5)) +
    scale_color_manual(values = c('#73C03F','grey80')) +
    scale_x_discrete(expand = c(0.05,0,0.15,0)) +
    scale_y_continuous(labels = scales::label_percent(),
                       limits = c(0,1)) +
    theme(
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 16, face='bold'),
        legend.text = element_text(size = 14),
        legend.position = 'top',
        legend.justification = 'left',
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face='bold'),
        panel.spacing.x = unit(0.75,'lines')
    )
# Exported 550 x 500
