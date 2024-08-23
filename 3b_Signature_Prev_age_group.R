# Author: Marcos Diaz-Gay
# Date: Jul 21, 2024
# RStudio

library(tidyverse)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')
met$country[met$country == 'Czech Republic'] = 'Czechia'

median_no_zeros = function(vector){
    vector[vector==0] = NA
    return(median(vector, na.rm = T))
}


df = met %>% filter(Status == 'MSS')

# Signatures
all_sigsint = grep('_c$', colnames(met), value = T)
sigsint = all_sigsint[colSums(apply(df[,all_sigsint],2,is.na))!=nrow(df)]
sigsint = sigsint[!str_detect(sigsint, 'Unknown')]
only_sigs = sigsint

# Drivers
d_mss = read.delim('../../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/genes_with_driver_muts_MSS.tsv')
mss_driver_genes = colnames(d_mss)[-1]
sigsint = c(sigsint, mss_driver_genes)



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
    mutate(vartype_group ='SB_DB_ID_DR')

dfplot4 = df %>%
    group_by(age_group) %>%
    summarise(Total_cases = sum(!is.na(CN1_c)))%>%
    mutate(vartype_group ='CN')

dfplot5 = df %>%
    group_by(age_group) %>%
    summarise(Total_cases = sum(!is.na(SV1_c)))%>%
    mutate(vartype_group ='SV')

dfplot3 = rbind(dfplot3, dfplot4)
dfplot3 = rbind(dfplot3, dfplot5)


dfplot = dfplot %>%
    mutate(vartype = ifelse(Signature %in% mss_driver_genes, 'DR', substr(Signature,1,2)),
           vartype_group = ifelse(vartype == 'CN', 'CN',
                                  ifelse(vartype == 'SV', 'SV',
                                         'SB_DB_ID_DR')))

dfplot = merge(dfplot, dfplot3, all.x = T)

dfplot = dfplot %>%
    mutate(Prevalence_perc = Prevalence / Total_cases,
           Signature_label = sapply(strsplit(Signature, '_c'), '[', 1),
    )

dfplot$Signature_label[dfplot$Signature_label == 'SBS288F_MSS'] = 'SBS_F'
dfplot$Signature_label[dfplot$Signature_label == 'SBS288H_MSS'] = 'SBS_H'
dfplot$Signature_label[dfplot$Signature_label == 'SBS288M_MSS'] = 'SBS_M'
dfplot$Signature_label[dfplot$Signature_label == 'SBS288O_MSS'] = 'SBS_O'
dfplot$Signature_label[dfplot$Signature_label == 'ID83J_MSS'] = 'ID_J'
dfplot$Signature_label[dfplot$Signature_label == 'CN68F_MSS'] = 'CN_F'
dfplot$Signature_label[dfplot$Signature_label == 'SV38B_MSS'] = 'SV_B'
dfplot$Signature_label[dfplot$Signature_label == 'SV38D_MSS'] = 'SV_D'
dfplot$Signature_label = factor(dfplot$Signature_label,
                                levels = gtools::mixedsort(unique(dfplot$Signature_label),
                                                           decreasing = T))

age_group_order = unique(dfplot$age_group)
dfplot$age_group = factor(dfplot$age_group, levels = age_group_order)

dfplot$vartype_plot = dfplot$vartype
dfplot$vartype_plot[dfplot$vartype=='SB'] = 'SBS'
dfplot$vartype_plot[dfplot$vartype=='DB'] = 'DBS'
dfplot$vartype_plot[dfplot$vartype=='ID'] = 'ID'
dfplot$vartype_plot[dfplot$vartype=='CN'] = 'CN'
dfplot$vartype_plot[dfplot$vartype=='SV'] = 'SV'
dfplot$vartype_plot[dfplot$vartype=='DR'] = 'Driver genes'
dfplot$vartype_plot = factor(dfplot$vartype_plot,
                             levels = c('SBS', 'DBS',
                                        'ID', 'CN',
                                        'SV', 'Driver genes'))

dfplot = dfplot %>%
mutate(tendency = ifelse(Signature_label %in% c('SBS1','SBS5', 'ID1', 'ID4','ID9','APC', 'ID2','ID10'),
                         'Enriched in older patients',
                         ifelse(Signature_label %in% c('SBS88','SBS_M','ID14','ID18'),
                                'Enriched in younger patients','No significant enrichment')))


dfplot %>%
    filter(Signature_label != 'Unknown') %>%
    filter(!vartype %in% c('DR')) %>%
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
    # geom_label_repel(aes(label = ifelse(age_group=='70+', #& tendency!='no_sign',
    #                                     as.character(Signature_label),NA)),
    #                  nudge_x = 0.35) +
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
         y='Signature prevalence',
         color = 'Age of onset enrichment'
    ) +
    guides(
        # color= guide_legend(override.aes = list(linewidth=1.5)),
        color = 'none',
        size = 'none'
        ) +
    scale_size_manual(values = c(1,1,0.5)) +
    # scale_size_manual(values = c(0.5)) +
    scale_color_manual(values = c('#73C03F','#8C3FC0','grey80')) +
    # scale_color_manual(values = c('grey80')) +
    scale_x_discrete(expand = c(0.05,0,0.15,0)) +
    scale_y_continuous(labels = scales::label_percent(),
                       limits = c(0,1)) +
    theme(#plot.margin = margin(5.5,5.5,b=45,5.5),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 14),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 16, face='bold'),
        legend.text = element_text(size = 14),
        legend.position = 'top',
        legend.justification = 'left',
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face='bold'),
        panel.spacing.x = unit(0.75,'lines')
    )
# Exported 1600 x 500
