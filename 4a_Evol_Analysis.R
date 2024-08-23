# Author: Marcos Diaz-Gay
# Date: Jul 21, 2024
# RStudio

library(tidyverse)
library(ggbeeswarm)
library(scales)

####################################
analysis_plot = 'Early clonal vs. late clonal'
# analysis_plot = 'Clonal vs. subclonal'
varianttype_plot = 'SBS'
# varianttype_plot = 'ID'
####################################


met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')

color_for_signatures = tibble(sig = c('SBS1',   'SBS2',   'SBS5',
                                            'SBS8',   'SBS17a', 'SBS17b',
                                            'SBS18',  'SBS34',  'SBS88',
                                            'SBS93',  'SBS_F',  'SBS_H',
                                            'SBS_M',
                                            'ID1',    'ID2',    'ID4',
                                            'ID6',    'ID8',    'ID9',
                                            'ID10',   'ID14',   'ID18',
                                            'ID_J'),
                              color_sig = c('grey50', 'grey50', 'grey50',
                                            'grey50', 'grey50', 'grey50',
                                            'grey50', 'grey50', 'skyblue2',
                                            'grey50', 'grey50', 'grey50',
                                            'grey50',
                                            'grey50', 'grey50', 'grey50',
                                            'grey50', 'grey50', 'grey50',
                                            'grey50', 'grey50', 'skyblue2',
                                            'grey50')
                              )



sbs_e_l = read.csv('../data_for_figures/Mariya_figure_data/evolutionary_plot_data_v2_2024-07-17/SBS_early_late.csv')
sbs_e_l$analysis = 'Early clonal vs. late clonal'
sbs_e_l$varianttype = 'SBS'
sbs_c_s = read.csv('../data_for_figures/Mariya_figure_data/evolutionary_plot_data_v2_2024-07-17/SBS_clonal_subclonal.csv')
sbs_c_s$analysis = 'Clonal vs. subclonal'
sbs_c_s$varianttype = 'SBS'
id_e_l = read.csv('../data_for_figures/Mariya_figure_data/evolutionary_plot_data_for_marcos_2024-05-10/INDEL_early_late.csv')
id_e_l$analysis = 'Early clonal vs. late clonal'
id_e_l$varianttype = 'ID'
id_c_s = read.csv('../data_for_figures/Mariya_figure_data/evolutionary_plot_data_for_marcos_2024-05-10/INDEL_clonal_subclonal.csv')
id_c_s$analysis = 'Clonal vs. subclonal'
id_c_s$varianttype = 'ID'

df = rbind(sbs_e_l, sbs_c_s, id_e_l, id_c_s)

df$sig[df$sig=='SBS288F'] = 'SBS_F'
df$sig[df$sig=='SBS288H'] = 'SBS_H'
df$sig[df$sig=='SBS288M'] = 'SBS_M'
df$sig[df$sig=='ID83J'] = 'ID_J'

signature_order = sbs_e_l %>% group_by(sig) %>% summarise(fc=median(fold.change, na.rm=T)) %>% arrange(fc) %>% pull(sig)
signature_order = c(signature_order, 'SBS_M') # Because SBS_M is not in sbs_e_l
signature_order_id = id_e_l %>% group_by(sig) %>% summarise(fc=median(fold.change, na.rm=T)) %>% arrange(fc) %>% pull(sig)
signature_order = c(signature_order, signature_order_id)

signature_order[signature_order=='SBS288F'] = 'SBS_F'
signature_order[signature_order=='SBS288H'] = 'SBS_H'
signature_order[signature_order=='SBS288M'] = 'SBS_M'
signature_order[signature_order=='ID83J'] = 'ID_J'

df = df %>%
    rename(donor_id = X) %>%
    left_join(color_for_signatures) %>%
    mutate(sig = factor(sig, levels = signature_order)) %>%
    left_join(met)

dfplot = df %>%
    filter(
        analysis == analysis_plot,
        varianttype == varianttype_plot
           )
p = ggplot(dfplot) +
    aes(x = sig, y =fold.change, fill = color_sig) +
    geom_quasirandom(aes(col = color_sig), size=1) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    scale_fill_identity()+
    scale_color_identity() +
    {if(unique(dfplot$varianttype)=='SBS')scale_y_continuous(limits = c(0.5-0.1,
                                                                        2+0.1))} +
    {if(unique(dfplot$varianttype)=='ID')scale_y_continuous(limits = c(0.6-0.1,
                                                                       1.6+0.1))} +
    theme_bw() +
    labs(title = ifelse(unique(dfplot$varianttype)=='ID',' ',
                        unique(dfplot$analysis)),
         subtitle = paste0(unique(dfplot$varianttype), ' signatures'),
         x = '',
         y = 'Fold change') +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 30, hjust=1),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = 'bold')
          ) +
    guides(fill = 'none',
           col = 'none'
           )

if(analysis_plot == 'Early clonal vs. late clonal'){
    p +
        geom_text(aes(label='Early clonal', x = -Inf, y=-Inf),
              check_overlap = TRUE, hjust = -0.2, vjust = -1.5,
              size = 5, fontface = 'bold') +
        geom_text(aes(label='Late clonal', x = -Inf, y=Inf),
                  check_overlap = TRUE, hjust = -0.2, vjust = 2.2,
                  size = 5, fontface = 'bold')
}
if(analysis_plot == 'Clonal vs. subclonal'){
    p+
        geom_text(aes(label='Clonal', x = -Inf, y=-Inf),
              check_overlap = TRUE, hjust = -0.2, vjust = -1.5,
              size = 5, fontface = 'bold') +
        geom_text(aes(label='Subclonal', x = -Inf, y=Inf),
              check_overlap = TRUE, hjust = -0.2, vjust = 2.2,
              size = 5, fontface = 'bold')
}
# 750 x 500



# EO vs. LO samples plot
p = ggplot(dfplot) +
    aes(x = sig, y =fold.change, fill = age_eo) +
    # facet_wrap(reorder(varianttype, desc(varianttype))~
    #                reorder(analysis,desc(analysis)), scales = 'free') +
    geom_hline(aes(yintercept = 1), linetype = 'dashed') +
    geom_quasirandom(aes(col = age_eo), size=1,
                     dodge.width = 0.8, width = 0.2) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c('#8C3FC0','#73C03F'))+
    scale_color_manual(values = c('#8C3FC0','#73C03F'))+
    {if(unique(dfplot$varianttype)=='SBS')scale_y_continuous(limits = c(0.5,2))} +
    {if(unique(dfplot$varianttype)=='ID')scale_y_continuous(limits = c(0.6,1.6))} +
    theme_bw() +
    labs(title = ifelse(unique(dfplot$varianttype)=='ID',' ',
                        unique(dfplot$analysis)),
         subtitle = paste0(unique(dfplot$varianttype), ' signatures'),
         x = '',
         y = 'Fold change') +
    {if(unique(dfplot$varianttype)=='SBS')labs(col = 'Age')} +
    {if(unique(dfplot$varianttype)=='ID')guides(col = 'none')} +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 30, hjust=1),
          legend.position = 'top',
          legend.justification = 'left',
          legend.text = element_text(size=14),
          legend.title = element_text(size=14),
          # panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = 'bold')
    ) +
    guides(
        fill = 'none',
    )
if(analysis_plot == 'Early clonal vs. late clonal'){
    p +
        geom_text(aes(label='Early clonal', x = -Inf, y=-Inf),
                  check_overlap = TRUE, hjust = -0.2, vjust = -1.5,
                  size = 5, fontface = 'bold') +
        geom_text(aes(label='Late clonal', x = -Inf, y=Inf),
                  check_overlap = TRUE, hjust = -0.2, vjust = 2.2,
                  size = 5, fontface = 'bold')
}
if(analysis_plot == 'Clonal vs. subclonal'){
    p+
        geom_text(aes(label='Clonal', x = -Inf, y=-Inf),
                  check_overlap = TRUE, hjust = -0.2, vjust = -1.5,
                  size = 5, fontface = 'bold') +
        geom_text(aes(label='Subclonal', x = -Inf, y=Inf),
                  check_overlap = TRUE, hjust = -0.2, vjust = 2.2,
                  size = 5, fontface = 'bold')
}
# 750 x 550 (SBS)
# 750 x 500 (ID)
