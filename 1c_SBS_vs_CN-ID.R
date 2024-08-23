# Author: Marcos Diaz-Gay
# Date: Jul 22, 2024
# RStudio

library(tidyverse)
library(scales)
library(patchwork)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')

met$Status[met$Status == 'POLD'] = 'POLD1'
met$Status[met$Status == 'BRCA'] = 'HRD'

all_status = met %>%
    group_by(Status) %>%
    summarise(n = n()) %>%
    mutate(Status_n = paste0(Status, ' (n=', n, ')'))
all_status$Status_n[all_status$Status=='Chemotherapy'] = "Chemotherapy\n(n=2)"

met = met %>%
    left_join(all_status)

met$Status_n = factor(met$Status_n,
                    levels = c("MSS (n=802)","MSI (n=153)",
                               "POLE (n=10)","POLD1 (n=3)",
                               "MUTYH (n=1)", "NTHL1 (n=2)",
                               "OGG1 (n=1)", "HRD (n=7)",
                               "Chemotherapy\n(n=2)"))

p1 = met %>%
        # filter(!Status %in% c('Chemotherapy')) %>%
    ggplot() +
    aes(x = PGA, y = TMB_SBS) +
    geom_point(aes(col = Status_n), alpha = 0.75, size = 2) +
    labs(y = 'Number of single base substitutions (SBS)',
         x = 'Percentage of genome aberrated (PGA)',
         title = 'Molecular subgroups',
         subtitle = 'SBS - PGA') +
    scale_x_continuous(label = label_percent()) +
    scale_y_log10(labels = label_comma()) +
    scale_color_manual(values = c('grey', '#e03f1f', '#1ecbe1', 'blue', 
                                  'darkgreen', 'green3', 'lightgreen', 'pink1',
                                  'purple')) +
    theme_bw() +
    theme(text = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'left',
          plot.title = element_text(face='bold', size=16),
          plot.subtitle = element_text(size = 14, face = 'bold', hjust = 0.5),
          legend.title = element_text(face='bold', size=16),
          legend.text = element_text(size=14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    guides(col = 'none')

p2 = met %>%
        # filter(!Status %in% c('Chemotherapy')) %>%
    ggplot() +
    aes(x = TMB_ID, y = TMB_SBS) +
    geom_point(aes(col = Status_n), alpha = 0.75, size = 2) +
    labs(y = '',
         x = 'Number of indels (ID)',
         col = 'Subgroup',
         subtitle = 'SBS - ID') +
    scale_x_log10(labels = label_comma()) +
    scale_y_log10(labels = label_comma()) +
    scale_color_manual(values = c('grey', '#e03f1f', '#1ecbe1', 'blue', 
                                  'darkgreen', 'green3', 'lightgreen','pink1',
                                  'purple')) +
    theme_bw() +
    theme(text = element_text(size = 14),
          legend.position = 'right',
          legend.justification = 'left',
          plot.title = element_text(face='bold', size=16),
          plot.subtitle = element_text(size = 14, face = 'bold', hjust = 0.5),
          legend.title = element_text(size=14),
          legend.text = element_text(size=14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))
    
p1+p2 
# Exported 1600 x 500
