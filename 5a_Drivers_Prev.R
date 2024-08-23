# Author: Marcos Diaz-Gay
# Date: Jul 21, 2024
# RStudio

library(tidyverse)
library(cowplot)
library(scales)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')

first_driver_gene = grep('APC',colnames(met))
last_driver_gene = grep('ACSL',colnames(met))
cols_driver_genes = first_driver_gene:last_driver_gene
driver_genes = colnames(met)[cols_driver_genes]

d_mss = read.delim('../../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/genes_with_driver_muts_MSS.tsv')
mss_driver_genes = colnames(d_mss)[-1]
driver_genes = mss_driver_genes

df = met %>%
    select(donor_id, Status, all_of(driver_genes)) %>%
    filter(Status %in% c('MSS', 'MSI')) %>%
    pivot_longer(cols = all_of(driver_genes),
                 names_to = 'Gene',
                 values_to = 'Mut_number') %>%
    mutate(Mut = Mut_number>0)

df$Driver_gene_Status = 'Known CRC driver gene'
df$Driver_gene_Status[df$Gene %in% c('MED12', 'NCOR1')] = 'Known cancer driver gene'
df$Driver_gene_Status[df$Gene %in% c('CCR4')] = 'Novel cancer driver gene'



total_cases_per_status = met %>%
    filter(Status %in% c('MSS', 'MSI')) %>%
    group_by(Status) %>%
    summarise(Total_samples = n())

dfplot = df %>%
    group_by(Gene, Status, Driver_gene_Status) %>%
    summarise(Mut_samples = sum(Mut)) %>%
    left_join(total_cases_per_status) %>%
    mutate(Mut_samples_prev = Mut_samples / Total_samples) %>%
    mutate(Mut_samples_prev = ifelse(Status == 'MSI', -Mut_samples_prev,
                                     Mut_samples_prev))

dfplot$Gene = factor(dfplot$Gene,
                     levels = unique(dfplot %>% arrange(desc(Mut_samples_prev)) %>% pull(Gene)))

dfplot$Driver_gene_Status = factor(dfplot$Driver_gene_Status,
                                   levels = c('Known CRC driver gene',
                                              'Known cancer driver gene',
                                              'Novel cancer driver gene'))

threshold = 0
mss_cases_above_th = dfplot %>%
    filter(Status == 'MSS') %>%
    filter(Mut_samples_prev>threshold) %>%
    pull(Gene)
msi_cases_above_th = dfplot %>%
    filter(Status == 'MSI') %>%
    filter(abs(Mut_samples_prev)>threshold) %>%
    pull(Gene)
cases_above_th = unique(c(mss_cases_above_th, msi_cases_above_th))

dfplot %>%
    filter(Gene %in% cases_above_th) %>%
    filter(Status == 'MSS') %>%
ggplot() +
    aes(x = Gene, y = Mut_samples_prev, fill = Driver_gene_Status) +
    geom_col() +
    # geom_hline(aes(yintercept = 0)) +
    # geom_vline(aes(xintercept = 20.5), linetype = 'dashed', col = 'grey40') +
    geom_text(aes(label = paste0(Mut_samples)),#,'\n(',abs(round(Mut_samples_prev,3)*100),'%)')),
              # nudge_y = ifelse(dfplot$Status == 'MSS', 0.05,-0.05),
              nudge_y = 0.05,
              size = 4) +
    coord_cartesian(ylim = c(0,1)) +
    theme_bw() +
    theme(plot.title = element_text(size = 18),
          plot.subtitle = element_text(size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 16, face ='bold'),
          legend.text = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'left',
          axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
          panel.grid.major.y = element_line(colour = 'grey90')) +
    scale_y_continuous(limits = c(0,1.1), breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
                       labels = c('100%', '75%', '50%', '25%', '0%',
                                  '25%', '50%', '75%', '100%'),) +
    scale_fill_manual(values = c('#d3de5f', '#277ad8', '#e03f1f')) +
    labs(x = '',
         y = 'Percentage of mutated samples',
         fill = 'Cancer driver genes')
# Exported 1600 x 500
