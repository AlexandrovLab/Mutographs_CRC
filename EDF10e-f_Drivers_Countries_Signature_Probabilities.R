# Author: Marcos Diaz-Gay
# Date: Dec 3, 2024
# RStudio

library(tidyverse)
library(cowplot)
library(scales)
library(ggstats)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')
all_met = met
met = met %>% filter(Status == 'MSS')

d_mss = read.delim('../../../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/genes_with_driver_muts_MSS.tsv')
mss_driver_genes = colnames(d_mss)[-1]

driver_genes = mss_driver_genes

df = met %>%
    select(donor_id, Status, all_of(driver_genes)) %>%
    filter(Status %in% c('MSS')) %>%
    pivot_longer(cols = all_of(driver_genes),
                 names_to = 'Gene',
                 values_to = 'Mut_number') %>%
    mutate(Mut = Mut_number>0)

df$Driver_gene_Status = 'MSS exclusive driver'

total_cases_per_status = met %>%
    filter(Status %in% c('MSS')) %>%
    group_by(Status) %>%
    summarise(Total_samples = n())

dfplot = df %>%
    group_by(Gene, Status, Driver_gene_Status) %>%
    summarise(Mut_samples = sum(Mut)) %>%
    left_join(total_cases_per_status) %>%
    mutate(Mut_samples_prev = Mut_samples / Total_samples) %>%
    mutate(Mut_samples_prev = ifelse(Status == 'MSI', Mut_samples_prev,
                                     Mut_samples_prev))


########################################################################################
# Generation of prob files

# library(data.table)
# dmuts = read.delim('../../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/driver_muts_MSS_all.tsv')
# probs = fread('../../../../A_PROJECTS/Mutographs/CRC/Attributions/CRC_Manuscript_2024FEB29/msa/CRC_Subset/802_Consensus_v2/cosmic_noCI_penalties/SigProfilerExtractor/SBS96/Probabilities_as_in_SPA/Probabilities_per_mutation_MSA_CRC_Manuscript_2024FEB29_MSS.csv',
#               sep = ',', check.names = F)
# 
# probs[is.na(probs)] = 0
# colnames(probs)[1] = 'Tumor_Sample_Barcode'
# colnames_sbs = grep('SBS', colnames(probs))
# colnames(probs)[colnames_sbs][16:19] = paste0(colnames(probs)[colnames_sbs][16:19], '_MSS')
# colnames(probs)[colnames_sbs] = paste0(colnames(probs)[colnames_sbs], '_prob')
# 
# dmuts$Chromosome = str_replace(dmuts$Chromosome, 'chr', '')
# 
# sig_drivers = merge(dmuts, probs, by = c('Tumor_Sample_Barcode', 'Chromosome', 'Start_Position'))
# 
# write.table(sig_drivers, '../data_for_figures/Prob_Drivers_Mutographs_CRC_MSS_SBS.tsv',
#             quote=F, row.names=F, sep='\t')
# rm(probs)
# gc()
# 
# # Generation of prob file (ID)
# dmuts = read.delim('../../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/driver_muts_MSS_all.tsv')
# probs = fread('../../../../A_PROJECTS/Mutographs/CRC/Attributions/CRC_Manuscript_2024FEB29/msa/CRC_Subset/802_Consensus_v2/cosmic_noCI_penalties/SigProfilerExtractor/ID83/Probabilities_as_in_SPA/Probabilities_per_mutation_MSA_CRC_Manuscript_2024FEB29_MSS_ID.csv',
#               sep = ',', check.names = F)
# 
# probs[is.na(probs)] = 0
# colnames(probs)[1] = 'Tumor_Sample_Barcode'
# colnames_sbs = grep('ID', colnames(probs))
# colnames(probs)[colnames_sbs][11] = paste0(colnames(probs)[colnames_sbs][11], '_MSS')
# colnames(probs)[colnames_sbs] = paste0(colnames(probs)[colnames_sbs], '_prob')
# 
# table(dmuts$Start_Position %in% unique(probs$Start_Position))
# 
# dmuts$Chromosome = str_replace(dmuts$Chromosome, 'chr', '')
# dmuts$Start_Position[dmuts$Variant_Type == 'DEL'] = dmuts$Start_Position[dmuts$Variant_Type == 'DEL'] - 1
# 
# sig_drivers = merge(dmuts, probs, by = c('Tumor_Sample_Barcode', 'Chromosome', 'Start_Position'))
# 
# write.table(sig_drivers, '../data_for_figures/Prob_Drivers_Mutographs_CRC_MSS_ID.tsv',
#             quote=F, row.names=F, sep='\t')
# rm(probs)
# gc()

########################################################################################

sig_drivers = tibble(read.delim('../data_for_figures/Prob_Drivers_Mutographs_CRC_MSS_SBS.tsv'))

probcols_sbs = colnames(sig_drivers)[str_detect(colnames(sig_drivers), '_prob')]

sig_drivers = sig_drivers %>%
    rename(donor_id = Tumor_Sample_Barcode) %>%
    left_join(met)

sig_drivers_id = tibble(read.delim('../data_for_figures/Prob_Drivers_Mutographs_CRC_MSS_ID.tsv'))

probcols_id = colnames(sig_drivers_id)[str_detect(colnames(sig_drivers_id), '_prob')]

sig_drivers_id = sig_drivers_id %>%
    rename(donor_id = Tumor_Sample_Barcode) %>%
    left_join(met)

sig_drivers[,probcols_id] = NA
sig_drivers_id[,probcols_sbs] = NA

sig_drivers = rbind(sig_drivers, sig_drivers_id)
probcols = c(probcols_sbs, probcols_id)


################################################################################
# COLOMBIA

dfplot = sig_drivers %>%
    filter(Variant_Type == 'SNP') %>%
    group_by(HGVSc, donor_id,country, Hugo_Symbol,age_eo) %>%
    select(all_of(probcols_sbs)) %>%
    pivot_longer(cols = all_of(probcols_sbs)) %>%
    mutate(name = sapply(strsplit(name, '_prob'),'[',1)) %>%
    group_by(HGVSc, donor_id,country, Hugo_Symbol,age_eo) %>%
    summarise(Signature = probcols_sbs[which.max(value)]) %>%
    mutate(Signature = sapply(strsplit(Signature, '_prob'),'[',1)) %>%
    mutate(category = ifelse(country == 'Colombia',
                             'Colombia', 'Other country'))
# dfplot$country = factor(dfplot$country, levels = c('Colombia', 'Other country'))
# dfplot$country[is.na(dfplot$country)] = 'Other country'

dfplot$Signature[which(dfplot$Signature=='SBS288F_MSS')] = 'SBS_F'
dfplot$Signature = factor(dfplot$Signature, levels = c('SBS94', 'SBS_F', 'SBS1', 'SBS5',
                                                       'SBS18','Others',
                                                       '<95% confidence')) 
dfplot$Signature[is.na(dfplot$Signature)] = 'Others'

table(dfplot$Signature, dfplot$category)

all = dfplot %>%
    group_by(category) %>%
    summarise(all=n())

dfplot_numbers = dfplot %>%
    group_by(category, Signature,.drop = FALSE) %>%
    summarise(prev=n()) %>%
    left_join(all) %>%
    filter(Signature %in% c('SBS94', 'SBS_F')) %>%
    mutate(prev_label = paste0(prev,'/',all),
           prev_perc = prev/all)

dfplot %>%
    left_join(dfplot_numbers) %>%
ggplot() +
    aes(y = reorder(category, desc(category)),
        fill = Signature) +
    geom_bar(position = position_fill(reverse = TRUE)) +
    geom_text(aes(x = ifelse(category == 'Colombia',
                             ifelse(Signature == 'SBS94',0.1, 0.27),
                             ifelse(Signature == 'SBS94',0.08, 0.14)),
                  y = ifelse(category == 'Colombia', 2,
                             ifelse(Signature == 'SBS94',1.2, 0.8)),
                  label = prev_label),
              check_overlap = T, size = 4.5) +
    scale_fill_manual(values = c('yellow3','red3', 'grey80', 'grey70',
                                 'grey60',  'grey50')) +
    theme_bw() +
    labs(x = 'Proportion of mutations probabilistically assigned',
         y = '', fill ='',
         title = 'Signatures assigned to individual driver mutations') +
    scale_x_continuous(labels = label_percent()) +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 16),
          axis.title = element_text(size = 14),
          axis.text.y = element_text(size = 14, hjust = 0.5),
          axis.text.x = element_text(size = 14,),
          legend.title = element_text(size = 16, face ='bold'),
          legend.text = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'left',
          strip.background = element_blank(),
          strip.text.x = element_text(size = 16, face = 'bold'),
          panel.spacing.x = unit(2,'lines')) +
    guides(fill = guide_legend(nrow=1))
# Exported 800 x 300

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# ARGENTINA

dfplot = sig_drivers %>%
    filter(Variant_Type == 'SNP') %>%
    group_by(HGVSc, donor_id,country, Hugo_Symbol,age_eo) %>%
    select(all_of(probcols_sbs)) %>%
    pivot_longer(cols = all_of(probcols_sbs)) %>%
    mutate(name = sapply(strsplit(name, '_prob'),'[',1)) %>%
    group_by(HGVSc, donor_id,country, Hugo_Symbol,age_eo) %>%
    summarise(Signature = probcols_sbs[which.max(value)]) %>%
    mutate(Signature = sapply(strsplit(Signature, '_prob'),'[',1)) %>%
    mutate(category = ifelse(country == 'Argentina',
                             'Argentina', 'Other country'))
# dfplot$country = factor(dfplot$country, levels = c('Colombia', 'Other country'))
# dfplot$country[is.na(dfplot$country)] = 'Other country'

dfplot$Signature = factor(dfplot$Signature, levels = c('SBS89', 'SBS1', 'SBS5',
                                                       'SBS18','Others',
                                                       '<95% confidence')) 
dfplot$Signature[is.na(dfplot$Signature)] = 'Others'

table(dfplot$Signature, dfplot$category)

all = dfplot %>%
    group_by(category) %>%
    summarise(all=n())

dfplot_numbers = dfplot %>%
    group_by(category, Signature,.drop = FALSE) %>%
    summarise(prev=n()) %>%
    left_join(all) %>%
    filter(Signature %in% c('SBS89')) %>%
    mutate(prev_label = paste0(prev,'/',all),
           prev_perc = prev/all)

dfplot %>%
    left_join(dfplot_numbers) %>%
    ggplot() +
    aes(y = reorder(category, desc(category)),
        fill = Signature) +
    geom_bar(position = position_fill(reverse = TRUE)) +
    geom_text(aes(x = 0.08,
                  label = prev_label),
              check_overlap = T, size = 4.5) +
    scale_fill_manual(values = c('blue4','grey90', 'grey80', 'grey70',
                                 'grey60',  'grey50')) +
    theme_bw() +
    labs(x = 'Proportion of mutations probabilistically assigned',
         y = '', fill ='',
         title = 'Signatures assigned to individual driver mutations') +
    scale_x_continuous(labels = label_percent()) +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 16),
          axis.title = element_text(size = 14),
          axis.text.y = element_text(size = 14, hjust = 0.5),
          axis.text.x = element_text(size = 14,),
          legend.title = element_text(size = 16, face ='bold'),
          legend.text = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'left',
          strip.background = element_blank(),
          strip.text.x = element_text(size = 16, face = 'bold'),
          panel.spacing.x = unit(2,'lines'))
# Exported 800 x 300


################################################################################

dfplot = sig_drivers %>%
    filter(!(Variant_Type == 'SNP')) %>%
    group_by(HGVSc, donor_id,country, Hugo_Symbol,age_eo) %>%
    select(all_of(probcols_id)) %>%
    pivot_longer(cols = all_of(probcols_id)) %>%
    mutate(name = sapply(strsplit(name, '_prob'),'[',1)) %>%
    group_by(HGVSc, donor_id,country, Hugo_Symbol,age_eo) %>%
    summarise(Signature = probcols_id[which.max(value)]) %>%
    mutate(Signature = sapply(strsplit(Signature, '_prob'),'[',1)) %>%
    mutate(category = ifelse(country == 'Argentina',
                             'Argentina', 'Other country'))
# dfplot$country = factor(dfplot$country, levels = c('Colombia', 'Other country'))
# dfplot$country[is.na(dfplot$country)] = 'Other country'
dfplot$Signature[which(dfplot$Signature=='ID83J_MSS')] = 'ID_J'

dfplot$Signature = factor(dfplot$Signature, levels = c('ID_J', 'ID1', 'ID2',
                                                       'ID14', 'Others',
                                                       '<95% confidence'))
dfplot$Signature[is.na(dfplot$Signature)] = 'Others'

table(dfplot$Signature, dfplot$category)

all = dfplot %>%
    group_by(category) %>%
    summarise(all=n())

dfplot_numbers = dfplot %>%
    group_by(category, Signature,.drop = FALSE) %>%
    summarise(prev=n()) %>%
    left_join(all) %>%
    filter(Signature %in% c('ID_J')) %>%
    mutate(prev_label = paste0(prev,'/',all),
           prev_perc = prev/all)

dfplot %>%
    left_join(dfplot_numbers) %>%
    ggplot() +
    aes(y = reorder(category, desc(category)),
        fill = Signature) +
    geom_bar(position = position_fill(reverse = TRUE)) +
    geom_text(aes(x = 0.08,
                  label = prev_label,
                  col = category),
              check_overlap = T, size = 4.5) +
    scale_fill_manual(values = c('blue4','grey90', 'grey80', 'grey70',
                                 'grey60',  'grey50')) +
    scale_color_manual(values = c('white', 'black')) +
    theme_bw() +
    labs(x = 'Proportion of mutations probabilistically assigned',
         y = '', fill ='',
         title = 'Signatures assigned to individual driver indels') +
    scale_x_continuous(labels = label_percent()) +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 16),
          axis.title = element_text(size = 14),
          axis.text.y = element_text(size = 14, hjust = 0.5),
          axis.text.x = element_text(size = 14,),
          legend.title = element_text(size = 16, face ='bold'),
          legend.text = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'left',
          legend.direction = 'horizontal',
          strip.background = element_blank(),
          strip.text.x = element_text(size = 16, face = 'bold'),
          panel.spacing.x = unit(2,'lines')) +
    guides(col = 'none')
# Exported 800 x 300
