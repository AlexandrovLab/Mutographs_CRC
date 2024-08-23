# Author: Marcos Diaz-Gay
# Date: Jul 21, 2024
# RStudio

library(tidyverse)
library(cowplot)
library(scales)
library(ggstats)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')
all_met = met
met = met %>% filter(Status == 'MSS')

d_mss = read.delim('../../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/genes_with_driver_muts_MSS.tsv')
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
# Generation of prob file
# library(data.table)
# dmuts = read.delim('../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/driver_muts_MSS_all.tsv')
# probs = fread('../../../A_PROJECTS/Mutographs/CRC/Attributions/CRC_Manuscript_2024FEB29/msa/CRC_Subset/802_Consensus_v2/cosmic_noCI_penalties/SigProfilerExtractor/SBS96/Probabilities_as_in_SPA/Probabilities_per_mutation_MSA_CRC_Manuscript_2024FEB29_MSS.csv',
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
# write.table(sig_drivers, 'Prob_Drivers_Mutographs_CRC_MSS_SBS.tsv',
#             quote=F, row.names=F, sep='\t')
# rm(probs)
# gc()
# 
# # Generation of prob file (ID)
# dmuts = read.delim('../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/driver_muts_MSS_all.tsv')
# probs = fread('../../../A_PROJECTS/Mutographs/CRC/Attributions/CRC_Manuscript_2024FEB29/msa/CRC_Subset/802_Consensus_v2/cosmic_noCI_penalties/SigProfilerExtractor/ID83/Probabilities_as_in_SPA/Probabilities_per_mutation_MSA_CRC_Manuscript_2024FEB29_MSS_ID.csv',
#               sep = ',', check.names = F)
# 
# probs[is.na(probs)] = 0
# colnames(probs)[1] = 'Tumor_Sample_Barcode'
# colnames_sbs = grep('ID', colnames(probs))
# colnames(probs)[colnames_sbs][11] = paste0(colnames(probs)[colnames_sbs][11], '_MSS')
# colnames(probs)[colnames_sbs] = paste0(colnames(probs)[colnames_sbs], '_prob')
# 
# 
# table(dmuts$Start_Position %in% unique(probs$Start_Position))
# 
# dmuts$Chromosome = str_replace(dmuts$Chromosome, 'chr', '')
# dmuts$Start_Position[dmuts$Variant_Type == 'DEL'] = dmuts$Start_Position[dmuts$Variant_Type == 'DEL'] - 1
# 
# sig_drivers = merge(dmuts, probs, by = c('Tumor_Sample_Barcode', 'Chromosome', 'Start_Position'))
# 
# write.table(sig_drivers, 'Prob_Drivers_Mutographs_CRC_MSS_ID.tsv',
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
# Fig. SBS

dfplot = sig_drivers %>%
    filter(Variant_Type == 'SNP') %>%
    group_by(HGVSc, donor_id,colibactin_signature_SBS88orID18, Hugo_Symbol,age_eo) %>%
    select(all_of(probcols_sbs)) %>%
    pivot_longer(cols = all_of(probcols_sbs)) %>%
    mutate(name = sapply(strsplit(name, '_prob'),'[',1)) %>%
    group_by(HGVSc, donor_id,colibactin_signature_SBS88orID18, Hugo_Symbol,age_eo) %>%
    summarise(Signature = probcols_sbs[which.max(value)]) %>%
    mutate(Signature = sapply(strsplit(Signature, '_prob'),'[',1)) %>%
    mutate(category = ifelse(!Hugo_Symbol %in% c('TP53', 'APC'),
                             'Other driver mutations',
                             ifelse(Hugo_Symbol == 'TP53',
                                    'TP53 driver mutations',
                                    ifelse(HGVSc != 'c.835-8A>G',
                                           'Other APC driver mutations',
                                    'APC c.835-8A>G'))),
           category = factor(category, levels = c('APC c.835-8A>G',
                                                  'Other APC driver mutations',
                                                  'TP53 driver mutations',
                                                  'Other driver mutations')))
dfplot$colibactin_signature_SBS88orID18 = factor(dfplot$colibactin_signature_SBS88orID18)
levels(dfplot$colibactin_signature_SBS88orID18) = c('Colibactin negative',
                                                    'Colibactin positive')

dfplot$Signature = factor(dfplot$Signature, levels = c('SBS88','SBS1', 'SBS5',
                                                       'SBS18','Others',
                                                       '<95% confidence'
                                                                   )) 
dfplot$Signature[is.na(dfplot$Signature)] = 'Others'

table(dfplot$Signature, dfplot$category, dfplot$colibactin_signature_SBS88orID18)

all = dfplot %>%
    filter(colibactin_signature_SBS88orID18=='Colibactin positive') %>%
    group_by(category,age_eo) %>%
    summarise(all=n())

dfplot_numbers = dfplot %>%
    filter(colibactin_signature_SBS88orID18=='Colibactin positive') %>%
    group_by(category, age_eo,Signature,.drop = FALSE) %>%
    summarise(prev=n()) %>%
    left_join(all) %>%
    filter(Signature %in% c('SBS88', 'ID18')) %>%
    mutate(prev_label = paste0(prev,'/',all),
           prev_perc = prev/all)

dfplot %>%
    filter(colibactin_signature_SBS88orID18=='Colibactin positive') %>%
    left_join(dfplot_numbers) %>%
ggplot() +
    aes(y = reorder(category, desc(category)),
        fill = Signature) +
    facet_grid(.~age_eo) +
    geom_bar(position = position_fill(reverse = TRUE)) +
    geom_text(aes(x = ifelse(category == 'APC c.835-8A>G', prev_perc/2,
                             ifelse(category == 'TP53 driver mutations',
                             0.07,0.12)),
                  label = ifelse(category %in% c('APC c.835-8A>G','Other APC driver mutations'),
                                 prev_label,
                                 ifelse(category =='TP53 driver mutations',
                                        ifelse(age_eo == '0-49', '0/31','0/79'),
                                        ifelse(age_eo == '0-49', '0/57', prev_label)))),
              check_overlap = T, size = 4.5) +
    scale_fill_manual(values = c('skyblue2','grey90', 'grey80', 'grey70',
                                 'grey60',  'grey50')) +
                                 # 'grey10', 'grey20', 'grey30',
                                 # 'grey40', 'grey50', 'grey60', 'grey70',
                                 # 'grey80', 'grey90','grey95')) +
    theme_bw() +
    labs(x = 'Proportion of mutations probabilistically assigned',
         y = '', fill ='',
         title = 'Signatures assigned to individual single base substitutions in colibactin positive cases') +
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
          # legend.box.background = element_rect(linewidth = 1),
          strip.background = element_blank(),
          strip.text = element_text(size = 16, face = 'bold',angle=0),
          panel.spacing.x = unit(2,'lines'))
# Exported 1600 x 500





# Stats SBS88
test =
    fisher.test(cbind(c(129-5, 5), c(14-9, 9)))
test

test =
    fisher.test(cbind(c(110-0, 0), c(14-9, 9)))
test

test =
    fisher.test(cbind(c(315-12, 12), c(14-9, 9)))
test


################################################################################
# Fig. ID
dfplot = sig_drivers %>%
    filter(Variant_Type != 'SNP') %>%
    group_by(HGVSc, donor_id,colibactin_signature_SBS88orID18, Hugo_Symbol,age_eo) %>%
    select(all_of(probcols_id)) %>%
    pivot_longer(cols = all_of(probcols_id)) %>%
    mutate(name = sapply(strsplit(name, '_prob'),'[',1)) %>%
    group_by(HGVSc, donor_id,colibactin_signature_SBS88orID18, Hugo_Symbol,age_eo) %>%
    summarise(Signature = probcols_id[which.max(value)]) %>%
    mutate(Signature = sapply(strsplit(Signature, '_prob'),'[',1)) %>%
    mutate(category = ifelse(!Hugo_Symbol %in% c('TP53', 'APC'),
                             '        Other driver indels        ',
                             ifelse(Hugo_Symbol == 'TP53',
                                    '        TP53 driver indels        ',
                                    '        APC driver indels        ')),
           category = factor(category, levels = c('        APC driver indels        ',
                                                  '        TP53 driver indels        ',
                                                  '        Other driver indels        ')))
dfplot$colibactin_signature_SBS88orID18 = factor(dfplot$colibactin_signature_SBS88orID18)
levels(dfplot$colibactin_signature_SBS88orID18) = c('Colibactin negative',
                                                    'Colibactin positive')

dfplot$Signature = factor(dfplot$Signature, levels = c('ID18','ID1', 'ID2',
                                                       'ID14','Others', '<95% confidence'))
dfplot$Signature[is.na(dfplot$Signature)] = 'Others'

table(dfplot$Signature, dfplot$category, dfplot$colibactin_signature_SBS88orID18)
table(dfplot$colibactin_signature_SBS88orID18, dfplot$category)


all = dfplot %>%
    filter(colibactin_signature_SBS88orID18=='Colibactin positive') %>%
    group_by(category,age_eo) %>%
    summarise(all=n())

dfplot_numbers = dfplot %>%
    filter(colibactin_signature_SBS88orID18=='Colibactin positive') %>%
    group_by(category,age_eo, Signature,.drop = FALSE) %>%
    summarise(prev=n()) %>%
    left_join(all) %>%
    filter(Signature %in% c('SBS88', 'ID18')) %>%
    mutate(prev_label = paste0(prev,'/',all),
           prev_perc = prev/all)


dfplot %>%
    filter(colibactin_signature_SBS88orID18=='Colibactin positive') %>%
    left_join(dfplot_numbers) %>%
ggplot() +
    aes(y = reorder(category, desc(category)),
        fill = Signature) +
    facet_grid(.~age_eo) +
    geom_bar(position = position_fill(reverse = TRUE)) +
    geom_text(aes(x = ifelse((category == '        TP53 driver indels        ' &
                                  age_eo == '0-49'),0.05,prev_perc/2),

                  label = ifelse((category == '        TP53 driver indels        ' &
                                      age_eo == '0-49'),'0/4',prev_label)),
              check_overlap = T, size = 4.5) +
    scale_fill_manual(values = c('skyblue2','grey90', 'grey80', 'grey70',
                                 'grey60',  'grey50')) +
    # 'grey10', 'grey20', 'grey30',
    # 'grey40', 'grey50', 'grey60', 'grey70',
    # 'grey80', 'grey90','grey95')) +
    theme_bw() +
    labs(x = 'Proportion of mutations probabilistically assigned',
         y = '', fill ='',
         title = 'Signatures assigned to individual indels in colibactin positive cases') +
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
          # legend.box.background = element_rect(linewidth = 1),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 16, face = 'bold'),
          panel.spacing.x = unit(2,'lines'))
# Exported 1600 x 500





# Stats ID18


test =
    fisher.test(cbind(c(15-2,2), c(83-21, 21)))
test

test =
    fisher.test(cbind(c(77-13,13), c(83-21, 21)))
test
