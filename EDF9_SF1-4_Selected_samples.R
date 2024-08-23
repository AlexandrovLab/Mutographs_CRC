# Author: Marcos Diaz-Gay
# Date: Jul 30, 2024
# RStudio

library(tidyverse)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')

# SBS288
mm = read.delim('Input_Data/CRC_Manuscript_v1.SBS288.all')

samples_to_filter_edf9 = c('PD54091a', 'PD54114a',
                      'PD54005a', 'PD57104a')

samples_dna_rep = met %>%
    filter(! Status %in% c('MSS', 'MSI')) %>%
    pull(donor_id)

samples_to_filter = c(samples_to_filter_edf9, samples_dna_rep)

mm_new = mm[, samples_to_filter]
mm_new = cbind(MutationType=mm[,1], mm_new)

write.table(mm_new, '../../../../A_PROJECTS/Mutographs/CRC/Signatures/CRC_Manuscript_2024FEB29/sigprofiler/matrix_generator/CRC_All/981_Consensus/plots_high_quality/CRC_Manuscript_v1.SBS288_EDF9_SF1-4.all',
            sep = '\t', quote = F, row.names=F)


mm = read.delim('Input_Data/CRC_Manuscript_v1.SBS288.all')

mm_new = cbind(MutationType=mm[,1], PD48980a=mm[, 'PD48980a'])
write.table(mm_new, '../../../../A_PROJECTS/Mutographs/CRC/Signatures/CRC_Manuscript_2024FEB29/sigprofiler/matrix_generator/CRC_All/981_Consensus/plots_high_quality/CRC_Manuscript_v1.SBS288_SF5.all',
            sep = '\t', quote = F, row.names=F)

# ID83
mm = read.delim('Input_Data/CRC_Manuscript_v1.ID83.all')
mm_new = mm[, samples_to_filter]
mm_new = cbind(MutationType=mm[,1], mm_new)
write.table(mm_new, '../../../../A_PROJECTS/Mutographs/CRC/Signatures/CRC_Manuscript_2024FEB29/sigprofiler/matrix_generator/CRC_All/981_Consensus/plots_high_quality/CRC_Manuscript_v1.ID83_EDF9_SF1-4.all',
            sep = '\t', quote = F, row.names=F)

mm = read.delim('Input_Data/CRC_Manuscript_v1.ID83.all')
mm_new = cbind(MutationType=mm[,1], PD48980a=mm[, 'PD48980a'])
write.table(mm_new, '../../../../A_PROJECTS/Mutographs/CRC/Signatures/CRC_Manuscript_2024FEB29/sigprofiler/matrix_generator/CRC_All/981_Consensus/plots_high_quality/CRC_Manuscript_v1.ID83_SF5.all',
            sep = '\t', quote = F, row.names=F)

# DBS78
mm = read.delim('Input_Data/CRC_Manuscript_v1.DBS78.all')
mm_new = mm[, samples_to_filter]
mm_new = cbind(MutationType=mm[,1], mm_new)
write.table(mm_new, '../../../../A_PROJECTS/Mutographs/CRC/Signatures/CRC_Manuscript_2024FEB29/sigprofiler/matrix_generator/CRC_All/981_Consensus/plots_high_quality/CRC_Manuscript_v1.DBS78_EDF9_SF1-4.all',
            sep = '\t', quote = F, row.names=F)

# Signature activities SBS
all_sigsint = grep('_c_rel$', colnames(met), value = T)
sigs = all_sigsint[c(1:19,31)]

dfplot = met %>%
    select(donor_id, all_of(sigs)) %>%
    filter(donor_id %in% samples_to_filter_edf9) %>%
    pivot_longer(-donor_id, names_to = 'signature') %>%
    mutate(signature = str_replace(signature, '_c_rel', '')) %>%
    mutate(signature = factor(signature, levels = c('SBS88','Other')))

dfplot$signature[is.na(dfplot$signature)] = 'Other'
dfplot = dfplot %>%
    group_by(donor_id, signature) %>%
    summarise(value = sum(value)) %>%
    mutate(label = ifelse(value == 0, '',
                          paste0(signature, ' (',
                         round(value*100,1), '%)')))
dfplot %>%
    # filter(donor_id == 'PD54091a') %>%
    # filter(donor_id == 'PD54114a') %>%
    filter(donor_id == 'PD54005a') %>%
    # filter(donor_id == 'PD57104a') %>%
ggplot() +
    aes(x = donor_id, y = value, fill = signature, label = label) +
    geom_col() +
    geom_text(position = position_fill(vjust = 0.5),
              size = 4.5, fontface = 'bold') +
    theme_void() +
    scale_fill_manual(values = c("skyblue2","grey70")) +
    guides(fill = 'none')
# Exported 200 x 200


# Signature activities ID
all_sigsint = grep('_c_rel$', colnames(met), value = T)
sigs = all_sigsint[c(32:43)]

dfplot = met %>%
    select(donor_id, all_of(sigs)) %>%
    filter(donor_id %in% samples_to_filter_edf9) %>%
    pivot_longer(-donor_id, names_to = 'signature') %>%
    mutate(signature = str_replace(signature, '_c_rel', '')) %>%
    mutate(signature = factor(signature, levels = c('ID18','Other')))

dfplot$signature[is.na(dfplot$signature)] = 'Other'
dfplot = dfplot %>%
    group_by(donor_id, signature) %>%
    summarise(value = sum(value)) %>%
    mutate(label = ifelse(value == 0, '',
                          paste0(signature, ' (',
                                 round(value*100,1), '%)')))
dfplot %>%
    # filter(donor_id == 'PD54091a') %>%
    # filter(donor_id == 'PD54114a') %>%
    # filter(donor_id == 'PD54005a') %>%
    filter(donor_id == 'PD57104a') %>%
    ggplot() +
    aes(x = donor_id, y = value, fill = signature, label = label) +
    geom_col() +
    geom_text(position = position_fill(vjust = 0.5),
              size = 4.5, fontface = 'bold') +
    theme_void() +
    scale_fill_manual(values = c("skyblue2","grey70")) +
    guides(fill = 'none')
# Exported 200 x 200
