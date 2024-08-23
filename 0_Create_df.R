# Creation of dataframe for figures of Diaz-Gay et al. 2024
# Author: Marcos Diaz-Gay
# Date: Jul 27, 2024
# RStudio

library(tidyverse)

################################################################################
# Metadata
met_core_data = read.csv('../../../../../../../../Restricted/Mutographs/CRC_Manuscript_2024FEB29/Somatic/metadata/MutWP1_CRC_core_data_Manuscript_v2.csv',
               check.names = F)
met_tumour_specific = read.csv('../../../../../../../../Restricted/Mutographs/CRC_Manuscript_2024FEB29/Somatic/metadata/MutWP1_CRC_tumour_specific_Manuscript_v2_2.csv',
                               check.names = F)
metadata = merge(met_core_data, met_tumour_specific, by = 'donor_id')
met_analysis_status = read.csv('../../../../../../../../Restricted/Mutographs/CRC_Manuscript_2024FEB29/Somatic/metadata/MutWP1_CRC_analysis_status_Manuscript_v2.csv',
                               check.names = F)
colnames(met_analysis_status)[1]='donor_id'
metadata = merge(metadata, met_analysis_status, by = 'donor_id')

## Recode tumorsite
metadata$tumorsite_group = NA
metadata$tumorsite_group[metadata$tumorsite == 'C18.0 - Cecum'] = 'Proximal colon'
metadata$tumorsite_group[metadata$tumorsite == 'C18.0 - Cecum; C18.2 - Right (ascending) colon'] = 'Proximal colon'
metadata$tumorsite_group[metadata$tumorsite == 'C18.2 - Right (ascending) colon'] = 'Proximal colon'
metadata$tumorsite_group[metadata$tumorsite == 'C18.2 - Right (ascending) colon; C18.4 - Transverse colon'] = 'Proximal colon'
metadata$tumorsite_group[metadata$tumorsite == 'C18.3 - Hepatic flexure'] = 'Proximal colon'
metadata$tumorsite_group[metadata$tumorsite == 'C18.4 - Transverse colon'] = 'Proximal colon'
metadata$tumorsite_group[metadata$tumorsite == 'C18.5 - Splenic flexure'] = 'Distal colon'
metadata$tumorsite_group[metadata$tumorsite == 'C18.6 - Left (descending) colon'] = 'Distal colon'
metadata$tumorsite_group[metadata$tumorsite == 'C18.6 - Left (descending) colon; C18.7 - Sigmoid'] = 'Distal colon'
metadata$tumorsite_group[metadata$tumorsite == 'C18.7 - Sigmoid'] = 'Distal colon'
metadata$tumorsite_group[metadata$tumorsite == 'C18.9 - Colon not otherwise specified'] = 'Other'
metadata$tumorsite_group[metadata$tumorsite == 'C18.7 - Sigmoid'] = 'Distal colon'
metadata$tumorsite_group[metadata$tumorsite == 'C19.9 - Rectosigmoid junction'] = 'Rectum'
metadata$tumorsite_group[metadata$tumorsite == 'C20 - Rectum'] = 'Rectum'
metadata$tumorsite_group[metadata$tumorsite == 'Other'] = 'Other'

metadata$tumorsite_group_v2 = as.character(metadata$tumorsite_group)
metadata$tumorsite_group_v2[metadata$tumorsite_group_v2 == 'Proximal colon'] = 'Colon'
metadata$tumorsite_group_v2[metadata$tumorsite_group_v2 == 'Distal colon'] = 'Colon'

metadata$tumorsite_group_v3 = as.character(metadata$tumorsite_group)
metadata$tumorsite_group_v3[metadata$tumorsite_group_v3 == 'Rectum'] = 'Left colon'
metadata$tumorsite_group_v3[metadata$tumorsite_group_v3 == 'Distal colon'] = 'Left colon'
metadata$tumorsite_group_v3[metadata$tumorsite_group_v3 == 'Proximal colon'] = 'Right colon'

## Recode age
metadata$age_group = 0
for (i in 1:nrow(metadata)){
    if (metadata$age_diag[i]<=39){
        metadata$age_group[i] = '0-39'
    }
    if (metadata$age_diag[i]<=49 & metadata$age_diag[i]>39){
        metadata$age_group[i] = '40-49'
    }
    if (metadata$age_diag[i]<=59 & metadata$age_diag[i]>49){
        metadata$age_group[i] = '50-59'
    }
    if (metadata$age_diag[i]<=69 & metadata$age_diag[i]>59){
        metadata$age_group[i] = '60-69'
    }
    if (metadata$age_diag[i]>69){
        metadata$age_group[i] = '70+'
    }
}

metadata$age_eo = 0
for (i in 1:nrow(metadata)){
    if (metadata$age_diag[i]<=49){
        metadata$age_eo[i] = '0-49'
    }
    if (metadata$age_diag[i]>49){
        metadata$age_eo[i] = '50+'
    }
}

## ASR
ASR = read.csv('Input_Data/ASR_countries_source_GLOBOCAN_2022.csv')
metadata = merge(metadata, ASR)
metadata = metadata[order(metadata$donor_id),]


################################################################################
# TMB
mm = read.delim('Input_Data/CRC_Manuscript_v1.SBS96.all')
tmb = colSums(mm[,-1])
tmb_df = data.frame(donor_id = colnames(mm[,-1]),
                    TMB_SBS = tmb)
metadata = merge(metadata, tmb_df)

mm = read.delim('Input_Data/CRC_Manuscript_v1.ID83.all')
tmb = colSums(mm[,-1])
tmb_df = data.frame(donor_id = colnames(mm[,-1]),
                    TMB_ID = tmb)
metadata = merge(metadata, tmb_df)

mm = read.delim('Input_Data/CRC_Manuscript_v1.DBS78.all')
tmb = colSums(mm[,-1])
tmb_df = data.frame(donor_id = colnames(mm[,-1]),
                    TMB_DBS = tmb)
metadata = merge(metadata, tmb_df)


################################################################################
# CN PGA
cn_stats = read.delim('Input_Data/cn_stats.txt')
cn_stats = cn_stats %>% dplyr::select(donor_id = sample, PGA = pAberrantFinal, WGD = GD,
                                      TMB_CN = nsegsRounded)
cn_stats$WGD = as.numeric(substr(cn_stats$WGD, 1, 1))
metadata = merge(metadata, cn_stats, all.x = T)


################################################################################
# SVs
sv38_mm = read.delim('Input_Data/SV-inMat-manualThresh-all.txt')
sv_burden = colSums(sv38_mm[,-1])
sv_burden = tibble(donor_id = names(sv_burden), TMB_SV = sv_burden)
metadata = merge(metadata, sv_burden, all.x = T)


################################################################################
# Signatures MSS (MSA activities)
## COSMIC
### SBS
sigs = read.xlsx('Input_Data/Mutographs_eoCRC_Manuscript_Tables_v16.xlsx',
                 sheet = 11, rows = 3:10000)

colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('SBS288', colnames(sigs))] = paste0(colnames(sigs)[grep('SBS288', colnames(sigs))],
                                                        '_MSS')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_c')

sigs_mss_and_msi = sigs

sigs = read.xlsx('Input_Data/Mutographs_eoCRC_Manuscript_Supplementary_Note_Tables_v16.xlsx',
                 sheet = 8, rows = 3:10000)
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('SBS288', colnames(sigs))] = paste0(colnames(sigs)[grep('SBS288', colnames(sigs))],
                                                        '_MSI')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_c')

sigs_mss_and_msi = merge(sigs_mss_and_msi, sigs, all = T)
metadata = merge(metadata, sigs_mss_and_msi, all.x = T)

metadata$SBS_Unknown_c = metadata$TMB_SBS - apply(metadata[,colnames(sigs_mss_and_msi)][,-1], 1, sum, na.rm = T)
sigsint = grep('_c', grep('SBS',colnames(metadata), value = T), value = T)
rel_mss_and_msi = metadata[,sigsint] / metadata$TMB_SBS
colnames(rel_mss_and_msi) = paste0(colnames(rel_mss_and_msi),'_rel')
rel_mss_and_msi = cbind(donor_id=metadata$donor_id, rel_mss_and_msi)
metadata = merge(metadata, rel_mss_and_msi, all.x = T)
metadata$SBS_Unknown_c[!metadata$Status %in% c('MSS','MSI')] = NA
metadata$SBS_Unknown_c_rel[!metadata$Status %in% c('MSS','MSI')] = NA

### ID
sigs = read.xlsx('Input_Data/Mutographs_eoCRC_Manuscript_Tables_v16.xlsx',
                 sheet = 12, rows = 3:10000)
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('ID83', colnames(sigs))] = paste0(colnames(sigs)[grep('ID83', colnames(sigs))],
                                                        '_MSS')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_c')

sigs_mss_and_msi = sigs

sigs = read.xlsx('Input_Data/Mutographs_eoCRC_Manuscript_Supplementary_Note_Tables_v16.xlsx',
                 sheet = 9, rows = 3:10000)
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('ID83', colnames(sigs))] = paste0(colnames(sigs)[grep('ID83', colnames(sigs))],
                                                      '_MSI')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_c')

sigs_mss_and_msi = merge(sigs_mss_and_msi, sigs, all = T)
metadata = merge(metadata, sigs_mss_and_msi, all.x = T)

metadata$ID_Unknown_c = metadata$TMB_ID - apply(metadata[,colnames(sigs_mss_and_msi)][,-1], 1, sum, na.rm = T)
sigsint = grep('_c', grep('ID',colnames(metadata), value = T), value = T)
rel_mss_and_msi = metadata[,sigsint] / metadata$TMB_ID
colnames(rel_mss_and_msi) = paste0(colnames(rel_mss_and_msi),'_rel')
rel_mss_and_msi = cbind(donor_id=metadata$donor_id, rel_mss_and_msi)
metadata = merge(metadata, rel_mss_and_msi, all.x = T)
metadata$ID_Unknown_c[!metadata$Status %in% c('MSS','MSI')] = NA
metadata$ID_Unknown_c_rel[!metadata$Status %in% c('MSS','MSI')] = NA

### Colibactin signatures
metadata$colibactin_signature_SBS88orID18 = ifelse(metadata$SBS88_c>0 | metadata$ID18_c>0, 'Yes', 'No')
metadata$colibactin_signature_SBS88orID18[is.na(metadata$colibactin_signature_SBS88orID18)] = 'No'

### DBS
sigs = read.xlsx('Input_Data/Mutographs_eoCRC_Manuscript_Tables_v16.xlsx',
                 sheet = 13, rows = 3:10000)
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('DBS78', colnames(sigs))] = paste0(colnames(sigs)[grep('DBS78', colnames(sigs))],
                                                      '_MSS')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_c')

sigs_mss_and_msi = sigs

sigs = read.xlsx('Input_Data/Mutographs_eoCRC_Manuscript_Supplementary_Note_Tables_v16.xlsx',
                 sheet = 10, rows = 3:10000)
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('DBS78', colnames(sigs))] = paste0(colnames(sigs)[grep('DBS78', colnames(sigs))],
                                                       '_MSI')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_c')

sigs_mss_and_msi = merge(sigs_mss_and_msi, sigs, all = T)
metadata = merge(metadata, sigs_mss_and_msi, all.x = T)

metadata$DBS_Unknown_c = metadata$TMB_DBS - apply(metadata[,colnames(sigs_mss_and_msi)][,-1], 1, sum, na.rm = T)
sigsint = grep('_c', grep('DBS',colnames(metadata), value = T), value = T)
rel_mss_and_msi = metadata[,sigsint] / metadata$TMB_DBS
colnames(rel_mss_and_msi) = paste0(colnames(rel_mss_and_msi),'_rel')
rel_mss_and_msi = cbind(donor_id=metadata$donor_id, rel_mss_and_msi)
metadata = merge(metadata, rel_mss_and_msi, all.x = T)
metadata$DBS_Unknown_c[!metadata$Status %in% c('MSS','MSI')] = NA
metadata$DBS_Unknown_c_rel[!metadata$Status %in% c('MSS','MSI')] = NA

### CN
sigs = read.xlsx('Input_Data/Mutographs_eoCRC_Manuscript_Tables_v16.xlsx',
                 sheet = 14, rows = 3:10000)
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('CNV48F',colnames(sigs))] = 'CN68F'
colnames(sigs)[grep('CN68', colnames(sigs))] = paste0(colnames(sigs)[grep('CN68', colnames(sigs))],
                                                       '_MSS')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_c')

sigs_mss_and_msi = sigs

sigs = read.xlsx('Input_Data/Mutographs_eoCRC_Manuscript_Supplementary_Note_Tables_v16.xlsx',
                 sheet = 11, rows = 3:10000)
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('CN68', colnames(sigs))] = paste0(colnames(sigs)[grep('CN68', colnames(sigs))],
                                                       '_MSI')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_c')

sigs_mss_and_msi = merge(sigs_mss_and_msi, sigs, all = T)
metadata = merge(metadata, sigs_mss_and_msi, all.x = T)

metadata$CN_Unknown_c = metadata$TMB_CN - apply(metadata[,colnames(sigs_mss_and_msi)][,-1], 1, sum, na.rm = T)
sigsint = grep('_c', grep('CN',colnames(metadata), value = T), value = T)
rel_mss_and_msi = metadata[,sigsint] / metadata$TMB_CN
colnames(rel_mss_and_msi) = paste0(colnames(rel_mss_and_msi),'_rel')
rel_mss_and_msi = cbind(donor_id=metadata$donor_id, rel_mss_and_msi)
metadata = merge(metadata, rel_mss_and_msi, all.x = T)
metadata$CN_Unknown_c[!metadata$Status %in% c('MSS', 'MSI')] = NA
metadata$CN_Unknown_c_rel[!metadata$Status %in% c('MSS', 'MSI')] = NA

### SV
sigs = read.xlsx('Input_Data/Mutographs_eoCRC_Manuscript_Tables_v16.xlsx',
                 sheet = 15, rows = 3:10000)
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('SV32B',colnames(sigs))] = 'SV38B'
colnames(sigs)[grep('SV32D',colnames(sigs))] = 'SV38D'
colnames(sigs)[grep('SV38', colnames(sigs))] = paste0(colnames(sigs)[grep('SV38', colnames(sigs))],
                                                       '_MSS')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_c')

sigs_mss_and_msi = sigs

sigs = read.xlsx('Input_Data/Mutographs_eoCRC_Manuscript_Supplementary_Note_Tables_v16.xlsx',
                 sheet = 12, rows = 3:10000)
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('SV38_MSS_D',colnames(sigs))] = 'SV38D_MSS'
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_c')

sigs_mss_and_msi = merge(sigs_mss_and_msi, sigs, all = T)
metadata = merge(metadata, sigs_mss_and_msi, all.x = T)

metadata$SV_Unknown_c = metadata$TMB_SV - apply(metadata[,colnames(sigs_mss_and_msi)][,-1], 1, sum, na.rm = T)
sigsint = grep('_c', grep('SV',colnames(metadata), value = T), value = T)
rel_mss_and_msi = metadata[,sigsint] / metadata$TMB_SV
colnames(rel_mss_and_msi) = paste0(colnames(rel_mss_and_msi),'_rel')
rel_mss_and_msi = cbind(donor_id=metadata$donor_id, rel_mss_and_msi)
metadata = merge(metadata, rel_mss_and_msi, all.x = T)
metadata$SV_Unknown_c[!metadata$Status %in% c('MSS', 'MSI')] = NA
metadata$SV_Unknown_c_rel[!metadata$Status %in% c('MSS', 'MSI')] = NA


## De Novo
### SBS
sigs = read.csv('../../../../A_PROJECTS/Mutographs/CRC/Attributions/CRC_Manuscript_2024FEB29/msa/CRC_Subset/802_Consensus_v2/denovo_noCI_penalties/SigProfilerExtractor/SBS288/output_tables/pruned_attribution_CRC_Manuscript_denovo_SBS288_abs_mutations.csv')
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('SBS288', colnames(sigs))] = paste0(colnames(sigs)[grep('SBS288', colnames(sigs))],
                                                        '_MSS')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_dn')

sigs_mss_and_msi = sigs

sigs = read.csv('../../../../A_PROJECTS/Mutographs/CRC/Attributions/CRC_Manuscript_2024FEB29/msa/CRC_MSI_Subset/153_Consensus/denovo_noCI_penalties/SBS288/output_tables/pruned_attribution_CRC_Manuscript_MSI_Subset_denovo_SBS288_abs_mutations.csv')
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('SBS288', colnames(sigs))] = paste0(colnames(sigs)[grep('SBS288', colnames(sigs))],
                                                        '_MSI')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_dn')

sigs_mss_and_msi = merge(sigs_mss_and_msi, sigs, all = T)
metadata = merge(metadata, sigs_mss_and_msi, all.x = T)

metadata$SBS_Unknown_dn = metadata$TMB_SBS - apply(metadata[,colnames(sigs_mss_and_msi)][,-1], 1, sum, na.rm = T)
sigsint = grep('_dn', grep('SBS',colnames(metadata), value = T), value = T)
rel_mss_and_msi = metadata[,sigsint] / metadata$TMB_SBS
colnames(rel_mss_and_msi) = paste0(colnames(rel_mss_and_msi),'_rel')
rel_mss_and_msi = cbind(donor_id=metadata$donor_id, rel_mss_and_msi)
metadata = merge(metadata, rel_mss_and_msi, all.x = T)
metadata$SBS_Unknown_dn[!metadata$Status %in% c('MSS','MSI')] = NA
metadata$SBS_Unknown_dn_rel[!metadata$Status %in% c('MSS','MSI')] = NA

### ID
sigs = read.csv('../../../../A_PROJECTS/Mutographs/CRC/Attributions/CRC_Manuscript_2024FEB29/msa/CRC_Subset/802_Consensus_v2/denovo_noCI_penalties/SigProfilerExtractor/ID83/output_tables/pruned_attribution_CRC_Manuscript_denovo_ID83_abs_mutations.csv')
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('ID83', colnames(sigs))] = paste0(colnames(sigs)[grep('ID83', colnames(sigs))],
                                                      '_MSS')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_dn')

sigs_mss_and_msi = sigs

sigs = read.csv('../../../../A_PROJECTS/Mutographs/CRC/Attributions/CRC_Manuscript_2024FEB29/msa/CRC_MSI_Subset/153_Consensus/denovo_noCI_penalties/ID83/output_tables/pruned_attribution_CRC_Manuscript_MSI_Subset_denovo_ID83_abs_mutations.csv')
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('ID83', colnames(sigs))] = paste0(colnames(sigs)[grep('ID83', colnames(sigs))],
                                                      '_MSI')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_dn')

sigs_mss_and_msi = merge(sigs_mss_and_msi, sigs, all = T)
metadata = merge(metadata, sigs_mss_and_msi, all.x = T)

metadata$ID_Unknown_dn = metadata$TMB_ID - apply(metadata[,colnames(sigs_mss_and_msi)][,-1], 1, sum, na.rm = T)
sigsint = grep('_dn', grep('ID',colnames(metadata), value = T), value = T)
rel_mss_and_msi = metadata[,sigsint] / metadata$TMB_ID
colnames(rel_mss_and_msi) = paste0(colnames(rel_mss_and_msi),'_rel')
rel_mss_and_msi = cbind(donor_id=metadata$donor_id, rel_mss_and_msi)
metadata = merge(metadata, rel_mss_and_msi, all.x = T)
metadata$ID_Unknown_dn[!metadata$Status %in% c('MSS','MSI')] = NA
metadata$ID_Unknown_dn_rel[!metadata$Status %in% c('MSS','MSI')] = NA

### DBS
sigs = read.csv('../../../../A_PROJECTS/Mutographs/CRC/Attributions/CRC_Manuscript_2024FEB29/msa/CRC_Subset/802_Consensus_v2/denovo_noCI_penalties/SigProfilerExtractor/DBS78/output_tables/pruned_attribution_CRC_Manuscript_denovo_DBS78_abs_mutations.csv')
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('DBS78', colnames(sigs))] = paste0(colnames(sigs)[grep('DBS78', colnames(sigs))],
                                                       '_MSS')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_dn')

sigs_mss_and_msi = sigs

sigs = read.csv('../../../../A_PROJECTS/Mutographs/CRC/Attributions/CRC_Manuscript_2024FEB29/msa/CRC_MSI_Subset/153_Consensus/denovo_noCI_penalties/DBS78/output_tables/pruned_attribution_CRC_Manuscript_MSI_Subset_denovo_DBS78_abs_mutations.csv')
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('DBS78', colnames(sigs))] = paste0(colnames(sigs)[grep('DBS78', colnames(sigs))],
                                                       '_MSI')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_dn')

sigs_mss_and_msi = merge(sigs_mss_and_msi, sigs, all = T)
metadata = merge(metadata, sigs_mss_and_msi, all.x = T)

metadata$DBS_Unknown_dn = metadata$TMB_DBS - apply(metadata[,colnames(sigs_mss_and_msi)][,-1], 1, sum, na.rm = T)
sigsint = grep('_dn', grep('DBS',colnames(metadata), value = T), value = T)
rel_mss_and_msi = metadata[,sigsint] / metadata$TMB_DBS
colnames(rel_mss_and_msi) = paste0(colnames(rel_mss_and_msi),'_rel')
rel_mss_and_msi = cbind(donor_id=metadata$donor_id, rel_mss_and_msi)
metadata = merge(metadata, rel_mss_and_msi, all.x = T)
metadata$DBS_Unknown_dn[!metadata$Status %in% c('MSS','MSI')] = NA
metadata$DBS_Unknown_dn_rel[!metadata$Status %in% c('MSS','MSI')] = NA

### CN
sigs = read.delim('../../../../A_PROJECTS/Mutographs/CRC/CN_and_SV/Chris_results/Signatures_MAY2024/CNsigs/extractor/cnsigs_input_matrix_mss/CH68/Suggested_Solution/CH68_De-Novo_Solution_cosmicfit/Assignment_Solution/Activities/Assignment_Solution_Activities.txt')
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('CH68', colnames(sigs))] = str_replace(colnames(sigs)[grep('CH68', colnames(sigs))],
                                                           'CH68', 'CN68')
colnames(sigs)[grep('CN68', colnames(sigs))] = paste0(colnames(sigs)[grep('CN68', colnames(sigs))],
                                                       '_MSS')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_dn')

sigs_mss_and_msi = sigs

sigs = read.delim('../../../../A_PROJECTS/Mutographs/CRC/CN_and_SV/Chris_results/Signatures_MAY2024/CNsigs/extractor/cnsigs_input_matrix_msi/CH68/Suggested_Solution/CH68_De-Novo_Solution_cosmicfit/Assignment_Solution/Activities/Assignment_Solution_Activities.txt')
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('CH68', colnames(sigs))] = str_replace(colnames(sigs)[grep('CH68', colnames(sigs))],
                                                           'CH68', 'CN68')
colnames(sigs)[grep('CN68', colnames(sigs))] = paste0(colnames(sigs)[grep('CN68', colnames(sigs))],
                                                       '_MSI')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_dn')

sigs_mss_and_msi = merge(sigs_mss_and_msi, sigs, all = T)
metadata = merge(metadata, sigs_mss_and_msi, all.x = T)

metadata$CN_Unknown_dn = metadata$TMB_CN - apply(metadata[,colnames(sigs_mss_and_msi)][,-1], 1, sum, na.rm = T)
sigsint = grep('_dn', grep('CN',colnames(metadata), value = T), value = T)
rel_mss_and_msi = metadata[,sigsint] / metadata$TMB_CN
colnames(rel_mss_and_msi) = paste0(colnames(rel_mss_and_msi),'_rel')
rel_mss_and_msi = cbind(donor_id=metadata$donor_id, rel_mss_and_msi)
metadata = merge(metadata, rel_mss_and_msi, all.x = T)
metadata$CN_Unknown_dn[!metadata$Status %in% c('MSS','MSI')] = NA
metadata$CN_Unknown_dn_rel[!metadata$Status %in% c('MSS','MSI')] = NA

### SV
sigs = read.delim('../../../../A_PROJECTS/Mutographs/CRC/CN_and_SV/Chris_results/Signatures_MAY2024/SVsigs/extractor/SV-inMat-manualThresh-mss/CH38/All_Solutions/CH38_6_Signatures/Activities_cosmicfit/Assignment_Solution/Activities/Assignment_Solution_Activities.txt')
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('CH38', colnames(sigs))] = str_replace(colnames(sigs)[grep('CH38', colnames(sigs))],
                                                           'CH38', 'SV38')
colnames(sigs)[grep('SV38', colnames(sigs))] = paste0(colnames(sigs)[grep('SV38', colnames(sigs))],
                                                      '_MSS')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_dn')

sigs_mss_and_msi = sigs

sigs = read.delim('../../../../A_PROJECTS/Mutographs/CRC/CN_and_SV/Chris_results/Signatures_MAY2024/SVsigs/extractor/SV-inMat-manualThresh-msi/CH38/Suggested_Solution/CH38_De-Novo_Solution_cosmicfit/Assignment_Solution/Activities/Assignment_Solution_Activities.txt')
colnames(sigs)[1] = 'donor_id'
colnames(sigs)[grep('CH38', colnames(sigs))] = str_replace(colnames(sigs)[grep('CH38', colnames(sigs))],
                                                           'CH38', 'SV38')
colnames(sigs)[grep('SV38', colnames(sigs))] = paste0(colnames(sigs)[grep('SV38', colnames(sigs))],
                                                      '_MSI')
colnames(sigs)[-1] = paste0(colnames(sigs)[-1], '_dn')

sigs_mss_and_msi = merge(sigs_mss_and_msi, sigs, all = T)
metadata = merge(metadata, sigs_mss_and_msi, all.x = T)

metadata$SV_Unknown_dn = metadata$TMB_SV - apply(metadata[,colnames(sigs_mss_and_msi)][,-1], 1, sum, na.rm = T)
sigsint = grep('_dn', grep('SV',colnames(metadata), value = T), value = T)
rel_mss_and_msi = metadata[,sigsint] / metadata$TMB_SV
colnames(rel_mss_and_msi) = paste0(colnames(rel_mss_and_msi),'_rel')
rel_mss_and_msi = cbind(donor_id=metadata$donor_id, rel_mss_and_msi)
metadata = merge(metadata, rel_mss_and_msi, all.x = T)
metadata$SV_Unknown_dn[!metadata$Status %in% c('MSS','MSI')] = NA
metadata$SV_Unknown_dn_rel[!metadata$Status %in% c('MSS','MSI')] = NA


################################################################################
# Drivers
d_mss = read.delim('../../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/genes_with_driver_muts_MSS.tsv')
colnames(d_mss)[1] = 'donor_id'

d_msi = read.delim('../../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/genes_with_driver_muts_MSI.tsv')
colnames(d_msi)[1] = 'donor_id'

d_all = merge(d_mss, d_msi, all = T)

d_all$total_driver_mutations = rowSums(d_all[,-1], na.rm = T)

metadata = merge(metadata, d_all, all.x = T)


################################################################################
# Recurrent driver mutations
d_mss_rdm = read.delim('../../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/genes_with_recurrent_driver_muts_MSS.tsv',
                   check.names = F)
colnames(d_mss_rdm)[1] = 'donor_id'

d_msi_rdm = read.delim('../../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/genes_with_recurrent_driver_muts_MSI.tsv',
                   check.names = F)
colnames(d_msi_rdm)[1] = 'donor_id'

d_all_rdm = merge(d_mss_rdm, d_msi_rdm, all = T)

metadata = merge(metadata, d_all_rdm, all.x = T)


################################################################################
# CN affecting driver genes
## MSS
path_cn_by_gene = '../../../../A_PROJECTS/Mutographs/CRC/CN_and_SV/Chris_results/CN_by_gene/'
amp = read.delim(paste0(path_cn_by_gene, '/mss-amp-matrix.txt'))
homdel = read.delim(paste0(path_cn_by_gene, '/mss-homdel-matrix.txt'))
loh = read.delim(paste0(path_cn_by_gene, '/mss-loh-matrix-cancers.txt'))

colnames(amp)[1] = 'Gene'
colnames(homdel)[1] = 'Gene'
colnames(loh)[1] = 'Gene'

amp = amp %>%
    pivot_longer(cols = all_of(starts_with('PD')),
                 names_to = 'donor_id',
                 values_to = 'CN_present') %>%
    mutate(Gene = paste0(Gene,'_amp'))

homdel = homdel %>%
    pivot_longer(cols = all_of(starts_with('PD')),
                 names_to = 'donor_id',
                 values_to = 'CN_present') %>%
    mutate(Gene = paste0(Gene,'_homdel'))
loh = loh %>%
    pivot_longer(cols = all_of(starts_with('PD')),
                 names_to = 'donor_id',
                 values_to = 'CN_present') %>%
    mutate(Gene = paste0(Gene,'_loh'))

cn = rbind(amp, homdel, loh)
cn = cn %>% pivot_wider(names_from = Gene,
                        values_from = CN_present)
driver_genes_mss = colnames(d_mss)[-1]
all_cols_to_keep = NULL
for (i in 1:length(driver_genes_mss)){
    gene = driver_genes_mss[i]
    grep_gene = grep(paste0('^',gene), colnames(cn))
    all_cols_to_keep = c(all_cols_to_keep, grep_gene)
}

all_cols_to_keep = sort(unique(all_cols_to_keep))
cn_mss = cn[,c(1,all_cols_to_keep)]

## MSI
amp = read.delim(paste0(path_cn_by_gene, '/msi-amp-matrix.txt'))
homdel = read.delim(paste0(path_cn_by_gene, '/msi-homdel-matrix.txt'))
loh = read.delim(paste0(path_cn_by_gene, '/msi-loh-matrix-cancers.txt'))

colnames(amp)[1] = 'Gene'
colnames(homdel)[1] = 'Gene'
colnames(loh)[1] = 'Gene'

amp = amp %>%
    pivot_longer(cols = all_of(starts_with('PD')),
                 names_to = 'donor_id',
                 values_to = 'CN_present') %>%
    mutate(Gene = paste0(Gene,'_amp'))

homdel = homdel %>%
    pivot_longer(cols = all_of(starts_with('PD')),
                 names_to = 'donor_id',
                 values_to = 'CN_present') %>%
    mutate(Gene = paste0(Gene,'_homdel'))
loh = loh %>%
    pivot_longer(cols = all_of(starts_with('PD')),
                 names_to = 'donor_id',
                 values_to = 'CN_present') %>%
    mutate(Gene = paste0(Gene,'_loh'))

cn = rbind(amp, homdel, loh)
cn = cn %>% pivot_wider(names_from = Gene,
                        values_from = CN_present)
driver_genes_msi = colnames(d_msi)[-1]
all_cols_to_keep = NULL
for (i in 1:length(driver_genes_msi)){
    gene = driver_genes_msi[i]
    grep_gene = grep(paste0('^',gene), colnames(cn))
    all_cols_to_keep = c(all_cols_to_keep, grep_gene)
}

all_cols_to_keep = sort(unique(all_cols_to_keep))
cn_msi = cn[,c(1,all_cols_to_keep)]

cn_mss_msi = merge(cn_mss, cn_msi, all = T)

metadata = merge(metadata, cn_mss_msi, all.x = T)


################################################################################
# Add Lynch cases
library(openxlsx)
germ = read.xlsx('Input_Data/Mutographs_eoCRC_Manuscript_Tables_v16.xlsx',
                 sheet = 2, rows = 3:1000)

colnames(germ) = paste0('germ_', colnames(germ))
colnames(germ)[2] = 'Status'
colnames(germ)[1] = 'donor_id'
germ$donor_id = paste0(str_sub(germ$donor_id,1,7),'a')

germ = germ %>%
    filter(germ_Gene %in% c('MLH1', 'MSH2', 'MSH6', 'PMS2'))

## Considering Lynch only Pathogenic / Likely_pathogenic variants
## Case reviewed manually in Clinvar (at May 30, 2024 is considered Pathogenic/Likely_pathogenic)
## https://www.ncbi.nlm.nih.gov/clinvar/variation/140847/?oq=c.2095G%3EC[varname]+PMS2&m=NM_000535.7(PMS2):c.2095G%3EC%20(p.Asp699His)
#germ[germ$`germ_Mutation.(cDNA)` == 'c.2095G>C','germ_ClinVar.Status'] = 'Pathogenic/Likely_pathogenic'

germ = germ %>%
    filter(germ_ClinVar.Status %in% c('Pathogenic', 'Likely pathogenic',
                                      'Pathogenic/Likely pathogenic')) %>%
    mutate(Lynch = 'Yes') %>%
    select(donor_id, Lynch)

metadata = merge(metadata, germ, all.x = T)
metadata$Lynch[is.na(metadata$Lynch)] = 'No'

################################################################################
# Purity
metadata$purity = as.numeric(metadata$`ASCAT %`)

## Estimate purity for CN_Normal_Pipeline
samples_to_check_purity_manually = metadata$donor_id[
    metadata$Analysis_Status=='CN Normal Pipeline']

additional_purity_based_on_vaf = read.csv('../data_for_figures/CN_Normal_Colon_Estimated_Purities.csv')

additional_purity_based_on_vaf = tibble(additional_purity_based_on_vaf) %>%
    select(PD_Number, purity_estimate) %>%
    rename(donor_id = PD_Number,
           purity = purity_estimate)

for (i in 1:nrow(metadata)){
    sample = metadata$donor_id[i]
    purity_row = metadata$purity[i]
    new_purity_row = ifelse(sample %in% samples_to_check_purity_manually,
                            additional_purity_based_on_vaf[additional_purity_based_on_vaf$donor_id == sample,
                                                           'purity'],
                            purity_row)
    metadata$purity[i] = new_purity_row
}
metadata$purity = as.numeric(metadata$purity)

################################################################################
# Print
write.table(metadata, '../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv',
            quote = F, row.names = F, sep = '\t')
