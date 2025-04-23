# Author: Marcos Diaz-Gay
# Date: Dec 2, 2024
# RStudio

library(tidyverse)
library(data.table)

# Get list of chromatin modifier genes (cmpgs)
## Source dbEM (https://doi.org/10.1038/srep19340)
dbem = readLines('Responses_to_Reviewers/epi_seq.fa')
dbem = dbem[grep('>',dbem)]
dbem = sapply(strsplit(dbem,'_'),'[',2)
dbem = sort(dbem)
dbem

# Annotate cmgs in mutographs samples

##### Oncokb annotation files ##################################################
all_maf = NULL
path_mafs = '../../../../../../../../Restricted/Mutographs/CRC_Manuscript_2023JUL26/Somatic/SNVs/oncokb_annotated/'
files_mafs = dir(path_mafs)
files_mafs = grep("_oncokb.maf", files_mafs, value = T)
for (i in 1:length(files_mafs)){
    print(i)
    print(files_mafs[i])
    maf = fread(paste0(path_mafs, '/', files_mafs[i]))
    maf = maf %>% filter(Hugo_Symbol %in% dbem)
    all_maf = rbind(all_maf, maf)
}

path_mafs_id = '../../../../../../../../Restricted/Mutographs/CRC_Manuscript_2023JUL26/Somatic/Indels/oncokb_annotated/'
files_mafs_id = dir(path_mafs_id)
files_mafs_id = grep("_oncokb.maf", files_mafs_id, value = T)
for (i in 1:length(files_mafs_id)){
    print(i)
    print(files_mafs_id[i])
    maf = fread(paste0(path_mafs_id, '/', files_mafs_id[i]))
    maf = maf %>% filter(Hugo_Symbol %in% dbem)
    all_maf = rbind(all_maf, maf)
}


##### Metadata #################################################################

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')
status = 'MSS'
df = met %>% filter(Status == status)

all_maf_802 = all_maf %>% filter(Tumor_Sample_Barcode %in% df$donor_id)

cmgs_unfiltered = all_maf_802

cmgs_filtered = cmgs_unfiltered %>%
    filter(Variant_Classification != 'Silent' &
               Variant_Classification != 'Splice_Region')
write.table(cmgs_filtered, 'Responses_to_Reviewers/cmgs_mutations_CRC_Mutographs_MSS_cases.tsv',
            sep = '\t', quote = F, row.names = F)

cmgs = read.delim('Responses_to_Reviewers/cmgs_mutations_CRC_Mutographs_MSS_cases.tsv')

bin_long = cmgs %>%
    group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
    summarise(Number_of_cmgs_muts = n())
bin_wide = bin_long %>%
    pivot_wider(names_from = Hugo_Symbol, values_from = Number_of_cmgs_muts)
bin_wide[is.na(bin_wide)] = 0
samples_to_merge = data.frame(Tumor_Sample_Barcode = df$donor_id)
bin_wide = merge(bin_wide, samples_to_merge, all.y = T)
bin_wide[is.na(bin_wide)] = 0

write.table(bin_wide, 'Responses_to_Reviewers/genes_with_cmgs_muts_MSS.tsv',
            row.names=F, quote=F, sep = '\t')
