# Author: Marcos Diaz-Gay
# Date: Jul 30, 2024
# RStudio

library(ggplot2)
library(dplyr)
library(ggrepel)

# devtools::install_github('marcos-diazg/utils.mdg')
library(utils.mdg)

# Average mutational profiles
met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')
samples_mss = met$donor_id[met$Status == 'MSS']
samples_msi = met$donor_id[met$Status == 'MSI']

contexts = c('SBS', 'ID', 'DBS')#, 'CN', 'SV', 'SV','CN')
subcontexts = c('SBS288', 'ID83', 'DBS78')#, 'CN48', 'SV32', 'SV38','CN68')

for (i in 1:length(contexts)){
    context = contexts[i]
    subcontext = subcontexts[i]
    
    sbs288 = read.delim(paste0('Input_Data/CRC_Manuscript_v1.', subcontext,'.all'),
                        row.names = 1)
    
    sbs288_mss = sbs288[,colnames(sbs288) %in% samples_mss]
    sbs288_msi = sbs288[,colnames(sbs288) %in% samples_msi]
    
    sbs288_mss = sapply(sbs288_mss,prop.table)
    sbs288_msi = sapply(sbs288_msi,prop.table)
    
    
    sbs288_mss = data.frame(rowMeans(sbs288_mss))
    colnames(sbs288_mss)[1] = paste0('Average MSS (n=802)')
    rownames(sbs288_mss) = rownames(sbs288)
    
    sbs288_msi = data.frame(rowMeans(sbs288_msi))
    colnames(sbs288_msi)[1] = paste0('Average MSI (n=153)')
    rownames(sbs288_msi) = rownames(sbs288)

    sbs288_mss_to_plot = cbind(rownames(sbs288_mss), sbs288_mss)
    colnames(sbs288_mss_to_plot)[1] = 'MutationType'
    dir.create('../data_for_figures/Avg_mut_profiles')
    write.table(sbs288_mss_to_plot, paste0('../data_for_figures/Avg_mut_profiles/Avg_MSS_',
                                           subcontext, '.all'),
                quote = F, row.names = F, sep = '\t')
    
    sbs288_msi_to_plot = cbind(rownames(sbs288_msi), sbs288_msi)
    colnames(sbs288_msi_to_plot)[1] = 'MutationType'
    write.table(sbs288_msi_to_plot, paste0('../data_for_figures/Avg_mut_profiles/Avg_MSI_',
                                           subcontext, '.all'),
                quote = F, row.names = F, sep = '\t')
}



##################################################
# Comparisons (Countries / Location / Age)

# comparisons = c("country", "tumorsite_group", "age_eo",'colibactin_signature_SBS88orID18','sex')
comparisons = c("country", "age_eo")

for (a in 1:length(comparisons)){
    comparison = comparisons[a]
    print(comparison)
    
    # groups = c('MSS', 'MSI')
    groups = c('MSS')
    countries = unique(met[order(met$ASR_CRC_both, decreasing = T), comparison])
    
    cs_contexts_groups_countries = NULL
    for (i in 1:length(contexts)){
        context = contexts[i]
        subcontext = subcontexts[i]
        print(subcontext)
        
        sbs288 = read.delim(paste0('Input_Data/CRC_Manuscript_v1.', subcontext,'.all'),
                            row.names = 1)
        cs_groups_countries = NULL
        profiles_countries = NULL
        for (j in 1:length(groups)){
            group = groups[j]
            print(group)
            
            numbers_per_country = met %>%
                filter(Status == group) %>%
                group_by(.data[[comparison]]) %>%
                summarise(all_samples = n()) %>%
                mutate(label = paste0(.data[[comparison]],
                                      ' (n=', all_samples,
                                      ')'))
            if (comparison == 'age_eo'){
                numbers_per_country$label = c(
                    'early-onset (0-49; n=97)',
                    'late-onset (50+; n=705)'
                )
            }
            
            cs_countries = NULL
            for (k in 1:length(countries)){
                country = countries[k]
                print(country)
                samples_to_filter = met$donor_id[met[,comparison] == country &
                                                     met$Status == group]
                
                avg_profile = sbs288[,colnames(sbs288) %in% samples_to_filter]
                avg_profile = sapply(avg_profile,prop.table)
                avg_profile = data.frame(rowMeans(avg_profile))
                colnames(avg_profile)[1] = numbers_per_country$label[numbers_per_country[,1] == country]
                rownames(avg_profile) = rownames(sbs288)
            
                if (j ==1 & k == 1){
                    profiles_countries = avg_profile
                } else {
                    profiles_countries = cbind(profiles_countries, avg_profile)
                }
            }
        }
        profiles_countries_to_plot = cbind(rownames(profiles_countries), profiles_countries)
        colnames(profiles_countries_to_plot)[1] = 'MutationType'
        colnames(profiles_countries_to_plot)[-1] = paste0('MSS ',
                                                          colnames(profiles_countries_to_plot)[-1])
        write.table(profiles_countries_to_plot,
                    paste0('../data_for_figures/Avg_mut_profiles/Matrix_',
                           subcontext, '_', comparison, '.all'),
                    quote = F, row.names = F, sep = '\t')
        
        cs_pc = cos_sim_matrix(profiles_countries,
                       profiles_countries)
        print(min(cs_pc))
    }
}







