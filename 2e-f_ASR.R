# Author: Marcos Diaz-Gay
# Date: Jul 21, 2024
# RStudio

library(tidyverse)
library(ggrepel)
library(ggpubr)
library(scales)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')
met$country[met$country == 'Czech Republic'] = 'Czechia'

###################
status = 'MSS'
# status = 'MSI'
###################

df = met %>% filter(Status == status)

all_sigsint = grep('_c$', colnames(met), value = T)
sigsint = all_sigsint[colSums(apply(df[,all_sigsint],2,is.na))!=nrow(df)]
sigsint = sigsint[!str_detect(sigsint, 'Unknown')]

total_samples = df %>%
    pivot_longer(cols = c('SBS1_c','ID1_c','DBS2_c','CN1_c','SV1_c')) %>%
    mutate(variant_type = str_sub(name,1,2)) %>%
    group_by(country, variant_type) %>%
    summarise(n = sum(!is.na(value))) 

dfplot_aux = df %>%
    group_by(country) %>%
    summarise(across(all_of(sigsint),~sum(. > 0, na.rm = T))) %>%
    pivot_longer(cols = sigsint, names_to = 'Signature', values_to = 'Prevalence') %>%
    mutate(variant_type = str_sub(Signature,1,2)) %>%
    left_join(total_samples) %>%
    mutate(Prevalence_perc = Prevalence / n,
           Signature_label = sapply(strsplit(Signature, '_c'), '[', 1))

dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'CN68F_MSS'] = 'CN_F'
dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'SBS288F_MSS'] = 'SBS_F'
dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'SBS288H_MSS'] = 'SBS_H'
dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'SBS288M_MSS'] = 'SBS_M'
dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'SBS288O_MSS'] = 'SBS_O'
dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'ID83J_MSS'] = 'ID_J'
dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'SV38B_MSS'] = 'SV_B'
dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'SV38D_MSS'] = 'SV_D'

dfplot = df %>%
    group_by(country, ASR_CRC_both) %>%
    summarise(across(all_of(sigsint), \(x) mean(x, na.rm = T))) %>%
    pivot_longer(cols = sigsint, names_to = 'Signature', values_to = 'mean_activity') %>%
    left_join(dfplot_aux)
    
df$tumorsite_group = factor(df$tumorsite_group)
df$tumorsite_group = relevel(df$tumorsite_group, 'Proximal colon')
df$sex = factor(df$sex)
df$sex = relevel(df$sex, 'Male')

# CRC
all_p_values = NULL
all_r_squares = NULL
all_betas = NULL
for (i in 1:length(sigsint)){
    signature = sigsint[i]
    print(signature)
    
    model = lm(get(signature) ~ ASR_CRC_both + age_diag + sex + tumorsite_group + purity,
               data = df)
    p_value = summary(model)$coefficients['ASR_CRC_both',4]
    r_square = summary(model)$r.squared
    beta = summary(model)$coefficients['ASR_CRC_both',1]

    all_p_values = c(all_p_values, p_value)
    all_r_squares = c(all_r_squares, r_square)
    all_betas = c(all_betas, beta)
}

df_sign = tibble(Signature = sigsint,
                 p_value = all_p_values,
                 q_value = c(p.adjust(all_p_values[which(str_detect(sigsint, 'SBS'))],
                                      method = 'fdr'),
                             p.adjust(all_p_values[which(str_detect(sigsint, 'ID'))],
                                      method = 'fdr'),
                             p.adjust(all_p_values[which(str_detect(sigsint, 'DBS'))],
                                      method = 'fdr'),
                             p.adjust(all_p_values[which(str_detect(sigsint, 'CN'))],
                                      method = 'fdr'),
                             p.adjust(all_p_values[which(str_detect(sigsint, 'SV'))],
                                      method = 'fdr')
                 ),
                 r_square = all_r_squares,
                 beta = all_betas)

# Plots
dfplot = dfplot %>%
    left_join(df_sign)
signature_of_interest = c('SBS88', 'ID18')

dfplot %>%
    filter(Signature_label %in% signature_of_interest) %>%
    ggplot() +
    aes(x = ASR_CRC_both, y = mean_activity) +
    geom_smooth(method = lm, fill = 'skyblue') +
    geom_point(aes(size = n, col = Prevalence_perc)) +
    facet_wrap(reorder(Signature_label,Prevalence_perc)~., scales = 'free', dir = 'v') +
    geom_text_repel(aes(label = country), size = 5) +
    geom_text(aes(x = -Inf,
                  y = Inf,
                  label=paste0('β = ',
                               formatC(beta,digits = 2, format = 'f'),
                               ', q = ',
                               formatC(q_value,digits = 3, format = 'f')
                               )),
              check_overlap = TRUE, size = 5, hjust = -0.1, vjust=1.4) +
    theme_bw() +
    scale_color_viridis_c(labels = scales::label_percent(),
                          direction = -1, breaks = c(0.1,0.3,0.5),
                          option = 'A') +
    scale_y_continuous(labels = label_comma()) +
    scale_size_continuous(breaks = c(25, 75, 125),
                          range = c(2,6)) +
    labs(x = 'ASR per 100,000 (colorectal cancer)',
         y = 'Average signature activity',
         title = 'Colorectal cancer incidence association',
         subtitle = 'Adjusted by age, sex, tumor subsite, and purity') +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.position = 'top',
          legend.direction = 'horizontal',
          legend.spacing.x = unit(2, 'lines'),
          legend.justification = 'left',
          strip.background = element_blank(),
          strip.text.x = element_text(size = 14, face = 'bold')) +
    guides(size = guide_bins(barwidth = 1.7, order = 1,
                             title = 'Total\ncases'),
           color = guide_colorbar(barwidth = 6, order = 2,
                                  title = 'Signature\nprevalence'))
# Exported 550 x 750


signature_of_interest = c('SBS1', 'SBS_H','CN_F')

dfplot %>%
    filter(Signature_label %in% signature_of_interest) %>%
    ggplot() +
    aes(x = ASR_CRC_both, y = mean_activity) +
    geom_smooth(method = lm, fill = 'skyblue') +
    geom_point(aes(size = n, col = Prevalence_perc)) +
    facet_wrap(reorder(Signature_label,-Prevalence_perc)~., scales = 'free', dir = 'v') +
    geom_text_repel(aes(label = country), size = 5) +
    geom_text(aes(x = -Inf,
                  y = Inf,
                  label=paste0('β = ',
                               formatC(beta,digits = 2, format = 'f'),
                               ', q = ',
                               formatC(q_value,digits = 3, format = 'f')
                  )),
              check_overlap = TRUE, size = 5, hjust = -0.1, vjust=1.4) +
    theme_bw() +
    scale_color_viridis_c(labels = scales::label_percent(),
                          direction = -1, breaks = c(0.2,0.6,1),
                          option = 'A') +
    scale_y_continuous(labels = label_comma(),
                       expand = expansion(mult = c(0.05, 0.2))
                       ) +
    scale_size_continuous(breaks = c(25, 75, 125),
                          range = c(2,6)) +
    labs(x = 'ASR per 100,000 (colorectal cancer)',
         y = 'Average signature activity',
         title = 'Colorectal cancer incidence association',
         subtitle = 'Adjusted by age, sex, tumor subsite, and purity') +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.position = 'top',
          legend.direction = 'horizontal',
          legend.spacing.x = unit(2, 'lines'),
          legend.justification = 'left',
          strip.background = element_blank(),
          strip.text.x = element_text(size = 14, face = 'bold')) +
    guides(size = guide_bins(barwidth = 1.7, order = 1,
                             title = 'Total\ncases'),
           color = guide_colorbar(barwidth = 6, order = 2,
                                  title = 'Signature\nprevalence'))
# Exported 550 x 1200



############################################################################################################################################################
# Separating CRC / Colon / Rectum

all_asrs = met %>% 
    group_by(country, ASR_CRC_both, ASR_colon_both, ASR_rectum_both) %>%
    summarise(ASR_CRC_both = unique(ASR_CRC_both)) %>%
    pivot_longer(cols = all_of(c('ASR_colon_both',
                                 'ASR_rectum_both')),
                 names_to = 'asr_type', values_to = 'asr_value') %>%
    mutate(asr_type_label = asr_type)
all_asrs$tumorsite_group_v2 = NA
all_asrs$tumorsite_group_v2[all_asrs$asr_type == 'ASR_colon_both'] = 'Colon'
all_asrs$tumorsite_group_v2[all_asrs$asr_type == 'ASR_rectum_both'] = 'Rectum'

all_asrs$asr_type_label[all_asrs$asr_type=='ASR_colon_both'] = 'Colon cancer'
all_asrs$asr_type_label[all_asrs$asr_type=='ASR_rectum_both'] = 'Rectal cancer'


# Colon
all_p_values = NULL
all_r_squares = NULL
all_betas = NULL
for (i in 1:length(sigsint)){
    signature = sigsint[i]
    print(signature)
    
    dfmodel = df %>% filter(tumorsite_group_v2 == 'Colon')
    
    model = lm(get(signature) ~ ASR_colon_both + age_diag + sex + purity,
               data = dfmodel)
    p_value = summary(model)$coefficients['ASR_colon_both',4]
    r_square = summary(model)$r.squared
    beta = summary(model)$coefficients['ASR_colon_both',1]
    
    all_p_values = c(all_p_values, p_value)
    all_r_squares = c(all_r_squares, r_square)
    all_betas = c(all_betas, beta)
}

df_sign_2 = tibble(Signature = sigsint,
                 p_value = all_p_values,
                 q_value = c(p.adjust(all_p_values[which(str_detect(sigsint, 'SBS'))],
                                      method = 'fdr'),
                             p.adjust(all_p_values[which(str_detect(sigsint, 'ID'))],
                                      method = 'fdr'),
                             p.adjust(all_p_values[which(str_detect(sigsint, 'DBS'))],
                                      method = 'fdr'),
                             p.adjust(all_p_values[which(str_detect(sigsint, 'CN'))],
                                      method = 'fdr'),
                             p.adjust(all_p_values[which(str_detect(sigsint, 'SV'))],
                                      method = 'fdr')
                 ),
                 r_square = all_r_squares,
                 beta = all_betas,
                 asr_type = 'ASR_colon_both')

# Rectum
all_p_values = NULL
all_r_squares = NULL
all_betas = NULL
for (i in 1:length(sigsint)){
    signature = sigsint[i]
    print(signature)
    
    dfmodel = df %>% filter(tumorsite_group_v2 == 'Rectum')
    
    model = lm(get(signature) ~ ASR_rectum_both + age_diag + sex + purity,
               data = dfmodel)
    p_value = summary(model)$coefficients['ASR_rectum_both',4]
    r_square = summary(model)$r.squared
    beta = summary(model)$coefficients['ASR_rectum_both',1]
    
    all_p_values = c(all_p_values, p_value)
    all_r_squares = c(all_r_squares, r_square)
    all_betas = c(all_betas, beta)
}

df_sign_3 = tibble(Signature = sigsint,
                   p_value = all_p_values,
                   q_value = c(p.adjust(all_p_values[which(str_detect(sigsint, 'SBS'))],
                                        method = 'fdr'),
                               p.adjust(all_p_values[which(str_detect(sigsint, 'ID'))],
                                        method = 'fdr'),
                               p.adjust(all_p_values[which(str_detect(sigsint, 'DBS'))],
                                        method = 'fdr'),
                               p.adjust(all_p_values[which(str_detect(sigsint, 'CN'))],
                                        method = 'fdr'),
                               p.adjust(all_p_values[which(str_detect(sigsint, 'SV'))],
                                        method = 'fdr')
                   ),
                   r_square = all_r_squares,
                   beta = all_betas,
                   asr_type = 'ASR_rectum_both')

df_sign = rbind(df_sign_2, df_sign_3)


total_samples = df %>%
    pivot_longer(cols = c('SBS1_c','ID1_c','DBS2_c','CN1_c','SV1_c')) %>%
    mutate(variant_type = str_sub(name,1,2)) %>%
    group_by(country, tumorsite_group_v2, variant_type) %>%
    summarise(n = sum(!is.na(value))) 

dfplot_aux = df %>%
    group_by(country, tumorsite_group_v2) %>%
    summarise(across(all_of(sigsint),~sum(. > 0, na.rm = T))) %>%
    pivot_longer(cols = sigsint, names_to = 'Signature', values_to = 'Prevalence') %>%
    mutate(variant_type = str_sub(Signature,1,2)) %>%
    left_join(total_samples) %>%
    mutate(Prevalence_perc = Prevalence / n,
           Signature_label = sapply(strsplit(Signature, '_c'), '[', 1))

dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'CN68F_MSS'] = 'CN_F'
dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'SBS288F_MSS'] = 'SBS_F'
dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'SBS288H_MSS'] = 'SBS_H'
dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'SBS288M_MSS'] = 'SBS_M'
dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'SBS288O_MSS'] = 'SBS_O'
dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'ID83J_MSS'] = 'ID_J'
dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'SV38B_MSS'] = 'SV_B'
dfplot_aux$Signature_label[dfplot_aux$Signature_label == 'SV38D_MSS'] = 'SV_D'

dfplot = df %>%
    group_by(country, tumorsite_group_v2) %>%
    summarise(across(all_of(sigsint), \(x) mean(x, na.rm = T))) %>%
    pivot_longer(cols = sigsint, names_to = 'Signature', values_to = 'mean_activity') %>%
    left_join(dfplot_aux)

dfplot = dfplot %>%
    left_join(all_asrs) %>%
    left_join(df_sign)

signature_of_interest = c('SBS88', 'ID18')

dfplot %>%
    filter(Signature_label %in% signature_of_interest) %>%
ggplot() +
    aes(x = asr_value, y = mean_activity) +
    geom_smooth(method = lm, fill = 'skyblue') +
    geom_point(aes(size = n, col = Prevalence_perc)) +
    facet_wrap(reorder(Signature_label, Prevalence_perc)~asr_type_label, scales = 'free') +
    # facet_grid(reorder(Signature_label, Prevalence_perc)~asr_type_label, scales = 'free') +
    geom_text_repel(aes(label = country), size = 5) +
    geom_text(aes(x = -Inf,
                  y = Inf,
                  label=paste0('β = ',
                               formatC(beta,digits = 2, format = 'f'),
                               ', q = ',
                               formatC(q_value,digits = 3, format = 'f')
                  )),
              check_overlap = TRUE, size = 5, hjust = -0.1, vjust=1.4) +
    theme_bw() +
    scale_color_viridis_c(labels = scales::label_percent(),
                          direction = -1,
                          option = 'A') +
    scale_y_continuous(labels = label_comma()) +
    scale_size_continuous(limits = c(0,80)) +
    labs(x = 'ASR per 100,000',
         y = 'Average signature activity',
         title = 'Colon and rectal cancer incidence association',
         subtitle = 'Adjusted by age, sex, and purity') +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.position = 'top',
          legend.direction = 'horizontal',
          legend.spacing.x = unit(2, 'lines'),
          legend.justification = 'left',
          strip.background = element_blank(),
          strip.text.x = element_text(size = 14, face = 'bold'),
          strip.text.y = element_text(size = 14, face = 'bold',
                                      angle = 0)) +
    guides(size = guide_bins(barwidth = 1.7, order = 1,
                             title = 'Total cases\t\t'),
           color = guide_colorbar(barwidth = 9, order = 2,
                                  title = 'Signature prevalence\t\t\t',
                                  title.hjust = -2))
# Exported 1000 x 750

signature_of_interest = c('SBS1', 'SBS_H','CN_F')

dfplot %>%
    # filter(asr_type == 'ASR_CRC_both') %>%
    filter(Signature_label %in% signature_of_interest) %>%
    ggplot() +
    aes(x = asr_value, y = mean_activity) +
    geom_smooth(method = lm, fill = 'skyblue') +
    geom_point(aes(size = n, col = Prevalence_perc)) +
    facet_wrap(reorder(Signature_label, -Prevalence_perc)~asr_type_label, scales = 'free',
               ncol = 2) +
    # facet_grid(reorder(Signature_label, -Prevalence_perc)~asr_type_label, scales = 'free') +
    geom_text_repel(aes(label = country), size = 5) +
    geom_text(aes(x = -Inf,
                  y = Inf,
                  label=paste0('β = ',
                               formatC(beta,digits = 2, format = 'f'),
                               ', q = ',
                               formatC(q_value,digits = 3, format = 'f')
                  )),
              check_overlap = TRUE, size = 5, hjust = -0.1, vjust=1.4) +
    theme_bw() +
    scale_color_viridis_c(labels = scales::label_percent(),
                          direction = -1,
                          option = 'A') +
    scale_y_continuous(labels = label_comma(),
                       expand = expansion(mult = c(0.05, 0.2))
    ) +
    scale_size_continuous(limits = c(0,80)) +
    labs(x = 'ASR per 100,000',
         y = 'Average signature activity',
         title = 'Colon and rectal cancer incidence association',
         subtitle = 'Adjusted by age, sex, and purity') +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.position = 'top',
          legend.direction = 'horizontal',
          legend.spacing.x = unit(2, 'lines'),
          legend.justification = 'left',
          strip.background = element_blank(),
          strip.text.x = element_text(size = 14, face = 'bold'),
          strip.text.y = element_text(size = 14, face = 'bold',
                                      angle = 0)) +
    guides(size = guide_bins(barwidth = 1.7, order = 1,
                             title = 'Total cases\t\t'),
           color = guide_colorbar(barwidth = 9, order = 2,
                                  title = 'Signature prevalence\t\t\t',
                                  title.hjust = -2))
# Exported 1000 x 1200


# Supplementary table
df_sign$Signature = str_replace(df_sign$Signature, '_c','')
df_sign$Signature[df_sign$Signature=='SBS288F_MSS'] = 'SBS_F'
df_sign$Signature[df_sign$Signature=='SBS288H_MSS'] = 'SBS_H'
df_sign$Signature[df_sign$Signature=='SBS288M_MSS'] = 'SBS_M'
df_sign$Signature[df_sign$Signature=='SBS288O_MSS'] = 'SBS_O'
df_sign$Signature[df_sign$Signature=='ID83J_MSS'] = 'ID_J'
df_sign$Signature[df_sign$Signature=='CN68F_MSS'] = 'CN_F'
df_sign$Signature[df_sign$Signature=='SV38B_MSS'] = 'SV_B'
df_sign$Signature[df_sign$Signature=='SV38D_MSS'] = 'SV_D'

signature_order = gtools::mixedsort(unique(df_sign$Signature))
signature_order = signature_order[c(which(str_detect(signature_order, 'SBS')),
                                    which(str_detect(signature_order, 'DBS')),
                                    which(str_detect(signature_order, 'ID')),
                                    which(str_detect(signature_order, 'CN')),
                                    which(str_detect(signature_order, 'SV')))]

testdata_to_st = df_sign %>%
    mutate(Signature = factor(Signature,levels = signature_order),
           Tumor_subsite = ifelse(asr_type=='ASR_colon_both', 'Colon',
                                  'Rectum')) %>%
    arrange(Signature,Tumor_subsite
    ) %>%
    select(Signature, Tumor_subsite, beta, p_value, q_value) %>%
    rename('Tumor subsite' = Tumor_subsite,
        'p-value' = p_value,
        'q-value' = q_value,
        'beta value' = beta)
write.csv(testdata_to_st, '../supplementary_tables_regressions/ST18.csv', quote=F, row.names=F)
