# Author: Marcos Diaz-Gay
# Date: Jul 21, 2024
# RStudio

library(tidyverse)
library(scales)
library(ggrepel)
library(broom)
library(ggh4x)

###################
status = 'MSS'
# status = 'MSI'
variable = 'ASR_CRC_both'
###################

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')
met$country[met$country == 'Czech Republic'] = 'Czechia'

df = met %>% filter(Status == status)
all_sigsint = grep('_c$', colnames(met), value = T)
sigsint = all_sigsint[colSums(apply(df[,all_sigsint],2,is.na))!=nrow(df)]
sigsint = sigsint[!str_detect(sigsint, 'Unknown')]

df$sex = factor(df$sex)
df$sex = relevel(df$sex, 'Male')
df$tumorsite_group = factor(df$tumorsite_group)
df$tumorsite_group = relevel(df$tumorsite_group, 'Proximal colon')
df$country = factor(df$country, levels = c('Brazil', 'Iran',
                                           'Colombia', 'Thailand',
                                           'Argentina', 'Russia',
                                           'Canada', 'Poland',
                                           'Czechia', 'Serbia',
                                           'Japan'))

df$ASR_group = ifelse(df$country %in%
                              c('Iran', 'Thailand', 'Colombia', 'Brazil'),
                          'Medium ASR', 'High ASR')
df$ASR_group = factor(df$ASR_group, levels = 
                              c('Medium ASR', 'High ASR'))

dfmodel = as_tibble(df[,c('donor_id',sigsint)])

covdata = df %>% select(donor_id, age_diag, sex, tumorsite_group, purity, country,
                        variable)

number_sigs_by_vartype = tibble(name = sigsint) %>%
    mutate(vartype = substr(name,1,2)) %>%
    group_by(vartype) %>%
    summarise(n_sigs_vartype = n())

testdata <- dfmodel %>%
    pivot_longer(-donor_id) %>%
    left_join(covdata) %>% 
    group_by(name) %>%
    # mutate(value=as.factor(value)) %>% 
    do(tresult = safely(stats::lm)(value ~ ASR_CRC_both + age_diag + sex + tumorsite_group + purity,
                                    data=.)) %>% 
    mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
    filter(!tresult_null) %>%
    mutate(fit = list(tidy(tresult[['result']],exponentiate=FALSE))) %>%
    select(name,fit) %>% 
    unnest(cols = c(fit)) %>% 
    ungroup() %>%
    mutate(vartype = substr(name,1,2)) %>%
    left_join(number_sigs_by_vartype) %>%
    filter(str_detect(term,variable)) %>%
    mutate(label = paste0(
                          str_replace(name,'_c','')
                          ),
           signature = str_replace(name,'_c',''),
           # country = country_to_test
           ) %>%
    rename(p_val = p.value,
           independent_vars = term,
           OR = estimate
    ) %>%
    select(independent_vars, signature, OR, p_val, vartype, label)


# Volcano plot LR
testdata = testdata %>%
    arrange(signature) %>%
    mutate(q_val = p.adjust(p_val, method='fdr')) %>%
    mutate(q_val_2 = c(p.adjust(p_val[which(str_detect(signature, 'CN'))],
                                method = 'fdr'),
                       p.adjust(p_val[which(str_detect(signature, 'DBS'))],
                                method = 'fdr'),
                       p.adjust(p_val[which(str_detect(signature, 'ID'))],
                                method = 'fdr'),
                       p.adjust(p_val[which(str_detect(signature, 'SBS'))],
                                method = 'fdr'),
                       p.adjust(p_val[which(str_detect(signature, 'SV'))],
                                method = 'fdr')
    )
    )

complete_testdata = testdata

# testdata$OR[testdata$OR>1e5] = Inf
# testdata$OR[testdata$OR<1e-5] = 0

# testdata$q_val[testdata$q_val<1e-6] = 1e-6
# testdata$q_val_2[testdata$q_val_2<1e-6] = 1e-6

if (status == 'MSS'){
    testdata$signature[testdata$signature=='SBS288F_MSS'] = 'SBS_F'
    testdata$signature[testdata$signature=='SBS288H_MSS'] = 'SBS_H'
    testdata$signature[testdata$signature=='SBS288M_MSS'] = 'SBS_M'
    testdata$signature[testdata$signature=='SBS288O_MSS'] = 'SBS_O'
    testdata$signature[testdata$signature=='ID83J_MSS'] = 'ID_J'
    testdata$signature[testdata$signature=='CN68F_MSS'] = 'CN_F'
    testdata$signature[testdata$signature=='SV38B_MSS'] = 'SV_B'
    testdata$signature[testdata$signature=='SV38D_MSS'] = 'SV_D'
} else {
    testdata$signature[testdata$signature=='SBS288I_MSI'] = 'SBS_MSI_I'
    testdata$signature[testdata$signature=='SBS288M_MSI'] = 'SBS_MSI_M'
    testdata$signature[testdata$signature=='SBS288N_MSI'] = 'SBS_MSI_N'
    testdata$signature[testdata$signature=='SBS288O_MSI'] = 'SBS_MSI_O'
    testdata$signature[testdata$signature=='DBS78B_MSI'] = 'DBS_MSI_B'
    testdata$signature[testdata$signature=='DBS78F_MSI'] = 'DBS_MSI_F'
}

testdata %>%
    mutate(
        col_or = ifelse(OR>1, 'Enriched in higher incidence countries',
                        'Enriched in lower incidence countries'),
        col_or = ifelse(q_val_2<0.05, col_or, NA)) %>%
    mutate(vartype = str_replace(vartype, 'SB', 'SBS'),
           vartype = str_replace(vartype, 'DB', 'DBS'),
           vartype = str_replace(vartype, 'ID', 'ID'),
           vartype = str_replace(vartype, 'CN', 'CN'),
           vartype = str_replace(vartype, 'SV', 'SV'),
           vartype = factor(vartype, levels = c('SBS',
                                                'DBS',
                                                'ID',
                                                'CN',
                                                'SV'))) %>%
    ggplot() +
    aes(x = OR, y = -log10(q_val_2), col = col_or) +
    facet_grid(.~vartype, scales = 'free') +
    geom_point(size = 4) +
    geom_hline(yintercept = -log10(0.05),
               linetype = 'dashed', col = 'orange') +
    geom_hline(yintercept = -log10(0.01),
               linetype = 'dashed', col = 'red') +
    geom_vline(xintercept = 0, col = 'grey') +
    theme_bw() +
    scale_x_continuous() +
    scale_y_continuous(breaks = c(0,2,4,6),
                       labels = c('0','2','4','>6')) +
    scale_color_manual(breaks = c('Enriched in higher incidence countries',
                                  'Enriched in lower incidence countries'),
                       values = c('#2D1256','#F4BD8B'),
                       na.value = 'grey80') +
    geom_text_repel(data = . %>% 
                        mutate(label = ifelse(q_val_2 < 0.05, signature, '')),
                    aes(label = label), size = 5, seed = 2, direction = 'y') +
    labs(#x = 'Mutations per 1 ASR per 100,000 (β value)',
         x = 'Number of mutations increasing with an increase of 1 ASR per 100,000 (β value)') +
    labs(y = 'q-value (-log10)',
         col = 'Colorectal cancer incidence association',
         subtitle = '  Adjusted by age, sex, tumor subsite, and purity') +
    theme(plot.title = element_text(size = 18),
          # plot.subtitle = element_text(size = 16),
          plot.subtitle = element_text(size = 14, vjust = -16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 16, face ='bold'),
          legend.text = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'left',
          legend.box.spacing = unit(25,'pt'),
          # legend.box.background = element_rect(linewidth = 1),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 14, face = 'bold'),
          panel.spacing.x = unit(1.5,'lines')
    ) +
    guides(col = guide_legend(override.aes = list(label = ''),))
# Exported 1600 x 500

# Supplementary table
testdata = complete_testdata
testdata$signature[testdata$signature=='SBS288F_MSS'] = 'SBS_F'
testdata$signature[testdata$signature=='SBS288H_MSS'] = 'SBS_H'
testdata$signature[testdata$signature=='SBS288M_MSS'] = 'SBS_M'
testdata$signature[testdata$signature=='SBS288O_MSS'] = 'SBS_O'
testdata$signature[testdata$signature=='ID83J_MSS'] = 'ID_J'
testdata$signature[testdata$signature=='CN68F_MSS'] = 'CN_F'
testdata$signature[testdata$signature=='SV38B_MSS'] = 'SV_B'
testdata$signature[testdata$signature=='SV38D_MSS'] = 'SV_D'

signature_order = gtools::mixedsort(unique(testdata$signature))
signature_order = signature_order[c(which(str_detect(signature_order, 'SBS')),
                                    which(str_detect(signature_order, 'DBS')),
                                    which(str_detect(signature_order, 'ID')),
                                    which(str_detect(signature_order, 'CN')),
                                    which(str_detect(signature_order, 'SV')))]
testdata_to_st = testdata %>%
    mutate(#country = factor(country, levels = c(country_order)),
           signature = factor(signature,levels = signature_order)) %>%
    arrange(signature,#country
            ) %>%
    select(signature, #country,
           OR, p_val, q_val_2) %>%
    rename(#Country = country,
           Signature = signature,
           'p-value' = p_val,
           'q-value' = q_val_2,
           'beta value' = OR)
write.csv(testdata_to_st, '../supplementary_tables_regressions/ST17.csv', quote=F, row.names=F)

