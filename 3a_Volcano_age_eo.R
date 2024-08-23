# Author: Marcos Diaz-Gay
# Date: Aug 11, 2024
# RStudio

library(tidyverse)
library(scales)
library(ggrepel)
library(broom)
library(logistf)

###################
status = 'MSS'
# status = 'MSI'
# variable = 'age_diag'
variable = 'age_eo'

###################

if (variable == 'age_diag') {
    title = 'Age of onset enrichment'
    subtitle = 'Adjusted by sex, country, and tumor subsite'
}
if (variable == 'age_eo') {
    title = 'Age of onset enrichment'
    subtitle = 'Adjusted by sex, country, and tumor subsite'
}

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')

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
                                           'Czech Republic', 'Serbia',
                                           'Japan'))
df$age_eo = factor(df$age_eo)
df$age_eo = relevel(df$age_eo, '0-49')

# Dichotomizing cases over 70% prevalence by above/below median
dfmodel = as_tibble(df[,c('donor_id',sigsint)])
threshold = 0.7
prevalence = (df %>%
    summarise(across(all_of(sigsint),~sum(. > 0))) / nrow(df))

sv_sigs = grep('SV', sigsint, value = T)
cn_sigs = grep('CN', sigsint, value = T)

prevalence[,cn_sigs] = df %>%
    summarise(across(all_of(cn_sigs),~sum(. > 0, na.rm = T))) / sum(!is.na(df$CN1_c))

prevalence[,sv_sigs] = df %>%
    summarise(across(all_of(sv_sigs),~sum(. > 0, na.rm = T))) / sum(!is.na(df$SV3_c))

sigs_high_prev = names(prevalence)[prevalence > threshold]
sigs_low_prev = names(prevalence)[prevalence <= threshold]

dfmodel_low_prev = dfmodel[,c('donor_id',sigs_low_prev)]
dfmodel_low_prev[,sigs_low_prev] = dfmodel_low_prev[,sigs_low_prev]>0

dfmodel_high_prev = dfmodel[,c('donor_id',sigs_high_prev)]
dfmodel_high_prev[,sigs_high_prev] = dfmodel_high_prev[,sigs_high_prev] %>%
    mutate(across(all_of(sigs_high_prev),~(. > median(., na.rm=T))))

dfmodel = cbind(dfmodel_high_prev, dfmodel_low_prev[,-1])


covdata = df %>% select(donor_id, age_eo, sex, tumorsite_group, country, purity)


number_sigs_by_vartype = tibble(name = sigsint) %>%
    mutate(vartype = substr(name,1,2)) %>%
    group_by(vartype) %>%
    summarise(n_sigs_vartype = n())
        

testdata <- dfmodel %>%
    pivot_longer(-donor_id) %>%
    left_join(covdata) %>% 
    group_by(name) %>%
    mutate(value=as.factor(value)) %>% 
    do(tresult = safely(stats::glm)(value ~ age_eo + sex + country + tumorsite_group + purity,family = binomial(),data=.)) %>% 
    mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
    filter(!tresult_null)
testdata = testdata %>%
        mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE)))
testdata = testdata %>%
    select(name,fit) %>% 
    unnest(cols = c(fit)) %>% 
    ungroup() %>%
    mutate(vartype = substr(name,1,2)) %>%
    left_join(number_sigs_by_vartype) %>%
    filter(str_detect(term,variable)) %>%
    mutate(label = paste0(str_replace(term,variable,''), ' ',
    # mutate(label = paste0(country_to_test, ' ',
                          str_replace(name,'_c','')),
           signature = str_replace(name,'_c',''),
           # country = country_to_test
    ) %>%
    rename(p_val = p.value,
           independent_vars = term,
           OR = estimate
           ) %>%
    # select(independent_vars, signature, OR, p_val, vartype, label, country)
    select(independent_vars, signature, OR, p_val, vartype, label)

testdata_GLM = testdata

# New testdata using the Firth method
testdata <- dfmodel %>%
    pivot_longer(-donor_id) %>%
    left_join(covdata) %>% 
    group_by(name) %>%
    mutate(value=as.factor(value)) %>% 
    do(tresult = safely(logistf::logistf)(value ~ age_eo + sex + country + tumorsite_group + purity,
                                          data=., control =  logistf.control(maxit = 1000))) %>% 
    mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
    filter(!tresult_null)

testdata$pre_fit = lapply(testdata$tresult,'[[','result')
list_to_testdata_fit = list()
for (j in 1:length(testdata$pre_fit)){
    list_to_testdata_fit[[j]] = tibble(
        term = testdata$pre_fit[[j]]$terms,
        estimate = exp(testdata$pre_fit[[j]]$coefficients),
        p.value = testdata$pre_fit[[j]]$prob
    )
}
testdata$fit = list_to_testdata_fit

testdata = testdata %>%
    select(name,fit) %>% 
    unnest(cols = c(fit)) %>% 
    ungroup() %>%
    mutate(vartype = substr(name,1,2)) %>%
    left_join(number_sigs_by_vartype) %>%
    filter(str_detect(term,variable)) %>%
    mutate(label = paste0(str_replace(term,variable,''), ' ',
                          # mutate(label = paste0(country_to_test, ' ',
                          str_replace(name,'_c','')),
           signature = str_replace(name,'_c',''),
           # country = country_to_test
    ) %>%
    rename(p_val = p.value,
           independent_vars = term,
           OR = estimate
    ) %>%
    # select(independent_vars, signature, OR, p_val, vartype, label, country)
    select(independent_vars, signature, OR, p_val, vartype, label)

testdata_FIRTH = testdata


# countries_ordered = unique(testdata_GLM$country)
signatures_ordered = unique(testdata_GLM$signature)

complete_separation = NULL
# for (i in 1:length(countries_ordered)){
#     country = countries_ordered[i]
#     print(country)
    for (j in 1:length(signatures_ordered)){
        signature = signatures_ordered[j]
        print(signature)
        
        if (paste0(signature,'_c') %in% c(sigs_high_prev)){
            is_complete_separation = as.logical(
                length(
                    which(table(df[,paste0(signature,'_c')]>median(df[,paste0(signature,'_c')], na.rm = T),
                                df$age_eo)<=0)
                ) +                length(
                    which(table(df[,paste0(signature,'_c')]>median(df[,paste0(signature,'_c')], na.rm = T),
                                df$country)<=0)
                ) +
                    length(
                        which(table(df[,paste0(signature,'_c')]>median(df[,paste0(signature,'_c')], na.rm = T),
                                    df$sex)<=0)
                    ) +
                    length(
                        which(table(df[,paste0(signature,'_c')]>median(df[,paste0(signature,'_c')], na.rm = T),
                                    df$tumorsite_group)<=0)
                    )
            )
        } else {
            is_complete_separation = as.logical(
                length(
                    which(table(df[,paste0(signature,'_c')]>0,
                                df$age_eo)<=0)
                ) +
                length(
                    which(table(df[,paste0(signature,'_c')]>0,
                                df$country)<=0)
                ) +
                    length(
                        which(table(df[,paste0(signature,'_c')]>0,
                                    df$sex)<=0)
                    ) +
                    length(
                        which(table(df[,paste0(signature,'_c')]>0,
                                    df$tumorsite_group)<=0)
                    )
            )
        }
        complete_separation = rbind(complete_separation,
                                    c(#as.character(country),
                                      signature, is_complete_separation))
    }
# }

complete_separation = tibble(data.frame(complete_separation))
colnames(complete_separation) = c(#'country',
                                  'signature', 'complete_separation')

testdata = NULL
model_used = NULL
for (i in 1:nrow(testdata_GLM)){
    # country_glm = testdata_GLM$country[i]
    # country_firth = testdata_FIRTH$country[i]
    # if (country_firth == country_glm) {
    #     country = as.character(country_glm)
    # } else {
    #     warning('The order of testdata data frames is different!')
    # }
    # print(country)
    
    sig_glm = testdata_GLM$signature[i]
    sig_firth = testdata_FIRTH$signature[i]
    if (sig_firth == sig_glm) {
        sig = sig_glm
    } else {
        warning('The order of testdata data frames is different!')
    }
    print(sig)
    
    is_cs = complete_separation$complete_separation[
        which(#complete_separation$country == country &
                  complete_separation$signature == sig)
    ]
    if (is_cs) {
        testdata = rbind(testdata, testdata_FIRTH[i,])
        model_used = c(model_used, "Firth's Bias-Reduced Logistic Regression")
    } else {
        testdata = rbind(testdata, testdata_GLM[i,])
        model_used = c(model_used, "Regular Logistic Regression")
    }
}
testdata$model_used = model_used




# Only for age (using OR for 10 years)
# if (variable == 'age_diag'){
#     testdata$OR = exp(testdata$OR * 10)
# }
# if (variable == 'age_eo'){
#     testdata$OR = exp(testdata$OR)
# }


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

testdata$OR[testdata$OR>1e5] = Inf
testdata$OR[testdata$OR<1e-5] = 0

testdata$q_val[testdata$q_val<1e-6] = 1e-6
testdata$q_val_2[testdata$q_val_2<1e-6] = 1e-6

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
           col_or = ifelse(OR<1, 'Enriched in early-onset patients',
                           'Enriched in late-onset patients'),
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
    aes(x = log2(OR), y = -log10(q_val_2), col = col_or) +
    facet_grid(.~vartype, ) +
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
    scale_color_manual(breaks = c('Enriched in early-onset patients',
                                  'Enriched in late-onset patients'
                                  ),
                       values = c('#8C3FC0','#73C03F'),
                       na.value = 'grey80') +
    geom_text_repel(data = . %>% 
                        mutate(signature = ifelse(q_val_2 < 0.05, signature, '')),
                    aes(label = signature), size = 5, seed = 2, nudge_y = 0.4) +
    {if (variable=='age_diag') labs(x = 'Odds ratio (10 years; log2)') else labs(x = 'Odds ratio (log2)')} +
    labs(y = 'q-value (-log10)',
         col = 'Age of onset enrichment',
         subtitle = '  Adjusted by sex, country, tumor subsite, and purity') +
    theme(plot.title = element_text(size = 18),
          plot.subtitle = element_text(size = 14, vjust = -16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 16, face ='bold'),
          legend.text = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'left',
          legend.box.spacing = unit(25,'pt'),
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
           OR, p_val, q_val_2, model_used) %>%
    rename(#Country = country,
           Signature = signature,
           'p-value' = p_val,
           'q-value' = q_val_2,
           'Model' = model_used)
write.csv(testdata_to_st, '../supplementary_tables_regressions/ST19.csv', quote=F, row.names=F)
