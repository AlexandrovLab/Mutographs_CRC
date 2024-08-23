# Author: Marcos Diaz-Gay
# Date: Aug 12, 2024
# RStudio

library(tidyverse)
library(scales)
library(ggrepel)
library(broom)
library(logistf)

###################
status = 'MSS'
# status = 'MSI'
variable = 'country'
###################

if (variable == 'country') {
    title = 'Country enrichment'
    subtitle = 'Adjusted by age, sex, tumor subsite, and purity'
}

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv',
                 check.names = F)
met$country[met$country=='Czech Republic'] = 'Czechia'

df = met %>% filter(Status == status)

d_mss = read.delim(paste0('../../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/genes_with_driver_muts_',
                   status, '.tsv'))
# threshold_driver_genes = 0.1 * nrow(df)
threshold_driver_genes = 0
driver_genes_to_consider = names(d_mss[,-1])[colSums((d_mss[,-1]>0)) > threshold_driver_genes]
mss_driver_genes = driver_genes_to_consider


d_rec_mss = read.delim(paste0('../../../../A_PROJECTS/Mutographs/CRC/Drivers/Driver_mutations/genes_with_recurrent_driver_muts_',
                              status,'.tsv'),
                       check.names = F)
# threshold_rec_muts = 19
threshold_rec_muts = 10
rec_muts_to_consider = names(colSums(d_rec_mss[,-1])[colSums(d_rec_mss[,-1]) > threshold_rec_muts])
mss_rec_driver_muts = rec_muts_to_consider

mss_cn_dg = c(grep('loh', colnames(df),value = T),
              grep('homdel',colnames(df),value = T),
              grep('amp',colnames(df),value = T))



sigsint = c(mss_driver_genes, mss_rec_driver_muts, mss_cn_dg)



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
df$age_eo = factor(df$age_eo)
df$age_eo = relevel(df$age_eo, '0-49')

# Dichotomizing cases over 70% prevalence by above/below median
dfmodel = as_tibble(df[,c('donor_id',sigsint)])
threshold = 1
prevalence = (df %>%
                  summarise(across(all_of(sigsint),~sum(. > 0))) / nrow(df))

prevalence[,mss_cn_dg] = df %>%
    summarise(across(all_of(mss_cn_dg),~sum(. > 0, na.rm = T))) / sum(!is.na(df$APC_loh))

sigs_high_prev = names(prevalence)[prevalence > threshold]
sigs_low_prev = names(prevalence)[prevalence <= threshold]

dfmodel_low_prev = dfmodel[,c('donor_id',sigs_low_prev)]
dfmodel_low_prev[,sigs_low_prev] = dfmodel_low_prev[,sigs_low_prev]>0

dfmodel_high_prev = dfmodel[,c('donor_id',sigs_high_prev)]
dfmodel_high_prev[,sigs_high_prev] = dfmodel_high_prev[,sigs_high_prev] %>%
    mutate(across(all_of(sigs_high_prev),~(. > median(., na.rm=T))))

dfmodel = cbind(dfmodel_high_prev, dfmodel_low_prev[,-1])


covdata = df %>% select(donor_id, age_diag, sex, tumorsite_group, purity, country)


number_sigs_by_vartype = tibble(name = sigsint) %>%
    mutate(vartype = c(rep('DR', length(mss_driver_genes)),
                       rep('RDM', length(mss_rec_driver_muts)),
                       rep('CND', length(mss_cn_dg)))) %>%
    group_by(vartype) %>%
    summarise(n_sigs_vartype = n())

countries = unique(df$country)
testdata_all_countries = NULL
for (i in 1:length(countries)){
    country_to_test = countries[i]
    print(country_to_test)
    
    testdata <- dfmodel %>%
        pivot_longer(-donor_id) %>%
        left_join(covdata) %>% 
        group_by(name) %>%
        mutate(value=as.factor(value)) %>% 
        do(tresult = safely(stats::glm)(value ~ age_diag + sex + tumorsite_group + purity + (country==country_to_test),family = binomial(),data=.)) %>% 
        mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
        filter(!tresult_null)
    if (variable == 'age_diag'){
        testdata = testdata %>%
            mutate(fit = list(tidy(tresult[['result']],exponentiate=FALSE)))
    } else {
        testdata = testdata %>%
            mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE)))
    }
    testdata = testdata %>%
        select(name,fit) %>% 
        unnest(cols = c(fit)) %>% 
        ungroup() %>%
        mutate(vartype = ifelse(str_detect(name,'_c'),'RDM',
                                ifelse(str_detect(name,'_'), 'CND','DR'))) %>%
        left_join(number_sigs_by_vartype) %>%
        filter(str_detect(term,variable)) %>%
        mutate(label = paste0(country_to_test, ' ',
                              str_replace(name,' ','')),
               signature = str_replace(name,' ',''),
               country = country_to_test) %>%
        rename(p_val = p.value,
               independent_vars = term,
               OR = estimate
        ) %>%
        select(independent_vars, signature, OR, p_val, vartype, label, country)
    
    testdata_all_countries = rbind(testdata_all_countries, testdata)
}
testdata_GLM = testdata_all_countries


# New testdata using the Firth method
testdata_all_countries = NULL
for (i in 1:length(countries)){
    country_to_test = countries[i]
    print(country_to_test)
    
    testdata <- dfmodel %>%
        pivot_longer(-donor_id) %>%
        left_join(covdata) %>% 
        group_by(name) %>%
        mutate(value=as.factor(value)) %>% 
        do(tresult = safely(logistf::logistf)(value ~ age_diag + sex + tumorsite_group + purity + (country==country_to_test),
                                              data=., control =  logistf.control(maxit = 1000),
                                              plcontrol = logistpl.control(maxit = 1000))) %>% 
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
        mutate(vartype = ifelse(str_detect(name,'_c'),'RDM',
                                ifelse(str_detect(name,'_'), 'CND','DR'))) %>%
        left_join(number_sigs_by_vartype) %>%
        filter(str_detect(term,variable)) %>%
        mutate(label = paste0(country_to_test, ' ',
                              str_replace(name,' ','')),
               signature = str_replace(name,' ',''),
               country = country_to_test) %>%
        rename(p_val = p.value,
               independent_vars = term,
               OR = estimate
        ) %>%
        select(independent_vars, signature, OR, p_val, vartype, label, country)
    
    testdata_all_countries = rbind(testdata_all_countries, testdata)
}
testdata_FIRTH = testdata_all_countries


countries_ordered = unique(testdata_GLM$country)
signatures_ordered = unique(testdata_GLM$signature)

complete_separation = NULL
for (i in 1:length(countries_ordered)){
    country = countries_ordered[i]
    print(country)
    for (j in 1:length(signatures_ordered)){
        signature = signatures_ordered[j]
        print(signature)
        
        if (paste0(signature) %in% c(sigs_high_prev)){
            is_complete_separation = as.logical(
                length(
                    which(table(df[,paste0(signature)]>median(df[,paste0(signature,'_c')], na.rm = T),
                                df$country == country)<=0)
                ) +
                    length(
                        which(table(df[,paste0(signature)]>median(df[,paste0(signature,'_c')], na.rm = T),
                                    df$sex)<=0)
                    ) +
                    length(
                        which(table(df[,paste0(signature)]>median(df[,paste0(signature,'_c')], na.rm = T),
                                    df$tumorsite_group)<=0)
                    )
            )
        } else {
            is_complete_separation = as.logical(
                length(
                    which(table(df[,paste0(signature)]>0,
                                df$country == country)<=0)
                ) +
                    length(
                        which(table(df[,paste0(signature)]>0,
                                    df$sex)<=0)
                    ) +
                    length(
                        which(table(df[,paste0(signature)]>0,
                                    df$tumorsite_group)<=0)
                    )
            )
        }
        
        complete_separation = rbind(complete_separation,
                                    c(as.character(country), signature, is_complete_separation))
    }
}

complete_separation = tibble(data.frame(complete_separation))
colnames(complete_separation) = c('country', 'signature', 'complete_separation')

testdata = NULL
model_used = NULL
for (i in 1:nrow(testdata_GLM)){
    country_glm = testdata_GLM$country[i]
    country_firth = testdata_FIRTH$country[i]
    if (country_firth == country_glm) {
        country = as.character(country_glm)
    } else {
        warning('The order of testdata data frames is different!')
    }
    print(country)
    
    sig_glm = testdata_GLM$signature[i]
    sig_firth = testdata_FIRTH$signature[i]
    if (sig_firth == sig_glm) {
        sig = sig_glm
    } else {
        warning('The order of testdata data frames is different!')
    }
    print(sig)
    
    is_cs = complete_separation$complete_separation[
        which(complete_separation$country == country &
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
if (variable == 'age_diag'){
    testdata$OR = exp(testdata$OR * 10)
}

# Volcano plot LR
testdata = testdata %>%
    arrange(vartype) %>%
    mutate(q_val = p.adjust(p_val, method='fdr')) %>%
    mutate(q_val_2 = c(p.adjust(p_val[which(str_detect(vartype, 'CND'))],
                                method = 'fdr'),
                       p.adjust(p_val[which(str_detect(vartype, 'DR'))],
                                method = 'fdr'),
                       p.adjust(p_val[which(str_detect(vartype, 'RDM'))],
                                method = 'fdr')
    )
    )

# testdata$q_val[testdata$q_val<1e-6] = 1e-6
# testdata$q_val_2[testdata$q_val_2<1e-6] = 1e-6
# testdata$OR[testdata$OR>2^4] = Inf
# testdata$OR[testdata$OR<2^-4] = 0


testdata %>%
    mutate(
        col_or = ifelse(OR>1, 'Enriched in late-onset patients',
                        'Enriched in early-onset patients'),
        col_or = ifelse(q_val_2<0.05, col_or, NA)) %>%
    mutate(
        vartype = str_replace(vartype, 'DR', 'Cancer genes'),
        vartype = str_replace(vartype, 'RDM', 'Recurrent driver mutations'),
        vartype = str_replace(vartype, 'CN signaturesD', 'CNA in driver genes'),
        vartype = factor(vartype, levels = c(
            'Cancer genes',
            'Recurrent driver mutations',
            'CNA in driver genes'))) %>%
    filter(vartype %in% c('Cancer genes',
                          'Recurrent driver mutations')) %>%
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
    scale_color_manual(breaks = c('Enriched in country',
                                  'Depleted in country'),
                       values = c('#8C3FC0','#73C03F'),
                       na.value = 'grey80') +
    geom_text_repel(data = . %>% 
                        mutate(label = ifelse(q_val_2 < 0.05, label, '')),
                    aes(label = label), size = 5, seed = 2) +
    {if (variable=='age_diag') labs(x = 'Odds ratio (10 years; log2)') else labs(x = 'Odds ratio (log2)')} +
    labs(y = 'q-value (-log10)',
         col = 'Country enrichment') +
    theme(plot.title = element_text(size = 18),
          plot.subtitle = element_text(size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 16, face ='bold'),
          legend.text = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'left',
          # legend.box.background = element_rect(linewidth = 1),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 14, face = 'bold'),
          panel.spacing.x = unit(1.5,'lines')
    ) +
    guides(col = guide_legend(override.aes = list(label = ''),))

country_order = met %>%
    group_by(country) %>%
    summarise(asr = mean(ASR_CRC_both)) %>%
    arrange(desc(asr)) %>%
    pull(country)

testdata_to_st = testdata %>%
    filter(vartype %in% c('DR')) %>%
    select(signature, country, OR, p_val, q_val_2, model_used) %>%
    mutate(country = factor(country, levels = country_order)) %>%
    arrange(signature, country) %>%
    rename(Country = country,
           'Cancer driver gene' = signature,
           'p-value' = p_val,
           'q-value' = q_val_2,
           'Model' = model_used)
write.csv(testdata_to_st, '../supplementary_tables_regressions/ST24.csv', quote=F, row.names=F)

testdata_to_st = testdata %>%
    filter(vartype %in% c('RDM')) %>%
    select(signature, country, OR, p_val, q_val_2, model_used) %>%
    mutate(country = factor(country, levels = country_order)) %>%
    arrange(signature, country) %>%
    rename(Country = country,
           'Recurrent driver mutation' = signature,
           'p-value' = p_val,
           'q-value' = q_val_2,
           'Model' = model_used)
write.csv(testdata_to_st, '../supplementary_tables_regressions/ST25.csv', quote=F, row.names=F)
