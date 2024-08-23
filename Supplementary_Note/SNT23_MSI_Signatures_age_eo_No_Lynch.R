# Author: Marcos Diaz-Gay
# Date: Aug 12, 2024
# RStudio

library(tidyverse)
library(scales)
library(ggrepel)
library(broom)
library(logistf)

###################
# status = 'MSS'
status = 'MSI'
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
met = met %>% filter(tumorsite_group != 'Other')

df = met %>% filter(Status == status)
df = df %>% filter(Lynch == 'No')

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


# # Only for age (using OR for 10 years)
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
    testdata$signature[testdata$signature=='SBS288I_MSI'] = 'SBS_I_MSI'
    testdata$signature[testdata$signature=='SBS288M_MSI'] = 'SBS_M_MSI'
    testdata$signature[testdata$signature=='SBS288N_MSI'] = 'SBS_N_MSI'
    testdata$signature[testdata$signature=='SBS288O_MSI'] = 'SBS_O_MSI'
    testdata$signature[testdata$signature=='DBS78B_MSI'] = 'DBS_B_MSI'
    testdata$signature[testdata$signature=='DBS78F_MSI'] = 'DBS_F_MSI'
    testdata$signature[testdata$signature=='SV38D_MSS'] = 'SV_D'
    
}

all_sigsint = str_replace(all_sigsint, '_c','')
all_sigsint[all_sigsint=='SBS288I_MSI'] = 'SBS_I_MSI'
all_sigsint[all_sigsint=='SBS288M_MSI'] = 'SBS_M_MSI'
all_sigsint[all_sigsint=='SBS288N_MSI'] = 'SBS_N_MSI'
all_sigsint[all_sigsint=='SBS288O_MSI'] = 'SBS_O_MSI'
all_sigsint[all_sigsint=='DBS78B_MSI'] = 'DBS_B_MSI'
all_sigsint[all_sigsint=='DBS78F_MSI'] = 'DBS_F_MSI'
all_sigsint[all_sigsint=='SV38D_MSS'] = 'SV_D'


testdata_to_st = testdata %>%
    mutate(signature = factor(signature, levels = c(
        all_sigsint))) %>%
    arrange(signature) %>%
    select(signature, OR, p_val, q_val_2, model_used) %>%
    rename(#Country = country,
           Signature = signature,
           'p-value' = p_val,
           'q-value' = q_val_2,
           'Model' = model_used)
write.csv(testdata_to_st, '../supplementary_tables_regressions/SNT23.csv', quote=F, row.names=F)
