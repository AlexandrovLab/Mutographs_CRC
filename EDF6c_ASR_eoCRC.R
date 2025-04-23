# Author: Marcos Diaz-Gay
# Date: Dec 27, 2024
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
variable = 'ASR_group'

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

df = met %>% filter(Status == status) %>% filter(age_eo=='0-49')
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

df$ASR_group = ifelse(df$ASR_CRC_both_EO_20_49<7,
                      'Low ASR', 'High ASR')
df$ASR_group = factor(df$ASR_group, levels =
                          c('Low ASR', 'High ASR'))

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


covdata = df %>% select(donor_id, ASR_group, age_diag, sex, tumorsite_group, purity)


number_sigs_by_vartype = tibble(name = sigsint) %>%
    mutate(vartype = substr(name,1,2)) %>%
    group_by(vartype) %>%
    summarise(n_sigs_vartype = n())
        

testdata <- dfmodel %>%
    pivot_longer(-donor_id) %>%
    left_join(covdata) %>% 
    group_by(name) %>%
    mutate(value=as.factor(value)) %>% 
    do(tresult = safely(stats::glm)(value ~ ASR_group + age_diag + sex + tumorsite_group + purity,family = binomial(),data=.)) %>% 
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
    do(tresult = safely(logistf::logistf)(value ~ ASR_group + age_diag + sex + tumorsite_group + purity,
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
                                df$ASR_group)<=0)
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
                                df$ASR_group)<=0)
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


########################################################################
########################################################################
# Plot
variable = 'ASR_group'
signatures_of_interest = c('SBS88_c', 'ID18_c')
                           
all_samples_by_group = df %>% 
    group_by(.data[[variable]]) %>%
    summarise(all_samples = n())

dfplot = NULL
for (i in 1:length(signatures_of_interest)){
    signature_of_interest = signatures_of_interest[i]
    print(signature_of_interest)
    
    dfplot_ind = df %>%
        filter(Status == status) %>%
        group_by(.data[[variable]]) %>%
        summarise(n_sig = sum(.data[[signature_of_interest]]>0)) %>%
        left_join(all_samples_by_group) %>%
        mutate(prev_sig = n_sig/all_samples,
               label_prev = paste0(n_sig,'/',all_samples),
               signature = str_replace(signature_of_interest,'_c',''))
    dfplot = rbind(dfplot, dfplot_ind)
}

dfplot2 = df %>%
    filter(Status == status) %>%
    group_by(.data[[variable]]) %>%
    summarise(n_sig = sum(.data[['colibactin_signature_SBS88orID18']]=='Yes')) %>%
    left_join(all_samples_by_group) %>%
    mutate(prev_sig = n_sig/all_samples,
           label_prev = paste0(n_sig,'/',all_samples),
           signature = str_replace(signature_of_interest,'_c','')) %>%
    mutate(signature = 'SBS88 or ID18')
dfplot = rbind(dfplot, dfplot2)

dfplot$signature = factor(dfplot$signature, levels = c('SBS88', 'ID18', 'SBS88 or ID18'))

q_values_age_group = testdata %>% filter(signature %in% c('SBS88', 'ID18')) %>%
    pull(p_val)

# model for joint colibactin prevalence
model_joint = summary(glm(factor(colibactin_signature_SBS88orID18) ~ ASR_group + age_diag + sex + tumorsite_group + purity,
                  data = df, family = binomial))
q_model_joint = model_joint$coefficients['ASR_groupHigh ASR','Pr(>|z|)']
#
q_values_age_group = c(q_values_age_group, q_model_joint)

q_values_age_group = paste0(expression('p = '),
                            ifelse(q_values_age_group < 0.001,
                                   formatC(q_values_age_group, format = 'e', digits = 1),
                                   formatC(q_values_age_group, format = 'f', digits = 3))
)

q_values_age_group = c(q_values_age_group[2], rep('',1),
                       q_values_age_group[1], rep('',1),
                       q_values_age_group[3], rep('',1))

####

                       
dfplot %>%
    arrange(signature) %>%
    ggplot() +
    aes(x = .data[[variable]], y = prev_sig, fill = .data[[variable]]) +
    facet_wrap(.~signature,
                 scales = 'free') +
    geom_col() +
    geom_text(aes(label = label_prev), 
              size = 4.5, vjust = -0.5,check_overlap = T) +
    geom_text(aes(y = +Inf, x = -Inf,
                  label = q_values_age_group),
              size = 4.5, vjust = 1.5, hjust = -0.5) +
    scale_y_continuous(
        labels = scales::label_percent(),
        expand = c(0,0,0.15,0)) +
    scale_fill_manual(values = c("#FFB37D", "#14061A")) +
    theme_bw() +
    labs(subtitle = 'Adjusted by age, sex, tumor subsite, and purity',
         x = 'ASR group',
         y = 'Signature prevalence',
         title = 'Early-onset colorectal cancer incidence association'
    ) +
    theme(
        plot.title = element_text(size = 16, face = 'bold'),
        plot.subtitle = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = 'bold')) +
    guides(fill = 'none')
# Exported 1600 x 500
