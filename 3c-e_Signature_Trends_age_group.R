# Author: Marcos Diaz-Gay
# Date: Aug 11, 2024
# RStudio

library(tidyverse)
library(scales)
library(ggbeeswarm)
library(broom) # For tidy
library(ggh4x) # For facet_manual
library(logistf)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')
met$country[met$country == 'Czech Republic'] = 'Czechia'

status = 'MSS'
variable = 'age_group'
signatures_of_interest = c('SBS88_c', 'ID18_c', 'SBS288M_MSS_c',
                           'ID14_c', 'SBS89_c')
# signatures_of_interest = c('SBS88_c','ID18_c')

df = met %>%
    filter(Status == status)

df$country = factor(df$country, levels = c('Brazil', 'Iran',
                                           'Colombia', 'Thailand',
                                           'Argentina', 'Russia',
                                           'Canada', 'Poland',
                                           'Czechia', 'Serbia',
                                           'Japan'))

all_samples_by_group = df %>% 
    group_by(.data[[variable]]) %>%
    summarise(all_samples = n())
    
dfplot = NULL
all_p_values_age = NULL
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

dfplot$signature[dfplot$signature=='SBS288M_MSS'] = 'SBS_M'

dfplot$signature = factor(dfplot$signature,
                          levels = c('SBS88', 'SBS_M', 'ID14', 'ID18', 'SBS89', 'SBS94'))

####################################################################################
# My models
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
df$age_group = factor(df$age_group, levels = c('0-39', '40-49', '50-59',
                                               '60-69', '70+'), ordered = T)

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


covdata = df %>% select(donor_id, age_group, sex, tumorsite_group, country, purity)


number_sigs_by_vartype = tibble(name = sigsint) %>%
    mutate(vartype = substr(name,1,2)) %>%
    group_by(vartype) %>%
    summarise(n_sigs_vartype = n())

testdata <- dfmodel %>%
    # select(donor_id,starts_with('SBS')) %>% 
    pivot_longer(-donor_id) %>%
    left_join(covdata) %>% 
    group_by(name) %>%
    mutate(value=as.factor(value)) %>% 
    do(tresult = safely(stats::glm)(value ~ age_group + sex + tumorsite_group + country + purity,family = binomial(),data=.)) %>% 
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
    mutate(q.value = ifelse((p.value * length(sigsint))>1, 1, (p.value * length(sigsint)))) %>%
    filter(str_detect(term,variable)) %>%
    mutate(label = paste0(str_replace(term,variable,''), ' ',
                          str_replace(name,'_c','')),
           signature = str_replace(name,'_c','')) %>%
    rename(p_val = p.value,
           independent_vars = term,
           OR = estimate,
           bonferroni = q.value) %>%
    select(independent_vars, signature, OR, p_val, bonferroni, vartype, label)

testdata_GLM = testdata

# New testdata using the Firth method
testdata <- dfmodel %>%
    # select(donor_id,starts_with('SBS')) %>% 
    pivot_longer(-donor_id) %>%
    left_join(covdata) %>% 
    group_by(name) %>%
    mutate(value=as.factor(value)) %>% 
    do(tresult = safely(logistf::logistf)(value ~ age_group + sex + tumorsite_group + country + purity,
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
    mutate(q.value = ifelse((p.value * length(sigsint))>1, 1, (p.value * length(sigsint)))) %>%
    filter(str_detect(term,variable)) %>%
    mutate(label = paste0(str_replace(term,variable,''), ' ',
                          str_replace(name,'_c','')),
           signature = str_replace(name,'_c','')) %>%
    rename(p_val = p.value,
           independent_vars = term,
           OR = estimate,
           bonferroni = q.value) %>%
    select(independent_vars, signature, OR, p_val, bonferroni, vartype, label)

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
                            df$age_group)<=0)
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
                            df$age_group)<=0)
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



testdata = testdata %>%
    arrange(signature)

testdata0=testdata %>%
    filter(independent_vars == 'age_group.L') %>%
    mutate(q_val_2 = c(p.adjust(p_val[which(str_detect(signature, 'CN'))],
                                method = 'fdr'),
                       p.adjust(p_val[which(str_detect(signature, 'DBS'))],
                                method = 'fdr'),
                       p.adjust(p_val[which(str_detect(signature, 'ID'))],
                                method = 'fdr'),
                       p.adjust(p_val[which(str_detect(signature, 'SBS'))],
                                method = 'fdr'),
                       p.adjust(p_val[which(str_detect(signature, 'SV'))],
                                method = 'fdr'))
    )

testdata = rbind(testdata0)

q_values_age_group = testdata %>%
    filter(signature %in% str_replace(signatures_of_interest,'_c','')) %>%
    mutate(signature = factor(signature, levels = c('SBS88', 'SBS288M_MSS',
                                                    'ID14', 'ID18', 'SBS89'))) %>%
    arrange(signature) %>%
    pull(p_val)

q_values_age_group = paste0(expression('p-trend = '),
                            ifelse(q_values_age_group < 0.001,
                            formatC(q_values_age_group, format = 'e', digits = 1),
                            formatC(q_values_age_group, format = 'f', digits = 3))
)

q_values_age_group = c(q_values_age_group[1], rep('',4),
                       q_values_age_group[2], rep('',4),
                       q_values_age_group[3], rep('',4),
                       q_values_age_group[4], rep('',4),
                       q_values_age_group[5], rep('',4)
                       )
####################################################################################


dfplot$signature = factor(dfplot$signature, levels = c('SBS88', 'SBS_M',
                                                       'ID14', 'ID18','SBS89'))

# Including p-values from multivariate LR model
design <- c(
"AAADDD
 EEBBCC"
)

# design <- c(
#     "AB"
# )

dfplot %>%
    arrange(signature) %>%
    # filter(signature %in% c('SBS88', 'ID18') )%>%
    ggplot() +
    aes(x = .data[[variable]], y = prev_sig, fill = .data[[variable]]) +
    facet_manual(.~signature, design = design,
                 scales = 'free') +
    geom_col() +
    geom_text(aes(label = label_prev), 
              size = 4, vjust = -0.5,check_overlap = T) +
    geom_text(aes(y = +Inf, 
                  label = q_values_age_group),
              size = 4.5, vjust = 1.5, hjust = 0.1) +
    scale_y_continuous(
        labels = scales::label_percent(),
                       expand = c(0,0,0.15,0)) +
    scale_fill_manual(values = scales::seq_gradient_pal("#8C3FC0", "#73C03F", "Lab")(seq(0,1,length.out=length(unique(met[,variable]))))) +
    theme_bw() +
    labs(subtitle = 'Adjusted by sex, country, tumor subsite, and purity',
        x = 'Age of onset',
        y = 'Signature prevalence',
        title = 'Age of onset trend enrichment'
    ) +
    theme(
        plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = 'bold')) +
    guides(fill = 'none')
# Exported 1100 x 1000

# Supplementary table
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
    mutate(signature = factor(signature,levels = signature_order)) %>%
    arrange(signature) %>%
    select(signature, p_val, model_used) %>%
    rename(Signature = signature,
        'p-trend' = p_val,
        'Model' = model_used)
write.csv(testdata_to_st, '../supplementary_tables_regressions/ST20.csv', quote=F, row.names=F)









#########################################################################################################
#########################################################################################################
# Colibactin positive vs. negative

# 3d
df_new = df %>%
    mutate(colibactin_signature_SBS88orID18_label = ifelse(
        colibactin_signature_SBS88orID18 == 'Yes', 'Positive', 'Negative'),
        colibactin_signature_SBS88orID18_label = factor(
            colibactin_signature_SBS88orID18_label, levels = c('Positive',
                                                               'Negative')))

model = (lm(age_diag ~ colibactin_signature_SBS88orID18 + sex +tumorsite_group+country+purity,
           data = df))
summary(model)
p_value_adj = summary(model)$coefficients['colibactin_signature_SBS88orID18Yes','Pr(>|t|)']

df_new %>%
    ggplot() +
    aes(x = colibactin_signature_SBS88orID18_label, y = age_diag,
        fill = colibactin_signature_SBS88orID18_label) +
    ggbeeswarm::geom_quasirandom(aes(col = colibactin_signature_SBS88orID18_label)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    scale_y_continuous(limits = c(18,102)) +
    scale_color_manual(values = c("skyblue2", "grey50")) +
    scale_fill_manual(values = c("skyblue2", "grey50")) +
    ggpubr::stat_pwc(
        aes(label = ifelse(
            after_stat(p) > 0,
            sprintf("p = %2.1e", p_value_adj),
        ))) +
    labs(x = '',
         y = 'Age of onset',
         title = 'Presence of colibactin\nsignatures (SBS88 or ID18)',
         subtitle = 'Adjusted by sex, country, tumor subsite,\nand purity') +
    theme_bw()+
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
          legend.title = element_text(size = 16, face ='bold'),
          legend.text = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'left',
          strip.background = element_blank(),
          strip.text.x = element_text(size = 14, face = 'bold'),
          panel.spacing.x = unit(2,'lines')) +
    guides(fill = 'none',
           col = 'none')
# Exported 425 x 500

# Fold changes
median(df[df$colibactin_signature_SBS88orID18=='Yes',]$age_diag)
median(df[df$colibactin_signature_SBS88orID18=='No',]$age_diag)

mean(df[df$colibactin_signature_SBS88orID18=='Yes',]$age_diag)
mean(df[df$colibactin_signature_SBS88orID18=='No',]$age_diag)


# 3e
df_new = df %>%
    mutate(colibactin_signature_SBS88orID18_label = ifelse(
        colibactin_signature_SBS88orID18 == 'Yes', 'Positive', 'Negative'),
        colibactin_signature_SBS88orID18_label = factor(
            colibactin_signature_SBS88orID18_label, levels = c('Positive',
                                                               'Negative')),
        tumorsite_group = factor(tumorsite_group, levels = c('Proximal colon',
                                                             'Distal colon',
                                                             'Rectum')))

all_p_values_adj = NULL
model = (lm(age_diag ~ colibactin_signature_SBS88orID18 + sex +country + purity,
            data = df %>% filter(tumorsite_group == 'Proximal colon')))
summary(model)
p_value_adj = summary(model)$coefficients['colibactin_signature_SBS88orID18Yes','Pr(>|t|)']
all_p_values_adj = c(all_p_values_adj, p_value_adj)

model = (lm(age_diag ~ colibactin_signature_SBS88orID18 + sex +country+ purity,
            data = df %>% filter(tumorsite_group == 'Distal colon')))
summary(model)
p_value_adj = summary(model)$coefficients['colibactin_signature_SBS88orID18Yes','Pr(>|t|)']
all_p_values_adj = c(all_p_values_adj, p_value_adj)

model = (lm(age_diag ~ colibactin_signature_SBS88orID18 + sex +country+ purity,
            data = df %>% filter(tumorsite_group == 'Rectum')))
summary(model)
p_value_adj = summary(model)$coefficients['colibactin_signature_SBS88orID18Yes','Pr(>|t|)']
all_p_values_adj = c(all_p_values_adj, p_value_adj)

all_q_values_adj = p.adjust(all_p_values_adj, method = 'fdr')

df_new %>%
    mutate(tumorsite_group_plot = ifelse(tumorsite_group == 'Proximal colon',
                                         'Proximal\ncolon\n',
                                         ifelse(tumorsite_group == 'Distal colon',
                                                'Distal\ncolon\n','Rectum\n')),
           tumorsite_group_plot = factor(tumorsite_group_plot,
                                         levels = c('Proximal\ncolon\n',
                                                    'Distal\ncolon\n','Rectum\n'))) %>%
ggplot() +
    aes(x = colibactin_signature_SBS88orID18_label, y = age_diag,
        fill = colibactin_signature_SBS88orID18_label) +
    facet_wrap(tumorsite_group_plot~., dir = 'h',nrow = 1)+
    ggbeeswarm::geom_quasirandom(aes(col = colibactin_signature_SBS88orID18_label)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    scale_y_continuous(limits = c(18,102)) +
    scale_color_manual(values = c("skyblue2", "grey50")) +
    scale_fill_manual(values = c("skyblue2", "grey50")) +
    # stat_pwc(label = 'p') +
    ggpubr::stat_pwc(
        aes(label = ifelse(
            after_stat(p) > 0.148,
            sprintf("q = %5.3f", all_q_values_adj[1]),
            ifelse(
                after_stat(p) < 0.0448,
                sprintf("q = %2.1e", all_q_values_adj[2]),
                sprintf("q = %5.3f", all_q_values_adj[3])
            )
    ))) +
    labs(x = '',
         y = 'Age of onset',
         title = NULL,
         subtitle = 'Adjusted by sex, country, and purity') +
    theme_bw()+
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14, vjust = -19.5),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
          legend.title = element_text(size = 16, face ='bold'),
          legend.text = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'left',
          strip.background = element_blank(),
          strip.text.x = element_text(size = 16, face = 'bold'),
          panel.spacing.x = unit(2,'lines')) +
    guides(fill = 'none',
           col = 'none')
# Exported 500 x 500

df %>%
    filter(colibactin_signature_SBS88orID18=='Yes',
           tumorsite_group == 'Distal colon'
           ) %>%
    pull(age_diag) %>%
    mean()
df %>%
    filter(colibactin_signature_SBS88orID18=='No',
           tumorsite_group == 'Distal colon'
    ) %>%
    pull(age_diag) %>%
    mean()
df %>%
    filter(colibactin_signature_SBS88orID18=='Yes',
           tumorsite_group == 'Rectum'
    ) %>%
    pull(age_diag) %>%
    mean()
df %>%
    filter(colibactin_signature_SBS88orID18=='No',
           tumorsite_group == 'Rectum'
    ) %>%
    pull(age_diag) %>%
    mean()
