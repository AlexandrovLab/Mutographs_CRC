# Author: Marcos Diaz-Gay
# Date: Jul 21, 2024
# RStudio

library(tidyverse)
library(scales)
library(ggbeeswarm)
library(broom) # For tidy
library(ggh4x) # For facet_manual

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')
met$country[met$country == 'Czech Republic'] = 'Czechia'

status = 'MSS'
variable = 'age_group'
signatures_of_interest = c('colibactin_signature_SBS88orID18')

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


dfplot_ind = df %>%
    filter(Status == status) %>%
    group_by(.data[[variable]]) %>%
    summarise(n_sig = sum(.data[[signatures_of_interest]]=='Yes')) %>%
    left_join(all_samples_by_group) %>%
    mutate(prev_sig = n_sig/all_samples,
           label_prev = paste0(n_sig,'/',all_samples),
           signature = signatures_of_interest)
dfplot = dfplot_ind


####################################################################################
# My models
sigsint = signatures_of_interest

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
# df$age_group = factor(df$age_group)
# df$age_group = relevel(df$age_group, '70+')
df$age_group = factor(df$age_group, levels = c('0-39', '40-49', '50-59',
                                               '60-69', '70+'), ordered = T)

# Dichotomizing cases over 70% prevalence by above/below median
dfmodel = as_tibble(df[,c('donor_id',sigsint)])
threshold = 0.7
prevalence = (df %>%
                  summarise(across(all_of(sigsint),~sum(. == 'Yes'))) / nrow(df))

sigs_high_prev = names(prevalence)[prevalence > threshold]
sigs_low_prev = names(prevalence)[prevalence <= threshold]

dfmodel_low_prev = dfmodel[,c('donor_id',sigs_low_prev)]
dfmodel_low_prev[,sigs_low_prev] = dfmodel_low_prev[,sigs_low_prev]=='Yes'

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
    mutate(vartype = substr(name,1,2)) %>%
    left_join(number_sigs_by_vartype) %>%
    # mutate(q.value = ifelse((p.value * n_sigs_vartype)>1, 1, (p.value * n_sigs_vartype))) %>%
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

testdata = testdata %>%
    arrange(signature)

testdata0=testdata %>%
    filter(independent_vars == 'age_group.L') %>%
    mutate(q_val_2 = c(p.adjust(p_val,
                                method = 'fdr'))
    )

testdata = rbind(testdata0)

q_values_age_group = testdata %>%
    # filter(signature %in% str_replace(signatures_of_interest,'_c','')) %>%
    # mutate(signature = factor(signature, levels = c('SBS88', 'SBS288M_MSS',
                                                    # 'ID14', 'ID18', 'SBS89'))) %>%
    arrange(signature) %>%
    pull(p_val)

q_values_age_group = paste0(expression('p-trend = '),
                            ifelse(q_values_age_group < 0.001,
                                   formatC(q_values_age_group, format = 'e', digits = 1),
                                   formatC(q_values_age_group, format = 'f', digits = 3))
)


q_values_age_group = c(q_values_age_group[1], rep('',4)
)
####################################################################################


dfplot %>%
    mutate(signature = 'Colibactin: SBS88 or ID18') %>%
    arrange(signature) %>%
    ggplot() +
    aes(x = .data[[variable]], y = prev_sig, fill = .data[[variable]]) +
    facet_wrap(.~signature,
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
         y = 'Prevalence',
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
# Exported 750 x 500


