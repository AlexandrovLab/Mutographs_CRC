# Author: Marcos Diaz-Gay
# Date: Aug 22, 2024
# RStudio

library(tidyverse)
library(scales)
library(ggrepel)
library(broom)
library(logistf)
library(openxlsx)

###################
# status = 'MSS'
status = 'MSI'
# variable = 'age_diag'
variable = 'age_eo'

###################


met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')

germline_vars = read.xlsx('Input_Data/Mutographs_eoCRC_Manuscript_Supplementary_Note_Tables_v16.xlsx',
                          sheet = 28, rows = 3:1000)
germline_vars = germline_vars %>%
    select(-Status) %>%
    mutate(Case.ID = str_replace(Case.ID,'b','a')) %>%
    rename(donor_id = Case.ID,
           Status = Molecular.subtype) %>%
    mutate(Status = ifelse(Status == 'HRD', 'BRCA',
                           Status)) %>%
    select(donor_id, Status, Gene) %>%
    group_by(donor_id, Status, Gene) %>%
    summarise(value = n()) %>%
    arrange(Gene) %>%
    mutate(Gene = paste0('germ_', Gene)) %>%
    pivot_wider(names_from = Gene)

met = merge(met, germline_vars, all.x = T)

colnames_met_germline = colnames(met)[str_detect(colnames(met), 'germ_')]
met[, colnames_met_germline][is.na(met[, colnames_met_germline])] = 0

met$germ_Lynch = met$germ_MLH1 + met$germ_MSH2 + met$germ_MSH6 + met$germ_PMS2
met$germ_Cowden = met$germ_SEC23B
met$germ_HRD = met$germ_BARD1 + met$germ_BLM + met$germ_BRCA1 + met$germ_BRCA2 +
    met$germ_BRIP1 + met$germ_MRE11 + met$germ_NBN + 
    met$germ_PALB2 + met$germ_RAD50 + met$germ_SLX4
met$germ_DDR = met$germ_ATM + met$germ_ATR + 
    met$germ_CHEK2 + met$germ_ERCC2 + met$germ_ERCC3 + met$germ_ERCC8 +
    met$germ_FANCC + met$germ_FANCI +met$germ_FANCL + met$germ_FANCM +
    met$germ_GTF2H5 + met$germ_MSH3 + met$germ_MUTYH + met$germ_NTHL1 +
    met$germ_PMS1 + met$germ_PNKP + met$germ_POLG + met$germ_POLH +
    met$germ_RECQL4 + met$germ_RRM2B + 
    met$germ_SPRTN + met$germ_TREX1 + met$germ_WRN +
    met$germ_XPA + met$germ_XRCC4 + met$germ_ZSWIM7
    

germline_syndromes = c('germ_Lynch', 'germ_Cowden',
                       'germ_HRD', 'germ_DDR')

df = met %>% filter(Status == status)

sigsint = c(colnames_met_germline, germline_syndromes)

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
threshold = 1
prevalence = (df %>%
                  summarise(across(all_of(sigsint),~sum(. > 0))) / nrow(df))
genes_to_not_consider = names(prevalence)[prevalence==0]

colnames_met_germline = colnames_met_germline[!colnames_met_germline %in% genes_to_not_consider]
sigsint = c(colnames_met_germline, germline_syndromes)
prevalence = prevalence[!names(prevalence) %in% genes_to_not_consider]

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
    mutate(vartype = c(rep('GV', length(colnames_met_germline)),
                       rep('GS', length(germline_syndromes))
                       )) %>%
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
    mutate(vartype = ifelse(name %in% c(germline_syndromes), 'GS', 'GV')) %>%
    left_join(number_sigs_by_vartype) %>%
    filter(str_detect(term,variable)) %>%
    mutate(label = paste0(str_replace(term,variable,''), ' ',
                          # mutate(label = paste0(country_to_test, ' ',
                          str_replace(name,' ','')),
           signature = str_replace(name,' ',''),
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
    mutate(vartype = ifelse(name %in% c(germline_syndromes), 'GS', 'GV')) %>%
    left_join(number_sigs_by_vartype) %>%
    filter(str_detect(term,variable)) %>%
    mutate(label = paste0(str_replace(term,variable,''), ' ',
                          # mutate(label = paste0(country_to_test, ' ',
                          str_replace(name,' ','')),
           signature = str_replace(name,' ',''),
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
    
    if (paste0(signature) %in% c(sigs_high_prev)){
        is_complete_separation = as.logical(
            length(
                which(table(df[,paste0(signature)]>median(df[,paste0(signature,'_c')], na.rm = T),
                            df$age_eo)<=0)
            ) +                length(
                which(table(df[,paste0(signature)]>median(df[,paste0(signature,'_c')], na.rm = T),
                            df$country)<=0)
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
                            df$age_eo)<=0)
            ) +
                length(
                    which(table(df[,paste0(signature)]>0,
                                df$country)<=0)
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
if (variable == 'age_diag'){
    testdata$OR = exp(testdata$OR * 10)
}


# Volcano plot LR
testdata = testdata %>%
    arrange(vartype) %>%
    mutate(q_val = p.adjust(p_val, method='fdr')) %>%
    mutate(q_val_2 = c(p.adjust(p_val[which(str_detect(vartype, 'GS'))],
                                method = 'fdr'),
                       p.adjust(p_val[which(str_detect(vartype, 'GV'))],
                                method = 'fdr')
    )
    )

complete_testdata = testdata

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
        vartype = str_replace(vartype, 'GV', 'Genes with P/LP variants'),
        vartype = str_replace(vartype, 'GS', 'P/LP variants in CRC predisposition\nsyndromes, HR, or DDR genes'),
        vartype = factor(vartype, levels = c(
            'Genes with P/LP variants',
            'P/LP variants in CRC predisposition\nsyndromes, HR, or DDR genes'))) %>%
    mutate(signature = str_replace(signature, 'germ_', ''),
           signature = paste0(signature,' syndrome')) %>%
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
    scale_y_continuous(
        # breaks = c(0,2,4,6),
        # labels = c('0','2','4','>6')
    ) +
    scale_color_manual(breaks = c('Enriched in early-onset patients',
                                  'Enriched in late-onset patients'),
                       values = c('#8C3FC0','#73C03F'),
                       na.value = 'grey80') +
    geom_text_repel(data = . %>% 
                        mutate(signature = ifelse(q_val_2 < 0.05, signature, '')),
                    aes(label = signature), size = 5, seed = 2, nudge_x = 0.08, nudge_y = 0.05) +
    {if (variable=='age_diag') labs(x = 'Odds ratio (10 years; log2)') else labs(x = 'Odds ratio (log2)')} +
    labs(y = 'q-value (-log10)',
         # col = 'Age of onset\nenrichment',
         col = '',
         title = 'Age of onset enrichment of P/LP germline variants in MSI cases (n=39/153)',
         subtitle = 'Adjusted by sex, country, tumor subsite, and purity'
    ) +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 16, face ='bold'),
          legend.text = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'left',
          # legend.box.spacing = unit(25,'pt'),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 14, face = 'bold'),
          panel.spacing.x = unit(1.5,'lines')
    ) +
    guides(col = guide_legend(override.aes = list(label = ''),))
# Exported 1100 x 500


# Supplementary table
testdata = complete_testdata
testdata_to_st = testdata %>%
    mutate(signature = str_replace(signature, 'germ_', '')) %>%
    mutate(Analysis = ifelse(vartype == 'GV', 'Individual genes',
                             'Gene classification')) %>%
    arrange(desc(Analysis)) %>%
    select(Analysis, signature, OR, p_val, q_val_2,model_used) %>%
    rename('Gene/Gene classification' = signature,
           'p-value' = p_val,
           'q-value' = q_val_2,
           'Model' = model_used)

write.csv(testdata_to_st, '../supplementary_tables_regressions/SNT30.csv', quote=F, row.names=F)
