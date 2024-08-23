# Author: Marcos Diaz-Gay
# Date: Aug 21, 2024
# RStudio

library(tidyverse)
library(logistf)

################################################################################
# Enrichment of eoCRC cases in countries, tumor subsites, and sexes

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')
met_no_Lynch = met %>%
    filter(Lynch == 'No')

met$age_eo = factor(met$age_eo,
                    levels = c('50+','0-49'))
met$sex = factor(met$sex,
                 levels = c('Male','Female'))
met$tumorsite_group = factor(met$tumorsite_group,
                             levels = c('Proximal colon', 'Distal colon', 'Rectum'))
met$tumorsite_group_v3 = factor(met$tumorsite_group_v3,
                             levels = c('Left colon', 'Right colon'))
met$country = factor(met$country)
met$country = relevel(met$country, ref = 'Brazil')

## Countries
countries = sort(unique(as.character(met$country)))

ps = NULL
ors = NULL
for (i in 1:length(countries)){
    country_model = countries[i]
    print(country_model)
    
    model = summary(glm(age_eo ~ (country==country_model) + sex + tumorsite_group_v3,
                data = met, family = binomial))
    p = model$coefficients['country == country_modelTRUE','Pr(>|z|)']
    or = exp(model$coefficients['country == country_modelTRUE','Estimate'])
    
    ps = c(ps, p)
    ors = c(ors, or)
}

sign = data.frame(country = countries,
                  OR = ors,
                  p_value = ps,
                  q_value = p.adjust(ps, method = 'fdr'))

apply(table(met$age_eo, met$country == 'Iran'),2,prop.table)

## Countries (no Lynch cases)
countries = sort(unique(as.character(met_no_Lynch$country)))

ps = NULL
ors = NULL
for (i in 1:length(countries)){
    country_model = countries[i]
    print(country_model)
    
    model = summary(glm(age_eo ~ (country==country_model) + sex + tumorsite_group_v3,
                        data = met_no_Lynch, family = binomial))
    p = model$coefficients['country == country_modelTRUE','Pr(>|z|)']
    or = exp(model$coefficients['country == country_modelTRUE','Estimate'])
    
    ps = c(ps, p)
    ors = c(ors, or)
}

sign = data.frame(country = countries,
                  OR = ors,
                  p_value = ps,
                  q_value = p.adjust(ps, method = 'fdr'))

apply(table(met_no_Lynch$age_eo, met_no_Lynch$country == 'Iran'),2,prop.table)

## Tumor subsites (distal and rectum together)
countries = sort(unique(as.character(met$tumorsite_group_v3)))

ps = NULL
ors = NULL
for (i in 1:length(countries)){
    country_model = countries[i]
    print(country_model)
    
    model = summary(glm(age_eo ~ (tumorsite_group_v3==country_model) + sex + country,
                        data = met, family = binomial))
    p = model$coefficients['tumorsite_group_v3 == country_modelTRUE','Pr(>|z|)']
    or = exp(model$coefficients['tumorsite_group_v3 == country_modelTRUE','Estimate'])
    
    ps = c(ps, p)
    ors = c(ors, or)
}

sign = data.frame(country = countries,
                  OR = ors,
                  p_value = ps,
                  q_value = p.adjust(ps, method = 'fdr'))

apply(table(met$age_eo, met$tumorsite_group_v3),2,prop.table)

## Sexes
countries = sort(unique(as.character(met$sex)))

ps = NULL
ors = NULL
for (i in 1:length(countries)){
    country_model = countries[i]
    print(country_model)
    
    model = summary(glm(age_eo ~ (sex==country_model) + tumorsite_group_v3 + country,
                        data = met, family = binomial))
    p = model$coefficients['sex == country_modelTRUE','Pr(>|z|)']
    or = exp(model$coefficients['sex == country_modelTRUE','Estimate'])
    
    ps = c(ps, p)
    ors = c(ors, or)
}

sign = data.frame(country = countries,
                  OR = ors,
                  p_value = ps,
                  q_value = p.adjust(ps, method = 'fdr'))

apply(table(met$age_eo, met$sex),2,prop.table)


################################################################################
# Enrichment of DNA repair deficiencies

## Require Firth correction due to complete separation
summary(logistf((Status=='MUTYH') ~ age_eo + sex + tumorsite_group_v3 + country,
            data = met, control =  logistf.control(maxit = 1000)))
summary(logistf((Status=='NTHL1') ~ age_eo + sex + tumorsite_group_v3 + country,
                data = met, control =  logistf.control(maxit = 1000)))
summary(logistf((Status=='OGG1') ~ age_eo + sex + tumorsite_group_v3 + country,
                data = met, control =  logistf.control(maxit = 1000)))
summary(logistf((Status=='POLD') ~ age_eo + sex + tumorsite_group_v3 + country,
                data = met, control =  logistf.control(maxit = 1000)))
summary(logistf((Status=='POLE') ~ age_eo + sex + tumorsite_group_v3 + country,
                data = met, control =  logistf.control(maxit = 1000)))
summary(logistf((Status=='BRCA') ~ age_eo + sex + tumorsite_group_v3 + country,
                data = met, control =  logistf.control(maxit = 1000)))

## Do not require Firth correction
summary(glm((Status=='MSI') ~ age_eo + sex + tumorsite_group_v3 + country,
            data = met, family = binomial))



################################################################################
# MSI cases without Lynch cases
met_no_Lynch = met %>%
    filter(Lynch == 'No')
met_no_Lynch$tumorsite_group_v3 = factor(met_no_Lynch$tumorsite_group_v3,
                                levels = c('Left colon', 'Right colon'))
summary(glm((Status=='MSI') ~ age_eo + sex + tumorsite_group_v3 + country,
            data = met_no_Lynch, family = binomial))

table(met$Lynch, met$Status)

################################################################################
## MSI - Countries
countries = sort(unique(as.character(met$country)))

ps = NULL
ors = NULL
for (i in 1:length(countries)){
    country_model = countries[i]
    print(country_model)
    
    model = summary(glm((Status=='MSI') ~ age_eo + (country==country_model) + sex + tumorsite_group_v3,
                        data = met, family = binomial))
    p = model$coefficients['country == country_modelTRUE','Pr(>|z|)']
    or = exp(model$coefficients['country == country_modelTRUE','Estimate'])
    
    ps = c(ps, p)
    ors = c(ors, or)
}

sign = data.frame(country = countries,
                  OR = ors,
                  p_value = ps,
                  q_value = p.adjust(ps, method = 'fdr'))

apply(table(met$Status=='MSI', met$country == 'Canada'),2,prop.table)


## MSI - Countries (no Lynch cases)
countries = sort(unique(as.character(met_no_Lynch$country)))

ps = NULL
ors = NULL
for (i in 1:length(countries)){
    country_model = countries[i]
    print(country_model)
    
    model = summary(glm((Status=='MSI') ~ age_eo + (country==country_model) + sex + tumorsite_group_v3,
                        data = met_no_Lynch, family = binomial))
    p = model$coefficients['country == country_modelTRUE','Pr(>|z|)']
    or = exp(model$coefficients['country == country_modelTRUE','Estimate'])
    
    ps = c(ps, p)
    ors = c(ors, or)
}

sign = data.frame(country = countries,
                  OR = ors,
                  p_value = ps,
                  q_value = p.adjust(ps, method = 'fdr'))

apply(table(met_no_Lynch$Status=='MSI', met_no_Lynch$country == 'Canada'),2,prop.table)
apply(table(met_no_Lynch$Status=='MSI', met_no_Lynch$country == 'Brazil'),2,prop.table)
