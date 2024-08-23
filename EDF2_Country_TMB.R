# Author: Marcos Diaz-Gay
# Date: Jul 21, 2024
# RStudio

library(tidyverse)
library(scales)
library(ggbeeswarm)
library(rstatix)
library(ggpubr)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')
met$country[met$country=='Czech Republic'] = 'Czechia'

met = met %>%
    filter(Status == 'MSS')
met$sex = factor(met$sex,
                    levels = c('Male','Female'))
met$tumorsite_group = factor(met$tumorsite_group,
                    levels = c('Proximal colon', 'Distal colon', 'Rectum'))

## TMB SBS/ID
tmb_types = c('TMB_SBS', 'TMB_ID', 'TMB_DBS',
              'TMB_CN', 'TMB_SV')
countries = unique(met$country)
all_p_values = NULL
for (h in 1:length(tmb_types)){
    tmb_to_test = tmb_types[h]
    print(tmb_to_test)
    
    for (i in 1:length(countries)){
        country_to_test = countries[i]
        print(country_to_test)
        
        model = lm(log10(get(tmb_to_test)) ~ (country==country_to_test) + age_diag + sex + tumorsite_group + purity,
                   data = met)
        p_value = summary(model)$coefficients['country == country_to_testTRUE', 'Pr(>|t|)']
        all_p_values = c(all_p_values, p_value)
    }
}
testdata = tibble(
    tmb_type = rep(tmb_types, each = length(countries)),
    country = rep(countries, length(tmb_types)),
    p_val = all_p_values,
    q_val = c(p.adjust(p_val[tmb_type == 'TMB_SBS'], method = 'fdr'),
              p.adjust(p_val[tmb_type == 'TMB_ID'], method = 'fdr'),
              p.adjust(p_val[tmb_type == 'TMB_DBS'], method = 'fdr'),
              p.adjust(p_val[tmb_type == 'TMB_CN'], method = 'fdr'),
              p.adjust(p_val[tmb_type == 'TMB_SV'], method = 'fdr')
              ),
    q_val_label = paste0('q=',
                         ifelse(q_val < 0.001,
                                formatC(q_val, format = 'e', digits = 0),
                                formatC(q_val, format = 'f', digits = 3))
    )
)

testdata = testdata %>%
    mutate(q_val_label = ifelse(q_val < 0.05, q_val_label, 'n.s.'))

met %>%
    pivot_longer(cols = all_of(tmb_types), names_to = 'tmb_type') %>%
    left_join(testdata) %>%
    mutate(tmb_type_label = str_replace(tmb_type,'TMB_',''),
           tmb_type_label = factor(tmb_type_label,
                                   levels = c('SBS','ID','DBS','CN','SV')))%>%
ggplot() +
    aes(x=reorder(country, ASR_CRC_both), y = value,
        fill = ASR_CRC_both) +
    facet_wrap(.~tmb_type_label, scales = 'free', ncol = 2) +
    geom_quasirandom(aes(col = ASR_CRC_both), size=1) +
    geom_hline(aes(yintercept = ifelse(tmb_type == 'TMB_SBS',median(met$TMB_SBS),
                                       ifelse(tmb_type == 'TMB_ID',median(met$TMB_ID),
                                              ifelse(tmb_type == 'TMB_DBS',median(met$TMB_DBS),
                                                     ifelse(tmb_type == 'TMB_CN',median(met$TMB_CN,na.rm = T),
                                                            median(met$TMB_SV,na.rm = T))))),
                   linetype = ''),
               col = 'royalblue', linewidth = 1) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_text(aes(y = +Inf, label = q_val_label), check_overlap = TRUE,
              vjust = 2, size = 4.2) +
    scale_y_continuous(transform = 'log10',
        labels = scales::label_comma(),
        expand = c(0.05,0,0.15,0)
        ) +
    scale_fill_viridis_c(option = 'A', direction = -1, na.value = 'grey80') +
    scale_color_viridis_c(option = 'A', direction = -1, na.value = 'grey80') +
    # scale_fill_steps(low = "#fde0b6", high = "#c17506",na.value = 'grey80',
    #                  # breaks = c(15 , 25, 35, 40), 
    # ) +
    # scale_color_steps(low = "#fde0b6", high = "#c17506",na.value = 'grey80',
    #                  # breaks = c(15 , 25, 35, 40), 
    # ) +
    scale_linetype_manual(name = "Median TMB across countries", values = 1) +
    theme_bw() +
    labs(title = 'MSS molecular subgroup - Country distribution of mutation burden',
         subtitle = 'Adjusted by age, sex, tumor subsite, and purity',
         y = 'Number of mutations') +
    theme(legend.position = "top",
          legend.title = element_text(face = 'bold', size = 15),
          legend.text = element_text(size=13),
          legend.box = 'horizontal', legend.spacing.x = unit(1,'cm'),
          legend.justification = 'left',
          plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 30, hjust = 1),
          strip.text = element_text(size=14, face = 'bold'),
          strip.background = element_blank(),
          legend.key.width = unit(2,'cm')) +
    guides(fill = guide_colorsteps(barwidth = 12, order = 2,
                                   title = 'ASR per 100,000   '),
           col = 'none')
# Exported 1650 x 1500
