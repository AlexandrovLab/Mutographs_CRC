# Author: Marcos Diaz-Gay
# Date: Jul 21, 2024
# RStudio

library(tidyverse)
library(ggbeeswarm)
library(ggrepel)
library(ggpubr)
library(scales)
library(logistf)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')

df = read.delim('../data_for_figures/Ammal_figure_data/CRC_MSS_PKS_ECOLI_METADATA.txt')
df = tibble(df)

df = df %>%
    mutate(gen_and_pks = ifelse((colibactin_signature_SBS88orID18 == 'Yes' &
                                     pks_status == 'Positive'),'Genomic+\nand pks+',
                                 ifelse((colibactin_signature_SBS88orID18 == 'No' &
                                            pks_status == 'Positive'),'Genomic-\nand pks+',
                                         ifelse((colibactin_signature_SBS88orID18 == 'Yes' &
                                                    pks_status == 'Negative'),'Genomic+\nand pks-',
                                                 ifelse((colibactin_signature_SBS88orID18 == 'No' &
                                                            pks_status == 'Negative'),'Genomic-\nand pks-','')
                                                )
                                 )
    ))
 
dfplot = df %>%
    pivot_longer(gen_and_pks) %>%
    mutate(value = factor(value))

df = merge(df, met)

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


model = (lm(age_diag ~ gen_and_pks + sex +tumorsite_group+country+purity,
            data = df))
summary(model)

p_values_model = summary(model)$coefficients[2:4,'Pr(>|t|)']
p_values_model = paste0('p=',
                            ifelse(p_values_model < 0.001,
                                   formatC(p_values_model, format = 'e', digits = 1),
                                   formatC(p_values_model, format = 'f', digits = 2))
)
p_values_model = c('Ref.', p_values_model[1:3])
                       
                       
ggplot(dfplot) +
    aes(x = value, y = age_diag, fill = value) +
    geom_quasirandom(aes(col = value)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_text(aes(x=1, y = +Inf, label = p_values_model[1]),
              size = 4.5, vjust = 1.5, check_overlap = T) +
    geom_text(aes(x=2, y = +Inf, label = p_values_model[2]),
              size = 4.5, vjust = 1.5, check_overlap = T) +
    geom_text(aes(x=3, y = +Inf, label = p_values_model[3]),
              size = 4.5, vjust = 1.5, check_overlap = T) +
    geom_text(aes(x=4, y = +Inf, label = p_values_model[4]),
              size = 4.5, vjust = 1.5, check_overlap = T) +
    scale_fill_manual(values = c('grey50','#1FE098',
                                 '#1FC7E0','#1F67E0'))+
    scale_color_manual(values = c('grey50','#1FE098',
                                  '#1FC7E0','#1F67E0'))+
    scale_y_continuous(limits = c(18,100)) +
    theme_bw() +
    labs(title = 'Colibactin exposure status according\nto genomics and microbiomics',
         subtitle = 'Adjusted by sex, country, tumor subsite, and purity',
         x = '',
         y = 'Age of onset') +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = 'bold')
    ) +
    guides(fill = 'none',
           col = 'none'
    )
# 500 x 500


################################################################################

## Colibactin genomics vs. microbiomics
fisher.test(table(df$pks_status, df$colibactin_signature_SBS88orID18))

df2=df

table(df$pks_status, df$country) # Considering no pks positive sample is found
                                 # in Colombia we use Firth correction

model = logistf(factor(pks_status) ~ colibactin_signature_SBS88orID18 + age_diag + sex + tumorsite_group + country + purity,
            data = df2)
summary(model)
p_value_adj = summary(model)$prob['colibactin_signature_SBS88orID18Yes']

all_samples_by_group = df2 %>%
    group_by(colibactin_signature_SBS88orID18) %>%
    summarise(all_samples = n())

dfplot3 = df2 %>%
    group_by(colibactin_signature_SBS88orID18) %>%
    summarise(n_sig = sum(pks_status == 'Positive')) %>%
    left_join(all_samples_by_group) %>%
    mutate(prev_sig = n_sig/all_samples,
           label_prev = paste0(n_sig,'/',all_samples))

dfplot3 %>%
    ggplot() +
    aes(x = colibactin_signature_SBS88orID18,
        y = prev_sig,
        fill = colibactin_signature_SBS88orID18) +
    geom_col() +
    geom_text(aes(label = label_prev), size = 4.5, vjust = -0.5) +
    geom_text(aes(x = 1.5,
                  y = Inf,
                  label=paste0('p = ',
                               formatC(p_value_adj,digits = 3, format = 'f')
                  )),
              check_overlap = TRUE, size = 4.5, hjust = 0.5, vjust=1.5) +
    
    scale_y_continuous(labels = scales::label_percent(),
                       limits = c(0,0.15),
                       # breaks = c(0,0.25,0.5,0.75)
                       ) +
    scale_fill_manual(values = c("grey50", "skyblue2")) +
    theme_bw() +
    labs(title = 'Detection of pks+ cases (colibactin\nexposure based on microbiomics)',
         x = 'Presence of colibactin signatures\n(SBS88 or ID18, colibactin\nexposure based on genomics)',
         y = 'Percentage of pks+ cases',
         subtitle = 'Adjusted by age, sex, country, tumor subsite,\nand purity') +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    guides(fill = 'none')
# Exported 475 x 500


################################################################################

## AGE GROUPS - Colibactin genomics vs. microbiomics

dfplot4 = dfplot %>%
    group_by(value, age_group) %>%
    count()

ggplot(dfplot4) +
    aes(x = age_group, y = n,
        fill = value) +
    geom_col(position = 'fill') +
    scale_fill_manual(values = c('grey50','#1FE098',
                                  '#1FC7E0','#1F67E0'))+
    scale_y_continuous(labels = scales::label_percent()) +
    theme_bw() +
    labs(title = 'Colibactin exposure status by age group',
         x = 'Age of onset',
         y = 'Percentage of cases',
         fill = '') +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'left',
          legend.title = element_text(size=14),
          legend.text = element_text(size=14)
          )
# Exported 600 x 500
