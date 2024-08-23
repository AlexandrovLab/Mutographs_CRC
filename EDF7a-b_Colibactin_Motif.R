# Author: Marcos Diaz-Gay
# Date: Jul 21, 2024
# RStudio

library(tidyverse)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')

mot = read.csv('../data_for_figures/Mariya_figure_data/motif_analysis/plots_for_marcos/data_for_motif_plots.csv')

df = merge(mot, met)
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
df$age_group = factor(df$age_group,ordered = T)

# Percentage motif mutations vs. age
model = lm(percent_with_motif_any ~ age_group + sex + country + tumorsite_group + purity,
           data = df)
summary(model)


p_values_model = summary(model)$coefficients[2,'Pr(>|t|)']
p_values_model = paste0('p-trend = ',
                        ifelse(p_values_model < 0.01,
                               formatC(p_values_model, format = 'e', digits = 1),
                               formatC(p_values_model, format = 'f', digits = 2))
)
p_values_model = c(p_values_model, rep('',4))


ggplot(mot) +
    aes(x=age_group, y = percent_with_motif_any, fill = age_group) +
    ggbeeswarm::geom_quasirandom(aes(col = age_group)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_text(aes(x=1, y = +Inf, label = p_values_model[1]),
              size = 4.5, vjust = 1.5, check_overlap = T, hjust = 0.1) +
    labs(x = 'Age of onset',
         title = 'W[T>N]W mutations with colibactin motif (WAWW[T>N]W)',
         y = 'Percentage of mutations',
         subtitle = 'Adjusted by sex, country, tumor subsite, and purity') +
    scale_fill_manual(values = scales::seq_gradient_pal("#8C3FC0", "#73C03F", "Lab")(seq(0,1,length.out=5))) +
    scale_color_manual(values = scales::seq_gradient_pal("#8C3FC0", "#73C03F", "Lab")(seq(0,1,length.out=5))) +
    scale_y_continuous(labels = scales::label_percent(),
                       limits = c(0.13, 0.7)) +
    theme_bw()+
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 16, face ='bold'),
          legend.text = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'left',
          # legend.box.background = element_rect(linewidth = 1),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 16, face = 'bold'),
          panel.spacing.x = unit(2,'lines')) +
    guides(col = 'none',
           fill = 'none')
# Exported 750 x 500


# Percentage motif mutations vs. colibactin sigs
model = lm(percent_with_motif_any ~ colibactin_signature_SBS88orID18 + age_group + sex + country + tumorsite_group + purity,
           data = df)
summary(model)


p_values_model = summary(model)$coefficients[2,'Pr(>|t|)']
p_values_model = paste0('p=',
                        ifelse(p_values_model < 0.01,
                               formatC(p_values_model, format = 'e', digits = 1),
                               formatC(p_values_model, format = 'f', digits = 2))
)

df$colibactin_signature_SBS88orID18[df$colibactin_signature_SBS88orID18=='Yes'] =
    'Positive'
df$colibactin_signature_SBS88orID18[df$colibactin_signature_SBS88orID18=='No'] =
    'Negative'
df$colibactin_signature_SBS88orID18 = factor(df$colibactin_signature_SBS88orID18)
df$colibactin_signature_SBS88orID18 = relevel(df$colibactin_signature_SBS88orID18,
                                              ref = 'Positive')
ggplot(df) +
    aes(x=colibactin_signature_SBS88orID18,
        y = percent_with_motif_any, fill = colibactin_signature_SBS88orID18) +
    ggbeeswarm::geom_quasirandom(aes(col = colibactin_signature_SBS88orID18)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    ggpubr::stat_pwc(
        aes(label = ifelse(
            after_stat(p) > 0,
            p_values_model,
        ))) +
    labs(x = 'Presence of colibactin signatures (SBS88 or ID18)',
         title = 'W[T>N]W mutations with colibactin motif (WAWW[T>N]W)',
         y = 'Percentage of mutations',
         subtitle = 'Adjusted by age, sex, country, tumor subsite, and purity') +
    scale_color_manual(values = c("skyblue2", "grey50")) +
    scale_fill_manual(values = c("skyblue2", "grey50")) +
    scale_y_continuous(labels = scales::label_percent(),
                       limits = c(0.13, 0.7)) +
    theme_bw()+
    theme(plot.title = element_text(size = 16, face = 'bold'),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 16, face ='bold'),
          legend.text = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'left',
          strip.background = element_blank(),
          strip.text.x = element_text(size = 16, face = 'bold'),
          panel.spacing.x = unit(2,'lines')) +
    guides(col = 'none',
           fill = 'none')
# Exported 750 x 500
