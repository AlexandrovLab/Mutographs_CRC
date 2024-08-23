# Author: Marcos Diaz-Gay
# Date: Jul 21, 2024
# RStudio

library(tidyverse)
library(scales)
library(sf)

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')

met$country_to_map = met$country
met$country_to_map[met$country=="Czech Republic"] = "Czechia"
met$country_to_map[met$country=="Serbia"] = "Republic of Serbia"

df = met %>%
    group_by(country_to_map) %>%
    summarise(ASR = mean(ASR_CRC_both),
              n = n(),
              perc_eo = sum(age_eo == '0-49')/n)
    
sdf <- rnaturalearthdata::countries50 %>% 
    st_as_sf() %>%
    left_join(df, by = c("admin" = "country_to_map"))

sdf$admin[sdf$admin=='Republic of Serbia'] = 'Serbia'
df$country_to_map[df$country_to_map=="Republic of Serbia"] = "Serbia"
met$country_to_map[met$country_to_map=="Republic of Serbia"] = "Serbia"

ranking <- st_geometry(sdf) %>% 
    st_point_on_surface() %>% 
    st_coordinates() %>% 
    as_tibble() %>% 
    bind_cols(tibble(country = sdf$admin, Total=sdf$n)) %>%
    filter(country %in% df$country_to_map)

sdf = sdf %>%
    mutate(labels=paste0(admin, '\n(n=', n, '; ',
                         round(perc_eo*100,1),'% eo)'))


df_tumorsite = met %>%
    group_by(country_to_map, tumorsite_group) %>%
    summarise(Total_by_tumorsite = n())
sdf_tumorsite = sdf %>%
    left_join(df_tumorsite, by = c("admin" = "country_to_map")) %>%
    select('label_x', 'label_y', 'admin', 'tumorsite_group', 'Total_by_tumorsite') %>%
    filter(admin %in% df$country_to_map)

# sdf_tumorsite = data.frame(sdf_tumorsite)[,-6]
sdf_tumorsite$Total_by_tumorsite = as.numeric(sdf_tumorsite$Total_by_tumorsite)

sdf_tumorsite = sdf_tumorsite %>%
    as_tibble() %>%
    pivot_wider(names_from = tumorsite_group, values_from = Total_by_tumorsite)

sdf_tumorsite[is.na(sdf_tumorsite)] = 0

ggplot() + 
    geom_sf(data = sdf,aes(fill=ASR), size = .2, color = 'white') +
    # geom_point(data = sdf %>% filter(admin %in% df$country_to_map),
    #            aes(x = label_x, y = label_y, col=perc_eo, size = n),
    #            pch=20)+
    # geom_scatterpie(data = sdf_tumorsite,
    #                 aes(x = label_x, y = label_y),
    #                 cols = c('Proximal colon','Distal colon', 'Rectum', 'Other'))
    ggrepel::geom_text_repel(data = sdf %>% filter(admin %in% df$country_to_map),
                             aes(x = label_x, y = label_y, label = labels),size=5,fontface=2,
                             box.padding = 1.25,min.segment.length = 1.5) +
    #geom_text(data = ranking, aes(x = X-.5, y = Y, label = country, color = SP_Group_New), hjust = 1, size = 2.5, nudge_y = .5) +
    # scale_fill_viridis_c(direction = -1,na.value = 'grey80')+
    # scale_fill_steps(low = "#fffaf4", high = "#DF8807",na.value = 'grey80',
    #                  # breaks = c(15 , 25, 35, 40), 
    #                  ) +
    scale_fill_viridis_c(option = 'A', direction = -1, na.value = 'grey80') +
    scale_color_steps(low = "#f5ebfc", high = "#8C3FC0",
                         labels = scales::label_percent()) +
    scale_size_continuous(breaks = c(50,100), range = c(3,8))+
    coord_sf(ylim = c(-50,90)) +
    #coord_sf(xlim=c(-150,150),ylim=c(5, 60),clip = "off") +
    theme_void() +
    # guides(fill = "none")+
    theme(legend.position = "top",
          legend.title = element_text(face = 'bold', size = 15),
          legend.text = element_text(size=13),
          legend.box = 'horizontal', legend.spacing.x = unit(1,'cm')) +
    guides(size = guide_bins(barwidth = 2, order = 1,
                               title = 'Number of\ncases'),
           fill = guide_colorsteps(barwidth = 12, order = 2,
                                 title = 'ASR per\n100,000'),
           color = guide_colorsteps(barwidth = 12, order = 3,
                                  title = '% of early onset\ncases (<50 y.o.)'))
# Exported 1200 x 720 (only for testing)
# Exported 12.5 x 7.5 inches (PDF) and later modified in Illustrator
