# Author: Marcos Diaz-Gay
# Date: Jul 21, 2024
# RStudio

library(tidyverse)
library(scales)
library(maps)
library(sf)

###################
# signature = 'SBS_F_c'
signature = 'ID_J_c'
###################

met = read.delim('../data_for_figures/Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21.tsv')

colnames(met)[colnames(met)=="SBS288F_MSS_c"] = 'SBS_F_c'
colnames(met)[colnames(met)=="SBS288H_MSS_c"] = 'SBS_H_c'
colnames(met)[colnames(met)=="SBS288M_MSS_c"] = 'SBS_M_c'
colnames(met)[colnames(met)=="SBS288O_MSS_c"] = 'SBS_O_c'
colnames(met)[colnames(met)=="ID83J_MSS_c"] = 'ID_J_c'

met$country_to_map = met$country
met$country_to_map[met$country=="Czech Republic"] = "Czechia"
met$country_to_map[met$country=="Serbia"] = "Republic of Serbia"

df = met %>%
    filter(Status == 'MSS') %>%
    group_by(country_to_map) %>%
    summarise(sig_tot = sum(.data[[signature]]>0),
              n = n(),
              sig_prev = sig_tot/n)

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
    mutate(labels=paste0(admin, ' (',sig_tot, '/', n, ')'))


df_tumorsite = met %>%
    group_by(country_to_map, tumorsite_group) %>%
    summarise(Total_by_tumorsite = n())
sdf_tumorsite = sdf %>%
    left_join(df_tumorsite, by = c("admin" = "country_to_map")) %>%
    select('label_x', 'label_y', 'admin', 'tumorsite_group', 'Total_by_tumorsite') %>%
    filter(admin %in% df$country_to_map)

sdf_tumorsite$Total_by_tumorsite = as.numeric(sdf_tumorsite$Total_by_tumorsite)

sdf_tumorsite = sdf_tumorsite %>%
    as_tibble() %>%
    pivot_wider(names_from = tumorsite_group, values_from = Total_by_tumorsite)

sdf_tumorsite[is.na(sdf_tumorsite)] = 0

ggplot() + 
    geom_sf(data = sdf,aes(fill=sig_prev), size = .2, color = 'white') +
    ggrepel::geom_text_repel(data = sdf %>% filter(admin %in% df$country_to_map),
                             aes(x = label_x, y = label_y, label = labels),size=4,
                             box.padding = 1.25,min.segment.length = 1.5) +
    scale_fill_viridis_c(direction=-1,na.value = 'grey80',labels = label_percent(),
                         breaks = pretty_breaks(n=4), option = 'A') +
    coord_sf(ylim = c(-50,90)) +
    theme_void() +
    theme(legend.position = "left",
          legend.title = element_text(size = 15,face = 'bold'),
          legend.text = element_text(size=13),
          legend.justification =  'top',
          legend.box.margin = margin(r=-20, t = 50),
          legend.box = 'horizontal',
          legend.spacing.x = unit(0.5,'cm'),
          legend.box.just = 'center',
          plot.margin = margin(1,1,1,1)) +
    guides(fill = guide_colorbar(title = paste0(str_replace(signature,'_c',''),
                                                '\nprevalence\n')))
# Exported 800 x 300
# Exported 8 x 3 inches in PDF to modify manually in Illustrator