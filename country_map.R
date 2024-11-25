library(tidyverse)
library(sf)
library(rnaturalearth)
library(tidygeocoder)
  

# Show Country Map of Sample Sizes ----------------------------------------

#setwd('~/Dropbox/UNC/Fall2024/BIOS784/Final_Project/BIOS-784-Project/')
sample_data <- read.csv('Pf7-samples.csv')
sample_data_qc <- sample_data %>% dplyr::filter(qc_pass==T)
  
# get sample counts by country
country_data <- sample_data_qc %>%
  group_by(country) %>%
  summarise(count = n())

# get sample counts/proportions by region
site_data <- sample_data_qc %>%
  group_by(country) %>%
  mutate(n = n()) %>%
  group_by(country, site) %>%
  summarise(siteCount = n(), siteProp = n()/unique(n)[1], .groups = "drop")
  
# get long/lat of each region
site_coords_city <- site_data %>%
  geocode(country = country, city = site, method = "osm")
site_coords_reg <- site_coords_city %>%
  dplyr::filter(is.na(lat)) %>%
  dplyr::select(-lat, -long) %>%
  geocode(country = country, state = site, method = "osm")
site_coords <- site_coords_city
site_coords[is.na(site_coords$lat),c("lat","long")] <- site_coords_reg[,c("lat","long")]

world <- ne_countries(scale = "medium", returnclass = "sf")
country_data_rename <- country_data %>%
  mutate(country = ifelse(country == "Democratic Republic of the Congo", "Dem. Rep. Congo", country))
world_data <- world %>%
  left_join(country_data_rename, by = c("name" = "country"))

ggplot()+
  geom_sf(data = world_data, aes(fill = count), color = "white")+
  scale_fill_viridis_c(option = "plasma", na.value = "gray90", name = "Sample Count per Country")+
  geom_point(data = site_coords, aes(x = long, y = lat, size = siteCount), color = "red", alpha = .7)+
  scale_size_continuous(name = "Sample Count per Site", range = c(1,5))+
  coord_sf(xlim = c(-120, 120), ylim = c(-30, 30)) +
  theme_minimal()+
  labs(title = "Global Map of MalariaGen Parasite Samples")
