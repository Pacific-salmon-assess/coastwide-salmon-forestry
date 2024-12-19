# Load libraries

library(tidyverse)
library(here)
library(ersst)
library(sf)
library(bcmaps)
library(ggplot2)
library(hues)

# Load data

data <- read.csv(here("origional-ecofish-data-models","Data","Processed","chum_SR_20_hat_yr_w_sstcobe_new.csv"))

data_pink_e <- read.csv(here("origional-ecofish-data-models","Data","Processed","pke_SR_10_hat_yr_w_sstcobe_new.csv"))

data_pink_o <- read.csv(here("origional-ecofish-data-models","Data","Processed","pko_SR_10_hat_yr_w_sstcobe_new.csv"))

min_max <- rbind(data %>% 
                   select(Y_LAT, X_LONG), data_pink_e %>% 
                   select(Y_LAT, X_LONG), data_pink_o %>% 
                   select(Y_LAT, X_LONG)) %>% 
  summarise_all(list(min = min, max = max))


years <- rbind(data %>% 
                 select(BroodYear), data_pink_e %>% 
                 select(BroodYear), data_pink_o %>% 
                 select(BroodYear)) %>% 
  distinct() %>% 
  arrange(BroodYear) %>% 
  mutate(sst_year = BroodYear + 1)


# sst_download(years = years$sst_year, months = 4:7, save.dir = here("data_processing", "sst_ersst"),
#              version = 5)

sst_download(years = 1996:2014, months = 4:7, save.dir = here("data_processing", "sst_ersst"),
                          version = 5)

sst <- sst_load(years$sst_year, 4:7, here("data_processing", "sst_ersst"), version = 5)


# subset data 

sst_subset <- sst_subset_space(sst, 
                               lat.min = min_max$Y_LAT_min-2, 
                               lat.max = min_max$Y_LAT_max,
                               lon.min = min_max$X_LONG_min -2 + 360,
                               lon.max = min_max$X_LONG_max + 360)

sst_df <- sst_dataframe(sst_subset) %>% 
  mutate(lon = ifelse(lon > 180, lon - 360, lon)) 

#save sst_df

write.csv(sst_df, here("data_processing", "sst_ersst", "sst_ersst_df.csv"), row.names = FALSE)


#plot

bc_boundary <- bc_bound() %>% st_transform(4326)

sst_df %>% filter(!is.na(sst), year == 1955 | year == 2012) %>% 
  ggplot() +
  geom_sf(data = bc_boundary, fill = "transparent", color = "slategray", alpha = 0.2) +
  facet_wrap(~year) +
  geom_raster(aes(x = lon, y = lat, fill = sst), alpha = 0.5) +
  #plot locations of sst data
  geom_point(aes(x = lon, y = lat), color = "slategray", alpha = 0.6) +
  scale_fill_viridis_c() +
  # geom_point(data = lighthouse_locations, aes(x = long, y = lat), color = "darkred", size = 3, alpha=0.8) +
  geom_point(data=salmon_data_location, aes(x = X_LONG, y = Y_LAT, color = "chum"), size = 2, alpha=0.2) +
  geom_point(data=pke_salmon_data_location, aes(x = X_LONG, y = Y_LAT, color = "pink-even"), size = 2, alpha=0.2) +
  geom_point(data=pko_salmon_data_location, aes(x = X_LONG, y = Y_LAT, color = "pink-odd"), size = 2, alpha=0.2) +
  # geom_text(data = lighthouse_locations, aes(x = long, y = lat, label = location), 
  #           nudge_x = -1.5, nudge_y = 0.2, size = 3) +
  # ggrepel::geom_label_repel(data = lighthouse_locations, aes(x = long, y = lat, label = location),
  #                           nudge_x = -1.5, nudge_y = 0.2, size = 3, background = "white", alpha = 0.5) +
  scale_color_manual(values = c("chum" = "#69C5C5",
                                "pink-even" = "#C76F6F",
                                "pink-odd" = "#9E70A1")) +
  scale_x_continuous(limits = c(-136, -122), breaks = seq(-136,-122,5)) +
  scale_y_continuous(limits = c(47, 58)) +
  guides(color = guide_legend(title = "Species"), override.aes = list(size = 4, alpha = 1),
         fill = guide_legend(title = "Spring SST (C)")) +
  labs(title = "Spring SST") + 
  xlab("Longitude") +
  ylab("Latitude") +
  theme_classic()+
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10))

ggsave(here("figures", "sstersst_1955_2012.png"), sst_1954_2012, width = 10, height = 5, dpi = 300)









haversine <- function(lat1, lon1, lat2, lon2) {
  ## This function computes the great circle distance between two points given
  ## their latitiude and longitude (in decimal degrees) using the haversine
  ## formula. The output is the distance between the two points in km.
  ##
  ## lat1 = latitude of first point
  ## lon1 = longitude of first point
  ## lat2 = latitude of second point
  ## lon2 = longitude of second point
  
  # Convert degrees to radians
  lat1 <- lat1 * pi / 180
  lon1 <- lon1 * pi / 180
  lat2 <- lat2 * pi / 180
  lon2 <- lon2 * pi / 180
  
  R <- 6371 # earth mean radius [km]
  delta.lon <- (lon2 - lon1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.lon/2)^2
  d <- 2 * R * asin(min(1, sqrt(a)))
  
  return(d) # distance in km
}

location_data_long_df_distinct <- sst_df %>% 
  # mutate(lon = ifelse(lon > 180, lon - 360, lon)) %>%
  filter(!is.na(sst)) %>%
  select(lat, lon) %>% 
  distinct()


salmon_data_location <- data %>% 
  select(CU,  Y_LAT, X_LONG, River) %>% 
  distinct()

distance_df <- tibble()
for(i in 1:nrow(salmon_data_location)){
  for(j in 1:nrow(location_data_long_df_distinct)){
    distance <- haversine(salmon_data_location$Y_LAT[i], salmon_data_location$X_LONG[i], 
                          location_data_long_df_distinct$lat[j], location_data_long_df_distinct$lon[j])
    distance_df <- distance_df %>% bind_rows(data.frame(CU = salmon_data_location$CU[i], 
                                                        River = salmon_data_location$River[i],
                                                        sst_ersst_lat = location_data_long_df_distinct$lat[j], 
                                                        sst_ersst_lon = location_data_long_df_distinct$lon[j],
                                                        distance = distance))
    
  }
  
}

min_distance_df <- distance_df %>% 
  group_by(CU, River) %>% 
  filter(distance == min(distance)) %>% 
  ungroup() 

sst_df_spring <- sst_df %>% 
  group_by(lat, lon ,year) %>% 
  summarise(spring_ersst = mean(sst)) %>%
  ungroup()

#check how many have NA

sst_df_spring %>% 
  filter(is.na(spring_ersst)) %>% 
  nrow()


salmon_data_distance_temp <- data %>% 
  left_join(min_distance_df %>% 
              select(CU, River, distance, sst_ersst_lat, sst_ersst_lon),
            by = c("CU" = "CU", "River" = "River")) %>% 
  left_join(sst_df_spring %>% 
              group_by(lat, lon , year) %>%
              mutate(BroodYear = year-1) %>% #sst fron year n will affect salmon whose BroodYear is n-1
              rename("sst_ersst_year" = "year"),
            by = c("BroodYear" = "BroodYear", "sst_ersst_lat" = "lat", "sst_ersst_lon" = "lon"))

#check how many rows have NA for spring_ersst

salmon_data_distance_temp %>% 
  filter(is.na(spring_ersst)) %>% 
  nrow()
#none

#save df

write.csv(salmon_data_distance_temp, here("origional-ecofish-data-models","Data","Processed","chum_SR_20_hat_yr_w_ersst.csv"), row.names = FALSE)

# look at the correlation between spring_sst and ersst
correlation = cor(salmon_data_distance_temp$spring_sst, salmon_data_distance_temp$spring_ersst)

salmon_data_distance_temp %>% 
  left_join(cu_names, by = c("CU" = "CU")) %>%
  ggplot() +
  geom_point(aes(x = spring_ersst, y = spring_sst, color = CU_name), , size = 2, alpha = 0.1) +
  labs(title = "Spring SST vs ERSST",
       x = "ERSST (C)",
       y = "COBE SST") +
  theme_classic() +
  theme(legend.position = "none")+
  scale_color_iwanthue(name = 'CU', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
                       lmin = 10, lmax = 95) +
  xlim(8, 14) +
  ylim(8, 14) +
  theme_classic()+
  #plot 1 to 1 line
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'black') +
  annotate("text", x = 12, y = 9, label = paste('Correlation:', round(correlation, 2)), size = 5)+
  theme(legend.key.height = unit(1, 'lines'),
        legend.position = 'right',
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14))+
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2), ncol = 1))

ggsave(here("figures", "spring_sst_vs_ersst.png"), width = 6, height = 6, dpi = 300)




