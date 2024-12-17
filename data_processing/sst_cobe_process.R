#read data from sst_cobe folder

# libraries

library(here)
library(tidyverse)
library(ggplot2)
library(ncdf4)
library(sf)
library(bcmaps)

# sst_cobe <- read_csv(here("data_processing", "sst_cobe", "sst_COBE2_full.csv"))

#split yr_mm to year and month
# 
# sst_cobe_long <- sst_cobe %>%
#   mutate(year = as.numeric(substr(yr_mn, 1, 4)),
#          month = as.numeric(substr(yr_mn, 6, 7))) %>%
#   #make longer
#   pivot_longer(cols = c(-yr_mn, -year, -month, -time_point), names_to = "River", values_to = "sst") %>%
#   group_by(year, River) %>%
#   filter(month %in% c(4,5,6,7)) %>%
#   summarise(mean_spring_sst = mean(sst, na.rm = TRUE))


# Look at Kingcome river by year

# sst_cobe_long %>%
#   filter(River == "KINGCOME RIVER", year >= 1954, year <= 2012) %>%
#   ggplot(aes(x = year, y = mean_spring_sst)) +
#   geom_point() +
#   geom_line(alpha = 0.5, linewidth = 1.5) +
#   scale_x_continuous(breaks = seq(1954, 2012, 5)) +
#   labs(title = "Kingcome River Spring SST",
#        x = "Year",
#        y = "SST (C)") +
#   theme_classic()

# 
# #read chum data
# 
# salmon_data <- read.csv(here("origional-ecofish-data-models","Data","Processed",
#                              "chum_SR_20_hat_yr_w_coord_w_SSTCOBE2.csv"))
# 
# #unique river names in salmon data
# 
# print(length(unique(salmon_data$River)))
# print(length(unique(salmon_data$GFE_ID)))
# 
# 
# 
# 
# #unique river names in sst data
# 
# print(length(unique(sst_cobe_long$River)))
# 
# # check if there are any 2 rivers (With different GFE_ID) with the same name
# 
# salmon_data %>%
#   group_by(River) %>%
#   summarise(n = n_distinct(GFE_ID)) %>%
#   filter(n > 1)
# 
# 
# #Check if SALMON RIVER, LAGOON CREEK has different X_LONG and Y_LAT
# salmon_data %>%
#   select(River, X_LONG, Y_LAT) %>%
#   filter(River == "LAGOON CREEK") %>%
#   group_by(X_LONG, Y_LAT) %>%
#   summarise(n = n_distinct(River)) 
# 
# #Coordinates are right but the names might be wrong
# 
# # need to get data again
# 
# sst_cobe_long %>%
#   group_by(River) %>%
#   summarise(n = n_distinct(year)) %>%
#   filter(n > 1)
# 
# 
# # get data from "http://psl.noaa.gov/thredds/dodsC/Datasets/COBE2/sst.mon.ltm.1981-2010.nc?lat[0:48:59],lon[0:200:300],time[0:1:11]"
# 
# # read data from the above link

# 
# 
# nc_file <- "http://psl.noaa.gov/thredds/dodsC/Datasets/COBE2/sst.mon.ltm.1981-2010.nc?lat[41.5:30.5],lon[5],time[3],sst[3][41.5:30.5][5]"
# 
# 
# nc <- nc_open(nc_file)
# 
# sst <- ncvar_get(nc, "sst")
# 
# year <- ncvar_get(nc, "time")
# 
# lat <- ncvar_get(nc, "lat")
# 
# lon <- ncvar_get(nc, "lon")
# 
# time_origin <- as.Date("1891-01-01")
# 
# time <- as.Date(time_origin + as.numeric(year), origin = time_origin)
# 
# 
# nc_file <- "http://psl.noaa.gov/thredds/dodsC/Datasets/COBE2/sst.mon.ltm.1981-2010.nc"
# 
# nc <- nc_open(nc_file)
# 
# sst <- ncvar_get(nc, "sst")
# 
# year <- ncvar_get(nc, "time")
# 
# lat <- ncvar_get(nc, "lat")
# 
# lon <- ncvar_get(nc, "lon")
# 
# time_origin <- as.Date("1891-01-01")
# 
# time_dates <- as.Date(year, origin = time_origin)
# 
# # Define your coordinates and years
# coordinates <- list(c(-123.5, 48.5), c(-125.5, 50.5))  # Example: (longitude, latitude)
# years <- 1954:2012
# 
# results <- list()
# for (coord in coordinates) {
#   lon_idx <- which.min(abs(lon - coord[1]))
#   lat_idx <- which.min(abs(lat - coord[2]))
#   for (year in years) {
#     start_date <- as.Date(paste0(year, "-01-01"))
#     end_date <- as.Date(paste0(year, "-12-31"))
#     time_idx <- which(time_dates >= start_date & time_dates <= end_date)
#     sst_values <- sst[lon_idx, lat_idx, time_idx]
#     monthly_avg <- tapply(sst_values, format(time_dates[time_idx], "%m"), mean, na.rm = TRUE)
#     results[[paste(coord, year)]] <- monthly_avg
#   }
# }
# 
# print(results)
# nc_close(nc)

nc_file <- nc_open(here("data_processing", "sst_cobe", "sst.mon.mean.nc"))

# Print the file summary to understand its structure
sst_data <- ncvar_get(nc_file, "sst")

# Get the dimensions of the data
lon <- ncvar_get(nc_file, "lon")
lat <- ncvar_get(nc_file, "lat")
time <- ncvar_get(nc_file, "time")


# Close the NetCDF file
nc_close(nc_file)

# time = as.POSIXct(time, origin = '1981-01-01 00:00:00', tz = 'UTC')

#time starts from 1981-01-01 00:00:00. monthly data


time <- as.Date(time, origin = '1891-01-01')

#make data frame with lat, lon, time and sst

sst_df <- expand.grid(lon = lon, lat = lat, time = time)


sst_df$sst <- as.vector(sst_data)


#subset sst data - lat 48-59, lon -122 to -123, time - 1954 to 2012


#plot the data using sf



# Define the coordinates of the bounding box

bc_boundary <- bc_bound() %>% st_transform(4326)



#plot sst for year 1954, all lat and lon

sst_1954 <- sst_df %>%
  filter(time == as.Date("1954-01-01"), lon > 220, lon < 240, lat > 48 , lat < 60) %>% 
  #convert lon to -180 to 180 (from 0 to 360)
  mutate(lon = ifelse(lon > 180, lon - 360, lon) %>% as.numeric()) %>% 
  ggplot() +
  geom_sf(data = bc_boundary) +
  geom_raster(aes(x = lon, y = lat, fill = sst), alpha = 0.5) +
  xlim(-140, -120) +
  ylim(48, 60) +
  scale_fill_viridis_c() +
  labs(title = "SST 1954",
       fill = "SST (C)") +
  theme_minimal()




sst_df %>%
  filter(time == as.Date("2012-01-01"), lon > 220, lon < 240, lat > 48 , lat < 60) %>% 
  #convert lon to -180 to 180 (from 0 to 360)
  mutate(lon = ifelse(lon > 180, lon - 360, lon) %>% as.numeric()) %>% filter(!is.na(sst)) %>% 
  ggplot() +
  geom_sf(data = bc_boundary) +
  geom_raster(aes(x = lon, y = lat, fill = sst), alpha = 0.5) +
  #plot locations of sst data
  geom_point(aes(x = lon, y = lat), color = "slategray") +
  xlim(-140, -120) +
  ylim(48, 60) +
  scale_fill_viridis_c() +
  labs(title = "SST 2012",
       fill = "SST (C)") +
  theme_minimal()


# plot salmon watersheds and sst data
pke_cu_data <- read_csv(here("data_processing","CU_data","PKE_CU_Sites_En.csv"))

pko_cu_data <- read_csv(here("data_processing","CU_data","PKO_CU_Sites_En.csv"))



pke_salmon_data <- read.csv(here("origional-ecofish-data-models","Data","Processed",
                                 "pke_SR_10_hat_yr_reduced_VRI90.csv"))

pko_salmon_data <- read.csv(here("origional-ecofish-data-models","Data","Processed",
                                 "PKO_SR_10_hat_yr_reduced_VRI90.csv"))


pke_cu_data <- pke_cu_data %>% 
  select(CU_NAME, FULL_CU_IN, SITE_NAME, Y_LAT, X_LONG, GFE_ID) %>% 
  distinct() 


pko_cu_data <- pko_cu_data %>%
  select(CU_NAME, FULL_CU_IN, SITE_NAME, Y_LAT, X_LONG, GFE_ID) %>% 
  distinct()

pko_salmon_data <- pko_salmon_data %>%
  left_join(pko_cu_data, by = c("CU" = "FULL_CU_IN", "GFE_ID" = "GFE_ID"))



pke_salmon_data <- pke_salmon_data %>% 
  left_join(pke_cu_data, by = c("CU" = "FULL_CU_IN", "GFE_ID" = "GFE_ID"))

pke_salmon_data_location <- pke_salmon_data %>% 
  select(CU,  Y_LAT, X_LONG, River, GFE_ID) %>% 
  distinct()

pko_salmon_data_location <- pko_salmon_data %>%
  select(CU,  Y_LAT, X_LONG, River, GFE_ID) %>% 
  distinct()

salmon_data <- read.csv(here("origional-ecofish-data-models","Data","Processed",
                             "chum_SR_20_hat_yr_w_coord_w_SSTCOBE2.csv"))


salmon_data_location <- salmon_data %>% 
  select(CU,  Y_LAT, X_LONG, River) %>% 
  distinct()


# make dataframe with average spring sst for each year

sst_spring <- sst_df %>% 
  filter(lon > 220, lon < 240, lat > 48 , lat < 60) %>% 
  mutate(lon = ifelse(lon > 180, lon - 360, lon) %>% as.numeric()) %>% 
  mutate(year = as.numeric(format(time, "%Y")), month = as.numeric(format(time, "%m"))) %>%
  filter(month %in% c(4,5,6,7), year >= 1954, year <= 2014) %>%
  group_by(year, lat, lon) %>% 
  summarize(spring_sst = mean(sst)) %>% 
  ungroup()
  

# plot the time series of average spring sst for each year

sst_spring %>% 
  group_by(year) %>% 
  summarize(mean_spatial_spring_sst = mean(spring_sst, na.rm = TRUE)) %>%
  ggplot(aes(x = year, y = mean_spatial_spring_sst)) +
  geom_point(color = "slategray") +
  geom_line(alpha = 0.5, linewidth = 1.5, color = "salmon") +
  ylim(8, 12) +
  scale_x_continuous(breaks = seq(1954, 2012, 5)) +
  labs(title = "Mean Spring SST",
       x = "Year",
       y = "SST (C)") +
  theme(legend.position = "none",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10)) +
  theme_classic()






sst_2012 <- sst_spring %>% filter(!is.na(spring_sst), year == 2012) %>% 
  mutate(lon = ifelse(lon > 180, lon - 360, lon) %>% as.numeric()) %>% filter(!is.na(spring_sst)) %>% 
  ggplot() +
  
  geom_raster(aes(x = lon, y = lat, fill = spring_sst), alpha = 0.5) +
  #plot locations of sst data
  geom_point(aes(x = lon, y = lat), color = "slategray", alpha = 0.6) +
  scale_fill_viridis_c() +
  geom_sf(data = bc_boundary, fill = "transparent", color = "slategray", alpha = 0.2) +
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
  scale_x_continuous(limits = c(-136, -122)) +
  scale_y_continuous(limits = c(48, 56)) +
  guides(color = guide_legend(title = "Species"), override.aes = list(size = 4, alpha = 1),
         fill = guide_legend(title = "Spring SST (C)")) +
  labs(title = "Spring SST 2012") + 
  xlab("Longitude") +
  ylab("Latitude") +
  theme_classic()+
  theme(legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10))



sst_1954_2012 <- sst_spring %>% filter(!is.na(spring_sst), year == 1954 | year == 2012) %>% filter(!is.na(spring_sst)) %>% 
  ggplot() +
  
  geom_raster(aes(x = lon, y = lat, fill = spring_sst), alpha = 0.5) +
  #plot locations of sst data
  geom_point(aes(x = lon, y = lat), color = "slategray", alpha = 0.6) +
  scale_fill_viridis_c() +
  geom_sf(data = bc_boundary, fill = "transparent", color = "slategray", alpha = 0.2) +
  # geom_point(data = lighthouse_locations, aes(x = long, y = lat), color = "darkred", size = 3, alpha=0.8) +
  geom_point(data=salmon_data_location, aes(x = X_LONG, y = Y_LAT, color = "chum"), size = 2, alpha=0.2) +
  geom_point(data=pke_salmon_data_location, aes(x = X_LONG, y = Y_LAT, color = "pink-even"), size = 2, alpha=0.2) +
  geom_point(data=pko_salmon_data_location, aes(x = X_LONG, y = Y_LAT, color = "pink-odd"), size = 2, alpha=0.2) +
  facet_wrap(~year) +
  # geom_text(data = lighthouse_locations, aes(x = long, y = lat, label = location), 
  #           nudge_x = -1.5, nudge_y = 0.2, size = 3) +
  # ggrepel::geom_label_repel(data = lighthouse_locations, aes(x = long, y = lat, label = location),
  #                           nudge_x = -1.5, nudge_y = 0.2, size = 3, background = "white", alpha = 0.5) +
  scale_color_manual(values = c("chum" = "#69C5C5",
                                "pink-even" = "#C76F6F",
                                "pink-odd" = "#9E70A1")) +
  scale_x_continuous(limits = c(-136, -122), breaks = seq(-136,-122,5)) +
  scale_y_continuous(limits = c(48, 56)) +
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


sst_1954_2012

ggsave(here("figures", "sstcobe2_1954_2012.png"), sst_1954_2012, width = 10, height = 5, dpi = 300)


# for each combination of salmon watershed location and lat, lon of the sst data,
# find the closest sst data point and get the sst value

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



location_data_long_df_distinct <- sst_spring %>% 
  filter(!is.na(spring_sst)) %>%
  select(lat, lon) %>% 
  distinct()

distance_df <- tibble()
for(i in 1:nrow(salmon_data_location)){
  for(j in 1:nrow(location_data_long_df_distinct)){
    distance <- haversine(salmon_data_location$Y_LAT[i], salmon_data_location$X_LONG[i], 
                          location_data_long_df_distinct$lat[j], location_data_long_df_distinct$lon[j])
    distance_df <- distance_df %>% bind_rows(data.frame(CU = salmon_data_location$CU[i], 
                                                        River = salmon_data_location$River[i],
                                                        sst_cobe_lat = location_data_long_df_distinct$lat[j], 
                                                        sst_cobe_lon = location_data_long_df_distinct$lon[j],
                                                        distance = distance))
    
  }
  
}


min_distance_df <- distance_df %>% 
  group_by(CU, River) %>% 
  filter(distance == min(distance)) %>% 
  ungroup() 


salmon_data_distance_temp <- salmon_data %>% 
  left_join(min_distance_df %>% 
              select(CU, River, distance, sst_cobe_lat, sst_cobe_lon),
            by = c("CU" = "CU", "River" = "River")) %>% 
  left_join(sst_spring %>% 
              group_by(lat, lon , year) %>%
              mutate(BroodYear = year-1) %>% #sst fron year n will affect salmon whose BroodYear is n-1
              rename("sst_cobe_year" = "year"),
            by = c("BroodYear" = "BroodYear", "sst_cobe_lat" = "lat", "sst_cobe_lon" = "lon"))

salmon_data_distance_temp %>% 
  write_csv(here("origional-ecofish-data-models","Data","Processed",
                 "chum_SR_20_hat_yr_w_sstcobe_new.csv"))

#check if any spring_sst is NA

salmon_data_distance_temp %>% 
  filter(is.na(spring_sst))


#do same for pinks even and pink odd



pke_distance_df <- tibble()

for(i in 1:nrow(pke_salmon_data_location)){
  for(j in 1:nrow(location_data_long_df_distinct)){
    distance <- haversine(pke_salmon_data_location$Y_LAT[i], pke_salmon_data_location$X_LONG[i], 
                          location_data_long_df_distinct$lat[j], location_data_long_df_distinct$lon[j])
    pke_distance_df <- pke_distance_df %>% bind_rows(data.frame(CU = pke_salmon_data_location$CU[i], 
                                                                River = pke_salmon_data_location$River[i],
                                                                sst_cobe_lat = location_data_long_df_distinct$lat[j], 
                                                                sst_cobe_lon = location_data_long_df_distinct$lon[j], 
                                                                distance = distance))
    
  }
  
}

pko_distance_df <- tibble()


for(i in 1:nrow(pko_salmon_data_location)){
  for(j in 1:nrow(location_data_long_df_distinct)){
    distance <- haversine(pko_salmon_data_location$Y_LAT[i], pko_salmon_data_location$X_LONG[i], 
                          location_data_long_df_distinct$lat[j], location_data_long_df_distinct$lon[j])
    pko_distance_df <- pko_distance_df %>% bind_rows(data.frame(CU = pko_salmon_data_location$CU[i], 
                                                                River = pko_salmon_data_location$River[i],
                                                                sst_cobe_lat = location_data_long_df_distinct$lat[j],
                                                                sst_cobe_lon = location_data_long_df_distinct$lon[j],
                                                                distance = distance))
    
  }
  
}


pke_min_distance_df <- pke_distance_df %>% 
  group_by(CU, River) %>% 
  # filter(lighthouse_location != "Nootka_Point", lighthouse_location != "Egg_Island") %>%
  filter(distance == min(distance)) %>% 
  ungroup()


pko_min_distance_df <- pko_distance_df %>%
  group_by(CU, River) %>% 
  # filter(lighthouse_location != "Nootka_Point", lighthouse_location != "Egg_Island") %>%
  filter(distance == min(distance)) %>% 
  ungroup()


pke_salmon_data_distance_temp <- pke_salmon_data %>%
  left_join(pke_min_distance_df %>% 
              select(CU, River, distance, sst_cobe_lat, sst_cobe_lon),
            by = c("CU" = "CU", "River" = "River")) %>% 
  left_join(sst_spring %>% 
              group_by(lat, lon , year) %>%
              mutate(BroodYear = year-1) %>% #sst fron year n will affect salmon whose BroodYear is n-1
              rename("sst_cobe_year" = "year"),
            by = c("BroodYear" = "BroodYear", "sst_cobe_lat" = "lat", "sst_cobe_lon" = "lon"))

pko_salmon_data_distance_temp <- pko_salmon_data %>%
  left_join(pko_min_distance_df %>% 
              select(CU, River, distance, sst_cobe_lat, sst_cobe_lon),
            by = c("CU" = "CU", "River" = "River")) %>% 
  left_join(sst_spring %>% 
              group_by(lat, lon , year) %>%
              mutate(BroodYear = year-1) %>% #sst fron year n will affect salmon whose BroodYear is n-1
              rename("sst_cobe_year" = "year"),
            by = c("BroodYear" = "BroodYear", "sst_cobe_lat" = "lat", "sst_cobe_lon" = "lon"))


pke_salmon_data_distance_temp %>%
  write_csv(here("origional-ecofish-data-models","Data","Processed",
                 "pke_SR_10_hat_yr_w_sstcobe_new.csv"))

pko_salmon_data_distance_temp %>%
  write_csv(here("origional-ecofish-data-models","Data","Processed",
                 "pko_SR_10_hat_yr_w_sstcobe_new.csv"))


pke_salmon_data_distance_temp %>% 
  filter(is.na(spring_sst))


pko_salmon_data_distance_temp %>% 
  filter(is.na(spring_sst)) %>% 
  select(BroodYear, sst_cobe_year)

pko_salmon_data_distance_temp %>% 
  select(BroodYear, sst_cobe_year) %>% 
  distinct()




