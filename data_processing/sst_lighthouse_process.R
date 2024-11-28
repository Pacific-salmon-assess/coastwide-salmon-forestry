# Goal - to read data from sst_lighthouse folder
# read folder name within sst_lighthouse folder


# load libraries

library(here)
library(tidyverse)
library(ggplot2)
library(bcmaps)
library(hues)


# Read data from sst_lighthouse folder

sst_lighthouse_folder <- here("data_processing","sst_lighthouse")

# read folder name within sst_lighthouse folder as location names

location_names <- list.files(sst_lighthouse_folder)

#read light house locations file

lighthouse_locations <- read_csv(here("data_processing","lighthouse_locations.csv"),skip = 2)
  # colnames(c("location","status","extent","lat","long", "comments"))

colnames(lighthouse_locations) <- c("location","status","extent","lat","long", "comments")

lighthouse_locations <- lighthouse_locations %>% 
  select(-c(comments)) %>% 
  filter(status == "ACTIVE")

#plot location map with bc boundary
bc_boundary <- bc_bound() %>% st_transform(4326)
ggplot() +
  geom_sf(data = bc_boundary, fill = "transparent", color = "slategray", alpha = 0.2) +
  geom_point(data = lighthouse_locations, aes(x = long, y = lat), color = "darkred", size = 2) +
  # geom_text(data = lighthouse_locations, aes(x = long, y = lat, label = location), 
  #           nudge_x = -1.5, nudge_y = 0.2, size = 3) +
  ggrepel::geom_label_repel(data = lighthouse_locations, aes(x = long, y = lat, label = location), 
                            nudge_x = -1.5, nudge_y = 0.2, size = 3) +
  scale_x_continuous(limits = c(-136, -122)) +
  scale_y_continuous(limits = c(48, 56)) +
  theme_minimal()

# read data from each location

# read file that ends with 2024.csv and has Monthly_Sea_Surface_Temperatures in the file name
location_data_long_df <- tibble()
for(i in 1:length(location_names)){
  location_name <- location_names[i]
  location_data_folder <- here(sst_lighthouse_folder,location_name)
  location_data_files <- list.files(location_data_folder)
  if(i == 3){
    location_data_files <- location_data_files[grepl("2023.csv",location_data_files)]
  }else{
    location_data_files <- location_data_files[grepl("2024.csv",location_data_files)]
  }
  
  location_data_files <- location_data_files[grepl("Monthly_Sea_Surface_Temperatures",location_data_files)]
  location_data <- read_csv(here(location_data_folder,location_data_files), skip = 1, na = c(999.99,"NA",999.90))
  location_data_long <- location_data %>% 
    #non capitalize column names
    rename_all(tolower) %>%
    pivot_longer(cols = -c("year"), names_to = "month", values_to = "temperature") %>% 
    mutate(location = location_name, location_name_subset = tolower(str_extract(location_name, "^[A-Za-z]+"))) %>% 
    #join with lighthouse_locations lat and long, join on first 2 words in location after making all but first letter lowercase
    left_join(lighthouse_locations %>% 
                mutate(location_name_subset = tolower(str_extract(location, "^[A-Za-z]+"))) %>% 
                select(location_name_subset, lat, long),
              by = "location_name_subset")
  location_data_long_df <- location_data_long_df %>% bind_rows(location_data_long)
  
}

#plot the time series of sea surface temperature for each location, facet wrap by location

location_data_long_df %>% 
  group_by(location,year) %>%
  filter(month %in% c("apr","may","jun","jul")) %>%
  summarize(spring_temperature = mean(temperature)) %>%
  ggplot(aes(x = year, y = spring_temperature, color ="darkred", alpha = 0.5)) +
  scale_x_continuous(limit = c(1954,2012), breaks = seq(1954,2012,20)) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~location) +
  theme_classic()+
  theme(legend.position = "none") 

#there are some missing values in sst

# need to read salmon dataset with coord for each CU, calculate the distance between each CU and each lighthouse location
# if the closest point is nootka_point, then look for the next closest point
# join by year and location

# read salmon dataset

salmon_data <- read.csv(here("origional-ecofish-data-models","Data","Processed",
                             "chum_SR_20_hat_yr_w_coord_w_SSTCOBE2.csv"))

#make function to calculate haversine distance between 2 points

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


#calculate distance between each combination of CU and lighthouse location

salmon_data_location <- salmon_data %>% 
  select(CU,  Y_LAT, X_LONG, River) %>% 
  distinct()

location_data_long_df_distinct <- location_data_long_df %>% 
  select(location, lat, long) %>% 
  distinct()

distance_df <- tibble()
for(i in 1:nrow(salmon_data_location)){
  for(j in 1:nrow(location_data_long_df_distinct)){
    distance <- haversine(salmon_data_location$Y_LAT[i], salmon_data_location$X_LONG[i], 
                          location_data_long_df_distinct$lat[j], location_data_long_df_distinct$long[j])
    distance_df <- distance_df %>% bind_rows(data.frame(CU = salmon_data_location$CU[i], 
                                           River = salmon_data_location$River[i],
                                           lighthouse_location = location_data_long_df_distinct$location[j], 
                                           distance = distance))
    
  }

}


#find the closest lighthouse location for each CU
# if it is Nootka, then find the next closest location

min_distance_df <- distance_df %>% 
  group_by(CU, River) %>% 
  filter(lighthouse_location != "Nootka_Point", lighthouse_location != "Egg_Island") %>%
  filter(distance == min(distance)) %>% 
  ungroup() 

median(min_distance_df$distance)
max(min_distance_df$distance)

ggplot() +
  geom_sf(data = bc_boundary, fill = "transparent", color = "slategray", alpha = 0.2) +
  geom_point(data = lighthouse_locations, aes(x = long, y = lat), color = "darkred", size = 2) +
  geom_point(data=salmon_data_location, aes(x = X_LONG, y = Y_LAT), color = "salmon", size = 1, alpha=0.2) +
  # geom_text(data = lighthouse_locations, aes(x = long, y = lat, label = location), 
  #           nudge_x = -1.5, nudge_y = 0.2, size = 3) +
  ggrepel::geom_label_repel(data = lighthouse_locations, aes(x = long, y = lat, label = location), 
                            nudge_x = -1.5, nudge_y = 0.2, size = 3) +
  scale_x_continuous(limits = c(-136, -122)) +
  scale_y_continuous(limits = c(48, 56)) +
  theme_minimal()

#join the spring_temperature data with salmon data,
#join by year, lighthouse location

# min_distance_temp <- min_distance_df %>% 
#   left_join(location_data_long_df %>% 
#               group_by(location,year) %>%
#               filter(month %in% c("apr","may","jun","jul")) %>%
#               summarize(spring_temperature = mean(temperature)),
#             by = c("lighthouse_location" = "location", "year")) %>% 
#   select(-c(lat, long))

salmon_data_distance_temp <- salmon_data %>% 
  left_join(min_distance_df %>% 
              select(CU, River, lighthouse_location, distance),
            by = c("CU" = "CU", "River" = "River")) %>% 
  left_join(location_data_long_df %>% 
              group_by(location,year) %>%
              filter(month %in% c("apr","may","jun","jul")) %>%
              summarize(spring_lighthouse_temperature = mean(temperature)) %>% 
              mutate(BroodYear = year-1) %>% #sst fron year n will affect salmon whose BroodYear is n-1
              rename("lighthouse_temp_year" = "year"),
            by = c("lighthouse_location" = "location", "BroodYear" = "BroodYear"))


#plot the relationship between spring temperature and salmon recruitment

salmon_data_distance_temp %>% 
  ggplot(aes(x = spring_lighthouse_temperature, y = ln_RS, color = CU)) +
  scale_color_iwanthue() +
  geom_point(alpha = 0.2) +
  # facet_wrap(~CU) +
  theme_minimal()


#save the data

salmon_data_distance_temp %>% 
  write_csv(here("origional-ecofish-data-models","Data","Processed",
                 "chum_SR_20_hat_yr_w_lighthouse_temp.csv"))


# Do the same for Pink data

# first get locations of CU and River

pke_cu_data <- read_csv(here("data_processing","CU_data","PKE_CU_Sites_En.csv"))

pko_cu_data <- read_csv(here("data_processing","CU_data","PKO_CU_Sites_En.csv"))


#get salmon data

pke_salmon_data <- read.csv(here("origional-ecofish-data-models","Data","Processed",
                             "pke_SR_10_hat_yr_reduced_VRI90.csv"))

pko_salmon_data <- read.csv(here("origional-ecofish-data-models","Data","Processed",
                             "PKO_SR_10_hat_yr_reduced_VRI90.csv"))


pke_cu_data <- pke_cu_data %>% 
  select(CU_NAME, FULL_CU_IN, SITE_NAME, Y_LAT, X_LONG, GFE_ID) %>% 
  distinct() 

#check if there are multiple rows for each combination of Y_LAT and X_LONG
pke_cu_data %>% 
  group_by(Y_LAT, X_LONG) %>% 
  summarize(n = n(), CU_NAME = unique(CU_NAME), SITE_NAME = unique(SITE_NAME), GFE_ID = unique(GFE_ID)) %>%
  filter(n > 1)



pke_salmon_data <- pke_salmon_data %>% 
  left_join(pke_cu_data, by = c("CU" = "FULL_CU_IN", "GFE_ID" = "GFE_ID"))

pko_cu_data <- pko_cu_data %>%
  select(CU_NAME, FULL_CU_IN, SITE_NAME, Y_LAT, X_LONG, GFE_ID) %>% 
  distinct()

pko_cu_data %>% 
  group_by(Y_LAT, X_LONG) %>% 
  summarize(n = n(), CU_NAME = unique(CU_NAME), SITE_NAME = unique(SITE_NAME), GFE_ID = unique(GFE_ID)) %>%
  filter(n > 1)


pko_salmon_data <- pko_salmon_data %>%
  left_join(pko_cu_data, by = c("CU" = "FULL_CU_IN", "GFE_ID" = "GFE_ID"))

#calculate distance between each combination of CU and lighthouse location

pke_salmon_data_location <- pke_salmon_data %>% 
  select(CU,  Y_LAT, X_LONG, River, GFE_ID) %>% 
  distinct()

pko_salmon_data_location <- pko_salmon_data %>%
  select(CU,  Y_LAT, X_LONG, River, GFE_ID) %>% 
  distinct()


pke_distance_df <- tibble()

for(i in 1:nrow(pke_salmon_data_location)){
  for(j in 1:nrow(location_data_long_df_distinct)){
    distance <- haversine(pke_salmon_data_location$Y_LAT[i], pke_salmon_data_location$X_LONG[i], 
                          location_data_long_df_distinct$lat[j], location_data_long_df_distinct$long[j])
    pke_distance_df <- pke_distance_df %>% bind_rows(data.frame(CU = pke_salmon_data_location$CU[i], 
                                           River = pke_salmon_data_location$River[i],
                                           lighthouse_location = location_data_long_df_distinct$location[j], 
                                           distance = distance))
    
  }

}

pko_distance_df <- tibble()


for(i in 1:nrow(pko_salmon_data_location)){
  for(j in 1:nrow(location_data_long_df_distinct)){
    distance <- haversine(pko_salmon_data_location$Y_LAT[i], pko_salmon_data_location$X_LONG[i], 
                          location_data_long_df_distinct$lat[j], location_data_long_df_distinct$long[j])
    pko_distance_df <- pko_distance_df %>% bind_rows(data.frame(CU = pko_salmon_data_location$CU[i], 
                                           River = pko_salmon_data_location$River[i],
                                           lighthouse_location = location_data_long_df_distinct$location[j], 
                                           distance = distance))
    
  }

}


pke_min_distance_df <- pke_distance_df %>% 
  group_by(CU, River) %>% 
  filter(lighthouse_location != "Nootka_Point", lighthouse_location != "Egg_Island") %>%
  filter(distance == min(distance)) %>% 
  ungroup()


pko_min_distance_df <- pko_distance_df %>%
  group_by(CU, River) %>% 
  filter(lighthouse_location != "Nootka_Point", lighthouse_location != "Egg_Island") %>%
  filter(distance == min(distance)) %>% 
  ungroup()


pke_salmon_data_distance_temp <- pke_salmon_data %>%
  left_join(pke_min_distance_df %>% 
              select(CU, River, lighthouse_location, distance),
            by = c("CU" = "CU", "River" = "River")) %>% 
  left_join(location_data_long_df %>% 
              group_by(location,year) %>%
              filter(month %in% c("apr","may","jun","jul")) %>%
              summarize(spring_lighthouse_temperature = mean(temperature)) %>% 
              mutate(BroodYear = year-1) %>% #sst fron year n will affect salmon whose BroodYear is n-1
              rename("lighthouse_temp_year" = "year"),
            by = c("lighthouse_location" = "location", "BroodYear" = "BroodYear"))

pko_salmon_data_distance_temp <- pko_salmon_data %>%
  left_join(pko_min_distance_df %>% 
              select(CU, River, lighthouse_location, distance),
            by = c("CU" = "CU", "River" = "River")) %>% 
  left_join(location_data_long_df %>% 
              group_by(location,year) %>%
              filter(month %in% c("apr","may","jun","jul")) %>%
              summarize(spring_lighthouse_temperature = mean(temperature)) %>% 
              mutate(BroodYear = year-1) %>% #sst fron year n will affect salmon whose BroodYear is n-1
              rename("lighthouse_temp_year" = "year"),
            by = c("lighthouse_location" = "location", "BroodYear" = "BroodYear"))


pke_salmon_data_distance_temp %>%
  write_csv(here("origional-ecofish-data-models","Data","Processed",
                 "pke_SR_10_hat_yr_w_lighthouse_temp.csv"))

pko_salmon_data_distance_temp %>%
  write_csv(here("origional-ecofish-data-models","Data","Processed",
                 "pko_SR_10_hat_yr_w_lighthouse_temp.csv"))



ggplot() +
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
  guides(color = guide_legend(title = "Species"), override.aes = list(size = 4, alpha = 1))+
  xlab("Longitude") +
  ylab("Latitude") +
  theme_classic()+
  theme(legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10))

ggsave(here("figures","chum_pink_watersheds_lighthouse_location_map.png"), width = 10, height = 10, dpi = 300)

ggsave(here("figures","chum_pink_watersheds_location_map.png"), width = 10, height = 10, dpi = 300)


#Check all rows in pko_salmon_data_distance_temp and pke_salmon_data_distance_temp where spring_lighthouse_temperature is NA

pko_na <- pko_salmon_data_distance_temp %>% 
  filter(is.na(spring_lighthouse_temperature)) #makes sense


location_data_long_df %>% 
  group_by(location,year) %>%
  filter(month %in% c("apr","may","jun","jul")) %>%
  summarize(spring_temperature = mean(temperature)) %>%
  ggplot(aes(x = year, y = spring_temperature, color ="darkred", alpha = 0.5)) +
  scale_x_continuous(limit = c(1954,2014), breaks = seq(1960,2000,20)) +
  geom_line(linewidth = 1.5) +
  xlab("Year") +
  ylab("Spring Sea Surface Temperature") +
  facet_wrap(~location) +
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10)) 

ggsave(here("figures","lighthouse_location_sst_time_series.png"), width = 10, height = 10, dpi = 300)



spring_temp <- location_data_long_df %>% 
  group_by(location,year) %>%
  filter(month %in% c("apr","may","jun","jul")) %>%
  summarize(spring_lighthouse_temperature = mean(temperature)) %>% 
  mutate(BroodYear = year-1) %>% #sst fron year n will affect salmon whose BroodYear is n-1
  rename("lighthouse_temp_year" = "year")
