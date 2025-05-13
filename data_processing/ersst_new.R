# write new csv data file with river coordinates and with ersst
# read file without duplicated data

data <- read.csv(here("origional-ecofish-data-models","Data","Processed","chum_SR_20_hat_yr.csv"))


glimpse(data)


chum_cu <- read.csv(here("data_processing","CU_data","CM_CU_SITES_En.csv"))

glimpse(chum_cu)

chum_cu_subset <- chum_cu %>% 
  select(FULL_CU_IN, X_LONG, Y_LAT, WS_CDE_20K,
         WS_CDE_50K, GFE_ID, CU_NAME)

glimpse(chum_cu_subset)

#check if all GFE ID in data are in chum_cu_subset

data %>% 
  pull(GFE_ID) %>% 
  {all (. %in% chum_cu_subset$GFE_ID)} #TRUE %in% chum_cu_subset$GFE_ID)} #TRUE}
  
# yes

data %>% 
  pull(WATERSHED_CDE) %>% 
  {all (. %in% chum_cu_subset$WS_CDE_50K)} #TRUE %in% chum_cu_subset$GFE_ID)} #TRUE}

#no 

chum_cu_subset <- chum_cu %>% 
  select(FULL_CU_IN, X_LONG, Y_LAT, GFE_ID, CU_NAME) %>% 
  unique()

#left join with lat long data
chum_data_w_coord <- data %>% 
  left_join(chum_cu_subset, by = c("CU" = "FULL_CU_IN", "GFE_ID" = "GFE_ID"),
       t(-CU, -GFE_ID) %>%     relationship = "many-to-one")

chum_data_w_coord %>% 
  select(River_GFE_ID) %>% 
  unique() %>% 
  nrow()

salmon_data_location <- chum_data_w_coord %>% 
  select(CU,  Y_LAT, X_LONG, River, GFE_ID) %>% 
  distinct()


# read sst df

sst_df <- read.csv(here("data_processing","sst_ersst","sst_ersst_df.csv"))

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

#distinct location of sst data

location_data_long_df_distinct <- sst_df %>% 
  # mutate(lon = ifelse(lon > 180, lon - 360, lon)) %>%
  filter(!is.na(sst)) %>%
  select(lat, lon) %>% 
  distinct()

distance_df <- tibble()

# looking at the distances between sst data points and salmon locations 
for(i in 1:nrow(salmon_data_location)){
  for(j in 1:nrow(location_data_long_df_distinct)){
    distance <- haversine(salmon_data_location$Y_LAT[i], salmon_data_location$X_LONG[i], 
                          location_data_long_df_distinct$lat[j], location_data_long_df_distinct$lon[j])
    distance_df <- distance_df %>% bind_rows(data.frame(CU = salmon_data_location$CU[i], 
                                                        River = salmon_data_location$River[i],
                                                        GFE_ID = salmon_data_location$GFE_ID[i],
                                                        sst_ersst_lat = location_data_long_df_distinct$lat[j], 
                                                        sst_ersst_lon = location_data_long_df_distinct$lon[j],
                                                        distance = distance))
    
  }
  
}

# looking the minimum of those distances for each river

min_distance_df <- distance_df %>% 
  group_by(CU, River, GFE_ID) %>% 
  filter(distance == min(distance)) %>% 
  ungroup() 

sst_df_spring <- sst_df %>% 
  group_by(lat, lon ,year) %>% 
  summarise(spring_ersst = mean(sst)) %>%
  ungroup()

salmon_data_distance_temp <- chum_data_w_coord %>% 
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

glimpse(salmon_data_distance_temp)

# look at all combinations of River_GFE_ID and BroodYears

salmon_data_distance_temp %>% 
  select(River_GFE_ID, BroodYear) %>% 
  group_by(River_GFE_ID, BroodYear) %>% 
  summarise(num_years = n()) %>% 
  # filter(num_years > 1) %>% 
  View()

salmon_data_distance_temp %>% 
  select(River_GFE_ID, BroodYear) %>% 
  group_by(River_GFE_ID) %>% 
  summarise(num_years = n()) %>% 
  # filter(num_years > 1) %>% 
  View()

url = "https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/index/ersst.v5.pdo.dat"

pdo = read.table(url, header = TRUE, skip = 1) #, col.names = c("year", "month", "pdo"))

# convert to long format

pdo_long <- pdo %>%
  pivot_longer(cols = -c(Year), names_to = "month", values_to = "pdo")

pdo_long_annual <- pdo_long %>%
  group_by(Year) %>%
  summarise(pdo = mean(pdo))


salmon_data_distance_temp_pdo <- salmon_data_distance_temp %>% 
  left_join(pdo_long_annual, by = c('sst_ersst_year' = 'Year')) %>% 
  mutate(pdo.std = (pdo - mean(pdo))/sd(pdo))

salmon_data_distance_temp_pdo %>% 
  select(River_GFE_ID, BroodYear) %>% 
  group_by(River_GFE_ID, BroodYear) %>% 
  summarise(num_years = n()) %>% 
  filter(num_years > 1) %>%
  View()

#save df

write.csv(salmon_data_distance_temp_pdo, here("origional-ecofish-data-models",
                                          "Data","Processed",
                                          "chum_SR_20_hat_yr_w_ocean_covariates.csv"), row.names = FALSE)

