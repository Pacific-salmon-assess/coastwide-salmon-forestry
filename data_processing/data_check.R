# read latest spawner recruit dataset and check for errors


library(ggplot2)
library(tidyverse)
library(here)


#read dataset

data1 <- read.csv(here("origional-ecofish-data-models","Data",
                       "Processed",
                       "chum_SR_20_hat_yr_w_ersst_pdo.csv"))


data1 %>%
  select(River, BroodYear) %>% 
  group_by(River) %>% 
  summarise(num_years = n()) %>% 
  filter(num_years > 58) %>% 
  View()

#look at all river, brood year combinations that are repeated

data1 %>%
  # select(River, BroodYear) %>% 
  group_by(River, BroodYear) %>% 
  summarise(num_years = n()) %>% 
  filter(num_years > 1) %>% 
  View()


data2 <- read.csv(here("origional-ecofish-data-models","Data",
                       "Processed",
                       "chum_SR_20_hat_yr_w_ersst.csv"))

data2 %>%
  # select(River, BroodYear) %>% 
  group_by(River, BroodYear) %>% 
  summarise(num_years = n()) %>% 
  filter(num_years > 1) %>% 
  View()

data3 <- read.csv(here("origional-ecofish-data-models","Data",
                       "Processed",
                       "chum_SR_20_hat_yr_reduced_VRI90.csv"))
nrow(data3)

data3 %>%
  # select(River, BroodYear) %>% 
  group_by(River, BroodYear) %>% 
  summarise(num_years = n()) %>% 
  filter(num_years > 1) %>% 
  View()



data4 <- read.csv(here("origional-ecofish-data-models","Data",
                       "Processed",
                       "chum_SR_20_hat_yr.csv"))
nrow(data4)


data4 %>%
  # select(River, BroodYear) %>% 
  group_by(River, BroodYear) %>% 
  summarise(num_years = n()) %>% 
  filter(num_years > 1) %>% 
  View()

#check number of rivers in data 4

data4 %>%
  select(River) %>% 
  unique()

data3 %>%
  select(River) %>% 
  unique()

data5 <- read.csv(here("origional-ecofish-data-models","Data","Processed","chum_SR_20_hat_yr_w_sstcobe_new.csv"))

data5 %>%
  # select(River, BroodYear) %>% 
  group_by(River, BroodYear) %>% 
  summarise(num_years = n()) %>% 
  filter(num_years > 1) %>% 
  View()


data6 <- read.csv(here("origional-ecofish-data-models","Data","Processed","chum_SR_20_hat_yr_w_sstcobe_new.csv"))

data6 %>%
  # select(River, BroodYear) %>% 
  group_by(River, BroodYear) %>% 
  summarise(num_years = n()) %>% 
  filter(num_years > 1) %>% 
  View()

data7 <- read.csv(here("origional-ecofish-data-models","Data","Processed",
                       "chum_SR_20_hat_yr_w_coord_w_SSTCOBE2.csv"))


data7 %>%
  # select(River, BroodYear) %>% 
  group_by(River, BroodYear) %>% 
  summarise(num_years = n()) %>% 
  filter(num_years > 1) %>% 
  View()




#check pink data

data8 <- read.csv(here("origional-ecofish-data-models","Data","Processed",
                       "pko_SR_10_hat_yr_w_ersst.csv"))

data8 %>%
  # select(River, BroodYear) %>% 
  group_by(River, BroodYear) %>% 
  summarise(num_years = n()) %>% 
  filter(num_years > 1) %>% 
  View()


data9 <- read.csv(here("origional-ecofish-data-models","Data","Processed",
                       "pko_SR_10_hat_yr_reduced_VRI90.csv"))

data9 %>%
  # select(River, BroodYear) %>% 
  group_by(River, BroodYear) %>% 
  summarise(num_years = n()) %>% 
  filter(num_years > 1) %>% 
  View()


data10 <- read.csv(here("origional-ecofish-data-models","Data","Processed",
                       "pko_SR_10.csv"))

data10 %>%
  # select(River, BroodYear) %>% 
  group_by(River, BroodYear) %>% 
  summarise(num_years = n()) %>% 
  filter(num_years > 1) %>% 
  View()
data10 %>% 
  filter(River=="KSI TS'OOHL TS'AP") %>% 
  View()


data11 <- read.csv(here("origional-ecofish-data-models","Data","Processed",
                        "PKE_SR_10.csv"))

data11 %>%
  # select(River, BroodYear) %>% 
  group_by(River, BroodYear) %>% 
  summarise(num_years = n()) %>% 
  filter(num_years > 1) %>% 
  View()

