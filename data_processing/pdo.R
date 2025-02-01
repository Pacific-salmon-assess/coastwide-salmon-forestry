# read pdo data from https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/index/ersst.v5.pdo.dat

# look at correlation between sst and pdo

# load libraries
library(tidyverse)
library(ggplot2)
library(here)
library(hues)

# read pdo data

url = "https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/index/ersst.v5.pdo.dat"

pdo = read.table(url, header = TRUE, skip = 1) #, col.names = c("year", "month", "pdo"))

# convert to long format

pdo_long <- pdo %>%
  pivot_longer(cols = -c(Year), names_to = "month", values_to = "pdo")

ch20rsc <- read.csv(here("origional-ecofish-data-models","Data","Processed","chum_SR_20_hat_yr_w_ersst.csv"))

ch20rsc$River=ifelse(ch20rsc$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',ch20rsc$River)
ch20rsc$River=ifelse(ch20rsc$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',ch20rsc$River)


ch20rsc=ch20rsc[order(factor(ch20rsc$River),ch20rsc$BroodYear),]

ch20rsc$River_n <- as.numeric(factor(ch20rsc$River))

#normalize ECA 2 - square root transformation (ie. sqrt(x))
ch20rsc$sqrt.ECA=sqrt(ch20rsc$ECA_age_proxy_forested_only)
ch20rsc$sqrt.ECA.std=(ch20rsc$sqrt.ECA-mean(ch20rsc$sqrt.ECA))/sd(ch20rsc$sqrt.ECA)

#normalize CPD 2 - square root transformation (ie. sqrt(x))
ch20rsc$sqrt.CPD=sqrt(ch20rsc$disturbedarea_prct_cs)
ch20rsc$sqrt.CPD.std=(ch20rsc$sqrt.CPD-mean(ch20rsc$sqrt.CPD))/sd(ch20rsc$sqrt.CPD)




# CU names
#####
cu_names <- data.frame(CU = c("CM-1","CM-2","CM-3","CM-4","CM-5","CM-6",
                              "CM-7","CM-8","CM-9","CM-10","CM-11","CM-12",
                              "CM-13","CM-14","CM-15","CM-16","CM-17","CM-18",
                              "CM-19","CM-20","CM-21","CM-22","CM-23","CM-24",
                              "CM-25","CM-26","CM-27","CM-28","CM-29","CM-30",
                              "CM-31","CM-32","CM-33","CM-34","CM-35","CM-36",
                              "CM-37","CM-38", "CM-39"),
                       CU_name = c("Fraser Canyon",
                                   "Lower Fraser",
                                   "Howe Sound-Burrard Inlet",
                                   "Georgia Strait",
                                   "East Vancouver Island",
                                   "Loughborough",
                                   "Bute Inlet",
                                   "Southern Coastal Streams",
                                   "Upper Knight",
                                   "Southwest Vancouver Island",
                                   "Northwest Vancouver Island",
                                   "Smith Inlet",
                                   "Rivers Inlet",
                                   "Wannock",
                                   "Spiller-Fitz Hugh-Burke",
                                   "Bella Coola - Dean Rivers",
                                   "Bella Coola River - Late",
                                   "Hecate Lowlands",
                                   "Mussel-Kynoch",
                                   "Douglas-Gardner",
                                   "East Haida Gwaii",
                                   "Skidegate",
                                   "West Haida Gwaii",
                                   "North Haida Gwaii",
                                   "North Haida Gwaii-Stanley Creek",
                                   "Skeena Estuary",
                                   "Lower Skeena",
                                   "Middle Skeena",
                                   "Upper Skeena",
                                   "Portland Inlet",
                                   "Lower Nass",
                                   "Portland Canal-Observatory",
                                   "Unuk",
                                   "Lower Stikine",
                                   "Whiting",
                                   "Taku",
                                   "Lynn Canal",
                                   "Teslin",
                                   "Lower Liard"))
#####

ch20rsc <- ch20rsc %>% left_join(cu_names, by = 'CU') 

pdo_long_annual <- pdo_long %>%
  group_by(Year) %>%
  summarise(pdo = mean(pdo))

#left join by year_SST and Year in ch20rsc and pdo_long_annual

ch20rsc <- ch20rsc %>% left_join(pdo_long_annual, by = c('year_SST' = 'Year'))

ch20rsc <- ch20rsc %>% mutate(pdo.std = (pdo - mean(pdo))/sd(pdo))

#plot correlation between spring_ersst and pdo

cor <- cor(ch20rsc$spring_ersst, ch20rsc$pdo, use = "complete.obs")

ch20rsc %>%
  select(CU_name, spring_ersst, pdo) %>%
  unique() %>%
  ggplot() +
  geom_point(aes(x = spring_ersst, y = pdo, color = CU_name), alpha = 0.2, size = 2) +
  geom_text(aes(label = paste("Correlation = ", round(cor, 2))),
            x = 12, y = -2, hjust = 0, vjust = 1) +
  scale_color_iwanthue(name = 'CU', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
                       lmin = 10, lmax = 95) + 
  labs(title = "Correlation between spring ERSST and PDO",
       x = "Spring ERSST",
       y = "PDO") +
  theme_classic()+
  theme(legend.position = "right",
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.4, "lines"),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.001, "cm"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5)
  )+
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 2), ncol = 1))

#save plot

ggsave(here("figures","pdo_sst_correlation.png"), width = 8, height = 6, dpi = 300)

#save data

write.csv(ch20rsc, here("origional-ecofish-data-models","Data","Processed","chum_SR_20_hat_yr_w_ersst_pdo.csv"), row.names = FALSE)


