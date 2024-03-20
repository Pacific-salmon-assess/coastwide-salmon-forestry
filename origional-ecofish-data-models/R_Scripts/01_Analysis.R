## Project Name:
##
## Project Code:
##
## Script name:
##
##Purpose of script:
##
##Author:Joao Braga
##
##Date Created:2023-12-19
##
## Ecofish Research(VictoriaOffice),2023
## Email:jbraga@ecofishresearch.com
##
##Datasources:
##
##
##DARrequest:
##
##
##Repositorylink:
##
##
##
##Drive path:
##
##
## Notes:
##  The objective of this script is to test the ln_RS models for Cushing and give more clarity to the Ricker model
## 
## Variable Names within the data sets:
##
##

##Requiredpackages
list.of.packages <- c("tidyverse", "extrafont", "lubridate", "ggplot2", "sf", "qs", "lme4", "lmerTest", "glmmTMB", "performance", "parallel", "openxlsx")

#Whatpackagesneedtobeinstalled?
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]

##Installmissingpackages
if(length(new.packages)){
  install.packages(new.packages, dependencies = TRUE , repos = c("http://cran.rstudio.com/","http://R-Forge.R-project.org"))
}

##Loadinglibraries
lapply(list.of.packages, library, character.only = TRUE)

cores <- detectCores()

# Loading Data ----
## Chum ----
#  The chum salmon and forestry data was compiled based on two data scenarions:
#' *1) Scenario 1 - SR_20_hat_yr_reduced_VRI90 - This scenario requires salmon populations to have at least 20 (chum); otherwise, they were excluded from the analysis. *
#' *The scenario also removes stock-recruit years that were enhanced by salmon hatcheries and requires complete forest age data to cover at least 90% of the total watershed area. For chum salmon, we also removed conversation unit 4.*
#' 
#' *2) Scenario 2 - SR_20 - This scenario requires salmon populations to have at least 20 (chum); otherwise, they were excluded. *
#' * This scenario also requires complete forest age data to cover at least 90% of the total watershed area. *

chum_SR_20_hat_yr_reduced_VRI90 <- read_csv("Data/Processed/chum_SR_20_hat_yr_reduced_VRI90.csv")
chum_SR_20 <- read_csv("Data/Processed/chum_SR_20.csv")

names(chum_SR_20)
# 	Species - Salmon species code (CM = Chum, PKE = Pink even years, PKO = Pink odd years);                             
# 	River - Name of the salmon population and denotes the spatial boundary of the population (i.e., watershed);
# 	WATERSHED_CDE - 50k watershed code;
# 	DFO - DFO area code;                   
# 	CU - Conservation Unit Index;
# 	BroodYear - Brood Year, or spawning year;
# 	hararea_m2_3005 - annual disturbed area;
# 	area_km2 - watershed area;
# 	ln_area_km2 - natural log of watershed area;
# 	Spawners - # of spawners;
#   	logS -    natural log of spawners;
# 	Recruits - # of estimated recruits;
#   	Returns - # of estimated returns;
#   	recruitPerSpawner - # of recruits per spawner. Obtained by dividing “Recruits” by “Spawners”;
#   	ln_RS - natural log of “recruitPerSpawner”;
# 	ECA_age_proxy - Equivalent clearcut area estimated for the entire watershed;
# 	ECA_age_proxy_forested_only - Equivalent clearcut area estimated for the productive areas only (areas that grow or could grow forests, if undisturbed)’;
# 	lake_cover_pct - proportion lake coverage within the watershed;
# 	vri_bec_zone - Bioclimatic Zone;
# 	npgo - North Pacific Gyro index;
# 	portion_reporting_vri_cover1_missing - Proportion of unreported VRI cover within the watershed;
# 	hatchery_enhancement - Logical variable. Indicates if the watershed contains hatchery enhancement data;
# 	year_of_enhacement - Logical variable. Indicates if the watershed and year contain hatchery enhancement data;
# 	ln_area_km2_std - Standardized “ln_area_km2”. Standardized variables refer to variables where they were centered around zero and divided by one standard deviation;
# 	ECA_age_proxy_forested_only_std - Standardized “ECA_age_proxy_forested_only”. Standardized variables refer to variables where they were centered around zero and divided by one standard deviation;

## Pink (even) ----
#  The pink (even) salmon and forestry data was compiled based on two data scenarions:
#' *1) Scenario 1 - SR_20_hat_yr_reduced_VRI90 - This scenario requires salmon populations to have at least 10 (pink); otherwise, they were excluded from the analysis. *
#' *The scenario also removes stock-recruit years that were enhanced by salmon hatcheries and requires complete forest age data to cover at least 90% of the total watershed area*
#' 
#' *2) Scenario 2 - SR_20 - This scenario requires salmon populations to have at least 10 (pink); otherwise, they were excluded. *
#' * This scenario also requires complete forest age data to cover at least 90% of the total watershed area. *


pke_SR_10_hat_yr_reduced_VRI90 <- read_csv("Data/Processed/pke_SR_10_hat_yr_reduced_VRI90.csv")
pke_SR_10 <- read_csv("Data/Processed/pke_SR_10.csv")

## Pink (odd) ----
#  The pink (odd) salmon and forestry data was compiled based on two data scenarions:
#' *1) Scenario 1 - SR_20_hat_yr_reduced_VRI90 - This scenario requires salmon populations to have at least 10 (pink); otherwise, they were excluded from the analysis. *
#' *The scenario also removes stock-recruit years that were enhanced by salmon hatcheries and requires complete forest age data to cover at least 90% of the total watershed area*
#' 
#' *2) Scenario 2 - SR_20 - This scenario requires salmon populations to have at least 10 (pink); otherwise, they were excluded. *
#' * This scenario also requires complete forest age data to cover at least 90% of the total watershed area. *

pko_SR_10_hat_yr_reduced_VRI90 <- read_csv("Data/Processed/PKO_SR_10_hat_yr_reduced_VRI90.csv")
pko_SR_10 <- read_csv("Data/Processed/PKO_SR_10.csv")

# Modelling ----
## Chum ----
# # of populations
length(unique(chum_SR_20_hat_yr_reduced_VRI90$River))

### Stock-Recruitment Comparison ----
#### Ricker ----
CM_mod_ricker_SR <- glmmTMB(ln_RS ~  River:Spawners + (1|BroodYear/CU) + (1|CU),
                            data= chum_SR_20_hat_yr_reduced_VRI90,
                            REML = FALSE,
                            na.action =  "na.fail",
                            control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max = 2e5), parallel = cores - 2))

#### Cushing ----
CM_mod_cush_SR <- glmmTMB(ln_RS ~  River:logS + (1|BroodYear/CU) + (1|CU),
                          data= chum_SR_20_hat_yr_reduced_VRI90, 
                          REML = FALSE,
                          na.action =  "na.fail", 
                          control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max = 2e5), parallel = cores - 2))

### Adding Watershed Area to CU models ----
##### Ricker ----
CM_mod_ricker_ECA_area_int <- glmmTMB(ln_RS ~  ECA_age_proxy_forested_only_std * ln_area_km2_std + River:Spawners + (1 | BroodYear/CU) + (ECA_age_proxy_forested_only_std | CU),
                                      data= chum_SR_20_hat_yr_reduced_VRI90, 
                                      REML = FALSE,
                                      na.action =  "na.fail", 
                                      control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max = 2e5), parallel = cores - 2))  
##### Cushing ----
CM_mod_cush_ECA_area_int <- glmmTMB(ln_RS ~  ECA_age_proxy_forested_only_std * ln_area_km2_std + River:logS + (1 | BroodYear/CU) + (ECA_age_proxy_forested_only_std | CU),
                                    data= chum_SR_20_hat_yr_reduced_VRI90, 
                                    REML = FALSE,
                                    na.action =  "na.fail", 
                                    control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max = 2e5), parallel = cores - 2))  

## Pink Salmon Even ----
### Stock-Recruitment Comparison ----
#### Ricker ----
PKE_mod_ricker_SR <- glmmTMB(ln_RS ~  River:Spawners + (1|BroodYear/CU) + (1|CU),
                             data= pke_SR_10_hat_yr_reduced_VRI90,
                             REML = FALSE,
                             na.action =  "na.fail",
                             control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max = 2e5), parallel = cores - 2))

#### Cushing ----
PKE_mod_cush_SR <- glmmTMB(ln_RS ~  River:logS + (1|BroodYear/CU) + (1|CU),
                           data= pke_SR_10_hat_yr_reduced_VRI90, 
                           REML = FALSE,
                           na.action =  "na.fail", 
                           control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max = 2e5), parallel = cores - 2))  

### Adding Watershed Area to CU models ----
#### Ricker ----
# This model has convergence issues
PKE_mod_ricker_ECA_area_int <- glmmTMB(ln_RS ~  ECA_age_proxy_forested_only_std * ln_area_km2_std + River:Spawners + (1|BroodYear/CU) + (ECA_age_proxy_forested_only_std | CU),
                                       data= pke_SR_10_hat_yr_reduced_VRI90,
                                       REML = FALSE,
                                       na.action =  "na.fail",
                                       control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max = 2e5), parallel = cores - 2))

#### Cushing ----
PKE_mod_cush_ECA_area_int <- glmmTMB(ln_RS ~  ECA_age_proxy_forested_only_std * ln_area_km2_std + River:logS +(1|BroodYear/CU) + (ECA_age_proxy_forested_only_std | CU),
                                     data= pke_SR_10_hat_yr_reduced_VRI90, 
                                     REML = FALSE,
                                     na.action =  "na.fail", 
                                     control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max = 2e5), parallel = cores - 2))  

## Pink Salmon Odd ----
### Stock-Recruitment Comparison ----
#### Ricker ----
# This model has convergence issues
PKO_mod_ricker_SR <- glmmTMB(ln_RS ~  River:Spawners + (1|BroodYear/CU) + (1|CU),
                             data= pko_SR_10_hat_yr_reduced_VRI90,
                             REML = FALSE,
                             na.action =  "na.fail",
                             control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max = 2e5), parallel = cores - 2))

#### Cushing ----
PKO_mod_cush_SR <- glmmTMB(ln_RS ~  River:logS + (1|BroodYear/CU) + (1|CU),
                           data= pko_SR_10_hat_yr_reduced_VRI90, 
                           REML = FALSE,
                           na.action =  "na.fail", 
                           control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max = 2e5), parallel = cores - 2))  

### Adding Watershed Area to CU models ----
#### Ricker ----
# This model has convergence issues
PKO_mod_ricker_ECA_area_int <- glmmTMB(ln_RS ~  ECA_age_proxy_forested_only_std * ln_area_km2_std + River:Spawners + (1|BroodYear/CU) + (ECA_age_proxy_forested_only_std | CU),
                                       data= pko_SR_10_hat_yr_reduced_VRI90,
                                       REML = FALSE,
                                       na.action =  "na.fail",
                                       control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max = 2e5), parallel = cores - 2))

#### Cushing ----
# This model has convergence issues
PKO_mod_cush_ECA_area_int <- glmmTMB(ln_RS ~  ECA_age_proxy_forested_only_std * ln_area_km2_std + River:logS +(1|BroodYear/CU) + (ECA_age_proxy_forested_only_std | CU),
                                     data= pko_SR_10_hat_yr_reduced_VRI90, 
                                     REML = FALSE,
                                     na.action =  "na.fail", 
                                     control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max = 2e5), parallel = cores - 2))  

# End