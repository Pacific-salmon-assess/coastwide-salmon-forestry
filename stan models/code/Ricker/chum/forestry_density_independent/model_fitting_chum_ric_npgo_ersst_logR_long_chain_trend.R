#local machine or server

if(Sys.info()[7] == "mariakur") {
  print("Running on local machine")
  library(cmdstanr)
  set_cmdstan_path("C:/Users/mariakur/.cmdstan/cmdstan-2.35.0")
} else {
  print("Running on server")
  .libPaths(new = "/home/mkuruvil/R_Packages")
  library(cmdstanr)
  set_cmdstan_path("/home/mkuruvil/R_Packages/cmdstan-2.35.0")
}

#for reproducibility
set.seed(1194)

library(here);library(dplyr)
library(rstan)
library(tidyverse)
library(tictoc)
library(future)
library(furrr)
library(sf)
source(here('stan models','code','funcs.R'))

# load Stan model sets####
file_ric=file.path(here('stan models', 'code','Ricker',
                        'chum', 'forestry_density_independent',
                        'ric_chm_static_npgo_sst_logR_new_trend.stan'))
mric=cmdstanr::cmdstan_model(file_ric) #compile stan code to C++

# load datasets####

# ch20r <- read.csv(here("origional-ecofish-data-models","Data","Processed","chum_SR_20_hat_yr_w_npgo.csv"))
ch20r <- read.csv(here('origional-ecofish-data-models','Data',
                       'Processed', 'chum_SR_20_hat_yr_w_ersst_npgo.csv'))

lookup <- read.csv(here('origional-ecofish-data-models','Data',"forestry_data", "salmon_watersheds_lookup.csv"))

# Load population sheds ----
pop_sheds <- st_read(dsn = here("origional-ecofish-data-models","Data","forestry_data","Max_ECA_sheds_Nov_2024.gpkg"))  # Note that there are 1746 polygons, but only 1745 have VRI information in 2022.




options(mc.cores=8)

## data formatting ####

#two rivers with duplicated names:
ch20r$River=ifelse(ch20r$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',ch20r$River)
ch20r$River=ifelse(ch20r$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',ch20r$River)


ch20r=ch20r[order(factor(ch20r$River),ch20r$BroodYear),]
rownames(ch20r)=seq(1:nrow(ch20r))

ch20r$River_n <- as.numeric(factor(ch20r$River))

#normalize ECA 2 - square root transformation (ie. sqrt(x))
ch20r$sqrt.ECA=sqrt(ch20r$ECA_age_proxy_forested_only)
ch20r$sqrt.ECA.std=(ch20r$sqrt.ECA-mean(ch20r$sqrt.ECA))/sd(ch20r$sqrt.ECA)

#normalize CPD 2 - square root transformation (ie. sqrt(x))
ch20r$sqrt.CPD=sqrt(ch20r$disturbedarea_prct_cs)
ch20r$sqrt.CPD.std=(ch20r$sqrt.CPD-mean(ch20r$sqrt.CPD))/sd(ch20r$sqrt.CPD)

#standardize npgo
ch20r$winter_npgo.std=(ch20r$winter_npgo-mean(ch20r$winter_npgo))/sd(ch20r$winter_npgo)

ch20r$sst.std = (ch20r$spring_ersst-mean(ch20r$spring_ersst))/sd(ch20r$spring_ersst)

ch20r$year = (ch20r$BroodYear-min(ch20r$BroodYear))

##average ECA by stock
#just to see an overview of ECA by river
eca_s=ch20r%>%group_by(River)%>%summarize(m=mean(ECA_age_proxy_forested_only*100),range=max(ECA_age_proxy_forested_only*100)-min(ECA_age_proxy_forested_only*100),cu=unique(CU))

#extract max S for priors on capacity & eq. recruitment
smax_prior=
  ch20r %>%
  group_by(River) %>%
  summarize(m.s=Spawners[which.max(Recruits)],m.r=max(Recruits))

#ragged start and end points for each SR series
N_s=rag_n(ch20r$River)

#cus by stock
cu=distinct(ch20r,River,.keep_all = T)
cu.nrv=summary(factor(cu$CU))

#time points for each series
L_i=ch20r%>%group_by(River)%>%summarize(l=n(),min=min(BroodYear),max=max(BroodYear),tmin=min(BroodYear)-1954+1,tmax=max(BroodYear)-1954+1)


lookup <- lookup %>% 
  left_join(ch20r %>% select(GFE_ID, River_n, BroodYear, sqrt.CPD.std, Species, disturbedarea_prct_cs) %>% 
              group_by(GFE_ID, River_n, Species) %>% #filter only max year
              filter(BroodYear == max(BroodYear)) %>% unique(), by = c("GFE_ID","Species"))

rivers <- lookup %>% select(River_n) %>%  filter(!is.na(River_n)) %>% pull(River_n) %>% unique()


ch20r_w_outlet_lfid <- ch20r %>% 
  left_join(lookup %>% select(CU, GFE_ID, Species, LINEAR_FEATURE_ID), by = c("CU" = "CU", "GFE_ID" = "GFE_ID", "Species" = "Species")) %>% 
  left_join(pop_sheds %>% select(outlet_lfid, Region) %>% mutate(outlet_lfid = as.integer(outlet_lfid)), 
            by = c("LINEAR_FEATURE_ID" = "outlet_lfid")) 

region=distinct(ch20r_w_outlet_lfid,River,.keep_all = T)

#data list for fits
dl_chm_eca_npgo_sst=list(N=nrow(ch20r),
                L=max(ch20r$BroodYear)-min(ch20r$BroodYear)+1,
                C=length(unique(ch20r$CU)),
                Region=length(unique(ch20r_w_outlet_lfid$Region)),
                J=length(unique(ch20r$River)),
                C_i=as.numeric(factor(cu$CU)), #CU index by stock
                Region_i=as.numeric(factor(region$Region)), #region index by stock
                year = ch20r$year,
                ii=as.numeric(factor(ch20r$BroodYear)), #brood year index
                R_S=ch20r$ln_RS,
                S=ch20r$Spawners, 
                logR=log(ch20r$Recruits),
                forest_loss=ch20r$sqrt.ECA.std, #design matrix for standardized ECA
                npgo=ch20r$winter_npgo.std, #design matrix for standardized npgo
                sst=ch20r$sst.std,
                start_y=N_s[,1],
                end_y=N_s[,2],
                start_t=L_i$tmin,
                end_t=L_i$tmax,
                pSmax_mean=smax_prior$m.s, #prior for Smax (spawners that maximize recruitment) based on max observed spawners
                pSmax_sig=3*smax_prior$m.s,
                pRk_mean=0.75*smax_prior$m.r, ##prior for Rk (recruitment capacity) based on max observed spawners
                pRk_sig=smax_prior$m.r)
                
dl_chm_cpd_npgo_sst=list(N=nrow(ch20r),
                     L=max(ch20r$BroodYear)-min(ch20r$BroodYear)+1,
                     C=length(unique(ch20r$CU)),
                     Region=length(unique(ch20r_w_outlet_lfid$Region)),
                     J=length(unique(ch20r$River)),
                     C_i=as.numeric(factor(cu$CU)), #CU index by stock
                     Region_i=as.numeric(factor(region$Region)), #region index by stock
                     year = ch20r$year,
                     ii=as.numeric(factor(ch20r$BroodYear)), #brood year index
                     R_S=ch20r$ln_RS,
                     S=ch20r$Spawners, 
                     logR=log(ch20r$Recruits),
                     forest_loss=ch20r$sqrt.CPD.std, #design matrix for standardized ECA
                     npgo=ch20r$winter_npgo.std, #design matrix for standardized npgo
                     sst=ch20r$sst.std,
                     start_y=N_s[,1],
                     end_y=N_s[,2],
                     start_t=L_i$tmin,
                     end_t=L_i$tmax,
                     pSmax_mean=smax_prior$m.s, #prior for Smax (spawners that maximize recruitment) based on max observed spawners
                     pSmax_sig=3*smax_prior$m.s,
                     pRk_mean=0.75*smax_prior$m.r, ##prior for Rk (recruitment capacity) based on max observed spawners
                     pRk_sig=smax_prior$m.r)


print("eca")


if(Sys.info()[7] == "mariakur") {
  print("Running on local machine")
  ric_chm_eca_npgo_sst <- mric$sample(data=dl_chm_eca_npgo_sst,
                            chains = 2, 
                            iter_warmup = 100,
                            iter_sampling = 200,
                            refresh = 10,
                            adapt_delta = 0.999,
                            max_treedepth = 20,
                            thin = 4)
  
  write.csv(ric_chm_eca_npgo_sst$summary(),
            here(#'salmon_forestry_data_analysis',
                 'stan models', 'outs', 'summary',
                 'ric_chm_eca_ocean_covariates_logR_long_chain_trend_trial.csv'))
  
  ric_chm_eca_npgo_sst$save_object(here(#'salmon_forestry_data_analysis',
                                        'stan models', 'outs', 
  'fits','ric_chm_eca_ocean_covariates_logR_long_chain_trend_trial.RDS'))
  
  post_ric_chm_eca_npgo_sst=ric_chm_eca_npgo_sst$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                           'b_npgo','b_npgo_cu','b_npgo_rv',
                                                           'b_sst','b_sst_cu','b_sst_rv',
                                                           'alpha_j','Smax','sigma'),format='draws_matrix')
  
  post_ric_chm_eca_npgo_sst_mu2=ric_chm_eca_npgo_sst$draws(variables=c('mu2'),format='draws_matrix')
  
  write.csv(post_ric_chm_eca_npgo_sst,here(#'salmon_forestry_data_analysis',
                                           'stan models','outs',
                                           'posterior','ric_chm_eca_ocean_covariates_logR_long_chain_trend_trial.csv'))
  write.csv(post_ric_chm_eca_npgo_sst_mu2,here(#'salmon_forestry_data_analysis',
                                               'stan models',
                                               'outs','posterior','ric_chm_eca_ocean_covariates_logR_long_chain_trend_trial_mu2.csv'))
  
} else {
  ric_chm_eca_npgo_sst <- mric$sample(data=dl_chm_eca_npgo_sst,
                            chains = 6, 
                            iter_warmup = 1000,
                            iter_sampling = 2000,
                            refresh = 100,
                            adapt_delta = 0.999,
                            max_treedepth = 20)
  
  write.csv(ric_chm_eca_npgo_sst$summary(),here(#'salmon_forestry_data_analysis',
                                                "stan models","outs","summary","ric_chm_eca_ocean_covariates_logR_long_chain_trend.csv"))
  ric_chm_eca_npgo_sst$save_object(here(#'salmon_forestry_data_analysis',
                                        "stan models","outs","fits","ric_chm_eca_ocean_covariates_logR_long_chain_trend.RDS"))
  
  post_ric_chm_eca_npgo_sst=ric_chm_eca_npgo_sst$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                           'b_npgo','b_npgo_cu','b_npgo_rv',
                                                           'b_sst','b_sst_cu','b_sst_rv',
                                                           'alpha_j','Smax','sigma'),format='draws_matrix')
  
  post_ric_chm_eca_npgo_sst_mu2=ric_chm_eca_npgo_sst$draws(variables=c('mu2'),format='draws_matrix')
  
  
  write.csv(post_ric_chm_eca_npgo_sst,here(#'salmon_forestry_data_analysis',
                                           'stan models','outs','posterior','ric_chm_eca_ocean_covariates_logR_long_chain_trend.csv'))
  write.csv(post_ric_chm_eca_npgo_sst_mu2,here(#'salmon_forestry_data_analysis',
                                               'stan models','outs','posterior','ric_chm_eca_ocean_covariates_logR_long_chain_trend_mu2.csv'))
  
}


print("cpd")


if(Sys.info()[7] == "mariakur") {
  print("Running on local machine")
  ric_chm_cpd_npgo_sst <- mric$sample(data=dl_chm_cpd_npgo_sst,
                                  chains = 2, 
                                  iter_warmup = 10,
                                  iter_sampling = 20,
                                  refresh = 10,
                                  adapt_delta = 0.999,
                                  max_treedepth = 20)
  
  write.csv(ric_chm_cpd_npgo_sst$summary(),
            here(#'salmon_forestry_data_analysis',
                 'stan models', 'outs', 'summary',
                 'ric_chm_cpd_ocean_covariates_logR_long_chain_trial.csv'))
  ric_chm_cpd_npgo_sst$save_object(here(#'salmon_forestry_data_analysis',
    'stan models','outs',
    'fits','ric_chm_cpd_ocean_covariates_logR_long_chain_trend_trial.RDS'))
  
  post_ric_chm_cpd_npgo_sst=ric_chm_cpd_npgo_sst$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                           'b_npgo','b_npgo_cu','b_npgo_rv',
                                                           'b_sst','b_sst_cu','b_sst_rv',
                                                           'alpha_j','Smax','sigma'),format='draws_matrix')
  
  post_ric_chm_cpd_npgo_sst_mu2=ric_chm_cpd_npgo_sst$draws(variables=c('mu2'),format='draws_matrix')
  
  write.csv(post_ric_chm_cpd_npgo_sst,here(#'salmon_forestry_data_analysis',
                                           'stan models','outs','posterior','ric_chm_cpd_ocean_covariates_logR_long_chain_trend_trial.csv'))
  write.csv(post_ric_chm_cpd_npgo_sst_mu2,here(#'salmon_forestry_data_analysis',
                                               'stan models','outs','posterior','ric_chm_cpd_ocean_covariates_logR_long_chain_trend_trial_mu2.csv'))
  
} else {
  ric_chm_cpd_npgo_sst <- mric$sample(data=dl_chm_cpd_npgo_sst,
                                  chains = 6, 
                                  iter_warmup = 1000,
                                  iter_sampling = 2000,
                                  refresh = 200,
                                  adapt_delta = 0.999,
                                  max_treedepth = 20)
  
  write.csv(ric_chm_cpd_npgo_sst$summary(),here(#'salmon_forestry_data_analysis',
                                                "stan models","outs","summary","ric_chm_cpd_ocean_covariates_logR_long_chain_trend.csv"))
  ric_chm_cpd_npgo_sst$save_object(here(#'salmon_forestry_data_analysis',
                                        "stan models","outs","fits","ric_chm_cpd_ocean_covariates_logR_long_chain_trend.RDS"))
  
  post_ric_chm_cpd_npgo_sst=ric_chm_cpd_npgo_sst$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                           'b_npgo','b_npgo_cu','b_npgo_rv',
                                                           'b_sst','b_sst_cu','b_sst_rv',
                                                           'alpha_j','Smax','sigma'),format='draws_matrix')
  
  post_ric_chm_cpd_npgo_sst_mu2=ric_chm_cpd_npgo_sst$draws(variables=c('mu2'),format='draws_matrix')
  
  write.csv(post_ric_chm_cpd_npgo_sst,here(#'salmon_forestry_data_analysis',
                                           'stan models','outs','posterior','ric_chm_cpd_ocean_covariates_logR_long_chain_trend.csv'))
  write.csv(post_ric_chm_cpd_npgo_sst_mu2,here(#'salmon_forestry_data_analysis',
                                               'stan models','outs','posterior','ric_chm_cpd_ocean_covariates_logR_long_chain_trend_mu2.csv'))
  
}

