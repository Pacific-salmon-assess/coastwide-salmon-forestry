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
source(here('stan models','code','funcs.R'))

# load Stan model sets####
file_ric=file.path(here('stan models', 'code',
                        'Ricker', 'chum', 'forestry_density_independent',
                        'ric_chm_static_npgo_sst_logR_new.stan'))
mric=cmdstanr::cmdstan_model(file_ric) #compile stan code to C++

# load datasets####

# ch20r <- read.csv(here("origional-ecofish-data-models","Data","Processed","chum_SR_20_hat_yr_w_npgo.csv"))
ch20r <- read.csv(here("origional-ecofish-data-models","Data","Processed","chum_SR_20_hat_yr_w_ocean_covariates.csv"))

options(mc.cores=8)

## data formatting ####

#two rivers with duplicated names:
ch20r$River=ifelse(ch20r$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',ch20r$River)
ch20r$River=ifelse(ch20r$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',ch20r$River)


ch20r=ch20r[order(factor(ch20r$River),ch20r$BroodYear),]
rownames(ch20r)=seq(1:nrow(ch20r))

#normalize ECA 2 - square root transformation (ie. sqrt(x))
ch20r$sqrt.ECA=sqrt(ch20r$ECA_age_proxy_forested_only)
ch20r$sqrt.ECA.std=(ch20r$sqrt.ECA-mean(ch20r$sqrt.ECA))/sd(ch20r$sqrt.ECA)

#normalize CPD 2 - square root transformation (ie. sqrt(x))
ch20r$sqrt.CPD=sqrt(ch20r$disturbedarea_prct_cs)
ch20r$sqrt.CPD.std=(ch20r$sqrt.CPD-mean(ch20r$sqrt.CPD))/sd(ch20r$sqrt.CPD)

#standardize npgo
ch20r$npgo.std=(ch20r$npgo-mean(ch20r$npgo))/sd(ch20r$npgo)

ch20r$sst.std = (ch20r$spring_ersst-mean(ch20r$spring_ersst))/sd(ch20r$spring_ersst)

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

#data list for fits
dl_chm_eca_npgo_sst=list(N=nrow(ch20r),
                L=max(ch20r$BroodYear)-min(ch20r$BroodYear)+1,
                C=length(unique(ch20r$CU)),
                J=length(unique(ch20r$River)),
                C_i=as.numeric(factor(cu$CU)), #CU index by stock
                ii=as.numeric(factor(ch20r$BroodYear)), #brood year index
                R_S=ch20r$ln_RS,
                S=ch20r$Spawners, 
                logR=log(ch20r$Recruits),
                forest_loss=ch20r$sqrt.ECA.std, #design matrix for standardized ECA
                npgo=ch20r$npgo.std, #design matrix for standardized npgo
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
                     J=length(unique(ch20r$River)),
                     C_i=as.numeric(factor(cu$CU)), #CU index by stock
                     ii=as.numeric(factor(ch20r$BroodYear)), #brood year index
                     R_S=ch20r$ln_RS,
                     S=ch20r$Spawners, 
                     logR=log(ch20r$Recruits),
                     forest_loss=ch20r$sqrt.CPD.std, #design matrix for standardized ECA
                     npgo=ch20r$npgo.std, #design matrix for standardized npgo
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
                            iter_warmup = 50,
                            iter_sampling =50,
                            refresh = 10,
                            adapt_delta = 0.999,
                            max_treedepth = 20)
  write.csv(ric_chm_eca_npgo_sst$summary(),'./stan models/outs/summary/ric_chm_eca_ocean_covariates_logR_new_trial.csv')
  ric_chm_eca_npgo_sst$save_object('./stan models/outs/fits/ric_chm_eca_ocean_covariates_logR_new_trial.RDS')
  
  post_ric_chm_eca_npgo_sst=ric_chm_eca_npgo_sst$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                           'b_npgo','b_npgo_cu','b_npgo_rv',
                                                           'b_sst','b_sst_cu','b_sst_rv',
                                                           'alpha_j','Smax','sigma','mu2'),format='draws_matrix')
  write.csv(post_ric_chm_eca_npgo_sst,here('stan models','outs','posterior','ric_chm_eca_ocean_covariates_logR_new_trial.csv'))
  
} else {
  ric_chm_eca_npgo_sst <- mric$sample(data=dl_chm_eca_npgo_sst,
                            chains = 6, 
                            iter_warmup = 200,
                            iter_sampling = 700,
                            refresh = 100,
                            adapt_delta = 0.999,
                            max_treedepth = 20)
  
  write.csv(ric_chm_eca_npgo_sst$summary(),here("stan models","outs","summary","ric_chm_eca_ocean_covariates_logR_new.csv"))
  ric_chm_eca_npgo_sst$save_object(here("stan models","outs","fits","ric_chm_eca_ocean_covariates_logR_new.RDS"))
  
  post_ric_chm_eca_npgo_sst=ric_chm_eca_npgo_sst$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                           'b_npgo','b_npgo_cu','b_npgo_rv',
                                                           'b_sst','b_sst_cu','b_sst_rv',
                                                           'alpha_j','Smax','sigma','mu2'),format='draws_matrix')
  write.csv(post_ric_chm_eca_npgo_sst,here('stan models','outs','posterior','ric_chm_eca_ocean_covariates_logR_new.csv'))
  
}


print("cpd")


if(Sys.info()[7] == "mariakur") {
  print("Running on local machine")
  ric_chm_cpd_npgo_sst <- mric$sample(data=dl_chm_cpd_npgo_sst,
                                  chains = 2, 
                                  iter_warmup = 50,
                                  iter_sampling =50,
                                  refresh = 10,
                                  adapt_delta = 0.999,
                                  max_treedepth = 20)
  write.csv(ric_chm_cpd_npgo_sst$summary(),'./stan models/outs/summary/ric_chm_cpd_ocean_covariates_logR_new_trial.csv')
  ric_chm_cpd_npgo_sst$save_object('./stan models/outs/fits/ric_chm_cpd_ocean_covariates_logR_new_trial.RDS')
  
  post_ric_chm_cpd_npgo_sst=ric_chm_cpd_npgo_sst$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                           'b_npgo','b_npgo_cu','b_npgo_rv',
                                                           'b_sst','b_sst_cu','b_sst_rv',
                                                           'alpha_j','Smax','sigma','mu2'),format='draws_matrix')
  write.csv(post_ric_chm_cpd_npgo_sst,here('stan models','outs','posterior','ric_chm_cpd_ocean_covariates_logR_new_trial.csv'))
  
} else {
  ric_chm_cpd_npgo_sst <- mric$sample(data=dl_chm_cpd_npgo_sst,
                                  chains = 6, 
                                  iter_warmup = 200,
                                  iter_sampling = 700,
                                  refresh = 100,
                                  adapt_delta = 0.999,
                                  max_treedepth = 20)
  
  write.csv(ric_chm_cpd_npgo_sst$summary(),here("stan models","outs","summary","ric_chm_cpd_ocean_covariates_logR_new.csv"))
  ric_chm_cpd_npgo_sst$save_object(here("stan models","outs","fits","ric_chm_cpd_ocean_covariates_logR_new.RDS"))
  
  post_ric_chm_cpd_npgo_sst=ric_chm_cpd_npgo_sst$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                           'b_npgo','b_npgo_cu','b_npgo_rv',
                                                           'b_sst','b_sst_cu','b_sst_rv',
                                                           'alpha_j','Smax','sigma','mu2'),format='draws_matrix')
  write.csv(post_ric_chm_cpd_npgo_sst,here('stan models','outs','posterior','ric_chm_cpd_ocean_covariates_logR_new.csv'))
  
}

