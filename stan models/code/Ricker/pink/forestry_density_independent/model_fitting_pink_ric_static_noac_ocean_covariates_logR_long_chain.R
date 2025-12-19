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
file_ric_st=file.path(here('stan models', 'code','Ricker','pink','forestry_density_independent',
                           'ric_pink_static_noac_ocean_covariates_logR_new.stan'))
mric_p_st=cmdstanr::cmdstan_model(file_ric_st) #compile stan code to C++


#even year pinks
pk10r_e <- read.csv(here("origional-ecofish-data-models","Data","Processed","pke_SR_10_hat_yr_w_ersst.csv"))

#odd year pinks
pk10r_o <- read.csv(here("origional-ecofish-data-models","Data","Processed","pko_SR_10_hat_yr_w_ersst.csv"))

options(mc.cores=8)

# Pink salmon - even/odd broodlines #####
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',pk10r_o$River)
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=='915-765500-18600-00000-0000-0000-000-000-000-000-000-000','HEAD CREEK 2',pk10r_o$River)
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=='915-488000-41400-00000-0000-0000-000-000-000-000-000-000','WINDY BAY CREEK 2',pk10r_o$River)
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',pk10r_o$River)
pk10r_o=pk10r_o[order(factor(pk10r_o$River),pk10r_o$BroodYear),]
rownames(pk10r_o)=seq(1:nrow(pk10r_o))

pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',pk10r_e$River)
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='915-765500-18600-00000-0000-0000-000-000-000-000-000-000','HEAD CREEK 2',pk10r_e$River)
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='915-488000-41400-00000-0000-0000-000-000-000-000-000-000','WINDY BAY CREEK 2',pk10r_e$River)
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',pk10r_e$River)
pk10r_e=pk10r_e[order(factor(pk10r_e$River),pk10r_e$BroodYear),]
rownames(pk10r_e)=seq(1:nrow(pk10r_e))


#normalize ECA 2 - square root transformation (ie. sqrt(x))
pk10r_o$sqrt.ECA=sqrt(pk10r_o$ECA_age_proxy_forested_only)
pk10r_o$sqrt.ECA.std=(pk10r_o$sqrt.ECA-mean(pk10r_o$sqrt.ECA))/sd(pk10r_o$sqrt.ECA)

#normalize CPD 2 - square root transformation (ie. sqrt(x))
pk10r_o$sqrt.CPD=sqrt(pk10r_o$disturbedarea_prct_cs)
pk10r_o$sqrt.CPD.std=(pk10r_o$sqrt.CPD-mean(pk10r_o$sqrt.CPD))/sd(pk10r_o$sqrt.CPD)

#normalize ECA 2 - square root transformation (ie. sqrt(x))
pk10r_e$sqrt.ECA=sqrt(pk10r_e$ECA_age_proxy_forested_only)
pk10r_e$sqrt.ECA.std=(pk10r_e$sqrt.ECA-mean(pk10r_e$sqrt.ECA))/sd(pk10r_e$sqrt.ECA)

#normalize CPD 2 - square root transformation (ie. sqrt(x))
pk10r_e$sqrt.CPD=sqrt(pk10r_e$disturbedarea_prct_cs)
pk10r_e$sqrt.CPD.std=(pk10r_e$sqrt.CPD-mean(pk10r_e$sqrt.CPD))/sd(pk10r_e$sqrt.CPD)


#standardize npgo
pk10r_o$npgo.std=(pk10r_o$npgo-mean(pk10r_o$npgo))/sd(pk10r_o$npgo)
pk10r_e$npgo.std=(pk10r_e$npgo-mean(pk10r_e$npgo))/sd(pk10r_e$npgo)

pk10r_o$sst.std = (pk10r_o$spring_ersst-mean(pk10r_o$spring_ersst))/sd(pk10r_o$spring_ersst)
pk10r_e$sst.std = (pk10r_e$spring_ersst-mean(pk10r_e$spring_ersst))/sd(pk10r_e$spring_ersst)



pk10r_o$escapement.t_1=pk10r_e$Spawners[match(paste(pk10r_o$WATERSHED_CDE,pk10r_o$BroodYear-1),paste(pk10r_e$WATERSHED_CDE,pk10r_e$BroodYear))]
pk10r_e$escapement.t_1=pk10r_o$Spawners[match(paste(pk10r_e$WATERSHED_CDE,pk10r_e$BroodYear-1),paste(pk10r_o$WATERSHED_CDE,pk10r_o$BroodYear))]

pk10r_o$Broodline='Odd'
pk10r_e$Broodline='Even'

L_o=pk10r_o%>%group_by(River)%>%summarize(l=n(),by=min(BroodYear),tmin=(min(BroodYear)-min(pk10r_o$BroodYear))/2+1,tmax=(max(BroodYear)-min(pk10r_o$BroodYear))/2)
L_e=pk10r_e%>%group_by(River)%>%summarize(l=n(),by=min(BroodYear),tmin=(min(BroodYear)-min(pk10r_e$BroodYear))/2+1,tmax=(max(BroodYear)-min(pk10r_e$BroodYear))/2)
L_o$River2=paste(L_o$River,'Odd',sep='_')
L_e$River2=paste(L_e$River,'Even',sep='_')
L_all=rbind(L_e,L_o)
L_all=L_all[order(factor(L_all$River2)),]

pk10r_o$ii=as.numeric(factor(pk10r_o$BroodYear))
pk10r_e$ii=as.numeric(factor(pk10r_e$BroodYear))

pk10r=rbind(pk10r_e,pk10r_o)
pk10r$River2=paste(pk10r$River,pk10r$Broodline,sep='_')
pk10r=pk10r[order(factor(pk10r$River2),pk10r$BroodYear),]

#extract max S for priors on capacity & eq. recruitment
smax_prior=pk10r%>%group_by(River2) %>%summarize(m.s=max(Spawners),m.r=max(Recruits))


smax_prior = pk10r %>% 
  group_by(River2) %>% 
  summarize(m.s=Spawners[which.max(Recruits)],m.r=max(Recruits))

#ragged start and end points for each SR series
N_s=rag_n(pk10r$River2)

#cus by stock
cu1=distinct(pk10r,CU,.keep_all = T)
cu2=distinct(pk10r,River,.keep_all = T)
cu3=distinct(pk10r,River2,.keep_all = T)


dl_pk_eca=list(N=nrow(pk10r),
               L=length(unique(pk10r$BroodYear))/2, #total time now /2 for each broodline
               C=length(unique(pk10r$CU)),
               R=length(unique(pk10r$River)),
               J=length(unique(pk10r$River2)),
               B_i=as.numeric(factor(cu1$Broodline)),
               C_r=as.numeric(factor(cu2$CU)), #CU index by stock
               C_i=as.numeric(factor(cu3$CU)), #CU index by stock
               R_i=as.numeric(factor(cu3$River)),
               BL=as.numeric(factor(cu3$Broodline)),
               ii=pk10r$ii, #brood year index
               R_S=pk10r$ln_RS,
               S=pk10r$Spawners,
               logR=log(pk10r$Recruits),
               forest_loss=pk10r$sqrt.ECA.std,
               npgo=pk10r$npgo.std,
               sst=pk10r$sst.std,
               start_y=N_s[,1],
               end_y=N_s[,2],
               start_t=L_all$tmin,
               pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
               pSmax_sig=smax_prior$m.s*3,
               pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
               pRk_sig=smax_prior$m.r)

dl_pk_cpd=list(N=nrow(pk10r),
               L=length(unique(pk10r$BroodYear))/2, #total time now /2 for each broodline
               C=length(unique(pk10r$CU)),
               R=length(unique(pk10r$River)),
               J=length(unique(pk10r$River2)),
               B_i=as.numeric(factor(cu1$Broodline)),
               C_r=as.numeric(factor(cu2$CU)), #CU index by stock
               C_i=as.numeric(factor(cu3$CU)), #CU index by stock
               R_i=as.numeric(factor(cu3$River)),
               BL=as.numeric(factor(cu3$Broodline)),
               ii=pk10r$ii, #brood year index
               R_S=pk10r$ln_RS,
               S=pk10r$Spawners,
               logR=log(pk10r$Recruits),
               forest_loss=pk10r$sqrt.CPD.std,
               npgo=pk10r$npgo.std,
               sst=pk10r$sst.std,
               start_y=N_s[,1],
               end_y=N_s[,2],
               start_t=L_all$tmin,
               pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
               pSmax_sig=smax_prior$m.s*3,
               pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
               pRk_sig=smax_prior$m.r)

print("eca")


if(Sys.info()[7] == "mariakur") {
  print("Running on local machine")
  ric_pk_eca_st <- mric_p_st$sample(data=dl_pk_eca,
                                      chains = 1,
                                      iter_warmup = 20,
                                      iter_sampling =50,
                                      refresh = 10,
                                      adapt_delta = 0.95,
                                      max_treedepth = 20)
  write.csv(ric_pk_eca_st$summary(),here("stan models","outs","summary","ric_pk_eca_st_noac_ocean_covariates_logR_long_chain_trial.csv"))
  ric_pk_eca_st$save_object(here("stan models","outs","fits","ric_pk_eca_st_noac_ocean_covariates_long_chain_logR_trial.RDS"))
  
  post_ric_pk_eca_st=ric_pk_eca_st$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                     'b_npgo','b_npgo_cu','b_npgo_rv',
                                                     'b_sst','b_sst_cu','b_sst_rv',
                                                     'alpha_j','Smax','sigma','sigma_for_cu','sigma_for_rv'),format='draws_matrix')
  
  post_ric_pk_eca_st_mu2 <- ric_pk_eca_st$draws(variables=c('mu2'),format='draws_matrix')
  
  write.csv(post_ric_pk_eca_st,here('stan models','outs','posterior','ric_pk_eca_st_noac_ocean_covariates_logR_long_chain_trial.csv'))
  
  write.csv(post_ric_pk_eca_st_mu2,here('stan models','outs','posterior','ric_pk_eca_st_noac_ocean_covariates_logR_long_chain_mu2_trial.csv'))
  
} else {
  ric_pk_eca_st <- mric_p_st$sample(data=dl_pk_eca,
                                    chains = 6,
                                    iter_warmup = 1000,
                                    iter_sampling = 2000,
                                    refresh = 100,
                                    adapt_delta = 0.95,
                                    max_treedepth = 20,
                                    thin = 4)
  
  write.csv(ric_pk_eca_st$summary(),here("stan models","outs","summary","ric_pk_eca_st_noac_ocean_covariates_logR_long_chain.csv"))
  ric_pk_eca_st$save_object(here("stan models","outs","fits","ric_pk_eca_st_noac_ocean_covariates_logR_long_chain.RDS"))
  
  post_ric_pk_eca_st=ric_pk_eca_st$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                     'b_npgo','b_npgo_cu','b_npgo_rv',
                                                     'b_sst','b_sst_cu','b_sst_rv',
                                                     'alpha_j','Smax','sigma','sigma_for_cu','sigma_for_rv'),format='draws_matrix')
  
  post_ric_pk_eca_st_mu2 <- ric_pk_eca_st$draws(variables=c('mu2'),format='draws_matrix')
  
  write.csv(post_ric_pk_eca_st,here('stan models','outs','posterior','ric_pk_eca_st_noac_ocean_covariates_logR_long_chain.csv'))
  
  write.csv(post_ric_pk_eca_st_mu2,here('stan models','outs','posterior','ric_pk_eca_st_noac_ocean_covariates_logR_long_chain_mu2.csv'))
  
}

print("cpd")
#
if(Sys.info()[7] == "mariakur") {
  print("Running on local machine")
  ric_pk_cpd_st <- mric_p_st$sample(data=dl_pk_cpd,
                                    chains = 1,
                                    iter_warmup = 20,
                                    iter_sampling =50,
                                    refresh = 10,
                                    adapt_delta = 0.95,
                                    max_treedepth = 20,
                                    thin =2)
  write.csv(ric_pk_cpd_st$summary(),here("stan models","outs","summary","ric_pk_cpd_st_noac_ocean_covariates_logR_long_chain_trial.csv"))
  ric_pk_cpd_st$save_object(here("stan models","outs","fits","ric_pk_cpd_st_noac_ocean_covariates_logR_long_chain_trial.RDS"))
  
  post_ric_pk_cpd_st=ric_pk_cpd_st$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                     'b_npgo','b_npgo_cu','b_npgo_rv',
                                                     'b_sst','b_sst_cu','b_sst_rv',
                                                     'alpha_j','Smax','sigma','sigma_for_cu','sigma_for_rv'),format='draws_matrix')
  
  post_ric_pk_cpd_st_mu2 <- ric_pk_cpd_st$draws(variables=c('mu2'),format='draws_matrix')
  
  write.csv(post_ric_pk_cpd_st,here('stan models','outs','posterior','ric_pk_cpd_st_noac_ocean_covariates_logR_long_chain_trial.csv'))
  
  write.csv(post_ric_pk_cpd_st_mu2,here('stan models','outs','posterior','ric_pk_cpd_st_noac_ocean_covariates_logR_long_chain_mu2_trial.csv'))
  
} else {
  ric_pk_cpd_st <- mric_p_st$sample(data=dl_pk_cpd,
                                    chains = 6,
                                    iter_warmup = 1000,
                                    iter_sampling = 2000,
                                    refresh = 100,
                                    adapt_delta = 0.95,
                                    max_treedepth = 20,
                                    thin = 4)
  
  write.csv(ric_pk_cpd_st$summary(),here("stan models","outs","summary","ric_pk_cpd_st_noac_ocean_covariates_logR_long_chain.csv"))
  ric_pk_cpd_st$save_object(here("stan models","outs","fits","ric_pk_cpd_st_noac_ocean_covariates_logR_long_chain.RDS"))
  
  post_ric_pk_cpd_st=ric_pk_cpd_st$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                     'b_npgo','b_npgo_cu','b_npgo_rv',
                                                     'b_sst','b_sst_cu','b_sst_rv',
                                                     'alpha_j','Smax','sigma','sigma_for_cu','sigma_for_rv'),format='draws_matrix')
  
  post_ric_pk_cpd_st_mu2 = ric_pk_cpd_st$draws(variables=c('mu2'),format='draws_matrix')
  
  write.csv(post_ric_pk_cpd_st,here('stan models','outs','posterior','ric_pk_cpd_st_noac_ocean_covariates_logR_long_chain.csv'))
  write.csv(post_ric_pk_cpd_st_mu2,here('stan models','outs','posterior','ric_pk_cpd_st_noac_ocean_covariates_logR_long_chain_mu2.csv'))
  
}

