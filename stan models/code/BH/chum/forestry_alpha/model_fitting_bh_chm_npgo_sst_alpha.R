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
file_bh=file.path(here('stan models', 'code','BH','chum','forestry_alpha','bh_chm_npgo_sst_alpha.stan'))
mbh=cmdstanr::cmdstan_model(file_bh) #compile stan code to C++

# load datasets####

# ch20r <- read.csv(here("origional-ecofish-data-models","Data","Processed","chum_SR_20_hat_yr_w_npgo.csv"))
ch20r <- read.csv(here("origional-ecofish-data-models","Data","Processed","chum_SR_20_hat_yr_w_ersst.csv"))


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

ch20r$npgo_std = (ch20r$npgo-mean(ch20r$npgo))/sd(ch20r$npgo)

ch20r$sst.std = (ch20r$spring_ersst-mean(ch20r$spring_ersst))/sd(ch20r$spring_ersst)

##average ECA by stock
#just to see an overview of ECA by river
eca_s=ch20r%>%group_by(River)%>%summarize(m=mean(ECA_age_proxy_forested_only*100),range=max(ECA_age_proxy_forested_only*100)-min(ECA_age_proxy_forested_only*100),cu=unique(CU))

#extract max S for priors on capacity & eq. recruitment
smax_prior=ch20r%>%group_by(River) %>%summarize(m.s=max(Spawners),m.r=max(Recruits))

#ragged start and end points for each SR series
N_s=rag_n(ch20r$River)

#cus by stock
cu=distinct(ch20r,River,.keep_all = T)
cu.nrv=summary(factor(cu$CU))

#time points for each series
L_i=ch20r%>%group_by(River)%>%summarize(l=n(),min=min(BroodYear),max=max(BroodYear),tmin=min(BroodYear)-1954+1,tmax=max(BroodYear)-1954+1)

#data list for fits
dl_chm_eca_npgo=list(N=nrow(ch20r),
                     L=max(ch20r$BroodYear)-min(ch20r$BroodYear)+1,
                     C=length(unique(ch20r$CU)),
                     J=length(unique(ch20r$River)),
                     C_i=as.numeric(factor(cu$CU)), #CU index by stock
                     ii=as.numeric(factor(ch20r$BroodYear)), #brood year index
                     R_S=ch20r$ln_RS,
                     S=ch20r$Spawners, 
                     forest_loss=ch20r$sqrt.ECA.std, #design matrix for standardized ECA
                     npgo=ch20r$npgo_std,
                     sst=ch20r$sst.std,
                     start_y=N_s[,1],
                     end_y=N_s[,2],
                     start_t=L_i$tmin,
                     end_t=L_i$tmax,
                     # pSmax_mean=0.5*smax_prior$m.s, #prior for Smax (spawners that maximize recruitment) based on max observed spawners
                     # pSmax_sig=smax_prior$m.s,
                     pRk_mean=0.75*smax_prior$m.r, ##prior for Rk (recruitment capacity) based on max observed spawners
                     pRk_sig=smax_prior$m.r)

dl_chm_cpd_npgo=list(N=nrow(ch20r),
                     L=max(ch20r$BroodYear)-min(ch20r$BroodYear)+1,
                     C=length(unique(ch20r$CU)),
                     J=length(unique(ch20r$River)),
                     C_i=as.numeric(factor(cu$CU)), #CU index by stock
                     ii=as.numeric(factor(ch20r$BroodYear)), #brood year index
                     R_S=ch20r$ln_RS,
                     S=ch20r$Spawners, 
                     forest_loss=as.vector(ch20r$sqrt.CPD.std), #design matrix for standardized ECA
                     npgo = as.vector(ch20r$npgo_std),
                     sst = as.vector(ch20r$sst.std),
                     ECA=as.vector(ch20r$sqrt.CPD.std), #design matrix for standardized ECA
                     start_y=N_s[,1],
                     end_y=N_s[,2],
                     start_t=L_i$tmin,
                     end_t=L_i$tmax,
                     # pSmax_mean=0.5*smax_prior$m.s, #prior for Smax (spawners that maximize recruitment) based on max observed spawners
                     # pSmax_sig=smax_prior$m.s,
                     pRk_mean=0.75*smax_prior$m.r, #prior for Rk (recruitment capacity) based on max observed spawners
                     pRk_sig=smax_prior$m.r)



print("eca")

if(Sys.info()[7] == "mariakur") {
  print("Running on local machine")
  bh_chm_eca_npgo <- mbh$sample(data=dl_chm_eca_npgo,
                                chains = 6,
                                iter_warmup = 20,
                                iter_sampling =50,
                                refresh = 10,
                                adapt_delta = 0.95,
                                max_treedepth = 20)
  write.csv(bh_chm_eca_npgo$summary(),'./stan models/outs/summary/bh_chm_eca_npgo_sst_alpha_trial.csv')
  bh_chm_eca_npgo$save_object('./stan models/outs/fits/bh_chm_eca_npgo_sst_alpha_trial.RDS')
  
  post_bh_chm_eca_npgo=bh_chm_eca_npgo$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                         'b_npgo','b_npgo_cu','b_npgo_rv',
                                                         'b_sst','b_sst_cu','b_sst_rv',
                                                         'alpha_j','Rk','sigma'),format='draws_matrix')
  write.csv(post_bh_chm_eca_npgo,here('stan models','outs','posterior','bh_chm_eca_npgo_sst_alpha_trial.csv'))
  
} else {
  bh_chm_eca_npgo <- mbh$sample(data=dl_chm_eca_npgo,
                                chains = 6,
                                iter_warmup = 200,
                                iter_sampling =500,
                                refresh = 100,
                                adapt_delta = 0.999,
                                max_treedepth = 20)
  
  write.csv(bh_chm_eca_npgo$summary(),here("stan models","outs","summary","bh_chm_eca_npgo_sst_alpha.csv"))
  bh_chm_eca_npgo$save_object(here("stan models","outs","fits","bh_chm_eca_npgo_sst_alpha.RDS"))
  
  post_bh_chm_eca_npgo=bh_chm_eca_npgo$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                         'b_npgo','b_npgo_cu','b_npgo_rv',
                                                         'b_sst','b_sst_cu','b_sst_rv',
                                                         'alpha_j','Rk','sigma'),format='draws_matrix')
  write.csv(post_bh_chm_eca_npgo,here('stan models','outs','posterior','bh_chm_eca_npgo_sst_alpha.csv'))
  
}

print("cpd")

if(Sys.info()[7] == "mariakur") {
  print("Running on local machine")
  bh_chm_cpd_npgo <- mbh$sample(data=dl_chm_cpd_npgo,
                                chains = 6, 
                                iter_warmup = 20,
                                iter_sampling =50,
                                refresh = 10,
                                adapt_delta = 0.95,
                                max_treedepth = 20)
  write.csv(bh_chm_cpd_npgo$summary(),'./stan models/outs/summary/bh_chm_cpd_npgo_sst_alpha_trial.csv')
  bh_chm_cpd_npgo$save_object('./stan models/outs/fits/bh_chm_cpd_npgo_sst_alpha_trial.RDS')
  
  post_bh_chm_cpd_npgo=bh_chm_cpd_npgo$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                         'b_npgo','b_npgo_cu','b_npgo_rv',
                                                         'b_sst','b_sst_cu','b_sst_rv',
                                                         'alpha_j','Rk','sigma'),format='draws_matrix')
  write.csv(post_bh_chm_cpd_npgo,here('stan models','outs','posterior','bh_chm_cpd_npgo_sst_alpha_trial.csv'))
  
} else {
  bh_chm_cpd_npgo <- mbh$sample(data=dl_chm_cpd_npgo,
                                chains = 6, 
                                iter_warmup = 200,
                                iter_sampling =500,
                                refresh = 100,
                                adapt_delta = 0.999,
                                max_treedepth = 20)
  
  write.csv(bh_chm_cpd_npgo$summary(),here("stan models","outs","summary","bh_chm_cpd_npgo_sst_alpha.csv"))
  bh_chm_cpd_npgo$save_object(here("stan models","outs","fits","bh_chm_cpd_npgo_sst_alpha.RDS"))
  
  post_bh_chm_cpd_npgo=bh_chm_cpd_npgo$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                                         'b_npgo','b_npgo_cu','b_npgo_rv',
                                                         'b_sst','b_sst_cu','b_sst_rv',
                                                         'alpha_j','Rk','sigma'),format='draws_matrix')
  write.csv(post_bh_chm_cpd_npgo,here('stan models','outs','posterior','bh_chm_cpd_npgo_sst_alpha.csv'))
  
}





