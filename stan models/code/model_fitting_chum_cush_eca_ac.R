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
file_cush=file.path(here('stan models', 'code', 'cush_chm_ac.stan'))
mcush=cmdstanr::cmdstan_model(file_cush) #compile stan code to C++

# load datasets####

ch20r <- read.csv(here("origional-ecofish-data-models","Data","Processed","chum_SR_20_hat_yr_reduced_VRI90.csv"))


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
dl_chm_eca=list(N=nrow(ch20r),
                L=max(ch20r$BroodYear)-min(ch20r$BroodYear)+1,
                C=length(unique(ch20r$CU)),
                J=length(unique(ch20r$River)),
                C_i=as.numeric(factor(cu$CU)), #CU index by stock
                ii=as.numeric(factor(ch20r$BroodYear)), #brood year index
                R_S=ch20r$ln_RS,
                S=ch20r$Spawners, 
                forest_loss=ch20r$sqrt.ECA.std, #design matrix for standardized ECA
                start_y=N_s[,1],
                end_y=N_s[,2],
                start_t=L_i$tmin,
                end_t=L_i$tmax,
                pSmax_mean=0.5*smax_prior$m.s, #prior for Smax (spawners that maximize recruitment) based on max observed spawners
                pSmax_sig=smax_prior$m.s,
                pRk_mean=0.75*smax_prior$m.r, ##prior for Rk (recruitment capacity) based on max observed spawners
                pRk_sig=smax_prior$m.r)
                

print("eca")


if(Sys.info()[7] == "mariakur") {
  print("Running on local machine")
  cush_chm_eca <- mcush$sample(data=dl_chm_eca,
                            chains = 1, 
                            iter_warmup = 20,
                            iter_sampling =50,
                            refresh = 10,
                            adapt_delta = 0.999,
                            max_treedepth = 20)
  write.csv(cush_chm_eca$summary(),'./stan models/outs/summary/cush_chm_eca_trial.csv')
  cush_chm_eca$save_object('./stan models/outs/fits/cush_chm_eca_trial.RDS')
  
  post_cush_chm_eca=cush_chm_eca$draws(variables=c('b_for','b_for_cu','b_for_rv','alpha_t','alpha_j','b','sigma'),format='draws_matrix')
  write.csv(post_cush_chm_eca,here('stan models','outs','posterior','cush_chm_eca_trial.csv'))
  
} else {
  cush_chm_eca <- mcush$sample(data=dl_chm_eca,
                            chains = 6, 
                            iter_warmup = 200,
                            iter_sampling =500,
                            refresh = 100,
                            adapt_delta = 0.999,
                            max_treedepth = 20)
  
  write.csv(cush_chm_eca$summary(),here("stan models","outs","summary","cush_chm_eca_ac.csv"))
  cush_chm_eca$save_object(here("stan models","outs","fits","cush_chm_eca_ac.RDS"))
  
  post_cush_chm_eca=cush_chm_eca$draws(variables=c('b_for','b_for_cu','b_for_rv','alpha_t','alpha_j','b','sigma'),format='draws_matrix')
  write.csv(post_cush_chm_eca,here('stan models','outs','posterior','cush_chm_eca_ac.csv'))
  
}
