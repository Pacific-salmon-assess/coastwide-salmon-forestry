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


if(Sys.info()[7] == "mariakur") {
  bh_pk_eca <- readRDS(here('stan models','outs','fits','bh_pk_eca_noac_w_npgo_trial.RDS'))
  post_bh_pk_eca=bh_pk_eca$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                             'b_npgo','b_npgo_cu','b_npgo_rv',
                                             'alpha_j','Rk','sigma_for_cu','sigma_for_rv'),format='draws_matrix')
  write.csv(post_bh_pk_eca,here('stan models','outs','posterior','bh_pk_eca_noac_w_npgo_trial.csv'))
  
  #repeat for cpd
  
  bh_pk_cpd <- readRDS(here('stan models','outs','fits','bh_pk_cpd_noac_w_npgo_trial.RDS'))
  post_bh_pk_cpd=bh_pk_cpd$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                             'b_npgo','b_npgo_cu','b_npgo_rv',
                                             'alpha_j','Rk','sigma_for_cu','sigma_for_rv'),format='draws_matrix')
  write.csv(post_bh_pk_cpd,here('stan models','outs','posterior','bh_pk_cpd_noac_w_npgo_trial.csv'))
  
  
  
}else{
  bh_pk_eca <- readRDS(here('stan models','outs','fits','bh_pk_eca_noac_w_npgo.RDS'))
  post_bh_pk_eca=bh_pk_eca$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                             'b_npgo','b_npgo_cu','b_npgo_rv',
                                             'alpha_j','Rk','sigma_for_cu','sigma_for_rv'),format='draws_matrix')
  write.csv(post_bh_pk_eca,here('stan models','outs','posterior','bh_pk_eca_noac_w_npgo.csv'))
  
  
  bh_pk_cpd <- readRDS(here('stan models','outs','fits','bh_pk_cpd_noac_w_npgo.RDS'))
  post_bh_pk_cpd=bh_pk_cpd$draws(variables=c('b_for','b_for_cu','b_for_rv',
                                             'b_npgo','b_npgo_cu','b_npgo_rv',
                                             'alpha_j','Rk','sigma_for_cu','sigma_for_rv'),format='draws_matrix')
  write.csv(post_bh_pk_cpd,here('stan models','outs','posterior','bh_pk_cpd_noac_w_npgo.csv'))
  
  
  
}