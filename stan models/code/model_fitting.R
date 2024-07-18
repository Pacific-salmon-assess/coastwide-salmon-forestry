library(here);library(cmdstanr);library(dplyr)
source('./stan models/code/funcs.R')

# load datasets####

ch20r <- read.csv("./origional-ecofish-data-models/Data/Processed/chum_SR_20_hat_yr_reduced_VRI90.csv")

#even year pinks
pk10r_e <- read.csv("./origional-ecofish-data-models/Data/Processed/pke_SR_10_hat_yr_reduced_VRI90.csv")

#odd year pinks
pk10r_o <- read.csv("./origional-ecofish-data-models/Data/Processed/pko_SR_10_hat_yr_reduced_VRI90.csv")

options(mc.cores=8)

#Note: you will need to set your cmdstanr path to a folder with hier_eca_mod.stan file, using the next line
cmdstanr::set_cmdstan_path(path='C:/Users/greenbergda/Documents/.cmdstan/cmdstan-2.29.2')

#basic model excluding watershed areas:

# load Stan model sets####
file_bh=file.path(cmdstanr::cmdstan_path(),'sr models', "bh_rw_prod_eca_mod_rv2.stan")
mbh=cmdstanr::cmdstan_model(file_bh) #compile stan code to C++

file_bh=file.path(cmdstanr::cmdstan_path(),'sr models', "bh_rw_prod_eca_mod_ac_rv_pink.stan")
mbh_p=cmdstanr::cmdstan_model(file_bh) #compile stan code to C++


file_ric=file.path(cmdstanr::cmdstan_path(),'sr models', "ric_rw_prod_eca_mod_ac_rv2.stan")
mric=cmdstanr::cmdstan_model(file_ric) #compile stan code to C++

#Chum salmon####

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

dl_chm_cpd=list(N=nrow(ch20r),
              L=max(ch20r$BroodYear)-min(ch20r$BroodYear)+1,
              C=length(unique(ch20r$CU)),
              J=length(unique(ch20r$River)),
              C_i=as.numeric(factor(cu$CU)), #CU index by stock
              ii=as.numeric(factor(ch20r$BroodYear)), #brood year index
              R_S=ch20r$ln_RS,
              S=ch20r$Spawners, 
              forest_loss=as.vector(ch20r$sqrt.CPD.std), #design matrix for standardized ECA
              ECA=as.vector(ch20r$sqrt.CPD.std), #design matrix for standardized ECA
              start_y=N_s[,1],
              end_y=N_s[,2],
              start_t=L_i$tmin,
              end_t=L_i$tmax,
              pSmax_mean=0.5*smax_prior$m.s, #prior for Smax (spawners that maximize recruitment) based on max observed spawners
              pSmax_sig=smax_prior$m.s,
              pRk_mean=0.75*smax_prior$m.r, #prior for Rk (recruitment capacity) based on max observed spawners
              pRk_sig=smax_prior$m.r)


## ECA predictor ####
bh_chm_eca <- mbh$sample(data=dl_chm_eca,
                                        chains = 6, 
                                        init=0,
                                        iter_warmup = 200,
                                        iter_sampling =500,
                                        refresh = 100,
                                        adapt_delta = 0.999,
                                        max_treedepth = 20)

write.csv(bh_chm_eca$summary(),'./stan models/outs/summary/bh_chm_eca.csv')
bh_chm_eca$save_object('./stan models/outs/fits/bh_chm_eca.RDS')

post_bh_chm_eca=bh_chm_eca$draws(variables=c('b_for','b_for_cu','b_for_rv','alpha_t','alpha_j','Rk'),format='draws_matrix')
write.csv(post_bh_chm_eca,here('stan models','outs','fits','posterior','bh_chm_eca.csv'))

ric_chm_eca <- mric$sample(data=dl_chm_eca,
                               chains = 6, 
                               init=0,
                               iter_warmup = 200,
                               iter_sampling =500,
                               refresh = 100,
                               adapt_delta = 0.999,
                               max_treedepth = 20)

write.csv(ric_chm_eca$summary(),'./stan models/outs/summary/ric_chm_eca.csv')
ric_chm_eca$save_object('./stan models/outs/fits/ric_chm_eca.RDS')

post_ric_chm_eca=fit4ricac_chm_eca$draws(variables=c('b_for','b_for_cu','b_for_rv','alpha_t','alpha_j','Rk'),format='draws_matrix')
write.csv(dfit4eca,here('stan models','outs','fits','posterior','ric_chm_eca.csv'))

## CPD predictor ####
bh_chm_cpd <- mbh$sample(data=dl_chm_cpd,
                               chains = 8, 
                               init=0,
                               iter_warmup = 200,
                               iter_sampling =500,
                               refresh = 100,
                               adapt_delta = 0.999,
                               max_treedepth = 20)

write.csv(bh_chm_cpd$summary(),'./stan models/outs/summary/bh_chm_cpd.csv')
bh_chm_cpd$save_object('./stan models/outs/fits/bh_chm_cpd.RDS')

post_chm_cpd=bh_chm_cpd$draws(variables=c('b_for','b_for_cu','b_for_rv','alpha_t','alpha_j','Rk'),format='draws_matrix')
write.csv(post_chm_cpd,here('stan models','outs','fits','posterior','bh_chm_cpd.csv'))


ric_chm_cpd <- mric$sample(data=dl_chm_cpd,
                           chains = 6, 
                           init=0,
                           iter_warmup = 200,
                           iter_sampling =500,
                           refresh = 100,
                           adapt_delta = 0.999,
                           max_treedepth = 20)

write.csv(ric_chm_cpd$summary(),'./stan models/outs/summary/ric_chm_cpd.csv')
ric_chm_cpd$save_object('./stan models/outs/fits/ric_chm_cpd.RDS')

post_ric_chm_cpd=fit4ricac_chm_cpd$draws(variables=c('b_for','b_for_cu','b_for_rv','alpha_t','alpha_j','Rk'),format='draws_matrix')
write.csv(dfit4cpd,here('stan models','outs','fits','posterior','ric_chm_cpd.csv'))

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
                forest_loss=pk10r$sqrt.ECA.std, 
                start_y=N_s[,1],
                end_y=N_s[,2],
                start_t=L_all$tmin,
                pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
                pSmax_sig=smax_prior$m.s,
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
               forest_loss=pk10r$sqrt.CPD.std, 
               start_y=N_s[,1],
               end_y=N_s[,2],
               start_t=L_all$tmin,
               pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
               pSmax_sig=smax_prior$m.s,
               pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
               pRk_sig=smax_prior$m.r)

bh_pk_eca <- mbh_p$sample(data=dl_pk_eca,
                         chains = 6, 
                         iter_warmup = 200,
                         iter_sampling =500,
                         refresh = 100,
                         adapt_delta = 0.999,
                         max_treedepth = 20)

write.csv(bh_pk_eca$summary(),'./stan models/outs/summary/bh_pk_eca.csv')
bh_pk_eca$save_object('./stan models/outs/fits/bh_pk_eca.RDS')

post_bh_pk_eca=bh_pk_eca$draws(variables=c('b_for','b_for_cu','b_for_rv','alpha_t','alpha_j','Rk'),format='draws_matrix')
write.csv(post_bh_pk_eca,here('stan models','outs','fits','posterior','bh_pke_eca.csv'))

bh_pk_cpd <- mbh_p$sample(data=dl_pk_cpd,
                          chains = 6, 
                          iter_warmup = 200,
                          iter_sampling =500,
                          refresh = 100,
                          adapt_delta = 0.999,
                          max_treedepth = 20)

write.csv(bh_pk_cpd$summary(),'./stan models/outs/summary/bh_pk_cpd.csv')
bh_pk_cpd$save_object('./stan models/outs/fits/bh_pk_cpd.RDS')

post_bh_pk_cpd=bh_pk_cpd$draws(variables=c('b_for','b_for_cu','b_for_rv','alpha_t','alpha_j','Rk'),format='draws_matrix')
write.csv(post_bh_pk_cpd,here('stan models','outs','fits','posterior','bh_pke_cpd.csv'))


# Pink salmon - even broodlines ####

##data formatting####
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',pk10r_e$River)
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='915-765500-18600-00000-0000-0000-000-000-000-000-000-000','HEAD CREEK 2',pk10r_e$River)
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='915-488000-41400-00000-0000-0000-000-000-000-000-000-000','WINDY BAY CREEK 2',pk10r_e$River)
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',pk10r_e$River)
pk10r_e=pk10r_e[order(factor(pk10r_e$River),pk10r_e$BroodYear),]
rownames(pk10r_e)=seq(1:nrow(pk10r_e))

#normalize ECA 2 - square root transformation (ie. sqrt(x))
pk10r_e$sqrt.ECA=sqrt(pk10r_e$ECA_age_proxy_forested_only)
pk10r_e$sqrt.ECA.std=(pk10r_e$sqrt.ECA-mean(pk10r_e$sqrt.ECA))/sd(pk10r_e$sqrt.ECA)

#normalize CPD 2 - square root transformation (ie. sqrt(x))
pk10r_e$sqrt.CPD=sqrt(pk10r_e$disturbedarea_prct_cs)
pk10r_e$sqrt.CPD.std=(pk10r_e$sqrt.CPD-mean(pk10r_e$sqrt.CPD))/sd(pk10r_e$sqrt.CPD)

#extract max S for priors on capacity & eq. recruitment
smax_prior=pk10r_e%>%group_by(River) %>%summarize(m.s=max(Spawners),m.r=max(Recruits))

#ragged start and end points for each SR series
N_s=rag_n(pk10r_e$River)

#cus by stock
cu=distinct(pk10r_e,River,.keep_all = T)
summary(factor(cu$CU))

#time points for each series
L_i=pk10r_e%>%group_by(River)%>%summarize(l=n(),by=min(BroodYear),tmin=(min(BroodYear)-min(pk10r_e$BroodYear))/2+1,tmax=(max(BroodYear)-min(pk10r_e$BroodYear))/2)


dl_pke_eca=list(N=nrow(pk10r_e),
             L=length(unique(pk10r_e$BroodYear)),
             C=length(unique(pk10r_e$CU)),
             J=length(unique(pk10r_e$River)),
             C_i=as.numeric(factor(cu$CU)), #CU index by stock
             ii=as.numeric(factor(pk10r_e$BroodYear)), #brood year index
             R_S=pk10r_e$ln_RS,
             S=pk10r_e$Spawners, 
             forest_loss=pk10r_e$sqrt.ECA.std, 
             start_y=N_s[,1],
             end_y=N_s[,2],
             start_t=L_i$tmin,
             pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
             pSmax_sig=smax_prior$m.s,
             pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
             pRk_sig=smax_prior$m.r)

dl_pke_cpd=list(N=nrow(pk10r_e),
                L=length(unique(pk10r_e$BroodYear)),
                C=length(unique(pk10r_e$CU)),
                J=length(unique(pk10r_e$River)),
                C_i=as.numeric(factor(cu$CU)), #CU index by stock
                ii=as.numeric(factor(pk10r_e$BroodYear)), #brood year index
                R_S=pk10r_e$ln_RS,
                S=pk10r_e$Spawners, 
                forest_loss=as.vector(pk10r_e$sqrt.CPD.std), #design matrix for standardized ECA
                start_y=N_s[,1],
                end_y=N_s[,2],
                start_t=L_i$tmin,
                end_t=L_i$tmax,
                pSmax_mean=0.5*smax_prior$m.s, #prior for Smax (spawners that maximize recruitment) based on max observed spawners
                pSmax_sig=smax_prior$m.s,
                pRk_mean=0.75*smax_prior$m.r, #prior for Rk (recruitment capacity) based on max observed spawners
                pRk_sig=smax_prior$m.r)

## ECA predictor ####
bh_pke_eca <- mbh$sample(data=dl_pke_eca,
                         chains = 6, 
                         init=0,
                         iter_warmup = 200,
                         iter_sampling =500,
                         refresh = 100,
                         adapt_delta = 0.999,
                         max_treedepth = 20)

write.csv(bh_pke_eca$summary(),'./stan models/outs/summary/bh_pke_eca.csv')
bh_pke_eca$save_object('./stan models/outs/fits/bh_pke_eca.RDS')

post_bh_pke_eca=bh_pke_eca$draws(variables=c('b_for','b_for_cu','b_for_rv','alpha_t','alpha_j','Rk'),format='draws_matrix')
write.csv(post_bh_pke_eca,here('stan models','outs','fits','posterior','bh_pke_eca.csv'))

ric_chm_eca <- mric$sample(data=dl_chm_eca,
                           chains = 6, 
                           init=0,
                           iter_warmup = 200,
                           iter_sampling =500,
                           refresh = 100,
                           adapt_delta = 0.999,
                           max_treedepth = 20)

write.csv(ric_chm_eca$summary(),'./stan models/outs/summary/ric_chm_eca.csv')
ric_chm_eca$save_object('./stan models/outs/fits/ric_chm_eca.RDS')

post_ric_chm_eca=fit4ricac_chm_eca$draws(variables=c('b_for','b_for_cu','b_for_rv','alpha_t','alpha_j','Rk'),format='draws_matrix')
write.csv(dfit4eca,here('stan models','outs','fits','posterior','ric_chm_eca.csv'))

## CPD predictor ####
bh_pke_cpd <- mbh$sample(data=dl_pke_cpd,
                         chains = 6, 
                         init=0,
                         iter_warmup = 200,
                         iter_sampling =500,
                         refresh = 100,
                         adapt_delta = 0.999,
                         max_treedepth = 20)

write.csv(bh_pke_cpd$summary(),'./stan models/outs/summary/bh_pke_cpd.csv')
bh_pke_cpd$save_object('./stan models/outs/fits/bh_pke_cpd.RDS')

post_pke_cpd=bh_pke_cpd$draws(variables=c('b_for','b_for_cu','b_for_rv','alpha_t','alpha_j','Rk'),format='draws_matrix')
write.csv(post_pke_cpd,here('stan models','outs','fits','posterior','bh_pke_cpd.csv'))

# Pink salmon - odd broodlines ####

##data formatting####
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',pk10r_o$River)
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=='915-765500-18600-00000-0000-0000-000-000-000-000-000-000','HEAD CREEK 2',pk10r_o$River)
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=='915-488000-41400-00000-0000-0000-000-000-000-000-000-000','WINDY BAY CREEK 2',pk10r_o$River)
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',pk10r_o$River)
pk10r_o=pk10r_o[order(factor(pk10r_o$River),pk10r_o$BroodYear),]
rownames(pk10r_o)=seq(1:nrow(pk10r_o))

#normalize ECA 2 - square root transformation (ie. sqrt(x))
pk10r_o$sqrt.ECA=sqrt(pk10r_o$ECA_age_proxy_forested_only)
pk10r_o$sqrt.ECA.std=(pk10r_o$sqrt.ECA-mean(pk10r_o$sqrt.ECA))/sd(pk10r_o$sqrt.ECA)

#normalize CPD 2 - square root transformation (ie. sqrt(x))
pk10r_o$sqrt.CPD=sqrt(pk10r_o$disturbedarea_prct_cs)
pk10r_o$sqrt.CPD.std=(pk10r_o$sqrt.CPD-mean(pk10r_o$sqrt.CPD))/sd(pk10r_o$sqrt.CPD)

#extract max S for priors on capacity & eq. recruitment
smax_prior=pk10r_o%>%group_by(River) %>%summarize(m.s=max(Spawners),m.r=max(Recruits))

#ragged start and end points for each SR series
N_s=rag_n(pk10r_o$River)

#cus by stock
cu=distinct(pk10r_o,River,.keep_all = T)
summary(factor(cu$CU))

#time points for each series
L_i=pk10r_o%>%group_by(River)%>%summarize(l=n(),by=min(BroodYear),tmin=(min(BroodYear)-min(pk10r_o$BroodYear))/2+1,tmax=(max(BroodYear)-min(pk10r_o$BroodYear))/2)

dl_pko_eca=list(N=nrow(pk10r_o),
                L=length(unique(pk10r_o$BroodYear)),
                C=length(unique(pk10r_o$CU)),
                J=length(unique(pk10r_o$River)),
                C_i=as.numeric(factor(cu$CU)), #CU index by stock
                ii=as.numeric(factor(pk10r_o$BroodYear)), #brood year index
                R_S=pk10r_o$ln_RS,
                S=pk10r_o$Spawners, 
                forest_loss=pk10r_o$sqrt.ECA.std, 
                start_y=N_s[,1],
                end_y=N_s[,2],
                start_t=L_i$tmin,
                pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
                pSmax_sig=smax_prior$m.s,
                pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
                pRk_sig=smax_prior$m.r)

dl_pko_cpd=list(N=nrow(pk10r_o),
                L=length(unique(pk10r_o$BroodYear)),
                C=length(unique(pk10r_o$CU)),
                J=length(unique(pk10r_o$River)),
                C_i=as.numeric(factor(cu$CU)), #CU index by stock
                ii=as.numeric(factor(pk10r_o$BroodYear)), #brood year index
                R_S=pk10r_o$ln_RS,
                S=pk10r_o$Spawners, 
                forest_loss=as.vector(pk10r_o$sqrt.CPD.std), #design matrix for standardized ECA
                start_y=N_s[,1],
                end_y=N_s[,2],
                start_t=L_i$tmin,
                end_t=L_i$tmax,
                pSmax_mean=0.5*smax_prior$m.s, #prior for Smax (spawners that maximize recruitment) based on max observed spawners
                pSmax_sig=smax_prior$m.s,
                pRk_mean=0.75*smax_prior$m.r, #prior for Rk (recruitment capacity) based on max observed spawners
                pRk_sig=smax_prior$m.r)

## ECA predictor ####
bh_pko_eca <- mbh$sample(data=dl_pko_eca,
                         chains = 6, 
                         init=0,
                         iter_warmup = 200,
                         iter_sampling =500,
                         refresh = 100,
                         adapt_delta = 0.999,
                         max_treedepth = 20)

write.csv(bh_pko_eca$summary(),'./stan models/outs/summary/bh_pko_eca.csv')
bh_pko_eca$save_object('./stan models/outs/fits/bh_pko_eca.RDS')

post_bh_pko_eca=bh_pko_eca$draws(variables=c('b_for','b_for_cu','b_for_rv','alpha_t','alpha_j','Rk'),format='draws_matrix')
write.csv(post_bh_pko_eca,here('stan models','outs','fits','posterior','bh_pko_eca.csv'))

  ric_pko_eca <- mric$sample(data=dl_pko_eca,
                           chains = 6, 
                           init=0,
                           iter_warmup = 200,
                           iter_sampling =500,
                           refresh = 100,
                           adapt_delta = 0.999,
                           max_treedepth = 20)

write.csv(ric_pko_eca$summary(),'./stan models/outs/summary/ric_pko_eca.csv')
ric_pko_eca$save_object('./stan models/outs/fits/ric_pko_eca.RDS')

post_ric_pko_eca=ric_pko_eca$draws(variables=c('b_for','b_for_cu','b_for_rv','alpha_t','alpha_j','b'),format='draws_matrix')
write.csv(ric_pko_eca,here('stan models','outs','fits','posterior','ric_pko_eca.csv'))

## CPD predictor ####
bh_pko_cpd <- mbh$sample(data=dl_pko_cpd,
                         chains = 6, 
                         init=0,
                         iter_warmup = 200,
                         iter_sampling =500,
                         refresh = 100,
                         adapt_delta = 0.999,
                         max_treedepth = 20)

write.csv(bh_pko_cpd$summary(),'./stan models/outs/summary/bh_pko_cpd.csv')
bh_pko_cpd$save_object('./stan models/outs/fits/bh_pko_cpd.RDS')

post_bh_pko_cpd=bh_pko_cpd$draws(variables=c('b_for','b_for_cu','b_for_rv','alpha_t','alpha_j','Rk'),format='draws_matrix')
write.csv(post_bh_pko_cpd,here('stan models','outs','fits','posterior','bh_pko_cpd.csv'))

ric_pko_cpd <- mric$sample(data=dl_pko_cpd,
                           chains = 6, 
                           init=0,
                           iter_warmup = 200,
                           iter_sampling =500,
                           refresh = 100,
                           adapt_delta = 0.999,
                           max_treedepth = 20)

write.csv(ric_pko_cpd$summary(),'./stan models/outs/summary/ric_pko_cpd.csv')
ric_pko_cpd$save_object('./stan models/outs/fits/ric_pko_cpd.RDS')

post_ric_pko_cpd=ric_pko_cpd$draws(variables=c('b_for','b_for_cu','b_for_rv','alpha_t','alpha_j','b'),format='draws_matrix')
write.csv(ric_pko_cpd,here('stan models','outs','fits','posterior','ric_pko_cpd.csv'))

#Chum salmon + south coast####

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

dl_chm_cpd=list(N=nrow(ch20r),
                L=max(ch20r$BroodYear)-min(ch20r$BroodYear)+1,
                C=length(unique(ch20r$CU)),
                J=length(unique(ch20r$River)),
                C_i=as.numeric(factor(cu$CU)), #CU index by stock
                ii=as.numeric(factor(ch20r$BroodYear)), #brood year index
                R_S=ch20r$ln_RS,
                S=ch20r$Spawners, 
                forest_loss=as.vector(ch20r$sqrt.CPD.std), #design matrix for standardized ECA
                start_y=N_s[,1],
                end_y=N_s[,2],
                start_t=L_i$tmin,
                end_t=L_i$tmax,
                pSmax_mean=0.5*smax_prior$m.s, #prior for Smax (spawners that maximize recruitment) based on max observed spawners
                pSmax_sig=smax_prior$m.s,
                pRk_mean=0.75*smax_prior$m.r, #prior for Rk (recruitment capacity) based on max observed spawners
                pRk_sig=smax_prior$m.r)


