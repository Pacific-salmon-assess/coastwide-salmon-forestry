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
cmdstanr::set_cmdstan_path(path = NULL)

#basic model excluding watershed areas:

# load Stan model sets####

#random effects for ECA at the CU but not River level (fit set 2) - main current results
file1ric=file.path(cmdstanr::cmdstan_path(),'sr models', "ric_hier_eca_mod_cu.stan")
m1ric=cmdstanr::cmdstan_model(file1ric) #compile stan code to C++

file1bh=file.path(cmdstanr::cmdstan_path(),'sr models', "bh_hier_eca_mod_cu.stan")
m1bh=cmdstanr::cmdstan_model(file1bh) #compile stan code to C++

file1cs=file.path(cmdstanr::cmdstan_path(),'sr models', "cush_hier_eca_mod_cu.stan")
m1cs=cmdstanr::cmdstan_model(file1cs) #compile stan code to C++

#random effects for ECA at River level - these probably shouldn't be used...
file2ric=file.path(cmdstanr::cmdstan_path(),'sr models', "ric_hier_eca_mod_rv.stan")
m2ric=cmdstanr::cmdstan_model(file2ric) #compile stan code to C++

file2bh=file.path(cmdstanr::cmdstan_path(),'sr models', "bh_hier_eca_mod_rv.stan")
m2bh=cmdstanr::cmdstan_model(file2bh) #compile stan code to C++

#file2cs=file.path(cmdstanr::cmdstan_path(),'sr models', "bh_hier_eca_mod_rv.stan")
#m2cs=cmdstanr::cmdstan_model(file2cs) #compile stan code to C++

#non linear effect of ECA using GAM splines
file2bh_gam=file.path(cmdstanr::cmdstan_path(),'sr models', "bh_spline_eca_mod_cu.stan")
m2bh_gam=cmdstanr::cmdstan_model(file2bh_gam) #compile stan code to C++

#watershed area interaction models (fit set 3)
file3=file.path(cmdstanr::cmdstan_path(),'sr models', "ric_hier_eca_mod_ExA.stan")
m3ric=cmdstanr::cmdstan_model(file3) #compile stan code to C++

file3bh=file.path(cmdstanr::cmdstan_path(),'sr models', "bh_hier_eca_mod_ExA.stan")
m3bh=cmdstanr::cmdstan_model(file3bh) #compile stan code to C++

file3bh2=file.path(cmdstanr::cmdstan_path(),'sr models', "bh_hier_eca_mod_ExA2.stan")
m3bh2=cmdstanr::cmdstan_model(file3bh2) #compile stan code to C++

file3bh3=file.path(cmdstanr::cmdstan_path(),'sr models', "bh_hier_eca_mod_ExA3.stan")
m3bh3=cmdstanr::cmdstan_model(file3bh3) #compile stan code to C++

#latent multivariate productivity trends (fit set 4)
file4ric=file.path(cmdstanr::cmdstan_path(),'sr models', "ric_rw_prod_eca_mod.stan")
m4ric=cmdstanr::cmdstan_model(file4ric) #compile stan code to C++

file4bh=file.path(cmdstanr::cmdstan_path(),'sr models', "bh_rw_prod_eca_mod2.stan")
m4bh=cmdstanr::cmdstan_model(file4bh) #compile stan code to C++

#null models - no ECA effect - to extract residuals (fit set 5)
file5bh=file.path(cmdstanr::cmdstan_path(),'sr models', "bh_hier_null_mod_cu.stan")
m5bh=cmdstanr::cmdstan_model(file5bh) #compile stan code to C++

file5ric=file.path(cmdstanr::cmdstan_path(),'sr models', "ric_hier_null_mod_cu.stan")
m5ric=cmdstanr::cmdstan_model(file5ric) #compile stan code to C++


#chum salmon####

##data formatting####

#censure implausible values (extremely high productivity)
ch20r<- subset(ch20r,exp(ln_RS)<=100)

#two rivers with duplicated names:
ch20r$River=ifelse(ch20r$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',ch20r$River)
ch20r$River=ifelse(ch20r$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',ch20r$River)


ch20r=ch20r[order(factor(ch20r$River),ch20r$BroodYear),]
rownames(ch20r)=seq(1:nrow(ch20r))

#normalize ECA - logit transformation (ie. log(x/(1-x)))
ch20r$logit.ECA=qlogis(ch20r$ECA_age_proxy_forested_only+0.005)
ch20r$logit.ECA.std=(ch20r$logit.ECA-mean(ch20r$logit.ECA))/sd(ch20r$logit.ECA)

#normalize cumulative % disturbed
ch20r$disturbedarea_prct_cs2=ifelse(ch20r$disturbedarea_prct_cs<1,1,ch20r$disturbedarea_prct_cs)
ch20r$disturbedarea_prct_cs2=ifelse(ch20r$disturbedarea_prct_cs2>99,99,ch20r$disturbedarea_prct_cs2)

ch20r$logit.pdisturb=qlogis(ch20r$disturbedarea_prct_cs2/100)
ch20r$logit.pdisturb.std=(ch20r$logit.pdisturb-mean(ch20r$logit.pdisturb))/sd(ch20r$logit.pdisturb)

#sparse matrix of spawners
S_mat=make_design_matrix(ch20r$Spawners,ch20r$River)

#sparse matrix of ECA

#average ECA by stock
#just to see an overview of ECA by river
eca_s=ch20r%>%group_by(River)%>%summarize(m=mean(ECA_age_proxy_forested_only*100),m.std=mean(logit.ECA.std),range=max(ECA_age_proxy_forested_only*100)-min(ECA_age_proxy_forested_only*100),cu=unique(CU))

ECA_mat=make_design_matrix(ch20r$logit.ECA.std,ch20r$River)

disturb_mat=make_design_matrix(ch20r$logit.pdisturb.std,ch20r$River)

#watershed area by river
area_mat=make_design_matrix(ch20r$ln_area_km2_std,ch20r$River)

#interaction matrix
ExA_mat=ECA_mat*area_mat

#extract max S for priors on capacity & eq. recruitment
smax_prior=ch20r%>%group_by(River) %>%summarize(m.s=max(Spawners),m.r=max(Recruits))

#ragged start and end points for each SR series
N_s=rag_n(ch20r$River)
  
#cus by stock
cu=distinct(ch20r,River,.keep_all = T)
summary(factor(cu$CU))

#time points for each series
L_i=ch20r%>%group_by(River)%>%summarize(l=n(),tmin=min(BroodYear)-1954+1,tmax=max(BroodYear)-1954+1)

CU_year=expand.grid(levels(factor(ch20r$CU)),seq(min(ch20r$BroodYear),max(ch20r$BroodYear)))
CU_year= CU_year[order(CU_year[,1]),]
CU_year[,3]=paste(CU_year[,1],CU_year[,2],sep="_")

ch20r$cu_yr=match(paste(ch20r$CU,ch20r$BroodYear,sep='_'),CU_year[,3])

#data list for fits 1-3 - ECA
dl_chm1=list(N=nrow(ch20r),
        L=max(ch20r$BroodYear)-min(ch20r$BroodYear)+1,
        C=length(unique(ch20r$CU)),
        J=length(unique(ch20r$River)),
        N_i=L_i$l,#series lengths
        C_i=as.numeric(factor(cu$CU)), #CU index by stock
        C_j=as.numeric(factor(ch20r$CU)), #CU index by observation
        C_ii=as.numeric(factor(ch20r$cu_yr)), #CU index by observation
        J_i=as.numeric(factor(ch20r$River)), #River index by observation
        J_ii=ch20r$cu_yr, #index specific to each unique CU-Year combination
        ii=as.numeric(factor(ch20r$BroodYear)), #brood year index
        R_S=ch20r$ln_RS,
        S=S_mat, #design matrix for spawner counts
        Sv=ch20r$Spawners,
        ECA=ECA_mat, #design matrix for standardized ECA
   #     ECA_vec=as.matrix(ch20r$logit.ECA.std), #vector of ECA ()
        Area=area_mat, #design matrix for watershed area
        ExA=ExA_mat, #design matrix for std ECA x watershed area
        start_y=N_s[,1],
        end_y=N_s[,2],
        start_t=L_i$tmin,
        end_t=L_i$tmax,
        pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
        pSmax_sig=smax_prior$m.s,
        pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
        pRk_sig=smax_prior$m.r)

#data list for fits 4 - ECA
dl_chm2=list(N=nrow(ch20r),
             L=max(ch20r$BroodYear)-min(ch20r$BroodYear)+1,
             C=length(unique(ch20r$CU)),
             J=length(unique(ch20r$River)),
             N_i=L_i$l,#series lengths
             C_i=as.numeric(factor(cu$CU)), #CU index by stock
             C_j=as.numeric(factor(ch20r$CU)), #CU index by observation
             C_ii=as.numeric(factor(ch20r$cu_yr)), #CU index by observation
             J_i=as.numeric(factor(ch20r$River)), #River index by observation
             J_ii=ch20r$cu_yr, #index specific to each unique CU-Year combination
             ii=as.numeric(factor(ch20r$BroodYear)), #brood year index
             R_S=ch20r$ln_RS,
             S=ch20r$Spawners, 
             ECA=ch20r$logit.ECA.std, #design matrix for standardized ECA
             #     ECA_vec=as.matrix(ch20r$logit.ECA.std), #vector of ECA ()
             #    Area=area_mat, #design matrix for watershed area
             ExA=ExA_mat, #design matrix for std ECA x watershed area
             start_y=N_s[,1],
             end_y=N_s[,2],
             start_t=L_i$tmin,
             end_t=L_i$tmax,
             pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
             pSmax_sig=smax_prior$m.s,
             pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
             pRk_sig=smax_prior$m.r)

#non-linear effects
bspline_ECA=t(splines::bs(ch20r$ECA_age_proxy_forested_only_std,knots=seq(-2,3,1),degree=3,intercept=F))
pred_ECA=t(splines::bs(seq(min(ch20r$ECA_age_proxy_forested_only_std),max(ch20r$ECA_age_proxy_forested_only_std),length.out=100),knots=seq(-2,3,1),degree=3,intercept=F))

dl_chm3=list(N=nrow(ch20r),
             L=max(ch20r$BroodYear)-min(ch20r$BroodYear)+1,
             C=length(unique(ch20r$CU)),
             J=length(unique(ch20r$River)),
             N_i=L_i$l,#series lengths
             C_i=as.numeric(factor(cu$CU)), #CU index
             ii=as.numeric(factor(ch20r$BroodYear)), #brood year index
             R_S=ch20r$ln_RS,
             S=S_mat, #design matrix for spawner counts
             ECA=ECA_mat, #design matrix for standardized ECA
             Area=area_mat, #design matrix for watershed area
             ExA=ExA_mat, #design matrix for std ECA x watershed area
             ECA_vec=ch20r$logit.ECA.std,
             B_ECA=bspline_ECA,
             Nb=nrow(bspline_ECA),
             start_y=N_s[,1],
             end_y=N_s[,2],
             start_t=L_i$tmin,
             end_t=L_i$tmax,
             pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
             pSmax_sig=smax_prior$m.s,
             pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
             pRk_sig=smax_prior$m.r,
             pred_ECA=pred_ECA,
             Np=ncol(pred_ECA))


##ECA predictor ####

### Beverton-Holt model sets ####


fit1bh_chm_eca <- m1bh$sample(data=dl_chm1,
                      seed = 12345,
                      chains = 8, 
                      iter_warmup = 200,
                      iter_sampling = 800,
                      refresh = 100,
                      adapt_delta = 0.995,
                      max_treedepth = 20)

write.csv(fit1bh_chm_eca$summary(),'./stan models/outs/summary/fit1bh_eca_summary.csv')
fit1bh_chm_eca$save_object('./stan models/outs/fits/fit1bh_chm_eca.RDS')




fit2bh_chm_eca <- m2bh$sample(data=dl_chm1,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit2bh_chm_eca$summary(),'./stan models/outs/summary/fit2bh_chm_eca_summary.csv')
fit2bh$save_object('stan models/outs/fits/fit2bh_chm_eca.RDS')


fit3bh_chm_eca <- m3bh$sample(data=dl_chm1,
                              seed = 1235,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit3bh$summary(),'./stan models/outs/summary/fit3bh_summary.csv')
fit3bh$save_object('./stan models/outs/fits/fit3bh.RDS')

fit4bh_chm_eca <- m4bh$sample(data=dl_chm2,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit4bh_eca$summary(),'./stan models/outs/summary/fit4bh_eca_summary.csv')
fit4bh_eca$save_object('./stan models/outs/fits/fit4bh_eca.RDS')


fit_bh_eca_gam <- m2bh_gam$sample(data=dl_chm3,
                                  seed = 1235,
                                  chains = 8, 
                                  iter_warmup = 200,
                                  iter_sampling = 800,
                                  refresh = 100,
                                  adapt_delta = 0.995,
                                  max_treedepth = 20)


write.csv(fit_bh_chm_eca_gam$summary(),'./stan models/outs/summary/fit_bh_eca_gam_summary.csv')
fit_bh_chm_eca_gam$save_object('./stan models/outs/fits/fit_bh_eca_gam.RDS')

### Ricker model sets####

fit1ric_chm_eca <- m1ric$sample(data=dl_chm1,
                                seed = 33456,
                                chains = 8, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.999,
                                max_treedepth = 20)

write.csv(fit1ric_chm_eca$summary(),'./stan models/outs/summary/fit1ric_chm_eca_summary.csv')
fit1ric_chm_eca$save_object('./stan models/outs/fits/fit1ric_chm_eca.RDS')
fit1ric_chm_eca=readRDS('stan models/outs/fits/fit1ric_chm_eca.RDS')


fit2ric_chm_eca <- m2ric$sample(data=dl_chm1,
                  seed = 12345,
                  chains = 8, 
                  iter_warmup = 200,
                  iter_sampling = 800,
                  refresh = 100,
                  adapt_delta = 0.995,
                  max_treedepth = 20)

write.csv(fit2ric_chm_eca$summary(),'./stan models/outs/summary/fit2ric_chm_eca_summary.csv')
fit2ric_chm_eca$save_object('stan models/outs/fits/fit2ric_chm_eca.RDS')

fit3ric_chm_eca <- m3ric$sample(data=dl_chm1,
                                seed = 33366,
                                chains = 8, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.995,
                                max_treedepth = 20)

write.csv(fit3ric$summary(),'./stan models/outs/summary/fit3ric_chm_eca_summary.csv')
fit3ric$save_object('./stan models/outs/fits/fit3ric_chm_eca.RDS')


fit4ric_chm_eca <- m4ric$sample(data=dl_chm2,
                                seed = 333,
                                chains = 4, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.995,
                                max_treedepth = 20)

write.csv(fit4ric_chm_eca$summary(),'./stan models/outs/summary/fit4ric_chm_eca_summary.csv')
fit4ric_chm_eca$save_object('./stan models/outs/fits/fit4ric_chm_eca.RDS')


###Cushing model ###

fit1cs_chm_eca <- m1cs$sample(data=dl_chm1,
                              seed = 1234,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.999,
                              max_treedepth = 20)

write.csv(fit1cs_chm_eca$summary(),'./stan models/outs/summary/fit1cs_summary.csv')
fit1cs_chm_eca$save_object('./stan models/outs/fits/fit1cs_chm_eca.RDS')


## Cumulative percent disturbed####

#Eca now set to cumulative % disturbed metric

dl_chm_4=list(N=nrow(ch20r),
        L=max(ch20r$BroodYear)-min(ch20r$BroodYear)+1,
        C=length(unique(ch20r$CU)),
        J=length(unique(ch20r$River)),
        N_i=L_i$l,#series lengths
        C_i=as.numeric(factor(cu$CU)), #CU index by stock
        C_j=as.numeric(factor(ch20r$CU)), #CU index by observation
        C_ii=as.numeric(factor(ch20r$cu_yr)), #CU index by observation
        J_i=as.numeric(factor(ch20r$River)), #River index by observation
        J_ii=ch20r$cu_yr, #index specific to each unique CU-Year combination
        ii=as.numeric(factor(ch20r$BroodYear)), #brood year index
        R_S=ch20r$ln_RS,
        S=S_mat, #design matrix for spawner counts
        Sv=ch20r$Spawners,
        ECA=disturb_mat, #design matrix for standardized ECA
        #     ECA_vec=as.matrix(ch20r$logit.ECA.std), #vector of ECA ()
        Area=area_mat, #design matrix for watershed area
        ExA=ExA_mat, #design matrix for std ECA x watershed area
        start_y=N_s[,1],
        end_y=N_s[,2],
        start_t=L_i$tmin,
        end_t=L_i$tmax,
        pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
        pSmax_sig=smax_prior$m.s,
        pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
        pRk_sig=smax_prior$m.r)

dl_chm_5=list(N=nrow(ch20r),
              L=max(ch20r$BroodYear)-min(ch20r$BroodYear)+1,
              C=length(unique(ch20r$CU)),
              J=length(unique(ch20r$River)),
              N_i=L_i$l,#series lengths
              C_i=as.numeric(factor(cu$CU)), #CU index by stock
              C_j=as.numeric(factor(ch20r$CU)), #CU index by observation
              C_ii=as.numeric(factor(ch20r$cu_yr)), #CU index by observation
              J_i=as.numeric(factor(ch20r$River)), #River index by observation
              J_ii=ch20r$cu_yr, #index specific to each unique CU-Year combination
              ii=as.numeric(factor(ch20r$BroodYear)), #brood year index
              R_S=ch20r$ln_RS,
              S=ch20r$Spawners, 
              ECA=ch20r$logit.ECA.std, #design matrix for standardized ECA
              #     ECA_vec=as.matrix(ch20r$logit.ECA.std), #vector of ECA ()
              #    Area=area_mat, #design matrix for watershed area
              ExA=ExA_mat, #design matrix for std ECA x watershed area
              start_y=N_s[,1],
              end_y=N_s[,2],
              start_t=L_i$tmin,
              end_t=L_i$tmax,
              pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
              pSmax_sig=smax_prior$m.s,
              pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
              pRk_sig=smax_prior$m.r)


### Beverton Holt model fits ####
fit1bh_chm_cpd <- m1bh$sample(data=dl_chm_4,
                      seed = 12345,
                      chains = 8, 
                      iter_warmup = 200,
                      iter_sampling = 800,
                      refresh = 100,
                      adapt_delta = 0.999,
                      max_treedepth = 20)

write.csv(fit1bh_chm_cpd$summary(),'./stan models/outs/summary/fit1bh_chm_cpd_summary.csv')
fit1bh_cpd$save_object('./stan models/outs/fits/fit1bh_chm_cpd.RDS')

fit2bh_chm_cpd <- m2bh$sample(data=dl_chm_4,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit2bh_cpd$summary(),'./stan models/outs/summary/fit2bh_chm_cpd_summary.csv')
fit2bh_cpd$save_object('stan models/outs/fits/fit2bh_chm_cpd.RDS')

fit3bh_chm_cpd <- m3bh$sample(data=dl_chm_4,
                              seed = 1235,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit3bh_chm_cpd$summary(),'./stan models/outs/summary/fit3bh_chm_cpd.csv')
fit3bh_chm_cpd$save_object('./stan models/outs/fits/fit3bh_chm_cpd.RDS')

fit4bh_chm_cpd <- m4bh$sample(data=dl_chm_5,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit4bh$summary(),'./stan models/outs/summary/fit4bh_chm_cpd_summary.csv')
fit4bh$save_object('./stan models/outs/fits/fit4bh_chm_cpd.RDS')


### Ricker model fits ####
fit1ric_chm_pd <- m1ric$sample(data=dl_chm_4,
                               seed = 33456,
                               chains = 8, 
                               iter_warmup = 200,
                               iter_sampling = 800,
                               refresh = 100,
                               adapt_delta = 0.995,
                               max_treedepth = 20)

write.csv(fit1ric$summary(),'./stan models/outs/summary/fit1ric_chm_cpd_summary.csv')
fit1ric_chm_pd$save_object('./stan models/outs/fits/fit1ric_chm_pd.RDS')

fit2ric_chm_cpd <- m2ric$sample(data=dl_chm_4,
                        seed = 12345,
                        chains = 8, 
                        iter_warmup = 200,
                        iter_sampling = 800,
                        refresh = 100,
                        adapt_delta = 0.995,
                        max_treedepth = 20)

write.csv(fit2ric$summary(),'./stan models/outs/summary/fit2ric_chm_cpd_summary.csv')
fit2ric$save_object('stan models/outs/fits/fit2ric_chm_cpd.RDS')

fit3ric_chm_cpd <- m3ric$sample(data=dl_chm_4,
                                seed = 33366,
                                chains = 8, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.999,
                                max_treedepth = 20)

write.csv(fit3ric$summary(),'./stan models/outs/summary/fit3ric_chm_cpd_summary.csv')
fit3ric_chm_cpd$save_object('./stan models/outs/fits/fit3ric_chm_cpd.RDS')

fit4ric_chm_cpd <- m4ric$sample(data=dl_chm_5,
                        seed = 333,
                        chains = 4, 
                        iter_warmup = 200,
                        iter_sampling = 800,
                        refresh = 100,
                        adapt_delta = 0.995,
                        max_treedepth = 20)

write.csv(fit4ric_chm_cpd$summary(),'./stan models/outs/summary/fit4ric_chm_cpd_summary.csv')
fit4ric_chm_cpd$save_object('./stan models/outs/fits/fit4ric_chm_cpd.RDS')


#Even Pinks####

##data formatting####

#censure implausible values (extremely high productivity)
pk10r_e<- subset(pk10r_e,exp(ln_RS)<=80)

#two rivers with duplicated names:
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',pk10r_e$River)
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',pk10r_e$River)


pk10r_e=pk10r_e[order(factor(pk10r_e$River),pk10r_e$BroodYear),]
rownames(pk10r_e)=seq(1:nrow(pk10r_e))

#normalize ECA - logit transformation (ie. log(x/(1-x)))
pk10r_e$logit.ECA=qlogis(pk10r_e$ECA_age_proxy_forested_only+0.005)
pk10r_e$logit.ECA.std=(pk10r_e$logit.ECA-mean(pk10r_e$logit.ECA))/sd(pk10r_e$logit.ECA)

#normalize cumulative % disturbed
pk10r_e$disturbedarea_prct_cs2=ifelse(pk10r_e$disturbedarea_prct_cs<1,1,pk10r_e$disturbedarea_prct_cs)
pk10r_e$disturbedarea_prct_cs2=ifelse(pk10r_e$disturbedarea_prct_cs2>99,99,pk10r_e$disturbedarea_prct_cs2)

pk10r_e$logit.pdisturb=qlogis(pk10r_e$disturbedarea_prct_cs2/100)
pk10r_e$logit.pdisturb.std=(pk10r_e$logit.pdisturb-mean(pk10r_e$logit.pdisturb))/sd(pk10r_e$logit.pdisturb)

#sparse matrix of spawners
S_mat=make_design_matrix(pk10r_e$Spawners,pk10r_e$River)

#sparse matrix of ECA

#average ECA by stock
#just to see an overview of ECA by river
eca_s=pk10r_e%>%group_by(River)%>%summarize(m=mean(ECA_age_proxy_forested_only*100),m.std=mean(logit.ECA.std),range=max(ECA_age_proxy_forested_only*100)-min(ECA_age_proxy_forested_only*100),cu=unique(CU))

ECA_mat=make_design_matrix(pk10r_e$logit.ECA.std,pk10r_e$River)

disturb_mat=make_design_matrix(pk10r_e$logit.pdisturb.std,pk10r_e$River)

#watershed area by river
area_mat=make_design_matrix(pk10r_e$ln_area_km2_std,pk10r_e$River)

#interaction matrix
ExA_mat=ECA_mat*area_mat

#extract max S for priors on capacity & eq. recruitment
smax_prior=pk10r_e%>%group_by(River) %>%summarize(m.s=max(Spawners),m.r=max(Recruits))

#ragged start and end points for each SR series
N_s=rag_n(pk10r_e$River)

#cus by stock
cu=distinct(pk10r_e,River,.keep_all = T)
summary(factor(cu$CU))

#time points for each series
L_i=pk10r_e%>%group_by(River)%>%summarize(l=n(),tmin=min(BroodYear)-1954+1,tmax=max(BroodYear)-1954+1)

CU_year=expand.grid(levels(factor(pk10r_e$CU)),seq(min(pk10r_e$BroodYear),max(pk10r_e$BroodYear)))
CU_year= CU_year[order(CU_year[,1]),]
CU_year[,3]=paste(CU_year[,1],CU_year[,2],sep="_")

pk10r_e$cu_yr=match(paste(pk10r_e$CU,pk10r_e$BroodYear,sep='_'),CU_year[,3])

#data list for fits 1-3 - ECA
dl_pke1=list(N=nrow(p10r_e),
             L=max(p10r_e$BroodYear)-min(p10r_e$BroodYear)+1,
             C=length(unique(p10r_e$CU)),
             J=length(unique(p10r_e$River)),
             N_i=L_i$l,#series lengths
             C_i=as.numeric(factor(cu$CU)), #CU index by stock
             C_j=as.numeric(factor(p10r_e$CU)), #CU index by observation
             C_ii=as.numeric(factor(p10r_e$cu_yr)), #CU index by observation
             J_i=as.numeric(factor(p10r_e$River)), #River index by observation
             J_ii=p10r_e$cu_yr, #index specific to each unique CU-Year combination
             ii=as.numeric(factor(p10r_e$BroodYear)), #brood year index
             R_S=p10r_e$ln_RS,
             S=S_mat, #design matrix for spawner counts
             Sv=p10r_e$Spawners,
             ECA=ECA_mat, #design matrix for standardized ECA
             #     ECA_vec=as.matrix(p10r_e$logit.ECA.std), #vector of ECA ()
             Area=area_mat, #design matrix for watershed area
             ExA=ExA_mat, #design matrix for std ECA x watershed area
             start_y=N_s[,1],
             end_y=N_s[,2],
             start_t=L_i$tmin,
             end_t=L_i$tmax,
             pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
             pSmax_sig=smax_prior$m.s,
             pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
             pRk_sig=smax_prior$m.r)

#data list for fits 4 - ECA
dl_pke2=list(N=nrow(p10r_e),
             L=max(p10r_e$BroodYear)-min(p10r_e$BroodYear)+1,
             C=length(unique(p10r_e$CU)),
             J=length(unique(p10r_e$River)),
             N_i=L_i$l,#series lengths
             C_i=as.numeric(factor(cu$CU)), #CU index by stock
             C_j=as.numeric(factor(p10r_e$CU)), #CU index by observation
             C_ii=as.numeric(factor(p10r_e$cu_yr)), #CU index by observation
             J_i=as.numeric(factor(p10r_e$River)), #River index by observation
             J_ii=p10r_e$cu_yr, #index specific to each unique CU-Year combination
             ii=as.numeric(factor(p10r_e$BroodYear)), #brood year index
             R_S=p10r_e$ln_RS,
             S=p10r_e$Spawners, 
             ECA=p10r_e$logit.ECA.std, #design matrix for standardized ECA
             #     ECA_vec=as.matrix(p10r_e$logit.ECA.std), #vector of ECA ()
             #    Area=area_mat, #design matrix for watershed area
             ExA=ExA_mat, #design matrix for std ECA x watershed area
             start_y=N_s[,1],
             end_y=N_s[,2],
             start_t=L_i$tmin,
             end_t=L_i$tmax,
             pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
             pSmax_sig=smax_prior$m.s,
             pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
             pRk_sig=smax_prior$m.r)

#non-linear effects
bspline_ECA=t(splines::bs(p10r_e$ECA_age_proxy_forested_only_std,knots=seq(-2,3,1),degree=3,intercept=F))
pred_ECA=t(splines::bs(seq(min(p10r_e$ECA_age_proxy_forested_only_std),max(p10r_e$ECA_age_proxy_forested_only_std),length.out=100),knots=seq(-2,3,1),degree=3,intercept=F))

dl_pke3=list(N=nrow(p10r_e),
             L=max(p10r_e$BroodYear)-min(p10r_e$BroodYear)+1,
             C=length(unique(p10r_e$CU)),
             J=length(unique(p10r_e$River)),
             N_i=L_i$l,#series lengths
             C_i=as.numeric(factor(cu$CU)), #CU index
             ii=as.numeric(factor(p10r_e$BroodYear)), #brood year index
             R_S=p10r_e$ln_RS,
             S=S_mat, #design matrix for spawner counts
             ECA=ECA_mat, #design matrix for standardized ECA
             Area=area_mat, #design matrix for watershed area
             ExA=ExA_mat, #design matrix for std ECA x watershed area
             ECA_vec=p10r_e$logit.ECA.std,
             B_ECA=bspline_ECA,
             Nb=nrow(bspline_ECA),
             start_y=N_s[,1],
             end_y=N_s[,2],
             start_t=L_i$tmin,
             end_t=L_i$tmax,
             pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
             pSmax_sig=smax_prior$m.s,
             pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
             pRk_sig=smax_prior$m.r,
             pred_ECA=pred_ECA,
             Np=ncol(pred_ECA))


##ECA predictor ####

### Beverton-Holt model sets ####


fit1bh_pke_eca <- m1bh$sample(data=dl_pke1,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit1bh_pke_eca$summary(),'./stan models/outs/summary/fit1bh_eca_summary.csv')
fit1bh_pke_eca$save_object('./stan models/outs/fits/fit1bh_pke_eca.RDS')




fit2bh_pke_eca <- m2bh$sample(data=dl_pke1,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit2bh_pke_eca$summary(),'./stan models/outs/summary/fit2bh_pke_eca_summary.csv')
fit2bh$save_object('stan models/outs/fits/fit2bh_pke_eca.RDS')


fit3bh_pke_eca <- m3bh$sample(data=dl_pke1,
                              seed = 1235,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit3bh$summary(),'./stan models/outs/summary/fit3bh_summary.csv')
fit3bh$save_object('./stan models/outs/fits/fit3bh.RDS')

fit4bh_pke_eca <- m4bh$sample(data=dl_pke2,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit4bh_eca$summary(),'./stan models/outs/summary/fit4bh_eca_summary.csv')
fit4bh_eca$save_object('./stan models/outs/fits/fit4bh_eca.RDS')


fit_bh_eca_gam <- m2bh_gam$sample(data=dl_pke3,
                                  seed = 1235,
                                  chains = 8, 
                                  iter_warmup = 200,
                                  iter_sampling = 800,
                                  refresh = 100,
                                  adapt_delta = 0.995,
                                  max_treedepth = 20)


write.csv(fit_bh_pke_eca_gam$summary(),'./stan models/outs/summary/fit_bh_eca_gam_summary.csv')
fit_bh_pke_eca_gam$save_object('./stan models/outs/fits/fit_bh_eca_gam.RDS')

### Ricker model sets####

fit1ric_pke_eca <- m1ric$sample(data=dl_pke1,
                                seed = 33456,
                                chains = 8, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.999,
                                max_treedepth = 20)

write.csv(fit1ric_pke_eca$summary(),'./stan models/outs/summary/fit1ric_pke_eca_summary.csv')
fit1ric_pke_eca$save_object('./stan models/outs/fits/fit1ric_pke_eca.RDS')
fit1ric_pke_eca=readRDS('stan models/outs/fits/fit1ric_pke_eca.RDS')


fit2ric_pke_eca <- m2ric$sample(data=dl_pke1,
                                seed = 12345,
                                chains = 8, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.995,
                                max_treedepth = 20)

write.csv(fit2ric_pke_eca$summary(),'./stan models/outs/summary/fit2ric_pke_eca_summary.csv')
fit2ric_pke_eca$save_object('stan models/outs/fits/fit2ric_pke_eca.RDS')

fit3ric_pke_eca <- m3ric$sample(data=dl_pke1,
                                seed = 33366,
                                chains = 8, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.995,
                                max_treedepth = 20)

write.csv(fit3ric$summary(),'./stan models/outs/summary/fit3ric_pke_eca_summary.csv')
fit3ric$save_object('./stan models/outs/fits/fit3ric_pke_eca.RDS')


fit4ric_pke_eca <- m4ric$sample(data=dl_pke2,
                                seed = 333,
                                chains = 4, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.995,
                                max_treedepth = 20)

write.csv(fit4ric_pke_eca$summary(),'./stan models/outs/summary/fit4ric_pke_eca_summary.csv')
fit4ric_pke_eca$save_object('./stan models/outs/fits/fit4ric_pke_eca.RDS')


###Cushing model ###

fit1cs_pke_eca <- m1cs$sample(data=dl_pke1,
                              seed = 1234,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.999,
                              max_treedepth = 20)

write.csv(fit1cs_pke_eca$summary(),'./stan models/outs/summary/fit1cs_summary.csv')
fit1cs_pke_eca$save_object('./stan models/outs/fits/fit1cs_pke_eca.RDS')


## Cumulative percent disturbed####

#Eca now set to cumulative % disturbed metric

dl_pke_4=list(N=nrow(p10r_e),
              L=max(p10r_e$BroodYear)-min(p10r_e$BroodYear)+1,
              C=length(unique(p10r_e$CU)),
              J=length(unique(p10r_e$River)),
              N_i=L_i$l,#series lengths
              C_i=as.numeric(factor(cu$CU)), #CU index by stock
              C_j=as.numeric(factor(p10r_e$CU)), #CU index by observation
              C_ii=as.numeric(factor(p10r_e$cu_yr)), #CU index by observation
              J_i=as.numeric(factor(p10r_e$River)), #River index by observation
              J_ii=p10r_e$cu_yr, #index specific to each unique CU-Year combination
              ii=as.numeric(factor(p10r_e$BroodYear)), #brood year index
              R_S=p10r_e$ln_RS,
              S=S_mat, #design matrix for spawner counts
              Sv=p10r_e$Spawners,
              ECA=disturb_mat, #design matrix for standardized ECA
              #     ECA_vec=as.matrix(p10r_e$logit.ECA.std), #vector of ECA ()
              Area=area_mat, #design matrix for watershed area
              ExA=ExA_mat, #design matrix for std ECA x watershed area
              start_y=N_s[,1],
              end_y=N_s[,2],
              start_t=L_i$tmin,
              end_t=L_i$tmax,
              pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
              pSmax_sig=smax_prior$m.s,
              pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
              pRk_sig=smax_prior$m.r)

dl_pke_5=list(N=nrow(p10r_e),
              L=max(p10r_e$BroodYear)-min(p10r_e$BroodYear)+1,
              C=length(unique(p10r_e$CU)),
              J=length(unique(p10r_e$River)),
              N_i=L_i$l,#series lengths
              C_i=as.numeric(factor(cu$CU)), #CU index by stock
              C_j=as.numeric(factor(p10r_e$CU)), #CU index by observation
              C_ii=as.numeric(factor(p10r_e$cu_yr)), #CU index by observation
              J_i=as.numeric(factor(p10r_e$River)), #River index by observation
              J_ii=p10r_e$cu_yr, #index specific to each unique CU-Year combination
              ii=as.numeric(factor(p10r_e$BroodYear)), #brood year index
              R_S=p10r_e$ln_RS,
              S=p10r_e$Spawners, 
              ECA=p10r_e$logit.ECA.std, #design matrix for standardized ECA
              #     ECA_vec=as.matrix(p10r_e$logit.ECA.std), #vector of ECA ()
              #    Area=area_mat, #design matrix for watershed area
              ExA=ExA_mat, #design matrix for std ECA x watershed area
              start_y=N_s[,1],
              end_y=N_s[,2],
              start_t=L_i$tmin,
              end_t=L_i$tmax,
              pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
              pSmax_sig=smax_prior$m.s,
              pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
              pRk_sig=smax_prior$m.r)


### Beverton Holt model fits ####
fit1bh_pke_cpd <- m1bh$sample(data=dl_pke_4,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.999,
                              max_treedepth = 20)

write.csv(fit1bh_pke_cpd$summary(),'./stan models/outs/summary/fit1bh_pke_cpd_summary.csv')
fit1bh_cpd$save_object('./stan models/outs/fits/fit1bh_pke_cpd.RDS')

fit2bh_pke_cpd <- m2bh$sample(data=dl_pke_4,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit2bh_cpd$summary(),'./stan models/outs/summary/fit2bh_pke_cpd_summary.csv')
fit2bh_cpd$save_object('stan models/outs/fits/fit2bh_pke_cpd.RDS')

fit3bh_pke_cpd <- m3bh$sample(data=dl_pke_4,
                              seed = 1235,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit3bh_pke_cpd$summary(),'./stan models/outs/summary/fit3bh_pke_cpd.csv')
fit3bh_pke_cpd$save_object('./stan models/outs/fits/fit3bh_pke_cpd.RDS')

fit4bh_pke_cpd <- m4bh$sample(data=dl_pke_5,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit4bh$summary(),'./stan models/outs/summary/fit4bh_pke_cpd_summary.csv')
fit4bh$save_object('./stan models/outs/fits/fit4bh_pke_cpd.RDS')


### Ricker model fits ####
fit1ric_pke_pd <- m1ric$sample(data=dl_pke_4,
                               seed = 33456,
                               chains = 8, 
                               iter_warmup = 200,
                               iter_sampling = 800,
                               refresh = 100,
                               adapt_delta = 0.995,
                               max_treedepth = 20)

write.csv(fit1ric$summary(),'./stan models/outs/summary/fit1ric_pke_cpd_summary.csv')
fit1ric_pke_pd$save_object('./stan models/outs/fits/fit1ric_pke_pd.RDS')

fit2ric_pke_cpd <- m2ric$sample(data=dl_pke_4,
                                seed = 12345,
                                chains = 8, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.995,
                                max_treedepth = 20)

write.csv(fit2ric$summary(),'./stan models/outs/summary/fit2ric_pke_cpd_summary.csv')
fit2ric$save_object('stan models/outs/fits/fit2ric_pke_cpd.RDS')

fit3ric_pke_cpd <- m3ric$sample(data=dl_pke_4,
                                seed = 33366,
                                chains = 8, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.999,
                                max_treedepth = 20)

write.csv(fit3ric$summary(),'./stan models/outs/summary/fit3ric_pke_cpd_summary.csv')
fit3ric_pke_cpd$save_object('./stan models/outs/fits/fit3ric_pke_cpd.RDS')

fit4ric_pke_cpd <- m4ric$sample(data=dl_pke_5,
                                seed = 333,
                                chains = 4, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.995,
                                max_treedepth = 20)

write.csv(fit4ric_pke_cpd$summary(),'./stan models/outs/summary/fit4ric_pke_cpd_summary.csv')
fit4ric_pke_cpd$save_object('./stan models/outs/fits/fit4ric_pke_cpd.RDS')

#Odd Pinks####

##data formatting####

#censure implausible values (extremely high productivity)
pk10r_o<- subset(pk10r_o,exp(ln_RS)<=100)

#two rivers with duplicated names:
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',pk10r_o$River)
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',pk10r_o$River)


pk10r_o=pk10r_o[order(factor(pk10r_o$River),pk10r_o$BroodYear),]
rownames(pk10r_o)=seq(1:nrow(pk10r_o))

#normalize ECA - logit transformation (ie. log(x/(1-x)))
pk10r_o$logit.ECA=qlogis(pk10r_o$ECA_age_proxy_forested_only+0.005)
pk10r_o$logit.ECA.std=(pk10r_o$logit.ECA-mean(pk10r_o$logit.ECA))/sd(pk10r_o$logit.ECA)

#normalize cumulative % disturbed
pk10r_o$disturbedarea_prct_cs2=ifelse(pk10r_o$disturbedarea_prct_cs<1,1,pk10r_o$disturbedarea_prct_cs)
pk10r_o$disturbedarea_prct_cs2=ifelse(pk10r_o$disturbedarea_prct_cs2>99,99,pk10r_o$disturbedarea_prct_cs2)

pk10r_o$logit.pdisturb=qlogis(pk10r_o$disturbedarea_prct_cs2/100)
pk10r_o$logit.pdisturb.std=(pk10r_o$logit.pdisturb-mean(pk10r_o$logit.pdisturb))/sd(pk10r_o$logit.pdisturb)

#sparse matrix of spawners
S_mat=make_design_matrix(pk10r_o$Spawners,pk10r_o$River)

#sparse matrix of ECA

#average ECA by stock
#just to see an overview of ECA by river
eca_s=pk10r_o%>%group_by(River)%>%summarize(m=mean(ECA_age_proxy_forested_only*100),m.std=mean(logit.ECA.std),range=max(ECA_age_proxy_forested_only*100)-min(ECA_age_proxy_forested_only*100),cu=unique(CU))

ECA_mat=make_design_matrix(pk10r_o$logit.ECA.std,pk10r_o$River)

disturb_mat=make_design_matrix(pk10r_o$logit.pdisturb.std,pk10r_o$River)

#watershed area by river
area_mat=make_design_matrix(pk10r_o$ln_area_km2_std,pk10r_o$River)

#interaction matrix
ExA_mat=ECA_mat*area_mat

#extract max S for priors on capacity & eq. recruitment
smax_prior=pk10r_o%>%group_by(River) %>%summarize(m.s=max(Spawners),m.r=max(Recruits))

#ragged start and end points for each SR series
N_s=rag_n(pk10r_o$River)

#cus by stock
cu=distinct(pk10r_o,River,.keep_all = T)
summary(factor(cu$CU))

#time points for each series
L_i=pk10r_o%>%group_by(River)%>%summarize(l=n(),tmin=min(BroodYear)-1954+1,tmax=max(BroodYear)-1954+1)

CU_year=expand.grid(levels(factor(pk10r_o$CU)),seq(min(pk10r_o$BroodYear),max(pk10r_o$BroodYear)))
CU_year= CU_year[order(CU_year[,1]),]
CU_year[,3]=paste(CU_year[,1],CU_year[,2],sep="_")

pk10r_o$cu_yr=match(paste(pk10r_o$CU,pk10r_o$BroodYear,sep='_'),CU_year[,3])

#data list for fits 1-3 - ECA
dl_pko1=list(N=nrow(pk10r_o),
             L=max(pk10r_o$BroodYear)-min(pk10r_o$BroodYear)+1,
             C=length(unique(pk10r_o$CU)),
             J=length(unique(pk10r_o$River)),
             N_i=L_i$l,#series lengths
             C_i=as.numeric(factor(cu$CU)), #CU index by stock
             C_j=as.numeric(factor(pk10r_o$CU)), #CU index by observation
             C_ii=as.numeric(factor(pk10r_o$cu_yr)), #CU index by observation
             J_i=as.numeric(factor(pk10r_o$River)), #River index by observation
             J_ii=pk10r_o$cu_yr, #index specific to each unique CU-Year combination
             ii=as.numeric(factor(pk10r_o$BroodYear)), #brood year index
             R_S=pk10r_o$ln_RS,
             S=S_mat, #design matrix for spawner counts
             Sv=pk10r_o$Spawners,
             ECA=ECA_mat, #design matrix for standardized ECA
             #     ECA_vec=as.matrix(pk10r_o$logit.ECA.std), #vector of ECA ()
             Area=area_mat, #design matrix for watershed area
             ExA=ExA_mat, #design matrix for std ECA x watershed area
             start_y=N_s[,1],
             end_y=N_s[,2],
             start_t=L_i$tmin,
             end_t=L_i$tmax,
             pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
             pSmax_sig=smax_prior$m.s,
             pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
             pRk_sig=smax_prior$m.r)

#data list for fits 4 - ECA
dl_pko2=list(N=nrow(pk10r_o),
             L=max(pk10r_o$BroodYear)-min(pk10r_o$BroodYear)+1,
             C=length(unique(pk10r_o$CU)),
             J=length(unique(pk10r_o$River)),
             N_i=L_i$l,#series lengths
             C_i=as.numeric(factor(cu$CU)), #CU index by stock
             C_j=as.numeric(factor(pk10r_o$CU)), #CU index by observation
             C_ii=as.numeric(factor(pk10r_o$cu_yr)), #CU index by observation
             J_i=as.numeric(factor(pk10r_o$River)), #River index by observation
             J_ii=pk10r_o$cu_yr, #index specific to each unique CU-Year combination
             ii=as.numeric(factor(pk10r_o$BroodYear)), #brood year index
             R_S=pk10r_o$ln_RS,
             S=pk10r_o$Spawners, 
             ECA=pk10r_o$logit.ECA.std, #design matrix for standardized ECA
             #     ECA_vec=as.matrix(pk10r_o$logit.ECA.std), #vector of ECA ()
             #    Area=area_mat, #design matrix for watershed area
             ExA=ExA_mat, #design matrix for std ECA x watershed area
             start_y=N_s[,1],
             end_y=N_s[,2],
             start_t=L_i$tmin,
             end_t=L_i$tmax,
             pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
             pSmax_sig=smax_prior$m.s,
             pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
             pRk_sig=smax_prior$m.r)

#non-linear effects
bspline_ECA=t(splines::bs(pk10r_o$ECA_age_proxy_forested_only_std,knots=seq(-2,3,1),degree=3,intercept=F))
pred_ECA=t(splines::bs(seq(min(pk10r_o$ECA_age_proxy_forested_only_std),max(pk10r_o$ECA_age_proxy_forested_only_std),length.out=100),knots=seq(-2,3,1),degree=3,intercept=F))

dl_pko3=list(N=nrow(pk10r_o),
             L=max(pk10r_o$BroodYear)-min(pk10r_o$BroodYear)+1,
             C=length(unique(pk10r_o$CU)),
             J=length(unique(pk10r_o$River)),
             N_i=L_i$l,#series lengths
             C_i=as.numeric(factor(cu$CU)), #CU index
             ii=as.numeric(factor(pk10r_o$BroodYear)), #brood year index
             R_S=pk10r_o$ln_RS,
             S=S_mat, #design matrix for spawner counts
             ECA=ECA_mat, #design matrix for standardized ECA
             Area=area_mat, #design matrix for watershed area
             ExA=ExA_mat, #design matrix for std ECA x watershed area
             ECA_vec=pk10r_o$logit.ECA.std,
             B_ECA=bspline_ECA,
             Nb=nrow(bspline_ECA),
             start_y=N_s[,1],
             end_y=N_s[,2],
             start_t=L_i$tmin,
             end_t=L_i$tmax,
             pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
             pSmax_sig=smax_prior$m.s,
             pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
             pRk_sig=smax_prior$m.r,
             pred_ECA=pred_ECA,
             Np=ncol(pred_ECA))


##ECA predictor ####

### Beverton-Holt model sets ####


fit1bh_pko_eca <- m1bh$sample(data=dl_pko1,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit1bh_pko_eca$summary(),'./stan models/outs/summary/fit1bh_eca_summary.csv')
fit1bh_pko_eca$save_object('./stan models/outs/fits/fit1bh_pko_eca.RDS')




fit2bh_pko_eca <- m2bh$sample(data=dl_pko1,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit2bh_pko_eca$summary(),'./stan models/outs/summary/fit2bh_pko_eca_summary.csv')
fit2bh$save_object('stan models/outs/fits/fit2bh_pko_eca.RDS')


fit3bh_pko_eca <- m3bh$sample(data=dl_pko1,
                              seed = 1235,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit3bh$summary(),'./stan models/outs/summary/fit3bh_summary.csv')
fit3bh$save_object('./stan models/outs/fits/fit3bh.RDS')

fit4bh_pko_eca <- m4bh$sample(data=dl_pko2,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit4bh_eca$summary(),'./stan models/outs/summary/fit4bh_eca_summary.csv')
fit4bh_eca$save_object('./stan models/outs/fits/fit4bh_eca.RDS')


fit_bh_eca_gam <- m2bh_gam$sample(data=dl_pko3,
                                  seed = 1235,
                                  chains = 8, 
                                  iter_warmup = 200,
                                  iter_sampling = 800,
                                  refresh = 100,
                                  adapt_delta = 0.995,
                                  max_treedepth = 20)


write.csv(fit_bh_pko_eca_gam$summary(),'./stan models/outs/summary/fit_bh_eca_gam_summary.csv')
fit_bh_pko_eca_gam$save_object('./stan models/outs/fits/fit_bh_eca_gam.RDS')

### Ricker model sets####

fit1ric_pko_eca <- m1ric$sample(data=dl_pko1,
                                seed = 33456,
                                chains = 8, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.999,
                                max_treedepth = 20)

write.csv(fit1ric_pko_eca$summary(),'./stan models/outs/summary/fit1ric_pko_eca_summary.csv')
fit1ric_pko_eca$save_object('./stan models/outs/fits/fit1ric_pko_eca.RDS')
fit1ric_pko_eca=readRDS('stan models/outs/fits/fit1ric_pko_eca.RDS')


fit2ric_pko_eca <- m2ric$sample(data=dl_pko1,
                                seed = 12345,
                                chains = 8, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.995,
                                max_treedepth = 20)

write.csv(fit2ric_pko_eca$summary(),'./stan models/outs/summary/fit2ric_pko_eca_summary.csv')
fit2ric_pko_eca$save_object('stan models/outs/fits/fit2ric_pko_eca.RDS')

fit3ric_pko_eca <- m3ric$sample(data=dl_pko1,
                                seed = 33366,
                                chains = 8, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.995,
                                max_treedepth = 20)

write.csv(fit3ric$summary(),'./stan models/outs/summary/fit3ric_pko_eca_summary.csv')
fit3ric$save_object('./stan models/outs/fits/fit3ric_pko_eca.RDS')


fit4ric_pko_eca <- m4ric$sample(data=dl_pko2,
                                seed = 333,
                                chains = 4, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.995,
                                max_treedepth = 20)

write.csv(fit4ric_pko_eca$summary(),'./stan models/outs/summary/fit4ric_pko_eca_summary.csv')
fit4ric_pko_eca$save_object('./stan models/outs/fits/fit4ric_pko_eca.RDS')


###Cushing model ###

fit1cs_pko_eca <- m1cs$sample(data=dl_pko1,
                              seed = 1234,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.999,
                              max_treedepth = 20)

write.csv(fit1cs_pko_eca$summary(),'./stan models/outs/summary/fit1cs_summary.csv')
fit1cs_pko_eca$save_object('./stan models/outs/fits/fit1cs_pko_eca.RDS')


## Cumulative percent disturbed####

#Eca now set to cumulative % disturbed metric

dl_pko_4=list(N=nrow(pk10r_o),
              L=max(pk10r_o$BroodYear)-min(pk10r_o$BroodYear)+1,
              C=length(unique(pk10r_o$CU)),
              J=length(unique(pk10r_o$River)),
              N_i=L_i$l,#series lengths
              C_i=as.numeric(factor(cu$CU)), #CU index by stock
              C_j=as.numeric(factor(pk10r_o$CU)), #CU index by observation
              C_ii=as.numeric(factor(pk10r_o$cu_yr)), #CU index by observation
              J_i=as.numeric(factor(pk10r_o$River)), #River index by observation
              J_ii=pk10r_o$cu_yr, #index specific to each unique CU-Year combination
              ii=as.numeric(factor(pk10r_o$BroodYear)), #brood year index
              R_S=pk10r_o$ln_RS,
              S=S_mat, #design matrix for spawner counts
              Sv=pk10r_o$Spawners,
              ECA=disturb_mat, #design matrix for standardized ECA
              #     ECA_vec=as.matrix(pk10r_o$logit.ECA.std), #vector of ECA ()
              Area=area_mat, #design matrix for watershed area
              ExA=ExA_mat, #design matrix for std ECA x watershed area
              start_y=N_s[,1],
              end_y=N_s[,2],
              start_t=L_i$tmin,
              end_t=L_i$tmax,
              pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
              pSmax_sig=smax_prior$m.s,
              pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
              pRk_sig=smax_prior$m.r)

dl_pko_5=list(N=nrow(pk10r_o),
              L=max(pk10r_o$BroodYear)-min(pk10r_o$BroodYear)+1,
              C=length(unique(pk10r_o$CU)),
              J=length(unique(pk10r_o$River)),
              N_i=L_i$l,#series lengths
              C_i=as.numeric(factor(cu$CU)), #CU index by stock
              C_j=as.numeric(factor(pk10r_o$CU)), #CU index by observation
              C_ii=as.numeric(factor(pk10r_o$cu_yr)), #CU index by observation
              J_i=as.numeric(factor(pk10r_o$River)), #River index by observation
              J_ii=pk10r_o$cu_yr, #index specific to each unique CU-Year combination
              ii=as.numeric(factor(pk10r_o$BroodYear)), #brood year index
              R_S=pk10r_o$ln_RS,
              S=pk10r_o$Spawners, 
              ECA=pk10r_o$logit.ECA.std, #design matrix for standardized ECA
              #     ECA_vec=as.matrix(pk10r_o$logit.ECA.std), #vector of ECA ()
              #    Area=area_mat, #design matrix for watershed area
              ExA=ExA_mat, #design matrix for std ECA x watershed area
              start_y=N_s[,1],
              end_y=N_s[,2],
              start_t=L_i$tmin,
              end_t=L_i$tmax,
              pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
              pSmax_sig=smax_prior$m.s,
              pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
              pRk_sig=smax_prior$m.r)


### Beverton Holt model fits ####
fit1bh_pko_cpd <- m1bh$sample(data=dl_pko_4,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.999,
                              max_treedepth = 20)

write.csv(fit1bh_pko_cpd$summary(),'./stan models/outs/summary/fit1bh_pko_cpd_summary.csv')
fit1bh_cpd$save_object('./stan models/outs/fits/fit1bh_pko_cpd.RDS')

fit2bh_pko_cpd <- m2bh$sample(data=dl_pko_4,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit2bh_cpd$summary(),'./stan models/outs/summary/fit2bh_pko_cpd_summary.csv')
fit2bh_cpd$save_object('stan models/outs/fits/fit2bh_pko_cpd.RDS')

fit3bh_pko_cpd <- m3bh$sample(data=dl_pko_4,
                              seed = 1235,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit3bh_pko_cpd$summary(),'./stan models/outs/summary/fit3bh_pko_cpd.csv')
fit3bh_pko_cpd$save_object('./stan models/outs/fits/fit3bh_pko_cpd.RDS')

fit4bh_pko_cpd <- m4bh$sample(data=dl_pko_5,
                              seed = 12345,
                              chains = 8, 
                              iter_warmup = 200,
                              iter_sampling = 800,
                              refresh = 100,
                              adapt_delta = 0.995,
                              max_treedepth = 20)

write.csv(fit4bh$summary(),'./stan models/outs/summary/fit4bh_pko_cpd_summary.csv')
fit4bh$save_object('./stan models/outs/fits/fit4bh_pko_cpd.RDS')


### Ricker model fits ####
fit1ric_pko_pd <- m1ric$sample(data=dl_pko_4,
                               seed = 33456,
                               chains = 8, 
                               iter_warmup = 200,
                               iter_sampling = 800,
                               refresh = 100,
                               adapt_delta = 0.995,
                               max_treedepth = 20)

write.csv(fit1ric$summary(),'./stan models/outs/summary/fit1ric_pko_cpd_summary.csv')
fit1ric_pko_pd$save_object('./stan models/outs/fits/fit1ric_pko_pd.RDS')

fit2ric_pko_cpd <- m2ric$sample(data=dl_pko_4,
                                seed = 12345,
                                chains = 8, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.995,
                                max_treedepth = 20)

write.csv(fit2ric$summary(),'./stan models/outs/summary/fit2ric_pko_cpd_summary.csv')
fit2ric$save_object('stan models/outs/fits/fit2ric_pko_cpd.RDS')

fit3ric_pko_cpd <- m3ric$sample(data=dl_pko_4,
                                seed = 33366,
                                chains = 8, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.999,
                                max_treedepth = 20)

write.csv(fit3ric$summary(),'./stan models/outs/summary/fit3ric_pko_cpd_summary.csv')
fit3ric_pko_cpd$save_object('./stan models/outs/fits/fit3ric_pko_cpd.RDS')

fit4ric_pko_cpd <- m4ric$sample(data=dl_pko_5,
                                seed = 333,
                                chains = 4, 
                                iter_warmup = 200,
                                iter_sampling = 800,
                                refresh = 100,
                                adapt_delta = 0.995,
                                max_treedepth = 20)

write.csv(fit4ric_pko_cpd$summary(),'./stan models/outs/summary/fit4ric_pko_cpd_summary.csv')
fit4ric_pko_cpd$save_object('./stan models/outs/fits/fit4ric_pko_cpd.RDS')
