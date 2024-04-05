library(here);library(cmdstanr);library(dplyr)
source('./stan models/code/funcs.R')

# load datasets####

ch20r <- read.csv("./origional-ecofish-data-models/Data/Processed/chum_SR_20_hat_yr_reduced_VRI90.csv")
ch20 <- read.csv("./origional-ecofish-data-models/Data/Processed/chum_SR_20.csv")

pk10r <- read.csv("./origional-ecofish-data-models/Data/Processed/pke_SR_10_hat_yr_reduced_VRI90.csv")
pk10 <- read.csv("./origional-ecofish-data-models/Data/Processed/pke_SR_10.csv")

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

file2cs=file.path(cmdstanr::cmdstan_path(),'sr models', "bh_hier_eca_mod_rv.stan")
m2cs=cmdstanr::cmdstan_model(file2cs) #compile stan code to C++

#Still needs some tweaks - non linear effect of ECA
#file2bh_gp=file.path(cmdstanr::cmdstan_path(),'sr models', "bh_gp_eca_mod.stan")
#m2bhgp=cmdstanr::cmdstan_model(file2bh_gp) #compile stan code to C++

#watershed area interaction models (fit set 3)
file3=file.path(cmdstanr::cmdstan_path(),'sr models', "ric_hier_eca_mod_ExA.stan")
m3ric=cmdstanr::cmdstan_model(file3) #compile stan code to C++

file3bh=file.path(cmdstanr::cmdstan_path(),'sr models', "bh_hier_eca_mod_ExA.stan")
m3bh=cmdstanr::cmdstan_model(file3bh) #compile stan code to C++

#latent multivariate productivity trends (fit set 4)
file4ric=file.path(cmdstanr::cmdstan_path(),'sr models', "ric_hier_eca_mod_rw_prod.stan")
m4ric=cmdstanr::cmdstan_model(file3) #compile stan code to C++

file3bh=file.path(cmdstanr::cmdstan_path(),'sr models', "bh_hier_eca_mod_rw_prod.stan")
m3bh=cmdstanr::cmdstan_model(file3bh) #compile stan code to C++



#data formatting####
ch20r=ch20r[order(factor(ch20r$River),ch20r$BroodYear),]
rownames(ch20r)=seq(1:nrow(ch20r))

#normalize ECA - logit transformation (ie. log(x/(1-x)))
ch20r$logit.ECA=qlogis(ch20r$ECA_age_proxy_forested_only+0.005)
ch20r$logit.ECA.std=(ch20r$logit.ECA-mean(ch20r$logit.ECA))/sd(ch20r$logit.ECA)

#sparse matrix of spawners
S_mat=make_design_matrix(ch20r$Spawners,ch20r$River)

#sparse matrix of ECA

#average ECA by stock
#just to see an overview of ECA by river
eca_s=ch20r%>%group_by(River)%>%summarize(m=mean(ECA_age_proxy_forested_only*100),m.std=mean(logit.ECA.std),range=max(ECA_age_proxy_forested_only*100)-min(ECA_age_proxy_forested_only*100),cu=unique(CU))

ECA_mat=make_design_matrix(ch20r$logit.ECA.std,ch20r$River)

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


dl=list(N=nrow(ch20r),
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
        start_y=N_s[,1],
        end_y=N_s[,2],
        start_t=L_i$tmin,
        end_t=L_i$tmax,
        pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
        pSmax_sig=smax_prior$m.s,
        pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
        pRk_sig=smax_prior$m.r)


# stan model fits####
## Fit set 1: CU-level varying effects for ECA effects ####


fit1ric <- m2$sample(data=dl,
                  seed = 123,
                  chains = 8, 
                  iter_warmup = 200,
                  iter_sampling = 800,
                  refresh = 100,
                  adapt_delta = 0.995,
                  max_treedepth = 20)

write.csv(fit1ric$summary(),'./stan models/outs/summary/fit1ric_summary.csv')
fit1ric$save_object('./stan models/outs/fits/fit1ric.RDS')
fit1ric=readRDS('stan models/outs/fits/fit1ric.RDS')
loo_f2=fit1ric$loo()

#residual plots to check
d2=fit2$draws(variables=c('e_t','mu1','mu2'),format='draws_matrix')
e_t=apply(d2[,grepl('e_t',colnames(d2))],2,median)
mu1=apply(d2[,grepl('mu1',colnames(d2))],2,median)
mu2=apply(d2[,grepl('mu2',colnames(d2))],2,median)
plot(e_t~mu1)
#plot(e_t~mu2)


fit1cs <- m2cs$sample(data=dl,
                      seed = 12345,
                      chains = 8, 
                      iter_warmup = 200,
                      iter_sampling = 800,
                      refresh = 100,
                      adapt_delta = 0.995,
                      max_treedepth = 20)

write.csv(fit1cs$summary(),'./stan models/outs/summary/fit1cs_summary.csv')
fit1cs$save_object('./stan models/outs/fits/fit1cs.RDS')
fit1cs <- readRDS('./stan models/outs/fits/fit1cs.RDS')
loo_f1cs=fit1cs$loo()


d2c=fit1cs$draws(variables=c('e_t','mu1','mu2'),format='draws_matrix')
e_t=apply(d2c[,grepl('e_t',colnames(d2))],2,median)
mu1=apply(d2c[,grepl('mu1',colnames(d2))],2,median)
mu2=apply(d2c[,grepl('mu2',colnames(d2))],2,median)
plot(e_t~mu1)
#plot(e_t~mu2)

fit1bh <- m1bh$sample(data=dl,
                  seed = 12345,
                  chains = 8, 
                  iter_warmup = 200,
                  iter_sampling = 800,
                  refresh = 100,
                  adapt_delta = 0.995,
                  max_treedepth = 20)

write.csv(fit1bh$summary(),'./stan models/outs/summary/fit2bh_summary.csv')
fit1bh$save_object('./stan models/outs/fits/fit1bh.RDS')
fit1bh <- readRDS('./stan models/outs/fits/fit1bh.RDS')
loo_f1bh=fit1h$loo()

#residual check - mu1 is unadjusted for autocorrelation
d2bh=fit2bh$draws(variables=c('e_t','mu1','mu2'),format='draws_matrix')
e_t=apply(d2bh[,grepl('e_t',colnames(d2))],2,median)
mu1=apply(d2bh[,grepl('mu1',colnames(d2))],2,median)
mu2=apply(d2bh[,grepl('mu2',colnames(d2))],2,median)
plot(e_t~mu1)
#plot(e_t~mu2)

loo_summary=loo::loo_compare(loo_f2,loo_f2bh,loo_f2cs)
saveRDS(loo_summary,'loo_results_fitset2.RDS')

loo::loo_compare(loo_f1bh,loo_f2bh)

## Fit set 2: River-level varying effects for ECA effects ####

##F it set 3: Watershed area as a mediator of clearcut effects ####

fit3ric <- m3ric$sample(data=dl,
                  seed = 1235,
                  chains = 8, 
                  iter_warmup = 200,
                  iter_sampling = 800,
                  refresh = 100,
                  adapt_delta = 0.995,
                  max_treedepth = 20)

write.csv(fit3$summary(),'./stan models/outs/summary/fit3ric_summary.csv')
fit3$save_object('./stan models/outs/fits/fit3ric.RDS')
fit3=readRDS('stan models/outs/fits/fit3ric.RDS')
loo_f3=fit3$loo()

fit3bh <- m3bh$sample(data=dl,
                  seed = 1235,
                  chains = 8, 
                  iter_warmup = 200,
                  iter_sampling = 800,
                  refresh = 100,
                  adapt_delta = 0.995,
                  max_treedepth = 20)

write.csv(fit3bh$summary(),'./stan models/outs/summary/fit3bh_summary.csv')
fit3bh$save_object('./stan models/outs/fits/fit3bh.RDS')
fit3bh=readRDS('stan models/outs/fits/fit3bh.RDS')
loo_f3bh=fit3bh$loo()



#pink salmon

#data formatting...
pk10r=pk10r[order(factor(pk10r$River),pk10r$BroodYear),]
rownames(pk10r)=seq(1:nrow(pk10r))

#normalize ECA - logit transformation (ie. log(x/(1-x)))
pk10r$logit.ECA=qlogis(pk10r$ECA_age_proxy_forested_only+0.005)
pk10r$logit.ECA.std=(pk10r$logit.ECA-mean(pk10r$logit.ECA))/sd(pk10r$logit.ECA)

#sparse matrix of spawners
S_mat=make_design_matrix(pk10r$Spawners,pk10r$River)

#sparse matrix of ECA

#average ECA by stock
#just to see an overview of ECA by river
eca_s=pk10r%>%group_by(River)%>%summarize(m=mean(ECA_age_proxy_forested_only*100),m.std=mean(logit.ECA.std),range=max(ECA_age_proxy_forested_only*100)-min(ECA_age_proxy_forested_only*100),cu=unique(CU))

ECA_mat=make_design_matrix(pk10r$logit.ECA.std,pk10r$River)

#watershed area by river
area_mat=make_design_matrix(pk10r$ln_area_km2_std,pk10r$River)

#interaction matrix
ExA_mat=ECA_mat*area_mat

#extract max S for priors on capacity & eq. recruitment
smax_prior=pk10r%>%group_by(River) %>%summarize(m.s=max(Spawners),m.r=max(Recruits))

#ragged start and end points for each SR series
N_s=rag_n(pk10r$River)

#cus by stock
cu=distinct(pk10r,River,.keep_all = T)
summary(factor(cu$CU))

#time points for each series
L_i=pk10r%>%group_by(River)%>%summarize(l=n(),tmin=min(BroodYear)-1954+1,tmax=max(BroodYear)-1954+1)


dl=list(N=nrow(ch20r),
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
        start_y=N_s[,1],
        end_y=N_s[,2],
        start_t=L_i$tmin,
        end_t=L_i$tmax,
        pSmax_mean=0.5*smax_prior$m.s, #prior for smax based on max observed spawners
        pSmax_sig=smax_prior$m.s,
        pRk_mean=0.75*smax_prior$m.r, #prior for smax based on max observed spawners
        pRk_sig=smax_prior$m.r)


#old models - varying effects to river level (too flexible, likely)

fit1 <- m1$sample(data=dl,
                  seed = 12345,
                  chains = 8, 
                  iter_warmup = 200,
                  iter_sampling = 800,
                  refresh = 100,
                  adapt_delta = 0.995,
                  max_treedepth = 20)

write.csv(fit1$summary(),'./stan models/outs/summary/fit1ric_summary.csv')
fit1$save_object('stan models/outs/fits/fit1ric.RDS')
fit1=readRDS('./stan models/outs/fits/fit1ric.RDS')
loo_f1=fit1$loo()

d1=fit1$draws(variables=c('e_t','mu1','mu2'),format='draws_matrix')
e_t=apply(d1[,grepl('e_t',colnames(d1))],2,median)
mu1=apply(d1[,grepl('mu1',colnames(d1))],2,median)
mu2=apply(d1[,grepl('mu2',colnames(d1))],2,median)
plot(e_t~mu1)
plot(e_t~mu2)

fit1$loo()

fit1bh <- m1bh$sample(data=dl,
                      seed = 12345,
                      chains = 8, 
                      iter_warmup = 200,
                      iter_sampling = 800,
                      refresh = 100,
                      adapt_delta = 0.995,
                      max_treedepth = 20)

write.csv(fit1bh$summary(),'./stan models/outs/summary/fit1bh_summary.csv')
fit1bh$save_object('stan models/outs/fits/fit1bh.RDS')
fit1bh=readRDS('stan models/outs/fits/fit1bh.RDS')
loo_f1bh=fit1bh$loo()

#
loo::loo_compare(loo_f1,loo_f1bh)

d1bh=fit1bh$draws(variables=c('e_t','mu1','mu2'),format='draws_matrix')
e_t=apply(d1bh[,grepl('e_t',colnames(d1bh))],2,median)
mu1=apply(d1bh[,grepl('mu1',colnames(d1bh))],2,median)
mu2=apply(d1bh[,grepl('mu2',colnames(d1bh))],2,median)
plot(e_t~mu1)
plot(e_t~mu2)

fit1bh$loo()

