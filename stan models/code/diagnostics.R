# Goal - to plot trace plots of chains in the model fitting process

# libraries
library(ggplot2)
library(tidyverse)
library(bayesplot)
library(patchwork)
library(here)


# load the model objects (RDS)

eca <- readRDS(here('stan models','outs','fits','bh_chm_eca_npgo_ersst_pdo_trial.RDS'))
cpd <- readRDS(here('stan models','outs','fits','bh_chm_cpd_npgo_ersst_pdo.RDS'))

# plot the trace plots

eca_posterior <- as.array(eca)
trace_eca <- mcmc_trace(eca, pars = c('b_for','b_for_cu','b_for_rv',
                                       'b_npgo','b_npgo_cu','b_npgo_rv',
                                       'b_sst','b_sst_cu','b_sst_rv','b_pdo_cu', 'b_pdo_rv',
                                       'alpha_j','Smax','sigma'))
mcmc_trace(post_ric_chm_cpd_npgo_sst, pars = c('b_for'))

mcmc_trace(post_ric_chm_cpd_npgo_sst, pars = c('b_for'))

cpd <- read.csv(here('stan models','outs','posterior','bh_chm_cpd_npgo_ersst_pdo.csv'))

mcmc_trace(cpd, pars = c('b_for'))

#read summary file

eca_summary <- read.csv(here('stan models','outs','summary','bh_chm_eca_npgo_ersst_pdo.csv'))

cpd_summary <- read.csv(here('stan models','outs','summary','bh_chm_cpd_npgo_ersst_pdo.csv'))

eca_summary %>% 
  select(variable, rhat) %>% 
  filter(rhat > 1.1)

cpd_summary %>% 
  select(variable, rhat) %>% 
  filter(rhat > 1.1)

# postrior - post_bh_chm_eca

mcmc_trace(post_bh_chm_eca, pars = c('b_for'))

posterior <- read.csv(here('stan models','outs','posterior','bh_chm_eca_npgo_ersst_corrected_trial.csv'))



mcmc_trace(posterior, pars = c('b_for'))


#read summary file

summary <- read.csv(here('stan models','outs','summary','bh_chm_eca_npgo_ersst_corrected.csv'))

glimpse(summary)

#check mean of variable "b_for"

summary %>% 
  select(variable, mean) %>% 
  filter(variable == "b_for")

#check if any rhat is greater than 1.1

summary %>% 
  select(variable, rhat) %>% 
  filter(rhat > 1.05)

#sigma_a_rv is more than 1.05

install.packages("DHARMa")

library(DHARMa)

#read rds

fit <- readRDS(here('stan models','outs','fits','bh_chm_cpd_static_trial_corrected.RDS'))

#simulate data

sim <- simulate


DHARMaRes = createDHARMa(simulatedResponse = sim, observedResponse = testData$observedResponse, 
                         fittedPredictedResponse = predict(fittedModel), integerResponse = T)
plot(DHARMaRes, quantreg = F)


cpd_summary <- read.csv(here('stan models','outs','summary','ric_chm_cpd_npgo_ersst_cc_trial.csv'))

cpd_summary %>% 
  select(variable, rhat) %>% 
  filter(rhat > 1.01)

posterior <- read.csv(here('stan models','outs','posterior','ric_chm_cpd_npgo_ersst_cc_trial.csv'))

color_scheme_set("viridis")
mcmc_trace_highlight(posterior, pars = c('b_for'),facet_args = list(nrow = 2, labeller = label_parsed))

cpd_summary %>% 
  filter(variable == "b_for")



bh_pk_eca_correct=read.csv(here('stan models',
                                'outs',
                                'posterior',
                                'bh_pk_eca_noac_corrected.csv'),check.names=F)


#plot all 6 (500 values) chains with different colors

mcmc_trace(bh_pk_eca_correct[1:500,], pars = c('b_for'))
mcmc_trace(bh_pk_eca_correct[501:1000,], pars = c('b_for'))
mcmc_trace(bh_pk_eca_correct[1001:1500,], pars = c('b_for'))
mcmc_trace(bh_pk_eca_correct[1501:2000,], pars = c('b_for'))
mcmc_trace(bh_pk_eca_correct[2001:2500,], pars = c('b_for'))
mcmc_trace(bh_pk_eca_correct[2501:3000,], pars = c('b_for'))

#look at summary

bh_pk_eca_summary <- read.csv(here('stan models','outs','summary','bh_pk_eca_noac_corrected.csv'))

bh_pk_eca_summary %>% 
  select(variable, rhat) %>% 
  filter(rhat > 1.02)

bh_pk_cpd_correct=read.csv(here('stan models',
                                'outs',
                                'posterior',
                                'bh_pk_cpd_noac_corrected.csv'),check.names=F)

ric_chm_eca_ersst_cc2 = read.csv(here('stan models','outs','posterior','ric_chm_eca_npgo_ersst_cc2.csv'),check.names=F)

ric_chm_cpd_ersst_cc2 = read.csv(here('stan models','outs','posterior','ric_chm_cpd_npgo_ersst_cc2.csv'),check.names=F)

# trace plot for b_for

mcmc_trace(ric_chm_eca_ersst_cc2, pars = c('b_for'))

mcmc_trace(ric_chm_cpd_ersst_cc2, pars = c('b_for'))
mcmc_trace(ric_chm_cpd_ersst_cc2[1:500,], pars = c('b_for'))
mcmc_trace(ric_chm_cpd_ersst_cc2[501:1000,], pars = c('b_for'))
mcmc_trace(ric_chm_cpd_ersst_cc2[1001:1500,], pars = c('b_for'))
mcmc_trace(ric_chm_cpd_ersst_cc2[1501:2000,], pars = c('b_for'))
mcmc_trace(ric_chm_cpd_ersst_cc2[2001:2500,], pars = c('b_for'))
mcmc_trace(ric_chm_cpd_ersst_cc2[2501:3000,], pars = c('b_for'))


#look at summary

ric_chm_eca_ersst_cc2_summary <- read.csv(here('stan models','outs','summary','ric_chm_eca_npgo_ersst_cc2.csv'))

ric_chm_eca_ersst_cc2_summary %>% 
  select(variable, rhat) %>% 
  filter(rhat > 1.01)

ric_chm_cpd_ersst_cc2_summary <- read.csv(here('stan models','outs','summary','ric_chm_cpd_npgo_ersst_cc2.csv'))

ric_chm_cpd_ersst_cc2_summary %>% 
  select(variable, rhat) %>% 
  filter(rhat > 1.04)

ric_chm_cpd_ersst_cc2_summary %>% 
  select(variable, rhat, mean, median) %>% 
  filter(startsWith(variable, "b_for_cu")) %>% 
  filter(rhat > 1.01)

# order by median value for b_for_cu

ric_chm_cpd_ersst_cc2_summary %>% 
  select(variable, rhat, mean, median) %>% 
  filter(startsWith(variable, "b_for_cu")) %>% 
  arrange(-median)

ric_chm_eca_ersst_cc2_summary %>% 
  select(variable, rhat, mean, median) %>% 
  filter(startsWith(variable, "b_for_cu")) %>% 
  arrange(-median)

# lowest b_for_cu is for CU 8 (southern coastal streams) 


ric_chm_eca_cc = read.csv(here('stan models','outs','posterior','ric_chm_eca_static_trial_cc.csv'),check.names=F)

# trace plot for b_for

mcmc_trace(ric_chm_eca_cc, pars = c('b_for'))


#summary

ric_chm_cpd_cc_summary <- read.csv(here('stan models','outs','summary','ric_chm_cpd_static_trial_cc.csv'))

ric_chm_cpd_cc_summary %>% 
  select(variable, mean) %>% 
  filter(variable == "alpha0")

mcmc_trace(post_bh_chm_cpd, pars = c('b_for'))
