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

ric_chm_cpd_cc_summary <- read.csv(here('stan models','outs','summary','ric_chm_cpd_npgo_ersst.csv'))

ric_chm_cpd_cc_summary %>% 
  select(variable, mean) %>% 
  filter(startsWith(variable, "alpha"))


bh_chm_cpd_summary <- read.csv(here('stan models','outs','summary','bh_chm_cpd_npgo_ersst.csv'))

bh_chm_cpd_summary %>% 
  select(variable, mean) %>% 
  filter(startsWith(variable, "alpha"))

mcmc_trace(post_bh_chm_cpd, pars = c('b_for'))

mcmc_trace(post_bh_chm_eca_npgo, c('b_for'))

#summary
bh_chm_eca_npgo$summary() %>% 
  select(variable, mean) %>% 
  filter(startsWith(variable, "alpha"))


# generic
file <- 'bh_chm_eca_npgo_ersst_fixed_alpha_trial.csv'

posterior <- read.csv(here('stan models','outs','posterior',file))

mcmc_trace(posterior[1:500,], pars = c('b_for'))
mcmc_trace(posterior[501:1000,], pars = c('b_for'))
mcmc_trace(posterior[1001:1500,], pars = c('b_for'))
mcmc_trace(posterior[1501:2000,], pars = c('b_for'))
mcmc_trace(posterior[2001:2500,], pars = c('b_for'))
mcmc_trace(posterior[2501:3000,], pars = c('b_for'))

mcmc_trace(posterior[1:500,], pars = c('b_npgo'))
mcmc_trace(posterior[501:1000,], pars = c('b_npgo'))
mcmc_trace(posterior[1001:1500,], pars = c('b_npgo'))
mcmc_trace(posterior[1501:2000,], pars = c('b_npgo'))
mcmc_trace(posterior[2001:2500,], pars = c('b_npgo'))
mcmc_trace(posterior[2501:3000,], pars = c('b_npgo'))

mcmc_trace(posterior[1:500,], pars = c('K.10.'))
mcmc_trace(posterior[1:500,], pars = c('alpha_j.52.'))

#plot all chains in one plot with different colors

summary <- read.csv(here('stan models','outs','summary',file))

summary %>% 
  select(variable, rhat, mean) %>% 
  # filter(startsWith(variable, "alpha")) %>% 
  filter(startsWith(variable, "sigma")) %>%
  filter(rhat > 1.011) %>% 
  #decresing rhat values
  arrange(-rhat)

#very high rhat - cnanot use this posterior ric_chm_eca_npgo_ersst_cc.csv,ric_chm_cpd_npgo_ersst_cc_ad95.csv,ric_chm_eca_npgo_ersst_cc_ad95.csv
#ric_chm_eca_npgo_sst_K has very high rhat and one chain was divergent


# look at the mean alpha values in the summary of the ricker models

summary <- read.csv(here('stan models','outs','summary','ric_chm_eca_static_cc.csv'))

summary %>% 
  select(variable, rhat, mean, median) %>% 
  filter(startsWith(variable, "sigma_a_rv"))

summary <- read.csv(here('stan models','outs','summary','bh_chm_eca_static_K_trial.csv'))

summary %>% 
  select(variable, rhat, mean, median) %>% 
  filter(startsWith(variable, "b_for"))


# plotting histogram of river level alphas

file <- 'bh_chm_cpd_npgo_ersst_fixed_alpha.csv'

posterior <- read.csv(here('stan models','outs','posterior',file))

summary <- read.csv(here('stan models','outs','summary',file))


summary %>% 
  select(variable, rhat, mean, median) %>% 
  filter(startsWith(variable, "alpha_j")) %>% 
  ggplot(aes(x = mean)) +
  geom_histogram(bins = 30, fill = "salmon", alpha = 0.5, color = "black")+
  theme_classic()+
  labs(x = "Mean river level alpha values", y = "Count")


summary %>% 
  select(variable, rhat, mean, median) %>% 
  filter(startsWith(variable, "alpha_cu")) %>% 
  ggplot(aes(x = mean)) +
  geom_histogram(bins = 30, fill = "salmon", alpha = 0.5, color = "black")+
  theme_classic()+
  labs(x = "Mean CU level alpha values", y = "Count")

summary %>% 
  select(variable, rhat, mean, median) %>% 
  filter(startsWith(variable, "sigma")) 


summary %>% 
  select(variable, rhat, mean, median) %>% 
  filter(startsWith(variable, "z_a_rv")) %>%
  ggplot(aes(x = mean)) +
  geom_histogram(bins = 30, fill = "salmon", alpha = 0.5, color = "black")+
  theme_classic()+
  labs(x = "Mean z score deviations", y = "Count")

summary %>% 
  select(variable, rhat, mean, median) %>% 
  filter(startsWith(variable, "z_a_cu")) %>%
  ggplot(aes(x = mean)) +
  geom_histogram(bins = 30, fill = "salmon", alpha = 0.5, color = "black")+
  theme_classic()+
  labs(x = "Mean z score deviations", y = "Count")











