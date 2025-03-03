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

posterior <- read.csv(here('stan models','outs','posterior','bh_chm_cpd_static_corrected.csv'))



mcmc_trace(posterior, pars = c('b_for'))


