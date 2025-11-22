#make pairs plot between fixed effects to see correlation
#between effect of forestry and effect of ocean covariates

library(here)
library(bayesplot)



#read posterior

posterior <- read.csv(here("stan models","outs","posterior","ric_chm_cpd_ocean_covariates_logR.csv"))


mcmc_pairs(
  posterior,
  pars = c("b_for","b_npgo","b_sst"),
  off_diag_args = list(size = 1, alpha = 0.1)
)

# do this for every CU

mcmc_pairs(
  posterior,
  pars = c("b_for_cu.1.","b_npgo_cu.1.","b_sst_cu.1."),
  off_diag_args = list(size = 1, alpha = 0.1)
)

mcmc_pairs(
  posterior,
  regex_pars = c("b_for_cu.[2,4].","b_npgo_cu.[2,4].","b_sst_cu.[2,4]."),
  off_diag_args = list(size = 1, alpha = 0.1)
)
