# Goal - make prior and posterior distribution plots for some parameters:
# alpha0, all fixed effects

library(tidyverse)
library(here)
library(ggplot2)


# Load the data:

ric_chm_eca_ocean_covariates_logR_long_chain=read.csv(here('stan models','outs','posterior',
                                                           'ric_chm_eca_ocean_covariates_logR_long_chain.csv'),check.names=F)
ric_chm_cpd_ocean_covariates_logR_long_chain=read.csv(here('stan models','outs','posterior',
                                                           'ric_chm_cpd_ocean_covariates_logR_long_chain.csv'),check.names=F)


# Make the prior distribution data frame:

prior_df=data.frame(parameter=c('alpha0','b_for','b_sst','b_npgo','mu_sigma'),
                    mean=c(1.2,0,0,0,1),
                    sd=c(1,1,1,1,1))

#draw from prior

set.seed(123)

prior_draws_df <- prior_df %>%
  rowwise() %>%
  do(data.frame(parameter = .$parameter,
                value = rnorm(3000, mean = .$mean, sd = .$sd))) %>% 
  mutate(type = "prior")

#remove all values <0 from mu_sigma prior
prior_draws_df <- prior_draws_df %>%
  filter(!(parameter == "mu_sigma" & value < 0))

post_ric_chm_cpd_npgo_sst=ric_chm_cpd_ocean_covariates_logR_long_chain$draws(variables=c('b_for',
                                                                 'b_npgo',
                                                                 'b_sst',
                                                                 'alpha0','mu_sigma'),format='draws_matrix')

# Plot the prior and posterior distribution of alpha0 in one figure:

ggplot() +
  geom_density(data = prior_draws_df %>% filter(parameter == 'alpha0'), aes(x = value, fill = "prior") , alpha = 0.5) +
  geom_density(data = as.data.frame(post_ric_chm_cpd_npgo_sst) %>% select("alpha0"), aes(x = alpha0, fill = "posterior"), alpha = 0.5) +
  labs(title = 'Intrinsic productivity',
       x = 'Value',
       y = 'Density') +
  theme_classic() +
  scale_fill_manual(values = c("prior" = "#B2A3B5", "posterior" = '#53474E'))+
  theme(legend.title=element_blank(),
        legend.position  = c(0.9,0.5))

#forestry
ggplot() +
  geom_density(data = prior_draws_df %>% filter(parameter == 'beta_for'), aes(x = value, fill = "prior") , alpha = 0.5) +
  geom_density(data = as.data.frame(post_ric_chm_cpd_npgo_sst) %>% select("b_for"), aes(x = b_for, fill = "posterior"), alpha = 0.5) +
  labs(title = 'Effect of Forestry',
       x = 'Value',
       y = 'Density') +
  theme_classic() +
  scale_fill_manual(values = c("prior" = "#B2A3B5", "posterior" = '#53474E'))+
  theme(legend.title=element_blank(),
        legend.position  = c(0.9,0.5))

#sst
ggplot() +
  geom_density(data = prior_draws_df %>% filter(parameter == 'beta_sst'), aes(x = value, fill = "prior") , alpha = 0.5) +
  geom_density(data = as.data.frame(post_ric_chm_cpd_npgo_sst) %>% select("b_sst"), aes(x = b_sst, fill = "posterior"), alpha = 0.5) +
  labs(title = 'Effect of SST',
       x = 'Value',
       y = 'Density') +
  theme_classic() +
  scale_fill_manual(values = c("prior" = "#B2A3B5", "posterior" = '#53474E'))+
  theme(legend.title=element_blank(),
        legend.position  = c(0.9,0.5))



#npgo
ggplot() +
  geom_density(data = prior_draws_df %>% filter(parameter == 'beta_npgo'), aes(x = value, fill = "prior") , alpha = 0.5) +
  geom_density(data = as.data.frame(post_ric_chm_cpd_npgo_sst) %>% select("b_npgo"), aes(x = b_npgo, fill = "posterior"), alpha = 0.5) +
  labs(title = 'Effect of NPGO',
       x = 'Value',
       y = 'Density') +
  theme_classic() +
  scale_fill_manual(values = c("prior" = "#B2A3B5", "posterior" = '#53474E'))+
  theme(legend.title=element_blank(),
        legend.position  = c(0.9,0.5))

#mu_sigma
ggplot() +
  geom_density(data = prior_draws_df %>% filter(parameter == 'mu_sigma'), aes(x = value, fill = "prior") , alpha = 0.5) +
  geom_density(data = as.data.frame(post_ric_chm_cpd_npgo_sst) %>% select("mu_sigma"), aes(x = mu_sigma, fill = "posterior"), alpha = 0.5) +
  labs(title = 'Standard deviation',
       x = 'Value',
       y = 'Density') +
  theme_classic() +
  scale_fill_manual(values = c("prior" = "#B2A3B5", "posterior" = '#53474E'))+
  theme(legend.title=element_blank(),
        legend.position  = c(0.9,0.5))



# make the posterior dataframe long and then use facet_wrap to make all the plots into one plot

post_ric_chm_cpd_npgo_sst_long=as.data.frame(post_ric_chm_cpd_npgo_sst) %>%
  select("alpha0","b_for","b_sst","b_npgo","mu_sigma") %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>% 
  mutate(type = "posterior")


post_prior_long=prior_draws_df %>%
  rbind(post_ric_chm_cpd_npgo_sst_long)


#convert parameter to factor
post_prior_long$parameter <- factor(post_prior_long$parameter, levels = c("b_for","b_sst","b_npgo","alpha0","mu_sigma"))


math_labels <- as_labeller(c(
  "alpha0" = "alpha[0]",
  "b_for" = "beta[0]^CDA",
  "b_sst" = "beta[0]^SST",
  "b_npgo" = "beta[0]^NPGO",
  "mu_sigma" = "sigma^mu"
  ), default = label_parsed)

math_labels2 <- as_labeller(c(
  "alpha0" = paste(expression("\u03B1")),
  "b_for" = paste(expression("\u03B2")),
  "b_sst" = paste(expression("\u03B2"),"^SST"),
  "b_npgo" = paste(expression("\u03B2"),expression("\u1D2C"),expression("\u1D2C")),
  "mu_sigma" = "mu"
))




ggplot(post_prior_long, aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.5, linewidth = 0.7) +
  #label parameters with math symbols
  facet_wrap(~parameter, scales = "free_y", ncol = 2, labeller = labeller(parameter = math_labels), dir = "v") +
  labs(#title = 'Prior and Posterior Distributions',
       x = 'Value',
       y = 'Density') +
  xlim(-1,3) +
  theme_classic() +
  scale_fill_manual(values = c("posterior" = '#72BDA3', "prior" = "#B2A3B5"))+
  scale_color_manual(values = c("posterior" = '#72BDA3', "prior" = "#B2A3B5"))+
  theme(legend.title=element_blank(),
        legend.position  = c(0.8,0.1),
        legend.text = element_text(size = 10),
        title = element_text(size = 12),
        #remove facet label box
        strip.background = element_rect(fill="gray90", color = "transparent"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))

#save
ggsave(here('figures','prior_posterior_dist.png'), width = 6, height = 6, dpi = 300)







