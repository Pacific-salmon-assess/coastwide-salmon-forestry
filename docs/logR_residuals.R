
# plot the residuals of logR

library(here)
library(ggplot2)
library(tidyverse)
library(latex2exp)



# data


ch20rsc <- read.csv(here("origional-ecofish-data-models","Data","Processed",
                         "chum_SR_20_hat_yr_w_ocean_covariates.csv"))


#two rivers with duplicated names:
ch20rsc$River=ifelse(ch20rsc$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',ch20rsc$River)
ch20rsc$River=ifelse(ch20rsc$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',ch20rsc$River)


ch20rsc=ch20rsc[order(factor(ch20rsc$River),ch20rsc$BroodYear),]

ch20rsc$River_n <- as.numeric(factor(ch20rsc$River))

#normalize ECA 2 - square root transformation (ie. sqrt(x))
ch20rsc$sqrt.ECA=sqrt(ch20rsc$ECA_age_proxy_forested_only)
ch20rsc$sqrt.ECA.std=(ch20rsc$sqrt.ECA-mean(ch20rsc$sqrt.ECA))/sd(ch20rsc$sqrt.ECA)

#normalize CPD 2 - square root transformation (ie. sqrt(x))
ch20rsc$sqrt.CPD=sqrt(ch20rsc$disturbedarea_prct_cs)
ch20rsc$sqrt.CPD.std=(ch20rsc$sqrt.CPD-mean(ch20rsc$sqrt.CPD))/sd(ch20rsc$sqrt.CPD)

ch20rsc$npgo.std=(ch20rsc$npgo-mean(ch20rsc$npgo))/sd(ch20rsc$npgo)
ch20rsc$sst.std=(ch20rsc$spring_ersst-mean(ch20rsc$spring_ersst))/sd(ch20rsc$spring_ersst)

cu = distinct(ch20rsc,.keep_all = T)

cu_n = as.numeric(factor(cu$CU))

ch20rsc$CU_n <- cu_n

# make CU_NAME values Title case instead of all caps
ch20rsc$CU_name <- str_to_title(ch20rsc$CU_NAME)

# glimpse(ch20rsc)

max_eca_df <- ch20rsc %>% 
  select(River, River_n,CU, CU_n, CU_name, ECA_age_proxy_forested_only) %>%
  group_by(River) %>%
  summarize(River_n = first(River_n),
            River = first(River),
            CU = first(CU),
            CU_n = first(CU_n),
            CU_name = first(CU_name),
            eca_max = 100*max(ECA_age_proxy_forested_only, na.rm =TRUE)) %>% 
  mutate(eca_level = case_when(eca_max < 12 ~ 'low',
                               eca_max >= 12 & eca_max < 24 ~ 'medium',
                               eca_max >= 24 ~ 'high'))

max_cpd_df <- ch20rsc %>% 
  select(River, River_n,CU, CU_n, CU_name, disturbedarea_prct_cs) %>%
  group_by(River) %>% 
  summarize(River_n = first(River_n),
            River = first(River),
            CU = first(CU),
            CU_n = first(CU_n),
            CU_name = first(CU_name),
            cpd_max = max(disturbedarea_prct_cs, na.rm = TRUE))





ric_chm_eca_ocean_covariates=read.csv(here('stan models','outs','posterior','ric_chm_eca_ocean_covariates.csv'),check.names=F)
ric_chm_cpd_ocean_covariates=read.csv(here('stan models','outs','posterior','ric_chm_cpd_ocean_covariates.csv'),check.names=F)

df <- ch20rsc

posterior <- ric_chm_cpd_ocean_covariates

# posterior_bh <- bh_chm_cpd_ocean_covariates



rivers <- unique(df$River_n)

# plotting residuals of ricker model vs predicted values of ricker model for all rivers combined

residual_logR_df_full_chum <- data.frame()

for(river in rivers){
  
  river_data <- df %>% filter(River_n == river)
  
  posterior_rv_b_for <- posterior %>% 
    select(starts_with('b_for_rv')) %>%
    select(ends_with(paste0("[",river,"]")))
  
  posterior_rv_b_npgo <- posterior %>% 
    select(starts_with('b_npgo_rv')) %>%
    select(ends_with(paste0("[",river,"]")))
  
  posterior_rv_b_sst <- posterior %>% 
    select(starts_with('b_sst_rv')) %>%
    select(ends_with(paste0("[",river,"]")))
  
  posterior_rv_alpha_j <- posterior %>% 
    select(starts_with('alpha_j')) %>%
    select(ends_with(paste0("[",river,"]")))
  
  posterior_rv_S_max <- posterior %>% 
    select(starts_with('Smax')) %>%
    select(ends_with(paste0("[",river,"]")))
  
  
  residual_logR_df <- data.frame(spawners = river_data$Spawners,
                            recruits = river_data$Recruits,
                            observed_log_R = log(river_data$Recruits),
                            predicted_log_R = apply(matrix(log(river_data$Spawners), ncol=length(river_data$Spawners), nrow = 3000, byrow = TRUE) + (matrix(posterior_rv_alpha_j[,1], ncol = length(river_data$Spawners), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%river_data$Spawners + as.matrix(posterior_rv_b_for)%*%river_data$sqrt.CPD.std +
                                                        as.matrix(posterior_rv_b_npgo)%*%river_data$npgo.std + 
                                                        as.matrix(posterior_rv_b_sst)%*%river_data$sst.std ), 2, median),
                            
                            cpd = river_data$disturbedarea_prct_cs,
                            sqrt.CPD.std = river_data$sqrt.CPD.std,
                            forestry_effect = apply(as.matrix(posterior_rv_b_for)%*%river_data$sqrt.CPD.std , 2, median),
                            CU_name = river_data$CU_name
  )
  
  residual_logR_df$residuals <- residual_logR_df$observed_log_R - residual_logR_df$predicted_log_R
  
  residual_logR_df$partial_residuals <- residual_logR_df$residuals + residual_logR_df$forestry_effect
  
  residual_logR_df_full_chum <- rbind(residual_logR_df_full_chum, residual_logR_df)
  
 
  
}

ggplot(residual_logR_df_full_chum)+
  geom_point(aes(x = predicted_log_R, y = residuals, color = cpd), alpha = 0.2, size = 2) +
  geom_hline(yintercept = 0, color = 'black', linetype = 'dashed') +
  labs(title = "Ricker model with CPD, NPGO, ERSST", x = TeX(r"(Predicted $\log (Recruits)$)"), y = "Residuals") +
  scale_color_gradient2(name = 'CPD (%)',
                        low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 50)+
  ylim(-6, 6) +
  theme_classic() +
  
  theme(legend.position = "right",
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(1, "lines"),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.001, "cm"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5)
  )

ggsave(here("figures", "residuals_logR_chum.png"),
       width = 8, height = 6, dpi = 300, units = "in")

# look at the correlation between min S and median alpha from the posterior for each river


alpha_minS_full <- data.frame()

posterior_bh <- bh_chm_eca_ocean_covariates

for(river in rivers){
  
  river_data <- df %>% filter(River_n == river)
  
  
  posterior_rv_alpha_j <- posterior_bh %>% 
    select(starts_with('alpha_j')) %>%
    select(ends_with(paste0("[",river,"]")))
  
  alpha_minS <- data.frame(alpha = median(posterior_rv_alpha_j[,1]),
                           min_spawners = min(river_data$Spawners))
  
  alpha_minS_full <- rbind(alpha_minS_full, alpha_minS)
}

ggplot(alpha_minS_full)+
  geom_point(aes(x = log(min_spawners), y = alpha), color = 'cyan', size = 2, alpha = 0.3) +
  # geom_smooth(aes(x = log(min_spawners), y = alpha), method = 'lm', color = 'black', se = FALSE) +
  labs(title = "Minimum Spawners vs Median Alpha", x = "Minimum Spawners", y = "Median Alpha") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5)
  )


#plotting the residuals of eca model with logR 

df <- ch20rsc %>% 
  mutate(eca_level = case_when(
    ECA_age_proxy_forested_only*100 < 12 ~ 'low',
    ECA_age_proxy_forested_only*100 >= 12 & ECA_age_proxy_forested_only*100 < 24 ~ 'medium',
    ECA_age_proxy_forested_only*100 >= 24 ~ 'high'
  ))



ric_chm_eca_ocean_covariates_logR=read.csv(here('stan models','outs','posterior','ric_chm_eca_ocean_covariates_logR.csv'),check.names=F)
posterior <- ric_chm_eca_ocean_covariates_logR


residual_logR_df_full_chum <- data.frame()

rivers <- unique(df$River_n)




for(river in rivers){
  
  river_data <- df %>% filter(River_n == river)
  
  posterior_rv_b_for <- posterior %>% 
    select(starts_with('b_for_rv')) %>%
    select(ends_with(paste0("[",river,"]")))
  
  posterior_rv_b_npgo <- posterior %>% 
    select(starts_with('b_npgo_rv')) %>%
    select(ends_with(paste0("[",river,"]")))
  
  posterior_rv_b_sst <- posterior %>% 
    select(starts_with('b_sst_rv')) %>%
    select(ends_with(paste0("[",river,"]")))
  
  posterior_rv_alpha_j <- posterior %>% 
    select(starts_with('alpha_j')) %>%
    select(ends_with(paste0("[",river,"]")))
  
  posterior_rv_S_max <- posterior %>% 
    select(starts_with('Smax')) %>%
    select(ends_with(paste0("[",river,"]")))
  
  
  residual_logR_df <- data.frame(spawners = river_data$Spawners,
                                 recruits = river_data$Recruits,
                                 observed_log_R = log(river_data$Recruits),
                                 predicted_log_R = apply(matrix(log(river_data$Spawners), 
                                                                ncol=length(river_data$Spawners), 
                                                                nrow = 3000, byrow = TRUE) + (matrix(posterior_rv_alpha_j[,1], 
                                                                                                     ncol = length(river_data$Spawners), 
                                                                                                     nrow = length(posterior_rv_alpha_j[,1])) - 
                                                                                                as.matrix(1/posterior_rv_S_max)%*%river_data$Spawners + 
                                                                                                as.matrix(posterior_rv_b_for)%*%river_data$sqrt.CPD.std +
                                                                                                                                                            
                                                                                                as.matrix(posterior_rv_b_npgo)%*%river_data$npgo.std + 
                                                                                                                                                            
                                                                                                as.matrix(posterior_rv_b_sst)%*%river_data$sst.std ), 2, median),
                                 
                                 eca = river_data$ECA_age_proxy_forested_only,
                                 eca_level = river_data$eca_level,
                                 sqrt.CPD.std = river_data$sqrt.ECA.std,
                                 forestry_effect = apply(as.matrix(posterior_rv_b_for)%*%river_data$sqrt.ECA.std , 2, median),
                                 CU_name = river_data$CU_name
  )
  
  residual_logR_df$residuals <- residual_logR_df$observed_log_R - residual_logR_df$predicted_log_R
  
  residual_logR_df$partial_residuals <- residual_logR_df$residuals + residual_logR_df$forestry_effect
  
  residual_logR_df_full_chum <- rbind(residual_logR_df_full_chum, residual_logR_df)
  
  
  
}

ggplot(residual_logR_df_full_chum)+
  geom_point(aes(x = predicted_log_R, y = residuals, color = eca_level), alpha = 0.2, size = 2) +
  geom_hline(yintercept = 0, color = 'black', linetype = 'dashed') +
  labs(title = "Ricker model with ECA, NPGO, ERSST", x = TeX(r"(Predicted $\log (Recruits)$)"), y = "Residuals") +
  scale_color_manual(name = "ECA level (%)",
                     values = c('low' = '#35978f', 'medium' = 'gray', 'high' = '#bf812d'))+
  ylim(-6, 6) +
  theme_classic() +
  
  theme(legend.position = "right",
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(1, "lines"),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.001, "cm"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5)
  )

#save

ggsave(here("figures", "residuals_logR_model_chum_eca.png"),
       width = 8, height = 6, dpi = 300, units = "in")


#update nov 2025 - plot resiudlas with autocorrelation
library(tidyverse)
library(here)
posterior_new <- read.csv(here('..','..','..','OneDrive - University of Victoria',
                               'Documents','coastwide_salmon_forestry_analysis','model_outs_too_big','ric_chm_cpd_ocean_covariates_logR_new.csv'),check.names=F)

posterior_new <- read.csv("~/coastwide_forestry_salmon_analysis/model_outs_too_big/ric_chm_cpd_ocean_covariates_logR_new.csv")

glimpse(posterior_new)

#calculate difference between mu2 and ln_RS from observed data

# posterior_new %>% 
#   select(starts_with("mu2")) %>% 
#   # mutate(residuals = apply(.,1,
#   #                          function(observed,predicted) predicted - observed, 
#   #                          observed = ch20rsc$ln_RS))
#   mutate(median_mu2 = apply(., 1, median)) %>%
#   select(starts_with("median"))

n_rows <- nrow(ch20rsc)

#make df

residual_df <- data.frame(observed = ch20rsc$ln_RS, residual = NA)

for(i in 1:n_rows){
  #subtract mu2 from observed
  mu_col <- posterior_new %>% 
    select(paste0("mu2",".",i,"."))
  residual <- mu_col - ch20rsc$ln_RS[i]
  
  residual_df$residual[i] <- median(residual[,1])
  
}


residual_df <- data.frame(observed = log(ch20rsc$Recruits), forestry = ch20rsc$disturbedarea_prct_cs) %>% 
  mutate(predicted = posterior_new %>% 
           select(starts_with("mu2")) %>%
           apply(., 2, median),
         residual = observed - predicted)

#plot residuals as a function of fitted

ggplot(residual_df)+
  geom_point(aes(x = predicted, y = residual, color = forestry), alpha = 0.2, size = 2) +
  geom_hline(yintercept = 0, color = 'black', linetype = 'dashed') +
  labs(title = "", x = TeX(r"(Predicted $\log (Recruits)$)"), y = "Residuals") +
  ylim(-6, 6) +
  theme_classic() +
  scale_color_gradient2(name = 'CPD (%)',
                        low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 50)+
  theme(legend.position = "right",
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(1, "lines"),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.001, "cm"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5)
  )

# save

ggsave(here("figures", "residuals_w_autocorrelation_logR_chum.png"),
       width = 6, height = 4, dpi = 300, units = "in")










