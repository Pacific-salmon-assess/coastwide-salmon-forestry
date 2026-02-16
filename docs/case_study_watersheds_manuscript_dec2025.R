# goal - make figures for case study watersheds part of the supplement

# figure 1 - predicted recruitment for different levels of CPD and ECA
# figure 2 - predicted recruitment as a time series


library(here);library(dplyr); library(stringr)
library(ggplot2)
library(tidyverse)
library(bayesplot)
library(patchwork)
library(bcmaps)
library(hues)
library(GGally)
library(latex2exp)
library(bayestestR)
library(HDInterval)

# Data wrangling

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

ric_chm_eca_ocean_covariates_logR=read.csv(here('stan models','outs','posterior','ric_chm_eca_ocean_covariates_logR.csv'),check.names=F)
ric_chm_cpd_ocean_covariates_logR=read.csv(here('stan models','outs','posterior','ric_chm_cpd_ocean_covariates_logR.csv'),check.names=F)
ric_chm_eca_ocean_covariates_logR_long_chain=read.csv(here('stan models','outs','posterior','ric_chm_eca_ocean_covariates_logR_long_chain.csv'),check.names=F)
ric_chm_cpd_ocean_covariates_logR_long_chain=read.csv(here('stan models','outs','posterior','ric_chm_cpd_ocean_covariates_logR_long_chain.csv'),check.names=F)




watersheds <- c("VINER SOUND CREEK","CARNATION CREEK", "PHILLIPS RIVER", "NIMPKISH RIVER", "DEENA CREEK", "NEEKAS CREEK")

case_study_watersheds_data <- ch20rsc %>% 
  filter(River %in% watersheds)

plot_predicted_recruits <- function(posterior1 = ric_chm_cpd_ocean_covariates_logR,
                                    posterior2 = ric_chm_eca_ocean_covariates_logR,
                                    river = "CARNATION CREEK",  
                                    effect1 = "cpd",
                                    effect2 = "eca",
                                    species = "chum", 
                                    model1 = "CPD",
                                    model2 = "ECA"){
  if(species == "chum"){
    df <- ch20rsc 
    river_data <- ch20rsc %>% filter(River == river)
    river <- river_data$River_n[1]
    # df$sst.std <- (ch20rsc$spring_ersst-mean(ch20rsc$spring_ersst))/sd(ch20rsc$spring_ersst)
    
  } else if(species == "pink"){
    df <- pk10r
    river_data <- pk10r %>% filter(River == river)
    river <- river_data$River_n2[1]
    # df$sst.std <- (pk10r$spring_ersst-mean(pk10r$spring_ersst))/sd(pk10r$spring_ersst)
  }
  
  posterior_df_b_for <- posterior1 %>%
    select(starts_with(paste0('b_for_rv[',as.character(river),']'))) %>%
    pivot_longer(cols = everything(), 
                 names_to = 'River', 
                 names_prefix = 'b_for_rv',
                 values_to = "forestry") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>% 
    mutate(River = river)
  
  posterior_df_alpha <- posterior1 %>%
    select(starts_with(paste0('alpha_j[',as.character(river),']'))) %>%
    pivot_longer(cols = everything(), 
                 names_to = 'River', 
                 names_prefix = 'alpha_j',
                 values_to = "alpha_j") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>% 
    mutate(River = river)
  
  posterior_df_Smax <- posterior1 %>%
    select(starts_with(paste0('Smax[',as.character(river),']'))) %>%
    pivot_longer(cols = everything(), 
                 names_to = 'River', 
                 names_prefix = 'Smax',
                 values_to = "Smax") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>% 
    mutate(River = river)
  
  posterior_df_npgo <- posterior1 %>%
    select(starts_with(paste0('b_npgo_rv[',as.character(river),']'))) %>%
    pivot_longer(cols = everything(), 
                 names_to = 'River', 
                 names_prefix = 'b_npgo_rv',
                 values_to = "npgo") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>% 
    mutate(River = river)
  
  posterior_df_sst <- posterior1 %>%
    select(starts_with(paste0('b_sst_rv[',as.character(river),']'))) %>%
    pivot_longer(cols = everything(), 
                 names_to = 'River', 
                 names_prefix = 'b_sst_rv',
                 values_to = "sst") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>% 
    mutate(River = river)
  
  posterior_df2_b_for <- posterior2 %>% 
    select(starts_with(paste0('b_for_rv[',as.character(river),']'))) %>%
    pivot_longer(cols = everything(), 
                 names_to = 'River', 
                 names_prefix = 'b_for_rv',
                 values_to = "forestry") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>% 
    mutate(River = river)
  
  #make a time series of median log Recruits with credible intervals using the spawner time series and forestry time series for each year
  
  predicted_recruits <- data.frame(BroodYear = river_data$BroodYear,
                                   Spawners = river_data$Spawners,
                                   Recruits = river_data$Recruits,
                                   predicted_log_recruits = apply(matrix(log(river_data$Spawners),ncol = length(river_data$Spawners), nrow = length(posterior_df_alpha$alpha_j) ) + 
                                                                    (matrix(posterior_df_alpha$alpha_j, ncol = length(river_data$Spawners), nrow = length(posterior_df_alpha$alpha_j)) - as.matrix(1/posterior_df_Smax[,1])%*%river_data$Spawners + as.matrix(posterior_df_b_for[,1])%*%river_data$sqrt.CPD.std +
                                                                       as.matrix(posterior_df_npgo[,1])%*%river_data$npgo.std + 
                                                                       as.matrix(posterior_df_sst[,1])%*%river_data$sst.std ), 2, median),
                                   predicted_log_recruits_upper <- apply(matrix(log(river_data$Spawners),ncol = length(river_data$Spawners), nrow = length(posterior_df_alpha$alpha_j) ) + 
                                                                           (matrix(posterior_df_alpha$alpha_j, ncol = length(river_data$Spawners), nrow = length(posterior_df_alpha$alpha_j)) - as.matrix(1/posterior_df_Smax[,1])%*%river_data$Spawners + as.matrix(posterior_df_b_for[,1])%*%river_data$sqrt.CPD.std +
                                                                              as.matrix(posterior_df_npgo[,1])%*%river_data$npgo.std + 
                                                                              as.matrix(posterior_df_sst[,1])%*%river_data$sst.std ), 2, quantile, 0.975),
                                   
                                   predicted_log_recruits_lower <- apply(matrix(log(river_data$Spawners),ncol = length(river_data$Spawners), nrow = length(posterior_df_alpha$alpha_j) ) + 
                                                                           (matrix(posterior_df_alpha$alpha_j, ncol = length(river_data$Spawners), nrow = length(posterior_df_alpha$alpha_j)) - as.matrix(1/posterior_df_Smax[,1])%*%river_data$Spawners + as.matrix(posterior_df_b_for[,1])%*%river_data$sqrt.CPD.std +
                                                                              as.matrix(posterior_df_npgo[,1])%*%river_data$npgo.std + 
                                                                              as.matrix(posterior_df_sst[,1])%*%river_data$sst.std ), 2, quantile, 0.025),
                                   
                                   forestry = river_data$sqrt.CPD.std,
                                   
                                   lowest_forestry = min(river_data$sqrt.CPD.std),
                                   
                                   predicted_log_recruits_low_forestry = apply(matrix(log(river_data$Spawners),ncol = length(river_data$Spawners), nrow = length(posterior_df_alpha$alpha_j) ) + 
                                                                                 (matrix(posterior_df_alpha$alpha_j, ncol = length(river_data$Spawners), nrow = length(posterior_df_alpha$alpha_j)) - as.matrix(1/posterior_df_Smax[,1])%*%river_data$Spawners + as.matrix(posterior_df_b_for[,1])%*%rep(min(river_data$sqrt.CPD.std),length(river_data$sqrt.CPD.std)) +
                                                                                    as.matrix(posterior_df_npgo[,1])%*%river_data$npgo.std + 
                                                                                    as.matrix(posterior_df_sst[,1])%*%river_data$sst.std ), 2, median)
                                   
                                   
                                   
  )
  
  
  
  
  
  
  #quick plot of log_recruits time series - observed and predicted
  plot1 <- ggplot() +
    geom_line(data = predicted_recruits, aes(x = BroodYear, y = predicted_log_recruits, color = "Predicted"), size = 1, alpha = 0.5) +
    geom_line(data = predicted_recruits, aes(x = BroodYear, y = predicted_log_recruits_low_forestry, color = "Predicted low forestry"), size = 1, alpha = 0.5) +
    geom_point(data = predicted_recruits, aes(x = BroodYear, y = log(Recruits),  fill = forestry), size = 2,shape = 21, color = "white",  alpha = 0.5) +
    # geom_line(data = predicted_recruits, aes(x = BroodYear, y = log(Recruits), color = "Observed"), size = 2, alpha = 0.5) +
    geom_ribbon(data = predicted_recruits, aes(x = BroodYear, ymin = predicted_log_recruits_lower, ymax = predicted_log_recruits_upper),
                fill = 'gray70', alpha = 0.4) +
    labs(x = "Brood Year", y = "log(Recruits)") +
    scale_color_manual(name = "Legend",
                       values = c("Predicted" = 'gray20',
                                  "Predicted low forestry" = '#A2C5AC')) +
    scale_fill_gradient2(name = 'CPD (standardized)',
                         low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 0)+
    theme_classic() +
    theme(legend.position = "right",
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 14, hjust = 0.5))
  
  plot2 <- ggplot() +
    geom_line(data = predicted_recruits, aes(x = BroodYear, y = exp(predicted_log_recruits), color = "Predicted"), size = 1, alpha = 0.5) +
    geom_line(data = predicted_recruits, aes(x = BroodYear, y = exp(predicted_log_recruits_low_forestry), color = "Predicted low forestry"), size = 1, alpha = 0.5) +
    geom_point(data = predicted_recruits, aes(x = BroodYear, y = Recruits,  fill = forestry), size = 2,shape = 21, color = "white",  alpha = 0.5) +
    # geom_line(data = predicted_recruits, aes(x = BroodYear, y = log(Recruits), color = "Observed"), size = 2, alpha = 0.5) +
    geom_ribbon(data = predicted_recruits, aes(x = BroodYear, ymin = exp(predicted_log_recruits_lower), 
                                               ymax = exp(predicted_log_recruits_upper)),
                fill = 'gray70', alpha = 0.4) +
    labs(x = "Brood Year", y = "log(Recruits)") +
    scale_color_manual(name = "Legend",
                       values = c("Predicted" = 'gray20',
                                  "Predicted low forestry" = '#A2C5AC')) +
    scale_fill_gradient2(name = 'CPD (standardized)',
                         low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 0)+
    theme_classic() +
    theme(legend.position = "right",
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 14, hjust = 0.5))
  
  return(plot1)
  
  
}

# plot_predicted_recruits(river = watersheds[1])


plot_recruit_spawner_river <- function(data, species = "chum", river_name,  posterior, posterior_a_t, posterior_a_t_bh){
  
  # river_data <- df %>% filter(River_n == river)
  
  if(species == "chum"){
    df <- ch20rsc 
    river_data <- ch20rsc %>% filter(River == river_name)
    river <- river_data$River_n[1]
    # df$sst.std <- (ch20rsc$spring_ersst-mean(ch20rsc$spring_ersst))/sd(ch20rsc$spring_ersst)
    
  } else if(species == "pink"){
    df <- pk10r
    river_data <- pk10r %>% filter(River == river_name)
    river <- river_data$River_n2[1]
    # df$sst.std <- (pk10r$spring_ersst-mean(pk10r$spring_ersst))/sd(pk10r$spring_ersst)
  }
  
  if(species == "chum"){
    posterior_rv_b_for <- posterior %>% 
      select(starts_with('b_for_rv')) %>%
      select(ends_with(paste0("[",river,"]")))
    
    posterior_rv_alpha_j <- posterior %>% 
      select(starts_with('alpha_j')) %>%
      select(ends_with(paste0("[",river,"]")))
    
    posterior_rv_S_max <- posterior %>% 
      select(starts_with('Smax')) %>%
      select(ends_with(paste0("[",river,"]")))
    
    
  } else if(species == "pink"){
    
    river_wo_broodline <- river_data$River_n
    river_w_broodline <- river_data$River_n2
    
    posterior_rv_b_for <- posterior %>% 
      select(starts_with('b_for_rv')) %>%
      select(ends_with(paste0("[",river_wo_broodline,"]")))
    
    posterior_rv_alpha_j <- posterior %>% 
      select(starts_with('alpha_j')) %>%
      select(ends_with(paste0("[",river_w_broodline,"]")))
    
    posterior_rv_S_max <- posterior %>% 
      select(starts_with('Smax')) %>%
      select(ends_with(paste0("[",river_w_broodline,"]")))
    
    #time varying productivity for ricker
    
    # if(river_data$Broodline[1] == "Even"){
    #   
    #   posterior_a_t_alpha_j <- posterior_a_t %>% 
    #   select(starts_with('alpha_j')) %>%
    #   select(ends_with(paste0("[",river_w_broodline,"]")))
    #   
    #   posterior_a_t_alpha_t <- posterior_a_t %>%
    #   select(starts_with('alpha_t[1,'))
    #   
    #   posterior_a_t_alpha_t_j_full <- matrix(posterior_a_t_alpha_j[,1], ncol = length(posterior_a_t_alpha_t), nrow = length(posterior_a_t_alpha_j[,1])) + posterior_a_t_alpha_t
    #   
    #   posterior_a_t_alpha_t_j <- posterior_a_t_alpha_t_j_full %>%
    #     pivot_longer(cols = everything(), names_to = "year", values_to = "alpha_t_odd") %>%
    #     mutate(year = as.numeric(substring(year, 11, ifelse(nchar(year) == 12, 11, 12))) + 1953)
    #   
    # } else if(river_data$Broodline[1] == "Odd"){
    #   
    #   
    #   posterior_a_t_alpha_j <- posterior_a_t %>% 
    #   select(starts_with('alpha_j')) %>%
    #   select(ends_with(paste0("[",river_w_broodline,"]")))
    #   
    #   posterior_a_t_alpha_t <- posterior_a_t %>%
    #   select(starts_with('alpha_t[2,'))
    #   
    #   posterior_a_t_alpha_t_j_full <- matrix(posterior_a_t_alpha_j[,1], ncol = length(posterior_a_t_alpha_t), nrow = length(posterior_a_t_alpha_j[,1])) + posterior_a_t_alpha_t
    #   
    #   posterior_a_t_alpha_t_j <- posterior_a_t_alpha_t_j_full %>%
    #     pivot_longer(cols = everything(), names_to = "year", values_to = "alpha_t_odd") %>%
    #     mutate(year = as.numeric(substring(year, 11, ifelse(nchar(year) == 12, 11, 12))) + 1953)
    #   
    #   
    #   
    #   
    #   
    #   
    # }
    
    
    
    
    
    
  }
  
  
  spawners_predicted <- seq(0, max(river_data$Spawners), length.out = 100)
  
  # calculate recruit prediction
  
  low_cpd <- min(river_data$sqrt.CPD.std)
  high_cpd <- max(river_data$sqrt.CPD.std)
  # avg_cpd <- mean(river_data$sqrt.CPD.std)
  mid_cpd <- min(river_data$sqrt.CPD.std) + (max(river_data$sqrt.CPD.std) - min(river_data$sqrt.CPD.std))/2
  
  mid_cpd_real <- min(river_data$disturbedarea_prct_cs) + (max(river_data$disturbedarea_prct_cs) - min(river_data$disturbedarea_prct_cs))/2
  
  recruits_predicted <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - 
                                     as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*mid_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, median))*spawners_predicted
  
  recruits_predicted_lower <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted), 2, quantile, c(0.025)))*spawners_predicted
  
  recruits_predicted_upper <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted), 2, quantile, c(0.975)))*spawners_predicted
  
  
  
  # if(low_cpd == high_cpd){
  #   scale_color <- scale_color_manual(name = 'CPD in River (%)', values = c("black"))
  #   
  # } else{
  #    scale_color <- scale_color_gradient2(name = 'CPD in River (%)',
  #                                     low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 50)
  # }
  
  
  recruits_predicted_low_cpd <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*low_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, median))*spawners_predicted
  
  recruits_predicted_high_cpd <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*high_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, median))*spawners_predicted
  
  
  recruits_predicted_low_cpd_lower <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*low_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, quantile, c(0.025)))*spawners_predicted
  
  recruits_predicted_low_cpd_upper <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*low_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, quantile, c(0.975)))*spawners_predicted
  
  recruits_predicted_high_cpd_lower <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*high_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, quantile, c(0.025)))*spawners_predicted
  
  recruits_predicted_high_cpd_upper <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*high_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, quantile, c(0.975)))*spawners_predicted
  
  
  
  #make dataframe
  
  prediction_df <- data.frame(spawners = spawners_predicted,
                              recruits = recruits_predicted,
                              recruits_lower = recruits_predicted_lower,
                              recruits_upper = recruits_predicted_upper,
                              recruits_low_cpd = recruits_predicted_low_cpd,
                              recruits_high_cpd = recruits_predicted_high_cpd,
                              recruits_low_cpd_lower = recruits_predicted_low_cpd_lower,
                              recruits_low_cpd_upper = recruits_predicted_low_cpd_upper,
                              recruits_high_cpd_lower = recruits_predicted_high_cpd_lower,
                              recruits_high_cpd_upper = recruits_predicted_high_cpd_upper,
                              log_RS = log(recruits_predicted/spawners_predicted))
  
  
  
  
  
  
  
  #plot the time varying productivity vs year, with log(R/S) data
  
  
  
  
  
  #plot recruit vs spawner as points
  
  p1 <- ggplot() +
    geom_point(data = river_data,aes(x = Spawners, y = Recruits, color = disturbedarea_prct_cs), alpha = 0.5, size = 2) +
    geom_line(data = prediction_df, aes(x = spawners, y = recruits_low_cpd), color = '#35978f', size = 1, alpha = 0.5) +
    geom_line(data = prediction_df, aes(x = spawners, y = recruits_high_cpd), color = '#bf812d', size = 1, alpha = 0.5) +
    geom_line(data = prediction_df, aes(x = spawners, y = recruits), color = "black", size = 1, alpha = 0.5) +
    geom_ribbon(data = prediction_df, aes(x = spawners, ymin = recruits_low_cpd_lower, ymax = recruits_low_cpd_upper), fill ='#35978f', alpha = 0.5) +
    geom_ribbon(data = prediction_df, aes(x = spawners, ymin = recruits_high_cpd_lower, ymax = recruits_high_cpd_upper), fill ='#bf812d', alpha = 0.5) +
    #make y log scale
    # scale_y_log10() +
    labs(title = "", x = "Spawners", y = "Recruits") +
    # scale_color_manual(name = "CPD", values = c("Low" = '#35978f', "High" = '#bf812d')) +
    scale_color_gradient2(name = 'CDA in River (%)',
                          low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 20)+
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
  
  
  # p3 log R/s vs spawners
  
  p3 <- ggplot(river_data) + 
    geom_point(aes(x = Spawners, y = log(Recruits/Spawners), color = disturbedarea_prct_cs), alpha = 0.5, size = 2) +
    geom_line(data = prediction_df, aes(x = spawners, y = log(recruits_predicted_low_cpd/spawners_predicted)), color = '#35978f', size = 1, alpha = 0.5) +
    geom_line(data = prediction_df, aes(x = spawners, y = log(recruits_predicted_high_cpd/spawners_predicted)), color = '#bf812d', size = 1, alpha = 0.5) +
    geom_line(data = prediction_df, aes(x = spawners, y = log_RS), color = "black", size = 1, alpha = 0.5) +
    geom_ribbon(data = prediction_df, aes(x = spawners, ymin = log(recruits_predicted_lower/spawners_predicted), ymax = log(recruits_predicted_upper/spawners_predicted)), fill = "gray", alpha = 0.5) +
    labs(#title = "Ricker model with CPD, NPGO, ERSST", 
      x = "Spawners", y = TeX(r"($\log \left(\frac{Recruits}{Spawners}\right)$)")) +
    scale_color_gradient2(name = 'CDA in River (%)',
                          low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 20)+
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
  
  

  
  
  
  return((p1+p3)+ plot_layout(guides = "collect"))
}

plot_recruit_spawner_river(data = case_study_watersheds_data,
                            species = "chum",
                            river_name = watersheds[1],
                            posterior = ric_chm_cpd_ocean_covariates_logR,
                            posterior_a_t = NULL,
                            posterior_a_t_bh = NULL)


plot_both_forestry_effects_river_together <- function(posterior1 = ric_chm_cpd_ocean_covariates_logR,
                                                      posterior2 = ric_chm_eca_ocean_covariates_logR,
                                                      river_name,
                                                      river,
                                                      effect1 = "cpd",
                                                      effect2 = "eca",
                                                      species = "chum", 
                                                      model1 = "Cumulative\nDisturbance",
                                                      model2 = "ECA",
                                                      xlim = c(-0.5, 0.5)){
  
  
  
  
  posterior_df <- posterior1 %>%
    select(starts_with('b_for_rv')) %>%
    pivot_longer(cols = everything(), 
                 names_to = 'River', 
                 names_prefix = 'b_for_rv',
                 values_to = "forestry") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) 
  
  posterior_df2 <- posterior2 %>% 
    select(starts_with('b_for_rv')) %>%
    pivot_longer(cols = everything(), 
                 names_to = 'River', 
                 names_prefix = 'b_for_rv',
                 values_to = "forestry") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River)
  
  
  posterior_df_river1 <- posterior_df %>% filter(River_n == river) %>% 
    mutate(model = model1)
  
  posterior_df_river2 <- posterior_df2 %>% filter(River_n == river) %>% 
    mutate(model = model2)
  
  posterior_df_river <- posterior_df_river1 %>% 
    bind_rows(posterior_df_river2)
  
  #color by river
  plot1 <- ggplot() +
    stat_density(data= posterior_df_river, aes(forestry,#!!sym(forestry), 
                                               
                                               group = model, 
                                               color = model,
                                               fill = model),
                 geom = 'area', position = 'identity', 
                 alpha = 0.4, linewidth = 0.8) +
    
    geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
    labs(x = "Standardized coefficients", y = "Posterior density", title = river_name) +
    xlim(xlim[1], xlim[2]) +
    scale_color_manual(name = "",
                       values = c("ECA" = '#7F6A93', 
                                  "Cumulative\nDisturbance" = '#A2C5AC')) +
    scale_fill_manual(name = "",
                      values = c("ECA" = '#7F6A93', 
                                 "Cumulative\nDisturbance" = '#A2C5AC')) +
    # scale_color +
    # geom_density(aes(posterior$b_for), color = 'black', linewidth = 1.2, alpha = 0.2)+
    #vline at the median value of the posterior
    theme_classic()+
    theme(legend.position = c(0.9,0.9),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 14, hjust = 0)
    )
  
  
  
  return(plot1)
}


plot_all_effects_river_together <- function(posterior1 = ric_chm_cpd_ocean_covariates_logR,
                                            river_name,
                                            river,
                                            effect1 = "cpd",
                                            effect2 = "sst",
                                            effect3 = "npgo",
                                            species = "chum", 
                                            xlim = c(-0.5, 0.5)){
  
  
  
  
  posterior_df <- posterior1 %>%
    select(starts_with('b_for_rv'),starts_with('b_sst_rv'),starts_with('b_npgo_rv')) %>%
    pivot_longer(cols = everything(),
                 names_to = c('Effect','River'),
                 names_pattern = 'b_(.*)_rv(.*)',
                 values_to = "coefficient") %>%
    
    
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) 
  
  
  if(effect1 == "cpd"){
    posterior_df_river <- posterior_df %>% filter(River_n == river) %>% 
      mutate(Effect = case_when(Effect == "for" ~ "CDA",
                                Effect == "sst" ~ "SST",
                                Effect == "npgo" ~ "NPGO"))
  } else if(effect1 == "eca"){
    posterior_df_river <- posterior_df %>% filter(River_n == river) %>% 
      mutate(Effect = case_when(Effect == "for" ~ "ECA",
                                Effect == "sst" ~ "SST",
                                Effect == "npgo" ~ "NPGO"))
  }
  
  

  
  #color by river
  plot1 <- ggplot() +
    stat_density(data= posterior_df_river, aes(coefficient,#!!sym(forestry), 
                                               
                                               group = Effect, 
                                               color = Effect,
                                               fill = Effect),
                 geom = 'area', position = 'identity', 
                 alpha = 0.4, linewidth = 0.8) +
    
    geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
    geom_vline(aes(xintercept = median(posterior_df_river$coefficient[posterior_df_river$Effect == "CDA"]), 
               color = "CDA"), linetype = "dashed", linewidth = 0.8) +
    geom_vline(aes(xintercept = median(posterior_df_river$coefficient[posterior_df_river$Effect == "SST"]), 
               color = "SST"), linetype = "dashed", linewidth = 0.8) +
    geom_vline(aes(xintercept = median(posterior_df_river$coefficient[posterior_df_river$Effect == "NPGO"]),
               color = "NPGO"), linetype = "dashed", linewidth = 0.8) +
    labs(x = "Standardized coefficients", y = "Posterior density", title = river_name) +
    xlim(xlim[1], xlim[2]) +
    scale_color_manual(name = "",
                       values = c("ECA" = "#ADcCA5", 
                                  "CDA" = "#ADcCA5",
                                  "SST" = "#C78c63",
                                  "NPGO" = "#829Dc6")) +
    scale_fill_manual(name = "",
                      values = c("ECA" = "#ADcCA5", 
                                 "CDA" = "#ADcCA5",
                                 "SST" = "#C78c63",
                                 "NPGO" = "#829Dc6")) +
    # scale_color +
    # geom_density(aes(posterior$b_for), color = 'black', linewidth = 1.2, alpha = 0.2)+
    #vline at the median value of the posterior
    theme_classic()+
    theme(legend.position = c(0.8,0.9),
          legend.background = element_rect(fill = alpha('white', 0.5)),
          legend.text = element_text(size = 7),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_blank()
    )
  
  
  
  return(plot1)
}


plot_productivity_change_river_together <- function(posterior1 = ric_chm_cpd_ocean_covariates_logR_long_chain,
                                                    posterior2 = ric_chm_eca_ocean_covariates_logR_long_chain,
                                                    river_name = "CARNATION CREEK",  
                                                    effect1 = "cpd",
                                                    effect2 = "eca",
                                                    species = "chum", 
                                                    model1 = "CPD",
                                                    model2 = "ECA"){
  
  
  if(species == "chum"){
    df <- ch20rsc 
    river_data <- ch20rsc %>% filter(River == river_name)
    river <- river_data$River_n[1]
    # df$sst.std <- (ch20rsc$spring_ersst-mean(ch20rsc$spring_ersst))/sd(ch20rsc$spring_ersst)
    
  } else if(species == "pink"){
    df <- pk10r
    river_data <- pk10r %>% filter(River == river_name)
    river <- river_data$River_n2[1]
    # df$sst.std <- (pk10r$spring_ersst-mean(pk10r$spring_ersst))/sd(pk10r$spring_ersst)
  }
  
  posterior_df <- posterior1 %>%
    select(starts_with(paste0('b_for_rv[',as.character(river),']'))) %>%
    pivot_longer(cols = everything(), 
                 names_to = 'River', 
                 names_prefix = 'b_for_rv',
                 values_to = "forestry") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>% 
    mutate(River = river)
  
  posterior_df2 <- posterior2 %>% 
    select(starts_with(paste0('b_for_rv[',as.character(river),']'))) %>%
    pivot_longer(cols = everything(), 
                 names_to = 'River', 
                 names_prefix = 'b_for_rv',
                 values_to = "forestry") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>% 
    mutate(River = river)
  
  
  
  
  eca <- seq(0,1,length.out=100)
  
  eca_sqrt <- sqrt(eca)
  
  eca_sqrt_std <- (eca_sqrt-mean(eca_sqrt))/sd(eca_sqrt)
  
  cpd <- seq(0,100,length.out=100)
  
  cpd_sqrt <- sqrt(cpd)
  
  cpd_sqrt_std <- (cpd_sqrt-mean(cpd_sqrt))/sd(cpd_sqrt)
  
  
  
  
  no_eca <- min(eca_sqrt_std)
  
  no_cpd <- min(cpd_sqrt_std)
  
  
  
  # need to change b_rv to posterior 1 and posterior 2 and then make figure
  
  productivity_cpd <- (exp(as.matrix(posterior_df[,1])%*%
                             (cpd_sqrt_std-no_cpd)))*100 - 100
  
  productivity_eca <- (exp(as.matrix(posterior_df2[,1])%*%
                             (eca_sqrt_std-no_eca)))*100 - 100
  
  
  productivity_cpd_df <- data.frame(productivity_cpd_median = apply(productivity_cpd,2,median),
                                    productivity_cpd_025 = apply(productivity_cpd,2,quantile,c(0.025), 
                                                                 row.names = c("q025")),
                                    productivity_cpd_975 = apply(productivity_cpd,2,quantile,c(0.975),
                                                                 row.names = c("q975")),
                                    productivity_cpd_100 = apply(productivity_cpd,2,quantile,c(0.1),
                                                                 row.names = c("q100")),
                                    productivity_cpd_900 = apply(productivity_cpd,2,quantile,c(0.9),
                                                                 row.names = c("q900")),
                                    
                                    productivity_cpd_025_hd = apply(productivity_cpd,2,HDInterval::hdi, credMass = 0.95)[1,],
                                    productivity_cpd_975_hd = apply(productivity_cpd,2,HDInterval::hdi, credMass = 0.95)[2,],
                                    productivity_cpd_100_hd = apply(productivity_cpd,2,HDInterval::hdi, credMass = 0.8)[1,],
                                    productivity_cpd_900_hd = apply(productivity_cpd,2,HDInterval::hdi, credMass = 0.8)[2,],
                                    productivity_cpd_50_lower_hd = apply(productivity_cpd,2,HDInterval::hdi, credMass = 0.5)[1,],
                                    productivity_cpd_50_upper_hd = apply(productivity_cpd,2,HDInterval::hdi, credMass = 0.5)[2,],
                                    
                                    
                                    
                                    cpd_sqrt_std = cpd_sqrt_std,
                                    cpd = cpd,
                                    max_cpd = max(river_data$disturbedarea_prct_cs),
                                    max_sqrt_cpd = max(river_data$sqrt.CPD),
                                    model = "CPD"
  )
  
  
  productivity_eca_df <- data.frame(productivity_eca_median = apply(productivity_eca,2,median),
                                    productivity_eca_025 = apply(productivity_eca,2,quantile,c(0.025), 
                                                                 row.names = c("q025")),
                                    productivity_eca_975 = apply(productivity_eca,2,quantile,c(0.975),
                                                                 row.names = c("q975")),
                                    productivity_eca_100 = apply(productivity_eca,2,quantile,c(0.1),
                                                                 row.names = c("q100")),
                                    productivity_eca_900 = apply(productivity_eca,2,quantile,c(0.9),
                                                                 row.names = c("q900")),
                                    
                                    productivity_eca_025_hd = apply(productivity_eca,2,HDInterval::hdi, credMass = 0.95)[1,],
                                    productivity_eca_975_hd = apply(productivity_eca,2,HDInterval::hdi, credMass = 0.95)[2,],
                                    productivity_eca_100_hd = apply(productivity_eca,2,HDInterval::hdi, credMass = 0.8)[1,],
                                    productivity_eca_900_hd = apply(productivity_eca,2,HDInterval::hdi, credMass = 0.8)[2,],
                                    productivity_eca_50_lower_hd = apply(productivity_eca,2,HDInterval::hdi, credMass = 0.5)[1,],
                                    productivity_eca_50_upper_hd = apply(productivity_eca,2,HDInterval::hdi, credMass = 0.5)[2,],
                                    
                                    eca_sqrt_std = eca_sqrt_std,
                                    eca = eca,
                                    max_eca = max(river_data$ECA_age_proxy_forested_only),
                                    max_sqrt_eca = max(river_data$sqrt.ECA),
                                    model = "ECA"
  )
  
  
  
  
  plot1 <- ggplot(productivity_cpd_df) +
    geom_line(aes(x = cpd, y = productivity_cpd_median, group = 1,
                  color = "median"),alpha=0.9, linewidth = 0.8) +
    geom_ribbon(aes(x = cpd, ymin = productivity_cpd_025_hd, ymax = productivity_cpd_975_hd, fill = "95% credible interval"),
                alpha = 0.2) +
    geom_ribbon(aes(x = cpd, ymin = productivity_cpd_100_hd, ymax = productivity_cpd_900_hd, fill = "80% credible interval"),
                alpha = 0.4) +
    
    geom_segment(aes(x = unique(productivity_cpd_df$max_cpd), 
                     xend = unique(productivity_cpd_df$max_cpd), 
                     y = -100, 
                     yend = productivity_cpd_median[which.min(abs(cpd - unique(productivity_cpd_df$max_cpd)))],
                     color = "gray"), 
                 linetype = "dashed", linewidth = 0.8) +
    # add horizontal line 
    geom_segment(aes(x = 0, 
                     xend = unique(productivity_cpd_df$max_cpd), 
                     y = productivity_cpd_median[which.min(abs(cpd - unique(productivity_cpd_df$max_cpd)))], 
                     yend = productivity_cpd_median[which.min(abs(cpd - unique(productivity_cpd_df$max_cpd)))],
                     color = "gray"), 
                 linetype = "dashed", linewidth = 0.8) +
    annotate("text",x = unique(productivity_cpd_df$max_cpd), 
             y = productivity_cpd_df$productivity_cpd_median[which.min(abs(cpd - unique(productivity_cpd_df$max_cpd)))],
             label = paste0("Max CPD: ",round(unique(productivity_cpd_df$max_cpd),1),"%\n",
                            "Median change: ",round(productivity_cpd_df$productivity_cpd_median[which.min(abs(cpd - unique(productivity_cpd_df$max_cpd)))],1),"%"),
             vjust = -0.5, color = "black", size = 2, hjust = 0
    ) +
    
    scale_color_manual(name = "Model",
                       values = c("ECA model" = '#A2C5AC', 
                                  "median" = '#A2C5AC')) +
    scale_fill_manual(name = "Model",
                      values = c("ECA 95% credible interval" = '#A2C5AC', 
                                 "ECA 80% credible interval" = '#A2C5AC',
                                 "95% credible interval" = '#A2C5AC', 
                                 "80% credible interval" = '#A2C5AC')) +
    ylim(-100,100) +
    scale_x_continuous(n.breaks = 5) +
    labs(x = "Cumulative Disturbance (%)",
         y = "Change in recruitment (%)") +
    theme_classic() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.5, "cm"),
          legend.spacing.y = unit(0.1, "cm"),
          legend.key.height = unit(0.5, "cm"),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(size = 10, hjust = 0.5))+
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
  
  plot2 <- ggplot(productivity_eca_df)+
    geom_line(aes(x = eca, y = productivity_eca_median, group = 1,
                  color = "median"),alpha=0.9, linewidth = 0.8) +
    geom_ribbon(aes(x = eca, ymin = productivity_eca_025_hd, ymax = productivity_eca_975_hd, fill = "95% credible interval"),
                alpha = 0.2) +
    geom_ribbon(aes(x = eca, ymin = productivity_eca_100_hd, ymax = productivity_eca_900_hd, fill = "80% credible interval"),
                alpha = 0.4) +
    #add vertical line at max eca - max height should be value of productivity eca median
    geom_segment(aes(x = unique(productivity_eca_df$max_eca), 
                     xend = unique(productivity_eca_df$max_eca), 
                     y = -100, 
                     yend = productivity_eca_median[which.min(abs(eca - unique(productivity_eca_df$max_eca)))],
                     color = "gray"), 
                 linetype = "dashed", linewidth = 0.8) +
    # add horizontal line 
    geom_segment(aes(x = 0, 
                     xend = unique(productivity_eca_df$max_eca), 
                     y = productivity_eca_median[which.min(abs(eca - unique(productivity_eca_df$max_eca)))], 
                     yend = productivity_eca_median[which.min(abs(eca - unique(productivity_eca_df$max_eca)))],
                     color = "gray"), 
                 linetype = "dashed", linewidth = 0.8) +
    annotate("text",x = unique(productivity_eca_df$max_eca), 
             y = productivity_eca_df$productivity_eca_median[which.min(abs(eca - unique(productivity_eca_df$max_eca)))],
             label = paste0("Max ECA: ",round(unique(productivity_eca_df$max_eca),2),"\n",
                            "Median change: ",round(productivity_eca_df$productivity_eca_median[which.min(abs(eca - unique(productivity_eca_df$max_eca)))],1),"%"),
             vjust = -0.5, color = "black", size = 2, hjust = 0
    ) +
    
    
    scale_color_manual(values = c("median" = '#A2C5AC', 
                                  "CPD model" = '#A2C5AC')) +
    scale_fill_manual(values = c("95% credible interval" = '#A2C5AC', 
                                 "80% credible interval" = '#A2C5AC',
                                 "CPD model 95% credible interval" = '#A2C5AC', 
                                 "CPD model 80% credible interval" = '#A2C5AC')) +
    ylim(-100,100) +
    scale_x_continuous(n.breaks = 5) +
    labs(x = "ECA",
         y = "Change in recruitment (%)") +
    theme_classic() +
    theme(legend.position = "inside",
          legend.justification = c("right", "top"),
          legend.byrow = FALSE,
          legend.title = element_blank(),
          # legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.5, "cm"),
          legend.spacing.y = unit(0, "cm"),
          legend.key.height = unit(0.2, "cm"),
          legend.text = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(size = 10, hjust = 0.5))+
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
  
  return(plot1+plot2 + plot_layout(axes = 'collect') +
           plot_annotation(title = str_to_title(river_name),tag_level = 'A')&
           theme(plot.tag.position = c(0.05, 1),
                 plot.tag = element_text(size = 10, hjust = 0, vjust = 0, face = "bold")))
  
}



plot_both_forestry_effects_river_together(
  river_name = str_to_title(watersheds[1]),
  river = unique(case_study_watersheds_data$River_n)[1]
)


plot_both_forestry_effects_river_together(
  river_name = str_to_title(watersheds[2]),
  river = unique(case_study_watersheds_data$River_n)[2]
)

plot_both_forestry_effects_river_together(
  river_name = str_to_title(watersheds[3]),
  river = unique(case_study_watersheds_data$River_n)[3]
)

for(i in watersheds){
  
  effects_plot <- plot_all_effects_river_together(
    posterior1 = ric_chm_cpd_ocean_covariates_logR_long_chain,
    river_name = str_to_title(i),
    river = unique(case_study_watersheds_data$River_n[case_study_watersheds_data$River == i])
  )
  
  change_plot <- plot_productivity_change_river_together(
    posterior1 = ric_chm_cpd_ocean_covariates_logR_long_chain,
    posterior2 = ric_chm_eca_ocean_covariates_logR_long_chain,
    river_name = i,
    effect1 = "cpd",
    effect2 = "eca",
    species = "chum",
    model1 = "CPD",
    model2 = "ECA"
  )
  
  
  ggsave(filename = here("figures",
                          paste0("case_study_",str_replace_all(str_to_lower(i), " ", "_"),".png")),
         plot = effects_plot+change_plot + plot_layout(widths = c(1,1.2)) +
           plot_annotation(tag_levels = 'A',title = str_to_title(i))&
           theme(plot.tag.position = c(0.0, 1.0),
                 plot.tag = element_text(size = 10, hjust = 0, vjust = 0, face = "bold"),
                 axis.title = element_text(size = 12),
                 axis.text = element_text(size = 10)),
         width = 8,
         height = 4,
         units = "in",
         dpi = 300)
  
  
}

plot_all_effects_river_together(
  posterior1 = ric_chm_cpd_ocean_covariates_logR_long_chain,
  river_name = str_to_title(watersheds[1]),
  river = unique(case_study_watersheds_data$River_n)[1]
)

plot_productivity_change_river_together(
  posterior1 = ric_chm_cpd_ocean_covariates_logR_long_chain,
  posterior2 = ric_chm_eca_ocean_covariates_logR_long_chain,
  river_name = "CARNATION CREEK",
  effect1 = "cpd",
  effect2 = "eca",
  species = "chum",
  model1 = "CPD",
  model2 = "ECA"
)


plot_productivity_change_river_together(
  posterior1 = ric_chm_cpd_ocean_covariates_logR,
  posterior2 = ric_chm_eca_ocean_covariates_logR,
  river_name = "CARNATION CREEK",
  effect1 = "cpd",
  effect2 = "eca",
  species = "chum",
  model1 = "CPD",
  model2 = "ECA"
)




ggsave(filename = here("figures",
                        paste0("case_study_effects_change_",str_replace_all(str_to_lower(watersheds[1]), " ", ""),".png")),
       plot = effects_plot+change_plot + plot_layout(widths = c(1,1.4),tag_level = 'new')&
         plot_annotation(tag_levels = 'A'),
       width = 8,
       height = 3,
       units = "in",
       dpi = 300)        
        
# making changes feb 2026 - changing to ET intervals and using long chain posteriors

ric_chm_eca_ocean_covariates_logR_long_chain=read.csv(here('stan models','outs','posterior','ric_chm_eca_ocean_covariates_logR_long_chain.csv'),check.names=F)
ric_chm_cpd_ocean_covariates_logR_long_chain=read.csv(here('stan models','outs','posterior','ric_chm_cpd_ocean_covariates_logR_long_chain.csv'),check.names=F)

watersheds <- c("VINER SOUND CREEK","CARNATION CREEK", "PHILLIPS RIVER", "NIMPKISH RIVER", "DEENA CREEK", "NEEKAS CREEK")

case_study_watersheds_data <- ch20rsc %>% 
  filter(River %in% watersheds)


plot_recruit_spawner_river(data = case_study_watersheds_data,
                           species = "chum",
                           river_name = watersheds[1],
                           posterior = ric_chm_cpd_ocean_covariates_logR_long_chain,
                           posterior_a_t = NULL,
                           posterior_a_t_bh = NULL)



plot_all_effects_river_together(
  posterior1 = ric_chm_cpd_ocean_covariates_logR_long_chain,
  river_name = str_to_title(watersheds[1]),
  river = unique(case_study_watersheds_data$River_n[case_study_watersheds_data$River == watersheds[1]])
)

plot_productivity_change_river_together(
  posterior1 = ric_chm_cpd_ocean_covariates_logR_long_chain,
  posterior2 = ric_chm_eca_ocean_covariates_logR_long_chain,
  river_name = "CARNATION CREEK",
  effect1 = "cpd",
  effect2 = "eca",
  species = "chum",
  model1 = "CPD",
  model2 = "ECA"
)

plot_recruitment_change_river_together <- function(posterior1 = ric_chm_cpd_ocean_covariates_logR_long_chain,
                                                   posterior2 = ric_chm_eca_ocean_covariates_logR_long_chain,
                                                   river_name = "CARNATION CREEK",  
                                                   effect1 = "cpd",
                                                   effect2 = "eca",
                                                   species = "chum", 
                                                   model1 = "CPD",
                                                   model2 = "ECA",
                                                   hd = FALSE){
  
  
  if(species == "chum"){
    df <- ch20rsc 
    river_data <- ch20rsc %>% filter(River == river_name)
    river <- river_data$River_n[1]
    # df$sst.std <- (ch20rsc$spring_ersst-mean(ch20rsc$spring_ersst))/sd(ch20rsc$spring_ersst)
    
  } else if(species == "pink"){
    df <- pk10r
    river_data <- pk10r %>% filter(River == river_name)
    river <- river_data$River_n[1]
    # df$sst.std <- (pk10r$spring_ersst-mean(pk10r$spring_ersst))/sd(pk10r$spring_ersst)
  }
  
  posterior_df <- posterior1 %>%
    select(starts_with(paste0('b_for_rv[',as.character(river),']'))) %>%
    pivot_longer(cols = everything(), 
                 names_to = 'River', 
                 names_prefix = 'b_for_rv',
                 values_to = "forestry") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>% 
    mutate(River = river)
  
  posterior_df2 <- posterior2 %>% 
    select(starts_with(paste0('b_for_rv[',as.character(river),']'))) %>%
    pivot_longer(cols = everything(), 
                 names_to = 'River', 
                 names_prefix = 'b_for_rv',
                 values_to = "forestry") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>% 
    mutate(River = river)
  
  
  
  
  eca <- seq(0,1,length.out=100)
  
  eca_sqrt <- sqrt(eca)
  
  eca_sqrt_std <- (eca_sqrt-mean(eca_sqrt))/sd(eca_sqrt)
  
  cpd <- seq(0,100,length.out=100)
  
  cpd_sqrt <- sqrt(cpd)
  
  cpd_sqrt_std <- (cpd_sqrt-mean(cpd_sqrt))/sd(cpd_sqrt)
  
  
  
  
  no_eca <- min(eca_sqrt_std)
  
  no_cpd <- min(cpd_sqrt_std)
  
  
  
  # need to change b_rv to posterior 1 and posterior 2 and then make figure
  
  productivity_cpd <- (exp(as.matrix(posterior_df[,1])%*%
                             (cpd_sqrt_std-no_cpd)))*100 - 100
  
  productivity_eca <- (exp(as.matrix(posterior_df2[,1])%*%
                             (eca_sqrt_std-no_eca)))*100 - 100
  
  
  productivity_cpd_df <- data.frame(productivity_cpd_median = apply(productivity_cpd,2,median),
                                    productivity_cpd_025 = apply(productivity_cpd,2,quantile,c(0.025), 
                                                                 row.names = c("q025")),
                                    productivity_cpd_975 = apply(productivity_cpd,2,quantile,c(0.975),
                                                                 row.names = c("q975")),
                                    productivity_cpd_100 = apply(productivity_cpd,2,quantile,c(0.1),
                                                                 row.names = c("q100")),
                                    productivity_cpd_900 = apply(productivity_cpd,2,quantile,c(0.9),
                                                                 row.names = c("q900")),
                                    productivity_cpd_250 = apply(productivity_cpd,2,quantile,c(0.25),
                                                                 row.names = c("q250")),
                                    productivity_cpd_750 = apply(productivity_cpd,2,quantile,c(0.75),
                                                                 row.names = c("q750")),
                                    
                                    
                                    productivity_cpd_025_hd = apply(productivity_cpd,2,HDInterval::hdi, credMass = 0.95)[1,],
                                    productivity_cpd_975_hd = apply(productivity_cpd,2,HDInterval::hdi, credMass = 0.95)[2,],
                                    productivity_cpd_100_hd = apply(productivity_cpd,2,HDInterval::hdi, credMass = 0.8)[1,],
                                    productivity_cpd_900_hd = apply(productivity_cpd,2,HDInterval::hdi, credMass = 0.8)[2,],
                                    productivity_cpd_50_lower_hd = apply(productivity_cpd,2,HDInterval::hdi, credMass = 0.5)[1,],
                                    productivity_cpd_50_upper_hd = apply(productivity_cpd,2,HDInterval::hdi, credMass = 0.5)[2,],
                                    
                                    
                                    
                                    cpd_sqrt_std = cpd_sqrt_std,
                                    cpd = cpd,
                                    max_cpd = max(river_data$disturbedarea_prct_cs),
                                    max_sqrt_cpd = max(river_data$sqrt.CPD),
                                    model = "CPD"
  )
  
  
  productivity_eca_df <- data.frame(productivity_eca_median = apply(productivity_eca,2,median),
                                    productivity_eca_025 = apply(productivity_eca,2,quantile,c(0.025), 
                                                                 row.names = c("q025")),
                                    productivity_eca_975 = apply(productivity_eca,2,quantile,c(0.975),
                                                                 row.names = c("q975")),
                                    productivity_eca_100 = apply(productivity_eca,2,quantile,c(0.1),
                                                                 row.names = c("q100")),
                                    productivity_eca_900 = apply(productivity_eca,2,quantile,c(0.9),
                                                                 row.names = c("q900")),
                                    productivity_eca_250 = apply(productivity_eca,2,quantile,c(0.25),
                                                                 row.names = c("q250")),
                                    productivity_eca_750 = apply(productivity_eca,2,quantile,c(0.75),
                                                                 row.names = c("q750")),
                                    
                                    productivity_eca_025_hd = apply(productivity_eca,2,HDInterval::hdi, credMass = 0.95)[1,],
                                    productivity_eca_975_hd = apply(productivity_eca,2,HDInterval::hdi, credMass = 0.95)[2,],
                                    productivity_eca_100_hd = apply(productivity_eca,2,HDInterval::hdi, credMass = 0.8)[1,],
                                    productivity_eca_900_hd = apply(productivity_eca,2,HDInterval::hdi, credMass = 0.8)[2,],
                                    productivity_eca_50_lower_hd = apply(productivity_eca,2,HDInterval::hdi, credMass = 0.5)[1,],
                                    productivity_eca_50_upper_hd = apply(productivity_eca,2,HDInterval::hdi, credMass = 0.5)[2,],
                                    
                                    eca_sqrt_std = eca_sqrt_std,
                                    eca = eca,
                                    max_eca = max(river_data$ECA_age_proxy_forested_only),
                                    max_sqrt_eca = max(river_data$sqrt.ECA),
                                    model = "ECA"
  )
  
  
  if(hd == TRUE){
    
    plot1 <- ggplot(productivity_cpd_df) +
      geom_line(aes(x = cpd, y = productivity_cpd_median, group = 1,
                    color = "median"),alpha=0.9, linewidth = 0.8) +
      geom_ribbon(aes(x = cpd, ymin = productivity_cpd_025_hd, ymax = productivity_cpd_975_hd, fill = "95% credible interval"),
                  alpha = 0.2) +
      geom_ribbon(aes(x = cpd, ymin = productivity_cpd_100_hd, ymax = productivity_cpd_900_hd, fill = "80% credible interval"),
                  alpha = 0.4) +
      
      geom_segment(aes(x = unique(productivity_cpd_df$max_cpd), 
                       xend = unique(productivity_cpd_df$max_cpd), 
                       y = -100, 
                       yend = productivity_cpd_median[which.min(abs(cpd - unique(productivity_cpd_df$max_cpd)))],
                       color = "gray"), 
                   linetype = "dashed", linewidth = 0.8) +
      # add horizontal line 
      geom_segment(aes(x = 0, 
                       xend = unique(productivity_cpd_df$max_cpd), 
                       y = productivity_cpd_median[which.min(abs(cpd - unique(productivity_cpd_df$max_cpd)))], 
                       yend = productivity_cpd_median[which.min(abs(cpd - unique(productivity_cpd_df$max_cpd)))],
                       color = "gray"), 
                   linetype = "dashed", linewidth = 0.8) +
      annotate("text",x = unique(productivity_cpd_df$max_cpd), 
               y = productivity_cpd_df$productivity_cpd_median[which.min(abs(cpd - unique(productivity_cpd_df$max_cpd)))],
               label = paste0("Max CPD: ",round(unique(productivity_cpd_df$max_cpd),1),"%\n",
                              "Median change: ",round(productivity_cpd_df$productivity_cpd_median[which.min(abs(cpd - unique(productivity_cpd_df$max_cpd)))],1),"%"),
               vjust = -0.5, color = "black", size = 2, hjust = 0
      ) +
      
      scale_color_manual(name = "Model",
                         values = c("ECA model" = '#A2C5AC', 
                                    "median" = '#A2C5AC')) +
      scale_fill_manual(name = "Model",
                        values = c("ECA 95% credible interval" = '#A2C5AC', 
                                   "ECA 80% credible interval" = '#A2C5AC',
                                   "95% credible interval" = '#A2C5AC', 
                                   "80% credible interval" = '#A2C5AC')) +
      ylim(-100,100) +
      scale_x_continuous(n.breaks = 5) +
      labs(x = "CDA (%)",
           y = "Change in recruitment (%)") +
      theme_classic() +
      theme(legend.position = "none",
            legend.title = element_blank(),
            legend.key.size = unit(0.5, "cm"),
            legend.key.width = unit(0.5, "cm"),
            legend.spacing.y = unit(0.1, "cm"),
            legend.key.height = unit(0.5, "cm"),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 10, hjust = 0.5))+
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    
    plot2 <- ggplot(productivity_eca_df)+
      geom_line(aes(x = eca, y = productivity_eca_median, group = 1,
                    color = "median"),alpha=0.9, linewidth = 0.8) +
      geom_ribbon(aes(x = eca, ymin = productivity_eca_025_hd, ymax = productivity_eca_975_hd, fill = "95% credible interval"),
                  alpha = 0.2) +
      geom_ribbon(aes(x = eca, ymin = productivity_eca_100_hd, ymax = productivity_eca_900_hd, fill = "80% credible interval"),
                  alpha = 0.4) +
      #add vertical line at max eca - max height should be value of productivity eca median
      geom_segment(aes(x = unique(productivity_eca_df$max_eca), 
                       xend = unique(productivity_eca_df$max_eca), 
                       y = -100, 
                       yend = productivity_eca_median[which.min(abs(eca - unique(productivity_eca_df$max_eca)))],
                       color = "gray"), 
                   linetype = "dashed", linewidth = 0.8) +
      # add horizontal line 
      geom_segment(aes(x = 0, 
                       xend = unique(productivity_eca_df$max_eca), 
                       y = productivity_eca_median[which.min(abs(eca - unique(productivity_eca_df$max_eca)))], 
                       yend = productivity_eca_median[which.min(abs(eca - unique(productivity_eca_df$max_eca)))],
                       color = "gray"), 
                   linetype = "dashed", linewidth = 0.8) +
      annotate("text",x = unique(productivity_eca_df$max_eca), 
               y = productivity_eca_df$productivity_eca_median[which.min(abs(eca - unique(productivity_eca_df$max_eca)))],
               label = paste0("Max ECA: ",round(unique(productivity_eca_df$max_eca),2),"\n",
                              "Median change: ",round(productivity_eca_df$productivity_eca_median[which.min(abs(eca - unique(productivity_eca_df$max_eca)))],1),"%"),
               vjust = -0.5, color = "black", size = 2, hjust = 0
      ) +
      
      
      scale_color_manual(values = c("median" = '#A2C5AC', 
                                    "CDA model" = '#A2C5AC')) +
      scale_fill_manual(values = c("95% credible interval" = '#A2C5AC', 
                                   "80% credible interval" = '#A2C5AC',
                                   "CPD model 95% credible interval" = '#A2C5AC', 
                                   "CPD model 80% credible interval" = '#A2C5AC')) +
      ylim(-100,100) +
      scale_x_continuous(n.breaks = 5) +
      labs(x = "ECA",
           y = "Change in recruitment (%)") +
      theme_classic() +
      theme(legend.position = "inside",
            legend.justification = c("right", "top"),
            legend.byrow = FALSE,
            legend.title = element_blank(),
            # legend.key.size = unit(0.5, "cm"),
            legend.key.width = unit(0.5, "cm"),
            legend.spacing.y = unit(0, "cm"),
            legend.key.height = unit(0.2, "cm"),
            legend.text = element_text(size = 6),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 10, hjust = 0.5))+
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    
  } else{
    plot1 <- ggplot(productivity_cpd_df) +
      geom_line(aes(x = cpd, y = productivity_cpd_median, group = 1,
                    color = "median"),alpha=0.9, linewidth = 0.8) +
      geom_ribbon(aes(x = cpd, ymin = productivity_cpd_025, ymax = productivity_cpd_975, fill = "95% credible interval"),
                  alpha = 0.2) +
      geom_ribbon(aes(x = cpd, ymin = productivity_cpd_100, ymax = productivity_cpd_900, fill = "80% credible interval"),
                  alpha = 0.4) +
      geom_ribbon(aes(x = cpd, ymin = productivity_cpd_250, ymax = productivity_cpd_750, fill = "50% credible interval"),
                  alpha = 0.6) +
      
      geom_segment(aes(x = unique(productivity_cpd_df$max_cpd), 
                       xend = unique(productivity_cpd_df$max_cpd), 
                       y = -100, 
                       yend = productivity_cpd_median[which.min(abs(cpd - unique(productivity_cpd_df$max_cpd)))],
                       color = "gray"), 
                   linetype = "dashed", linewidth = 0.8) +
      # add horizontal line 
      geom_segment(aes(x = 0, 
                       xend = unique(productivity_cpd_df$max_cpd), 
                       y = productivity_cpd_median[which.min(abs(cpd - unique(productivity_cpd_df$max_cpd)))], 
                       yend = productivity_cpd_median[which.min(abs(cpd - unique(productivity_cpd_df$max_cpd)))],
                       color = "gray"), 
                   linetype = "dashed", linewidth = 0.8) +
      annotate("text",x = unique(productivity_cpd_df$max_cpd), 
               y = productivity_cpd_df$productivity_cpd_median[which.min(abs(cpd - unique(productivity_cpd_df$max_cpd)))],
               label = paste0("Max CDA: ",round(unique(productivity_cpd_df$max_cpd),1),"%\n",
                              "Median change: ",round(productivity_cpd_df$productivity_cpd_median[which.min(abs(cpd - unique(productivity_cpd_df$max_cpd)))],1),"%"),
               vjust = -0.5, color = "black", size = 2, hjust = 0
      ) +
      
      scale_color_manual(name = "Model",
                         values = c("ECA model" = '#A2C5AC', 
                                    "median" = '#A2C5AC')) +
      scale_fill_manual(name = "Model",
                        values = c("ECA 95% credible interval" = '#A2C5AC', 
                                   "ECA 80% credible interval" = '#A2C5AC',
                                   "95% credible interval" = '#A2C5AC', 
                                   "80% credible interval" = '#A2C5AC',
                                   "50% credible interval" = '#A2C5AC')) +
      ylim(-100,100) +
      scale_x_continuous(n.breaks = 5) +
      labs(x = "CDA (%)",
           y = "Change in recruitment (%)") +
      theme_classic() +
      theme(legend.position = "none",
            legend.title = element_blank(),
            legend.key.size = unit(0.5, "cm"),
            legend.key.width = unit(0.5, "cm"),
            legend.spacing.y = unit(0.1, "cm"),
            legend.key.height = unit(0.5, "cm"),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 10, hjust = 0.5))+
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    
    plot2 <- ggplot(productivity_eca_df)+
      geom_line(aes(x = eca, y = productivity_eca_median, group = 1,
                    color = "median"),alpha=0.9, linewidth = 0.8) +
      geom_ribbon(aes(x = eca, ymin = productivity_eca_025, ymax = productivity_eca_975, fill = "95% credible interval"),
                  alpha = 0.2) +
      geom_ribbon(aes(x = eca, ymin = productivity_eca_100, ymax = productivity_eca_900, fill = "80% credible interval"),
                  alpha = 0.4) +
      geom_ribbon(aes(x = eca, ymin = productivity_eca_250, ymax = productivity_eca_750, fill = "50% credible interval"),
                  alpha = 0.6) +
      #add vertical line at max eca - max height should be value of productivity eca median
      geom_segment(aes(x = unique(productivity_eca_df$max_eca), 
                       xend = unique(productivity_eca_df$max_eca), 
                       y = -100, 
                       yend = productivity_eca_median[which.min(abs(eca - unique(productivity_eca_df$max_eca)))],
                       color = "gray"), 
                   linetype = "dashed", linewidth = 0.8) +
      # add horizontal line 
      geom_segment(aes(x = 0, 
                       xend = unique(productivity_eca_df$max_eca), 
                       y = productivity_eca_median[which.min(abs(eca - unique(productivity_eca_df$max_eca)))], 
                       yend = productivity_eca_median[which.min(abs(eca - unique(productivity_eca_df$max_eca)))],
                       color = "gray"), 
                   linetype = "dashed", linewidth = 0.8) +
      annotate("text",x = unique(productivity_eca_df$max_eca), 
               y = productivity_eca_df$productivity_eca_median[which.min(abs(eca - unique(productivity_eca_df$max_eca)))],
               label = paste0("Max ECA: ",round(unique(productivity_eca_df$max_eca),2),"\n",
                              "Median change: ",round(productivity_eca_df$productivity_eca_median[which.min(abs(eca - unique(productivity_eca_df$max_eca)))],1),"%"),
               vjust = -0.5, color = "black", size = 2, hjust = 0
      ) +
      
      
      scale_color_manual(values = c("median" = '#A2C5AC', 
                                    "CPD model" = '#A2C5AC')) +
      scale_fill_manual(values = c("95% credible interval" = '#A2C5AC', 
                                   "80% credible interval" = '#A2C5AC',
                                   "50% credible interval" = '#A2C5AC',
                                   "CPD model 95% credible interval" = '#A2C5AC', 
                                   "CPD model 80% credible interval" = '#A2C5AC')) +
      ylim(-100,100) +
      scale_x_continuous(n.breaks = 5) +
      labs(x = "ECA",
           y = "Change in recruitment (%)") +
      theme_classic() +
      theme(legend.position = "inside",
            legend.justification = c("right", "top"),
            legend.byrow = FALSE,
            legend.title = element_blank(),
            # legend.key.size = unit(0.5, "cm"),
            legend.key.width = unit(0.5, "cm"),
            legend.spacing.y = unit(0, "cm"),
            legend.key.height = unit(0.2, "cm"),
            legend.text = element_text(size = 6),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 10, hjust = 0.5))+
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
  }
  
 
  
  return(plot1+plot2 + plot_layout(axes = 'collect') +
           plot_annotation(title = str_to_title(river_name),tag_level = 'A')&
           theme(plot.tag.position = c(0.05, 1),
                 plot.tag = element_text(size = 10, hjust = 0, vjust = 0, face = "bold")))
  
}


plot_recruitment_change_river_together(
  posterior1 = ric_chm_cpd_ocean_covariates_logR_long_chain,
  posterior2 = ric_chm_eca_ocean_covariates_logR_long_chain,
  river_name = "CARNATION CREEK",
  effect1 = "cpd",
  effect2 = "eca",
  species = "chum",
  model1 = "CPD",
  model2 = "ECA",
  hd = FALSE
)


plot_recruit_spawner_river(data = case_study_watersheds_data,
                           species = "chum",
                           river_name = watersheds[6],
                           posterior = ric_chm_cpd_ocean_covariates_logR_long_chain,
                           posterior_a_t = NULL,
                           posterior_a_t_bh = NULL)

plot_recruit_spawner_river_new <- function(species = "chum", river_name,  posterior){
  
  # river_data <- df %>% filter(River_n == river)
  
  if(species == "chum"){
    df <- ch20rsc 
    river_data <- ch20rsc %>% filter(River == river_name)
    river <- river_data$River_n[1]
    # df$sst.std <- (ch20rsc$spring_ersst-mean(ch20rsc$spring_ersst))/sd(ch20rsc$spring_ersst)
    
  } else if(species == "pink"){
    df <- pk10r
    river_data <- pk10r %>% filter(River == river_name)
    river <- river_data$River_n2[1]
    # df$sst.std <- (pk10r$spring_ersst-mean(pk10r$spring_ersst))/sd(pk10r$spring_ersst)
  }
  
  if(species == "chum"){
    posterior_rv_b_for <- posterior %>% 
      select(starts_with('b_for_rv')) %>%
      select(ends_with(paste0("[",river,"]")))
    
    posterior_rv_alpha_j <- posterior %>% 
      select(starts_with('alpha_j')) %>%
      select(ends_with(paste0("[",river,"]")))
    
    posterior_rv_S_max <- posterior %>% 
      select(starts_with('Smax')) %>%
      select(ends_with(paste0("[",river,"]")))
    
    
  } else if(species == "pink"){
    
    river_wo_broodline <- river_data$River_n
    river_w_broodline <- river_data$River_n2
    
    posterior_rv_b_for <- posterior %>% 
      select(starts_with('b_for_rv')) %>%
      select(ends_with(paste0("[",river_wo_broodline,"]")))
    
    posterior_rv_alpha_j <- posterior %>% 
      select(starts_with('alpha_j')) %>%
      select(ends_with(paste0("[",river_w_broodline,"]")))
    
    posterior_rv_S_max <- posterior %>% 
      select(starts_with('Smax')) %>%
      select(ends_with(paste0("[",river_w_broodline,"]")))
    
    
    
  }
  
  
  spawners_predicted <- seq(0, max(river_data$Spawners), length.out = 100)
  
  # calculate recruit prediction
  
  low_cpd <- min(river_data$sqrt.CPD.std)
  high_cpd <- max(river_data$sqrt.CPD.std)
  # avg_cpd <- mean(river_data$sqrt.CPD.std)
  mid_cpd <- min(river_data$sqrt.CPD.std) + (max(river_data$sqrt.CPD.std) - min(river_data$sqrt.CPD.std))/2
  
  mid_cpd_real <- min(river_data$disturbedarea_prct_cs) + (max(river_data$disturbedarea_prct_cs) - min(river_data$disturbedarea_prct_cs))/2
  
  
  recruits_predicted_low_cpd <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*low_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, median))*spawners_predicted
  
  recruits_predicted_high_cpd <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*high_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, median))*spawners_predicted
  
  
  recruits_predicted_low_cpd_lower <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*low_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, quantile, c(0.025)))*spawners_predicted
  
  recruits_predicted_low_cpd_upper <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*low_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, quantile, c(0.975)))*spawners_predicted
  
  recruits_predicted_high_cpd_lower <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*high_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, quantile, c(0.025)))*spawners_predicted
  
  recruits_predicted_high_cpd_upper <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*high_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, quantile, c(0.975)))*spawners_predicted
  
  
  
  #make dataframe
  
  prediction_df <- data.frame(spawners = spawners_predicted,
                              
                              recruits_low_cpd = recruits_predicted_low_cpd,
                              recruits_high_cpd = recruits_predicted_high_cpd,
                              recruits_low_cpd_lower = recruits_predicted_low_cpd_lower,
                              recruits_low_cpd_upper = recruits_predicted_low_cpd_upper,
                              recruits_high_cpd_lower = recruits_predicted_high_cpd_lower,
                              recruits_high_cpd_upper = recruits_predicted_high_cpd_upper)
  
  
  
  
  
  
  
  #plot the time varying productivity vs year, with log(R/S) data
  
  
  
  
  
  #plot recruit vs spawner as points
  
  p1 <- ggplot() +
    geom_point(data = river_data,aes(x = Spawners, y = Recruits, color = disturbedarea_prct_cs), alpha = 0.5, size = 2) +
    geom_line(data = prediction_df, aes(x = spawners, y = recruits_low_cpd), color = '#35978f', size = 1, alpha = 0.5) +
    geom_line(data = prediction_df, aes(x = spawners, y = recruits_high_cpd), color = '#bf812d', size = 1, alpha = 0.5) +
    # geom_line(data = prediction_df, aes(x = spawners, y = recruits), color = "black", size = 1, alpha = 0.5) +
    geom_ribbon(data = prediction_df, aes(x = spawners, ymin = recruits_low_cpd_lower, ymax = recruits_low_cpd_upper), fill ='#35978f', alpha = 0.1) +
    geom_ribbon(data = prediction_df, aes(x = spawners, ymin = recruits_high_cpd_lower, ymax = recruits_high_cpd_upper), fill ='#bf812d', alpha = 0.1) +
    #make y log scale
    # scale_y_log10() +
    labs(title = "", x = "Spawners", y = "Recruits") +
    # scale_color_manual(name = "CPD", values = c("Low" = '#35978f', "High" = '#bf812d')) +
    scale_color_gradient2(name = 'CDA (%)', guide = guide_colorbar(barwidth = 3, barheight = 0.5),
                          low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 20, n.breaks = 4) +
    theme_classic() +
    theme(legend.position = c(0.8,0.9),
          legend.background = element_rect(fill = alpha('white', 0.5)),
          legend.direction = "horizontal",
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8, vjust = 1, hjust = 1),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_blank()
          # legend.key.width = unit(0.5, "cm"),
          # legend.key.height = unit(1, "lines"),
          # legend.spacing.y = unit(0.001, "cm")
    )
  
  
  # p3 log R/s vs spawners
  
  p2 <- ggplot(river_data) + 
    geom_point(aes(x = Spawners, y = log(Recruits/Spawners), color = disturbedarea_prct_cs), alpha = 0.5, size = 2) +
    geom_line(data = prediction_df, aes(x = spawners, y = log(recruits_predicted_low_cpd/spawners_predicted)), color = '#35978f', size = 1, alpha = 0.5) +
    geom_line(data = prediction_df, aes(x = spawners, y = log(recruits_predicted_high_cpd/spawners_predicted)), color = '#bf812d', size = 1, alpha = 0.5) +
    # geom_line(data = prediction_df, aes(x = spawners, y = log_RS), color = "black", size = 1, alpha = 0.5) +
    geom_ribbon(data = prediction_df, aes(x = spawners, ymin = log(recruits_predicted_low_cpd_lower/spawners_predicted), ymax = log(recruits_predicted_low_cpd_upper/spawners_predicted)), fill = "#35978f", alpha = 0.1) +
    geom_ribbon(data = prediction_df, aes(x = spawners, ymin = log(recruits_predicted_high_cpd_lower/spawners_predicted), ymax = log(recruits_predicted_high_cpd_upper/spawners_predicted)), fill = "#bf812d", alpha = 0.1) +
    labs(#title = "Ricker model with CPD, NPGO, ERSST", 
      x = "Spawners", y = TeX(r"($\log \left(\frac{Recruits}{Spawners}\right)$)")) +
    scale_color_gradient2(name = 'CDA (%)',
                          low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 20, guide = guide_colorbar(barwidth = 3, barheight = 0.5)) +
    theme_classic() +
    theme(legend.position = c(0.8,0.9),
          legend.background = element_rect(fill = alpha('white', 0.5)),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_blank(),
          legend.key.width = unit(0.5, "cm"),
          legend.key.height = unit(1, "lines"),
          legend.spacing.y = unit(0.001, "cm")
    )
  
  p3 <- ggplot(river_data) + 
    geom_point(aes(x = Spawners, y = log(Recruits), color = disturbedarea_prct_cs), alpha = 0.5, size = 2) +
    geom_line(data = prediction_df, aes(x = spawners, y = log(recruits_predicted_low_cpd)), color = '#35978f', size = 1, alpha = 0.5) +
    geom_line(data = prediction_df, aes(x = spawners, y = log(recruits_predicted_high_cpd)), color = '#bf812d', size = 1, alpha = 0.5) +
    # geom_line(data = prediction_df, aes(x = spawners, y = log_RS), color = "black", size = 1, alpha = 0.5) +
    geom_ribbon(data = prediction_df, aes(x = spawners, ymin = log(recruits_predicted_low_cpd_lower), ymax = log(recruits_predicted_low_cpd_upper)), fill = "#35978f", alpha = 0.1) +
    geom_ribbon(data = prediction_df, aes(x = spawners, ymin = log(recruits_predicted_high_cpd_lower), ymax = log(recruits_predicted_high_cpd_upper)), fill = "#bf812d", alpha = 0.1) +
    labs(#title = "Ricker model with CPD, NPGO, ERSST", 
      x = "Spawners", y = TeX(r"($\log \left(Recruits\right)$)")) +
    scale_color_gradient2(name = 'CDA (%)',
                          low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 20)+
    theme_classic() +
    theme(legend.position = c(0.8,0.9),
          legend.background = element_rect(fill = alpha('white', 0.5)),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_blank(),
          legend.key.width = unit(0.5, "cm"),
          legend.key.height = unit(1, "lines"),
          legend.spacing.y = unit(0.001, "cm")
    )
  
  
  
  
  
  
  return((p1))
}



for(i in watersheds){
  
  effects_plot <- plot_all_effects_river_together(
    posterior1 = ric_chm_cpd_ocean_covariates_logR_long_chain,
    river_name = str_to_title(i),
    river = unique(case_study_watersheds_data$River_n[case_study_watersheds_data$River == i])
  )
  
  change_plot <- plot_recruitment_change_river_together(
    posterior1 = ric_chm_cpd_ocean_covariates_logR_long_chain,
    posterior2 = ric_chm_eca_ocean_covariates_logR_long_chain,
    river_name = i,
    effect1 = "cpd",
    effect2 = "eca",
    species = "chum",
    model1 = "CPD",
    model2 = "ECA",
    hd = FALSE
  )
  
  spawner_recruit_plot <- plot_recruit_spawner_river_new(species = "chum",
                                 river_name = i,
                                 posterior = ric_chm_cpd_ocean_covariates_logR_long_chain)
  
  ggsave(filename = here("figures",
                         paste0("case_study_et_w_sr_curves_",str_replace_all(str_to_lower(i), " ", "_"),".png")),
         plot = ((effects_plot+spawner_recruit_plot + plot_layout(widths = c(1,1.2)))/change_plot) +
           plot_annotation(tag_levels = 'A',title = paste(str_to_title(i), "- Chum Salmon"))&
           theme(plot.tag.position = c(0.0, 1.0),
                 plot.tag = element_text(size = 10, hjust = 0, vjust = 0, face = "bold")),
         width = 7,
         height = 6,
         units = "in",
         dpi = 300)
  
  
}



# Pink --------------------------------------------------------------------



pk10r_e <- read.csv(here("origional-ecofish-data-models","Data","Processed","pke_SR_10_hat_yr_w_ersst.csv"))

#odd year pinks
pk10r_o <- read.csv(here("origional-ecofish-data-models","Data","Processed","pko_SR_10_hat_yr_w_ersst.csv"))

options(mc.cores=8)

# Pink salmon - even/odd croodlines #####
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',pk10r_o$River)
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=='915-765500-18600-00000-0000-0000-000-000-000-000-000-000','HEAD CREEK 2',pk10r_o$River)
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=='915-488000-41400-00000-0000-0000-000-000-000-000-000-000','WINDY cAY CREEK 2',pk10r_o$River)
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',pk10r_o$River)
pk10r_o=pk10r_o[order(factor(pk10r_o$River),pk10r_o$BroodYear),]
rownames(pk10r_o)=seq(1:nrow(pk10r_o))

pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',pk10r_e$River)
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='915-765500-18600-00000-0000-0000-000-000-000-000-000-000','HEAD CREEK 2',pk10r_e$River)
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='915-488000-41400-00000-0000-0000-000-000-000-000-000-000','WINDY cAY CREEK 2',pk10r_e$River)
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


#standardize npgo
pk10r_o$npgo.std=(pk10r_o$npgo-mean(pk10r_o$npgo))/sd(pk10r_o$npgo)
pk10r_e$npgo.std=(pk10r_e$npgo-mean(pk10r_e$npgo))/sd(pk10r_e$npgo)

# pk10r_o$sst.std = (pk10r_o$spring_lighthouse_temperature-mean(pk10r_o$spring_lighthouse_temperature))/sd(pk10r_o$spring_lighthouse_temperature)
# pk10r_e$sst.std = (pk10r_e$spring_lighthouse_temperature-mean(pk10r_e$spring_lighthouse_temperature))/sd(pk10r_e$spring_lighthouse_temperature)

pk10r_o$sst.std = (pk10r_o$spring_ersst-mean(pk10r_o$spring_ersst))/sd(pk10r_o$spring_ersst)
pk10r_e$sst.std = (pk10r_e$spring_ersst-mean(pk10r_e$spring_ersst))/sd(pk10r_e$spring_ersst)


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
# N_s=rag_n(pk10r$River2)

#cus by stock
cu1=distinct(pk10r,CU,.keep_all = T)
cu2=distinct(pk10r,River,.keep_all = T)
cu3=distinct(pk10r,River2,.keep_all = T)

pk10r$River_n <- as.numeric(factor(pk10r$River))
pk10r$River_n2 <- as.numeric(factor(pk10r$River2))
pk10r$CU_name <- pk10r$CU_NAME

max_eca_pink_df <- pk10r %>% 
  select(River2, River, River_n, CU, ECA_age_proxy_forested_only) %>%
  group_by(River) %>%
  summarize(River = first(River),
            River_n = first(River_n),
            CU = first(CU),
            eca_max = 100*max(ECA_age_proxy_forested_only, na.rm =TRUE)) %>% 
  mutate(eca_level = case_when(eca_max < 12 ~ 'low',
                               eca_max >= 12 & eca_max < 24 ~ 'medium',
                               eca_max >= 24 ~ 'high'))

max_cpd_pink_df <- pk10r %>%
  select(River2, River, River_n,  CU, disturbedarea_prct_cs) %>%
  group_by(River) %>% 
  summarize(River = first(River),
            River_n = first(River_n),
            CU = first(CU),
            cpd_max = max(disturbedarea_prct_cs, na.rm = TRUE))


ric_pk_cpd_ersst_long_chain <- read.csv(here('stan models','outs','posterior',"ric_pk_cpd_st_noac_ocean_covariates_logR_long_chain.csv"),
                                        check.names = FALSE)

ric_pk_eca_ersst_long_chain <- read.csv(here('stan models','outs','posterior',"ric_pk_eca_st_noac_ocean_covariates_logR_long_chain.csv"),
                                        check.names = FALSE)


# do same for pink but only for some watersheds - Phillips River, Deena Creek, Neekas
pink_watersheds_w_broodline = c("PHILLIPS RIVER_Even", 
                    "PHILLIPS RIVER_Odd",
                    "DEENA CREEK_Even",
                    "DEENA CREEK_Odd", 
                    "NEEKAS CREEK_Even",
                    "NEEKAS CREEK_Odd")

pink_watersheds <- c("PHILLIPS RIVER", 
                    "DEENA CREEK", 
                    "NEEKAS CREEK")


case_study_watersheds_data_pk <- pk10r %>% 
  filter(River2 %in% pink_watersheds)





plot_recruit_spawner_river_pink <- function(species = "pink", river_name,  posterior){
  
  # river_data <- df %>% filter(River_n == river)
  even_or_odd <- str_extract(river_name, "Even|Odd")
  
  df <- pk10r
  river_data <- pk10r %>% filter(River2 == river_name)
  
  # df$sst.std <- (pk10r$spring_ersst-mean(pk10r$spring_ersst))/sd(pk10r$spring_ersst)
  
  
  river_wo_broodline <- unique(river_data$River_n)
  river_w_broodline <- unique(river_data$River_n2)
  
  posterior_rv_b_for <- posterior %>% 
    select(starts_with('b_for_rv')) %>%
    select(ends_with(paste0("[",river_wo_broodline,"]")))
  
  posterior_rv_alpha_j <- posterior %>% 
    select(starts_with('alpha_j')) %>%
    select(ends_with(paste0("[",river_w_broodline,"]")))
  
  posterior_rv_S_max <- posterior %>% 
    select(starts_with('Smax')) %>%
    select(ends_with(paste0("[",river_w_broodline,"]")))
  
  # for(river in river_w_broodline_both){
  #   
  #   posterior_rv_alpha_j <- posterior %>% 
  #     select(starts_with('alpha_j')) %>%
  #     select(ends_with(paste0("[",river_w_broodline,"]")))
  #   
  #   posterior_rv_S_max <- posterior %>% 
  #     select(starts_with('Smax')) %>%
  #     select(ends_with(paste0("[",river_w_broodline,"]")))
  #   
  # }
  
  
  
  
  
  
  
  
  spawners_predicted <- seq(0, max(river_data$Spawners), length.out = 100)
  
  # calculate recruit prediction
  
  low_cpd <- min(river_data$sqrt.CPD.std)
  high_cpd <- max(river_data$sqrt.CPD.std)
  # avg_cpd <- mean(river_data$sqrt.CPD.std)
  mid_cpd <- min(river_data$sqrt.CPD.std) + (max(river_data$sqrt.CPD.std) - min(river_data$sqrt.CPD.std))/2
  
  mid_cpd_real <- min(river_data$disturbedarea_prct_cs) + (max(river_data$disturbedarea_prct_cs) - min(river_data$disturbedarea_prct_cs))/2
  
  
  recruits_predicted_low_cpd <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*low_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, median))*spawners_predicted
  
  recruits_predicted_high_cpd <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*high_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, median))*spawners_predicted
  
  
  recruits_predicted_low_cpd_lower <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*low_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, quantile, c(0.025)))*spawners_predicted
  
  recruits_predicted_low_cpd_upper <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*low_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, quantile, c(0.975)))*spawners_predicted
  
  recruits_predicted_high_cpd_lower <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*high_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, quantile, c(0.025)))*spawners_predicted
  
  recruits_predicted_high_cpd_upper <- exp(apply((matrix(posterior_rv_alpha_j[,1], ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1])) - as.matrix(1/posterior_rv_S_max)%*%spawners_predicted  + matrix(posterior_rv_b_for[,1]*high_cpd, ncol = length(spawners_predicted), nrow = length(posterior_rv_alpha_j[,1]))), 2, quantile, c(0.975)))*spawners_predicted
  
  
  
  #make dataframe
  
  prediction_df <- data.frame(spawners = spawners_predicted,
                              
                              recruits_low_cpd = recruits_predicted_low_cpd,
                              recruits_high_cpd = recruits_predicted_high_cpd,
                              recruits_low_cpd_lower = recruits_predicted_low_cpd_lower,
                              recruits_low_cpd_upper = recruits_predicted_low_cpd_upper,
                              recruits_high_cpd_lower = recruits_predicted_high_cpd_lower,
                              recruits_high_cpd_upper = recruits_predicted_high_cpd_upper)
  
  
  
  
  
  
  
  #plot the time varying productivity vs year, with log(R/S) data
  
  
  
  
  
  #plot recruit vs spawner as points
  
  p1 <- ggplot() +
    geom_point(data = river_data,aes(x = Spawners, y = Recruits, color = disturbedarea_prct_cs), alpha = 0.5, size = 2) +
    geom_line(data = prediction_df, aes(x = spawners, y = recruits_low_cpd), color = '#35978f', size = 1, alpha = 0.5) +
    geom_line(data = prediction_df, aes(x = spawners, y = recruits_high_cpd), color = '#bf812d', size = 1, alpha = 0.5) +
    # geom_line(data = prediction_df, aes(x = spawners, y = recruits), color = "black", size = 1, alpha = 0.5) +
    geom_ribbon(data = prediction_df, aes(x = spawners, ymin = recruits_low_cpd_lower, ymax = recruits_low_cpd_upper), fill ='#35978f', alpha = 0.1) +
    geom_ribbon(data = prediction_df, aes(x = spawners, ymin = recruits_high_cpd_lower, ymax = recruits_high_cpd_upper), fill ='#bf812d', alpha = 0.1) +
    #make y log scale
    # scale_y_log10() +
    labs(title = even_or_odd, x = "Spawners", y = "Recruits") +
    # scale_color_manual(name = "CPD", values = c("Low" = '#35978f', "High" = '#bf812d')) +
    scale_color_gradient2(name = 'CDA (%)', guide = guide_colorbar(barwidth = 3, barheight = 0.5),
                          low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 20, n.breaks = 4) +
    theme_classic() +
    theme(legend.position = c(0.8,0.9),
          legend.background = element_rect(fill = alpha('white', 0.5)),
          legend.direction = "horizontal",
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8, vjust = 1, hjust = 1),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(size = 9)
          # legend.key.width = unit(0.5, "cm"),
          # legend.key.height = unit(1, "lines"),
          # legend.spacing.y = unit(0.001, "cm")
    )
  
  
  # p3 log R/s vs spawners
  
  p2 <- ggplot(river_data) + 
    geom_point(aes(x = Spawners, y = log(Recruits/Spawners), color = disturbedarea_prct_cs), alpha = 0.5, size = 2) +
    geom_line(data = prediction_df, aes(x = spawners, y = log(recruits_predicted_low_cpd/spawners_predicted)), color = '#35978f', size = 1, alpha = 0.5) +
    geom_line(data = prediction_df, aes(x = spawners, y = log(recruits_predicted_high_cpd/spawners_predicted)), color = '#bf812d', size = 1, alpha = 0.5) +
    # geom_line(data = prediction_df, aes(x = spawners, y = log_RS), color = "black", size = 1, alpha = 0.5) +
    geom_ribbon(data = prediction_df, aes(x = spawners, ymin = log(recruits_predicted_low_cpd_lower/spawners_predicted), ymax = log(recruits_predicted_low_cpd_upper/spawners_predicted)), fill = "#35978f", alpha = 0.1) +
    geom_ribbon(data = prediction_df, aes(x = spawners, ymin = log(recruits_predicted_high_cpd_lower/spawners_predicted), ymax = log(recruits_predicted_high_cpd_upper/spawners_predicted)), fill = "#bf812d", alpha = 0.1) +
    labs(#title = "Ricker model with CPD, NPGO, ERSST", 
      x = "Spawners", y = TeX(r"($\log \left(\frac{Recruits}{Spawners}\right)$)")) +
    scale_color_gradient2(name = 'CDA (%)',
                          low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 20, guide = guide_colorbar(barwidth = 3, barheight = 0.5)) +
    theme_classic() +
    theme(legend.position = c(0.8,0.9),
          legend.background = element_rect(fill = alpha('white', 0.5)),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_blank(),
          legend.key.width = unit(0.5, "cm"),
          legend.key.height = unit(1, "lines"),
          legend.spacing.y = unit(0.001, "cm")
    )
  
  p3 <- ggplot(river_data) + 
    geom_point(aes(x = Spawners, y = log(Recruits), color = disturbedarea_prct_cs), alpha = 0.5, size = 2) +
    geom_line(data = prediction_df, aes(x = spawners, y = log(recruits_predicted_low_cpd)), color = '#35978f', size = 1, alpha = 0.5) +
    geom_line(data = prediction_df, aes(x = spawners, y = log(recruits_predicted_high_cpd)), color = '#bf812d', size = 1, alpha = 0.5) +
    # geom_line(data = prediction_df, aes(x = spawners, y = log_RS), color = "black", size = 1, alpha = 0.5) +
    geom_ribbon(data = prediction_df, aes(x = spawners, ymin = log(recruits_predicted_low_cpd_lower), ymax = log(recruits_predicted_low_cpd_upper)), fill = "#35978f", alpha = 0.1) +
    geom_ribbon(data = prediction_df, aes(x = spawners, ymin = log(recruits_predicted_high_cpd_lower), ymax = log(recruits_predicted_high_cpd_upper)), fill = "#bf812d", alpha = 0.1) +
    labs(#title = "Ricker model with CPD, NPGO, ERSST", 
      x = "Spawners", y = TeX(r"($\log \left(Recruits\right)$)")) +
    scale_color_gradient2(name = 'CDA (%)',
                          low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 20)+
    theme_classic() +
    theme(legend.position = c(0.8,0.9),
          legend.background = element_rect(fill = alpha('white', 0.5)),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_blank(),
          legend.key.width = unit(0.5, "cm"),
          legend.key.height = unit(1, "lines"),
          legend.spacing.y = unit(0.001, "cm")
    )
  
  
  
  
  
  
  return((p1))
}




for(i in pink_watersheds){
  
  effects_plot <- plot_all_effects_river_together(
    posterior1 = ric_pk_cpd_ersst_long_chain,
    river_name = str_to_title(i),
    river = unique(case_study_watersheds_data_pk$River_n[case_study_watersheds_data_pk$River == i])
  )
  
  change_plot <- plot_recruitment_change_river_together(
    posterior1 = ric_pk_cpd_ersst_long_chain,
    posterior2 = ric_pk_eca_ersst_long_chain,
    river_name = i,
    effect1 = "cpd",
    effect2 = "eca",
    species = "pink",
    model1 = "CPD",
    model2 = "ECA",
    hd = FALSE
  )
  
  spawner_recruit_plot_even <- plot_recruit_spawner_river_pink(species = "pink",
                                                         river_name = paste0(i,"_Even"),
                                                         posterior = ric_pk_cpd_ersst_long_chain)
  spawner_recruit_plot_odd <- plot_recruit_spawner_river_pink(species = "pink",
                                                        river_name = paste0(i,"_Odd"),
                                                        posterior = ric_pk_cpd_ersst_long_chain)
  
  ggsave(filename = here("figures",
                         paste0("case_study_pink_et_w_sr_curves_",str_replace_all(str_to_lower(i), " ", "_"),".png")),
         plot = ((effects_plot+spawner_recruit_plot_even/spawner_recruit_plot_odd + 
                    plot_layout(widths = c(1,1.2)))/change_plot) + plot_layout(heights = c(1.5,1)) +
           plot_annotation(tag_levels = 'A',title = paste(str_to_title(i), "- Pink Salmon"))&
           theme(plot.tag.position = c(0.0, 1.0),
                 plot.tag = element_text(size = 10, hjust = 0, vjust = 0, face = "bold")
                 ),
         width = 7,
         height = 6,
         units = "in",
         dpi = 300)
  
  
}

#print effect sizes and CIs

all_river_df <- data.frame()
for(i in watersheds){
  
  river_data <- case_study_watersheds_data %>% filter(River == i)
  river <- river_data$River_n[1]
  
  posterior_cpd_df <- ric_chm_cpd_ocean_covariates_logR_long_chain %>%
    select(starts_with('b_for_rv'),starts_with('b_sst_rv'),starts_with('b_npgo_rv')) %>%
    pivot_longer(cols = everything(),
                 names_to = c('Effect','River'),
                 names_pattern = 'b_(.*)_rv(.*)',
                 values_to = "coefficient") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>% 
    filter(River_n == river)
  
  posterior_eca_df <- ric_chm_eca_ocean_covariates_logR_long_chain %>%
    select(starts_with('b_for_rv'),starts_with('b_sst_rv'),starts_with('b_npgo_rv')) %>%
    pivot_longer(cols = everything(),
                 names_to = c('Effect','River'),
                 names_pattern = 'b_(.*)_rv(.*)',
                 values_to = "coefficient") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>% 
    filter(River_n == river)
  
  effect_sizes_cpd <- posterior_cpd_df %>%
    group_by(Effect) %>%
    summarise(median = round(median(coefficient),2),
              lower_95 = round(quantile(coefficient, 0.025),2),
              upper_95 = round(quantile(coefficient, 0.975),2))
  
  effect_sizes_eca <- posterior_eca_df %>%
    group_by(Effect) %>%
    summarise(median = round(median(coefficient),2),
              lower_95 = round(quantile(coefficient, 0.025),2),
              upper_95 = round(quantile(coefficient, 0.975),2))
  
  # print(paste0("River: ", i))
  # print("CPD effect sizes:")
  # print(effect_sizes_cpd)
  # print("ECA effect sizes:")
  # print(effect_sizes_eca)
  
  #print predicted recruitment change (%) at most recent CPD levels and most recent ECA levels
  
  eca <- seq(0,1,length.out=100)
  
  eca_sqrt <- sqrt(eca)
  
  eca_sqrt_std <- (eca_sqrt-mean(eca_sqrt))/sd(eca_sqrt)
  
  eca_df <- data.frame(eca = eca,
                        eca_sqrt = eca_sqrt,
                        eca_sqrt_std = eca_sqrt_std)
  
  cpd <- seq(0,100,length.out=100)
  
  cpd_sqrt <- sqrt(cpd)
  
  cpd_sqrt_std <- (cpd_sqrt-mean(cpd_sqrt))/sd(cpd_sqrt)
  
  cpd_df <- data.frame(cpd = cpd,
                        cpd_sqrt = cpd_sqrt,
                        cpd_sqrt_std = cpd_sqrt_std)
  
  
  
  
  no_eca <- min(eca_sqrt_std)
  
  no_cpd <- min(cpd_sqrt_std)
  
  max_eca <- eca_df$eca_sqrt_std[which.min(abs(eca_df$eca - max(river_data$ECA_age_proxy_forested_only)))]
  
  max_cpd <- cpd_df$cpd_sqrt_std[which.min(abs(cpd_df$cpd - max(river_data$disturbedarea_prct_cs)))]
  
  
  
  # need to change b_rv to posterior 1 and posterior 2 and then make figure
  
  recruitment_cpd <- (exp(as.matrix(posterior_cpd_df %>% filter(Effect == "for") %>% select(coefficient))%*%
                             (max_cpd-no_cpd)))*100 - 100
  
  recruitment_eca <- (exp(as.matrix(posterior_eca_df %>% filter(Effect == "for") %>% select(coefficient))%*%
                             (max_eca-no_eca)))*100 - 100
  
  recruitment_df <- data.frame(
    Effect = c("Recruitment change at max CPD level", "Recruitment change at max ECA level"),
    model = c("CDA", "ECA"),
    median = c(round(median(recruitment_cpd),1), round(median(recruitment_eca),1)),
    lower_95 = c(round(quantile(recruitment_cpd, 0.025),1), round(quantile(recruitment_eca, 0.025),2)),
    upper_95 = c(round(quantile(recruitment_cpd, 0.975),1), round(quantile(recruitment_eca, 0.975),1))
  )
  
  #put the effect sizes and recuitment df together
  
  river_df <- effect_sizes_cpd %>% 
    mutate(model = "CDA") %>% 
    rbind(effect_sizes_eca %>% mutate(model = "ECA")) %>% 
    rbind(recruitment_df) %>% 
    mutate(River = i, Species = "chum")
  
  all_river_df <- rbind(all_river_df, river_df)
  
  
}

  
# pink watersheds

all_river_df_pk <- data.frame()

for(i in pink_watersheds){
  
  river_data <- case_study_watersheds_data_pk %>% filter(River == i)
  river <- river_data$River_n[1]
  
  posterior_cpd_df <- ric_pk_cpd_ersst_long_chain %>%
    select(starts_with('b_for_rv'),starts_with('b_sst_rv'),starts_with('b_npgo_rv')) %>%
    pivot_longer(cols = everything(),
                 names_to = c('Effect','River'),
                 names_pattern = 'b_(.*)_rv(.*)',
                 values_to = "coefficient") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>% 
    filter(River_n == river)
  
  posterior_eca_df <- ric_pk_eca_ersst_long_chain %>%
    select(starts_with('b_for_rv'),starts_with('b_sst_rv'),starts_with('b_npgo_rv')) %>%
    pivot_longer(cols = everything(),
                 names_to = c('Effect','River'),
                 names_pattern = 'b_(.*)_rv(.*)',
                 values_to = "coefficient") %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>% 
    filter(River_n == river)
  
  effect_sizes_cpd <- posterior_cpd_df %>%
    group_by(Effect) %>%
    summarise(median = round(median(coefficient),2),
              lower_95 = round(quantile(coefficient, 0.025),2),
              upper_95 = round(quantile(coefficient, 0.975),2))
  
  effect_sizes_eca <- posterior_eca_df %>%
    group_by(Effect) %>%
    summarise(median = round(median(coefficient),2),
              lower_95 = round(quantile(coefficient, 0.025),2),
              upper_95 = round(quantile(coefficient, 0.975),2))
  
  
  eca <- seq(0,1,length.out=100)
  
  eca_sqrt <- sqrt(eca)
  
  eca_sqrt_std <- (eca_sqrt-mean(eca_sqrt))/sd(eca_sqrt)
  
  eca_df <- data.frame(eca = eca,
                       eca_sqrt = eca_sqrt,
                       eca_sqrt_std = eca_sqrt_std)
  
  cpd <- seq(0,100,length.out=100)
  
  cpd_sqrt <- sqrt(cpd)
  
  cpd_sqrt_std <- (cpd_sqrt-mean(cpd_sqrt))/sd(cpd_sqrt)
  
  cpd_df <- data.frame(cpd = cpd,
                       cpd_sqrt = cpd_sqrt,
                       cpd_sqrt_std = cpd_sqrt_std)
  
  
  
  
  no_eca <- min(eca_sqrt_std)
  
  no_cpd <- min(cpd_sqrt_std)
  
  max_eca <- eca_df$eca_sqrt_std[which.min(abs(eca_df$eca - max(river_data$ECA_age_proxy_forested_only)))]
  
  max_cpd <- cpd_df$cpd_sqrt_std[which.min(abs(cpd_df$cpd - max(river_data$disturbedarea_prct_cs)))]
  
  
  
  # need to change b_rv to posterior 1 and posterior 2 and then make figure
  
  recruitment_cpd <- (exp(as.matrix(posterior_cpd_df %>% filter(Effect == "for") %>% select(coefficient))%*%
                            (max_cpd-no_cpd)))*100 - 100
  
  recruitment_eca <- (exp(as.matrix(posterior_eca_df %>% filter(Effect == "for") %>% select(coefficient))%*%
                            (max_eca-no_eca)))*100 - 100
  
  recruitment_df <- data.frame(
    Effect = c("Recruitment change at max CPD level", "Recruitment change at max ECA level"),
    model = c("CDA", "ECA"),
    median = c(round(median(recruitment_cpd),1), round(median(recruitment_eca),1)),
    lower_95 = c(round(quantile(recruitment_cpd, 0.025),1), round(quantile(recruitment_eca, 0.025),2)),
    upper_95 = c(round(quantile(recruitment_cpd, 0.975),1), round(quantile(recruitment_eca, 0.975),1))
  )
  
  
  
  river_df <- effect_sizes_cpd %>% 
    mutate(model = "CDA") %>% 
    rbind(effect_sizes_eca %>% mutate(model = "ECA"))  %>% 
    rbind(recruitment_df) %>% 
    mutate(River = i, Species = "pink")
  
  all_river_df_pk <- rbind(all_river_df_pk, river_df)
  
}

# format table

# make a wide table with rivers-species combination in the columns and "median (lower_95, upper_95)"values


all_river_df_chum_pink <- all_river_df %>% 
  rbind(all_river_df_pk) %>% 
  mutate(Estimate = paste(median, " (", lower_95, ", ", upper_95, ")", sep = "")) %>%
  mutate(Effect = case_when(Effect == "for" ~ "Forestry effect size",
                        Effect == "sst" ~ "SST effect size",
                        Effect == "npgo" ~ "NPGO effect size",
                        Effect == "Recruitment change at max CPD level" ~ "Recruitment change",
                        Effect == "Recruitment change at max ECA level" ~ "Recruitment change")) %>%
  rename(Model = model) %>%
  select(River, Species, Effect, Estimate, Model) %>%
  pivot_wider(#id = c(River, Species, model, Effect),
              names_from = c(River, Species),
              values_from = Estimate,
              names_vary = "slowest"
              ) %>% 
  arrange(Model)

# save table

write.csv(all_river_df_chum_pink, here("tables","case_study_effect_sizes_recruitment_chum_pink.csv"), row.names = FALSE)








































