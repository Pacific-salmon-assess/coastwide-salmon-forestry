# make figure of the results of the models without cm-18 and with cm-18

library(here);library(dplyr); library(stringr)
library(ggplot2)
library(tidyverse)
library(bayesplot)
library(patchwork)
library(hues)
library(GGally)
library(latex2exp)



ric_chm_eca_ocean_covariates_logR=read.csv(here('stan models','outs','posterior',
                                                'ric_chm_eca_ocean_covariates_logR.csv'),
                                           check.names=F)
ric_chm_cpd_ocean_covariates_logR=read.csv(here('stan models','outs','posterior',
                                                'ric_chm_cpd_ocean_covariates_logR.csv'),
                                           check.names=F)
ric_chm_eca_ocean_covariates_logR_wo_cm18=read.csv(here('stan models','outs','posterior',
                                                'ric_chm_eca_ocean_covariates_logR_wo_CM18.csv'),
                                           check.names=F)
ric_chm_cpd_ocean_covariates_logR_wo_cm18=read.csv(here('stan models','outs','posterior',
                                                'ric_chm_cpd_ocean_covariates_logR_wo_CM18.csv'),
                                           check.names=F)

bh_chm_eca_ocean_covariates=read.csv(here('stan models','outs','posterior',
                                                'bh_chm_eca_ocean_covariates.csv'),
                                           check.names=F)
bh_chm_cpd_ocean_covariates=read.csv(here('stan models','outs','posterior',
                                                'bh_chm_cpd_ocean_covariates.csv'))

bh_chm_cpd_ocean_covariates_wo_cm18=read.csv(here('stan models','outs','posterior',
                                                'bh_chm_cpd_ocean_covariates_wo_CM18.csv'))

bh_chm_eca_ocean_covariates_wo_cm18=read.csv(here('stan models','outs','posterior',
                                                'bh_chm_eca_ocean_covariates_wo_CM18.csv'))



x_name = "Standardized coefficients of CPD"
y_name = "Posterior density\nof CPD effect"
model = "Ricker models with CPD"
xlim = c(-1,1)

# combine two posteriors with label of the model

posterior1 <- data.frame(b_for = ric_chm_cpd_ocean_covariates_logR$b_for)
posterior1$model <- "with CM18"
posterior2 <- data.frame(b_for = ric_chm_cpd_ocean_covariates_logR_wo_cm18$b_for)
posterior2$model <- "without CM18"


posteriors <- rbind(posterior1, posterior2)



ggplot(posteriors) +
  geom_density(aes(b_for, color = model, group = model, fill = model), 
               linewidth = 0.5, alpha = 0.5)+
  
  geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
  labs(title = "Ricker models with cpd", x = x_name, y = y_name) +
  xlim(xlim[1], xlim[2]) +
  # scale_color +
  # geom_density(aes(posterior$b_for), color = 'black', linewidth = 1.2, alpha = 0.2)+
  #vline at the median value of the posterior
  theme_classic()+
  theme(legend.position = 'right',
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
  )

#save
ggsave(here("figures","chum_ricker_cpd_wo_cm18.png"), width = 8, height = 6, dpi = 300)




posterior1_bh <- data.frame(b_for = bh_chm_cpd_ocean_covariates$b_for)
posterior1_bh$model <- "with CM18"
posterior2_bh <- data.frame(b_for = bh_chm_cpd_ocean_covariates_wo_cm18$b_for)
posterior2_bh$model <- "without CM18"

posteriors_bh <- rbind(posterior1_bh, posterior2_bh)

ggplot(posteriors_bh) +
  geom_density(aes(b_for, color = model, group = model, fill = model), 
               linewidth = 0.5, alpha = 0.5)+
  
  geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
  labs(title = "Beverton-Holt models with cpd", x = x_name, y = y_name) +
  xlim(xlim[1], xlim[2]) +
  # scale_color +
  # geom_density(aes(posterior$b_for), color = 'black', linewidth = 1.2, alpha = 0.2)+
  #vline at the median value of the posterior
  theme_classic()+
  theme(legend.position = 'right',
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
  )

#save
ggsave(here("figures","chum_bh_cpd_wo_cm18.png"), width = 8, height = 6, dpi = 300)ggsave(here("figures","chum_bh_cpd_wo_cm18.png"), width = 8, height = 6, dpi = 300)


# Data wrangling

ch20rsc <- read.csv(here("origional-ecofish-data-models","Data","Processed",
                         "chum_SR_20_hat_yr_w_ocean_covariates.csv"))

#remove CM-18
ch20rsc_wo_cm18 <- ch20rsc %>%
  filter(!(CU == "CM-18"))

#two rivers with duplicated names:
ch20rsc$River=ifelse(ch20rsc$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',ch20rsc$River)
ch20rsc$River=ifelse(ch20rsc$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',ch20rsc$River)

ch20rsc_wo_cm18$River=ifelse(ch20rsc_wo_cm18$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',ch20rsc_wo_cm18$River)
ch20rsc_wo_cm18$River=ifelse(ch20rsc_wo_cm18$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000",'LAGOON CREEK 2',ch20rsc_wo_cm18$River)

ch20rsc=ch20rsc[order(factor(ch20rsc$River),ch20rsc$BroodYear),]

ch20rsc$River_n <- as.numeric(factor(ch20rsc$River))

ch20rsc_wo_cm18=ch20rsc_wo_cm18[order(factor(ch20rsc_wo_cm18$River),ch20rsc_wo_cm18$BroodYear),]
ch20rsc_wo_cm18$River_n <- as.numeric(factor(ch20rsc_wo_cm18$River))

#normalize ECA 2 - square root transformation (ie. sqrt(x))
ch20rsc$sqrt.ECA=sqrt(ch20rsc$ECA_age_proxy_forested_only)
ch20rsc$sqrt.ECA.std=(ch20rsc$sqrt.ECA-mean(ch20rsc$sqrt.ECA))/sd(ch20rsc$sqrt.ECA)

#normalize CPD 2 - square root transformation (ie. sqrt(x))
ch20rsc$sqrt.CPD=sqrt(ch20rsc$disturbedarea_prct_cs)
ch20rsc$sqrt.CPD.std=(ch20rsc$sqrt.CPD-mean(ch20rsc$sqrt.CPD))/sd(ch20rsc$sqrt.CPD)

ch20rsc$npgo.std=(ch20rsc$npgo-mean(ch20rsc$npgo))/sd(ch20rsc$npgo)
ch20rsc$sst.std=(ch20rsc$spring_ersst-mean(ch20rsc$spring_ersst))/sd(ch20rsc$spring_ersst)

ch20rsc_wo_cm18$sqrt.ECA=sqrt(ch20rsc_wo_cm18$ECA_age_proxy_forested_only)
ch20rsc_wo_cm18$sqrt.ECA.std=(ch20rsc_wo_cm18$sqrt.ECA-mean(ch20rsc_wo_cm18$sqrt.ECA))/sd(ch20rsc_wo_cm18$sqrt.ECA)

ch20rsc_wo_cm18$sqrt.CPD=sqrt(ch20rsc_wo_cm18$disturbedarea_prct_cs)
ch20rsc_wo_cm18$sqrt.CPD.std=(ch20rsc_wo_cm18$sqrt.CPD-mean(ch20rsc_wo_cm18$sqrt.CPD))/sd(ch20rsc_wo_cm18$sqrt.CPD)


cu = distinct(ch20rsc,.keep_all = T)
cu_wo_cm18 = distinct(ch20rsc_wo_cm18,.keep_all = T)

cu_n = as.numeric(factor(cu$CU))
cu_n_wo_cm18 = as.numeric(factor(cu_wo_cm18$CU))

ch20rsc$CU_n <- cu_n
ch20rsc_wo_cm18$CU_n <- cu_n_wo_cm18

# make CU_NAME values Title case instead of all caps
ch20rsc$CU_name <- str_to_title(ch20rsc$CU_NAME)

ch20rsc_wo_cm18$CU_name <- str_to_title(ch20rsc_wo_cm18$CU_NAME)

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

max_eca_df_wo_cm18 <- ch20rsc_wo_cm18 %>% 
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


max_cpd_df_wo_cm18 <- ch20rsc_wo_cm18 %>% 
  select(River, River_n,CU, CU_n, CU_name, disturbedarea_prct_cs) %>%
  group_by(River) %>% 
  summarize(River_n = first(River_n),
            River = first(River),
            CU = first(CU),
            CU_n = first(CU_n),
            CU_name = first(CU_name),
            cpd_max = max(disturbedarea_prct_cs, na.rm = TRUE))


plot_forestry_effect <- function(posterior = bh_chm_eca, 
                                 effect = "eca", 
                                 species = "chum", 
                                 model = "BH model with NPGO", 
                                 wo_cm_18 = F,
                                 xlim = c(-2, 2), by_cu = FALSE, no_color = T,
                                 color_by_cu = FALSE, global_only = FALSE){
  
  if(effect == "eca"){
    x_name = "Standardized coefficients of ECA"
    y_name = "Posterior density\nof ECA effect"
    if(!no_color){
      color_var = "eca_level"
      color_name = "Max ECA level\nin River (%)"
      scale_color = scale_color_manual(name = "Max ECA level\nin River (%)",
                                       values = c('low' = '#35978f', 'medium' = 'gray', 'high' = '#bf812d'))
    }
    
    if(species == "chum"){
      max_df = max_eca_df
      if(wo_cm_18){
        max_df = max_eca_df_wo_cm18
      }
    } else if(species == "pink"){
      max_df = max_eca_pink_df
      
    }
    
  } else if(effect == "cpd"){
    x_name = "Standardized coefficients of CPD"
    y_name = "Posterior density\nof CPD effect"
    if(!no_color){
      color_var = "cpd_max"
      color_name = "Max CPD\nin River (%)"
      scale_color = scale_color_gradient2(name = 'Max CPD\nin River (%)',
                                          low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 50)
      
    }
    if(species == "chum"){
      max_df = max_cpd_df
      if(wo_cm_18){
        max_df = max_cpd_df_wo_cm18
      }
    } else if(species == "pink"){
      max_df = max_cpd_pink_df
    }
  }
  
  if(color_by_cu){
    color_var = "CU_name"
    color_name = "CU"
    scale_color = scale_color_iwanthue(name = 'CU', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
                                       lmin = 10, lmax = 95)
  }
       
       
       
  if(by_cu){
    if(wo_cm_18){
      posterior_df <- posterior %>%
        select(starts_with('b_for_cu')) %>%
        pivot_longer(cols = everything(), 
                     names_to = 'CU', 
                     names_prefix = 'b_for_cu',
                     values_to = effect) %>%
        mutate(CU_n = as.numeric(str_extract(CU, '\\d+'))) %>% 
        select(-CU) %>%
        left_join(ch20rsc_wo_cm18 %>% select(CU, CU_n, CU_name) %>% distinct(), by = 'CU_n') 
      
    } else{
      posterior_df <- posterior %>%
        select(starts_with('b_for_cu')) %>%
        pivot_longer(cols = everything(), 
                     names_to = 'CU', 
                     names_prefix = 'b_for_cu',
                     values_to = effect) %>%
        mutate(CU_n = as.numeric(str_extract(CU, '\\d+'))) %>% 
        select(-CU) %>%
        left_join(ch20rsc %>% select(CU, CU_n, CU_name) %>% distinct(), by = 'CU_n') 
    }
    
    scale_color = scale_color_iwanthue(name = 'CU', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
                                       lmin = 10, lmax = 95)
    color_var = "CU_name"
    
    #which CU has minimum and maximum median effect
    min_cu <- posterior_df %>% 
      group_by(CU_n, CU_name) %>% 
      summarize(median_effect = median(!!sym(effect)),
                CU_name = first(CU_name)) %>% 
      ungroup() %>%
      filter(median_effect == min(median_effect))
    
    max_cu <- posterior_df %>%
      group_by(CU_n, CU_name) %>% 
      summarize(median_effect = median(!!sym(effect)),
                CU_name = first(CU_name)) %>% 
      ungroup() %>%
      filter(median_effect == max(median_effect))
    
    #color by CU and label the CU which has minimum and maximum median effect
  plot1 <- ggplot() +
    stat_density(data= posterior_df, aes(!!sym(effect), color = !!sym(color_var), group = CU_n),
                 geom = 'line', position = 'identity', 
                 alpha = 0.5, linewidth = 1.5) +
    geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
    labs(title = model, x = x_name, y = y_name) +
    xlim(xlim[1], xlim[2]) +
    scale_color +
    geom_density(aes(posterior$b_for), color = 'black', linewidth = 1.2, alpha = 0.2)+
    geom_text(data = min_cu, aes(x = median_effect, y = 1, label = CU_name, color = !!sym(color_var)), size = 2,
              vjust = -0.5, hjust = 0.5)+
    geom_text(data = max_cu, aes(x = median_effect, y = 1, label = CU_name, color = !!sym(color_var)), size = 2,
              vjust = -0.5, hjust = 0.5)+
    #vline at the median value of the posterior
    geom_vline(xintercept = median(posterior$b_for), color = 'black', linetype = 'dashed')+
    theme_classic()+
    theme(legend.position = "right",
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.8, "lines"),
        #1 column for legend
        legend.box = "vertical",
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.001, "cm"),
        axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 14, hjust = 0.5)
        )+
    guides(color = guide_legend(override.aes = list(alpha = 1), ncol =1))
  
  } else {
    
    
    posterior_df <- posterior %>%
    select(starts_with('b_for_rv')) %>%
    pivot_longer(cols = everything(), 
                 names_to = 'River', 
                 names_prefix = 'b_for_rv',
                 values_to = effect) %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>%
    left_join(max_df, by = 'River_n')
    
    #color by river
  plot1 <- ggplot() +
    stat_density(data= posterior_df, aes(!!sym(effect), color = !!sym(color_var), group = River),
                 geom = 'line', position = 'identity', 
                 alpha = 0.1, linewidth = 1) + 
    geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
    labs(title = model, x = x_name, y = y_name) +
    xlim(xlim[1], xlim[2]) +
    scale_color +
    geom_density(aes(posterior$b_for), color = 'black', linewidth = 1.2, alpha = 0.2)+
    #vline at the median value of the posterior
    geom_vline(xintercept = median(posterior$b_for), color = 'black', linetype = 'dashed')+
    theme_classic()+
    theme(legend.position = 'right',
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.text = element_text(size = 7),
          legend.spacing.y = unit(0.001, "cm"),
          legend.key.width = unit(0.5, "cm"),
          legend.key.height = unit(0.8, "lines"),
          plot.title = element_text(size = 14, hjust = 0.5)
          )+
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.6), ncol = 1))
    
    if(global_only){
      
      
      plot1 <- ggplot() +
        geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
        labs(title = model, x = x_name, y = y_name) +
        xlim(xlim[1], xlim[2]) +
        scale_color +
        geom_density(aes(posterior$b_for), color = 'black', 
                     linewidth = 1.2, alpha = 0.1)+
        geom_text(x = median(posterior$b_for), y = 1, 
                  label = paste("Median:", round(median(posterior$b_for), 2)),
                  color = 'black', size = 4, hjust = 0.5)+
        #vline at the median value of the posterior
        geom_vline(xintercept = median(posterior$b_for), color = 'black', linetype = 'dashed')+
        theme_classic()+
        theme(legend.position = c(0.8,0.5),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.x = element_text(size = 10),
              axis.text.y = element_text(size = 10),
              legend.text = element_text(size = 7),
              legend.spacing.y = unit(0.001, "cm"),
              legend.key.width = unit(0.5, "cm"),
              legend.key.height = unit(0.8, "lines"),
              plot.title = element_text(size = 14, hjust = 0.5)
        )+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.6), ncol = 1))
      
    } else{
      
      
      posterior_df <- posterior %>%
        select(starts_with('b_for_rv')) %>%
        pivot_longer(cols = everything(), 
                     names_to = 'River', 
                     names_prefix = 'b_for_rv',
                     values_to = effect) %>%
        mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
        select(-River) %>%
        left_join(max_df, by = 'River_n')
      
      #color by river
      plot1 <- ggplot() +
        stat_density(data= posterior_df, aes(!!sym(effect), color = !!sym(color_var), group = River),
                     geom = 'line', position = 'identity', 
                     alpha = 0.1, linewidth = 1) +
        geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
        labs(title = model, x = x_name, y = y_name) +
        xlim(xlim[1], xlim[2]) +
        scale_color +
        geom_density(aes(posterior$b_for), color = 'black', linewidth = 1.2, alpha = 0.2)+
        #vline at the median value of the posterior
        geom_vline(xintercept = median(posterior$b_for), color = 'black', linetype = 'dashed')+
        theme_classic()+
        theme(legend.position = c(0.8,0.5),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.x = element_text(size = 10),
              axis.text.y = element_text(size = 10),
              legend.text = element_text(size = 7),
              legend.spacing.y = unit(0.001, "cm"),
              legend.key.width = unit(0.5, "cm"),
              legend.key.height = unit(0.8, "lines"),
              plot.title = element_text(size = 14, hjust = 0.5)
        )+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.6), ncol = 1))
      
    }
  
  }
  
  
  
  
  
  return(plot1)
}


(plot_forestry_effect(posterior = ric_chm_eca_ocean_covariates_logR, 
                      effect = "eca", species = "chum", 
                      model = "Ricker model with ECA", color_by_cu = F, no_color = F,
                      by_cu = F, xlim = c(-0.5,0.5), global_only = FALSE)+ 
    plot_forestry_effect(posterior = ric_chm_eca_ocean_covariates_logR_wo_cm18, 
                         effect = "eca", species = "chum", color_by_cu = F, no_color = F,
                         wo_cm_18 = T,
                         model = "Ricker model with ECA, without CM-18", 
                         xlim = c(-0.5,0.5), global_only = F)) + 
  plot_layout(guides = "collect")

(plot_forestry_effect(posterior = bh_chm_cpd_ocean_covariates, 
                      effect = "cpd", species = "chum", 
                      model = "BH model with CPD", color_by_cu = F, no_color = F,
                      by_cu = T,
                      xlim = c(-2,2), global_only = FALSE)+ 
    plot_forestry_effect(posterior = bh_chm_cpd_ocean_covariates_wo_cm18, 
                         effect = "cpd", species = "chum", color_by_cu = F, no_color = F,
                         wo_cm_18 = T,by_cu = T,
                         model = "BH model with CPD, without CM-18", 
                         xlim = c(-2,2), global_only = F)) + 
  plot_layout(guides = "collect")

#save
ggsave(here("figures","chum_bh_cu_cpd_wo_cm18.png"), width = 12, height = 6, dpi = 300)


#save
ggsave(here("figures","chum_bh_river_cpd_wo_cm18.png"), width = 10, height = 6, dpi = 300)
