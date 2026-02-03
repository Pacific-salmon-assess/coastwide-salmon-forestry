

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
library(posterior)

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






plot_productivity_decline_manuscript1 <- function(posterior = ric_chm_cpd_ocean_covariates_logR_long_chain, 
                                                  effect = "cpd", species = "chum", 
                                                  model = "Ricker", 
                                                  by_river = FALSE){
  
  
  if(species == "chum"){
    df <- ch20rsc 
    df$sst.std <- (ch20rsc$spring_ersst-mean(ch20rsc$spring_ersst))/sd(ch20rsc$spring_ersst)
    
  } else if(species == "pink"){
    df <- pk10r
    df$sst.std <- (pk10r$spring_ersst-mean(pk10r$spring_ersst))/sd(pk10r$spring_ersst)
  }
  
  if(effect == "eca"){
    forestry <- seq(0,1,length.out=100)
    b <- posterior %>% select(ends_with("b_for"))
    
    forestry_sqrt <- sqrt(forestry)
    
    forestry_sqrt_std = (forestry_sqrt-mean(forestry_sqrt))/sd(forestry_sqrt)
    
  }else if(effect == "cpd"){
    forestry <- seq(0,100,length.out=100)
    b <- posterior %>% select(ends_with("b_for"))
    forestry_sqrt <- sqrt(forestry)
    
    forestry_sqrt_std = (forestry_sqrt-mean(forestry_sqrt))/sd(forestry_sqrt)
    
  }else if(effect == "sst"){
    forestry <- seq(10,15,length.out=100)
    b <- posterior %>% select(ends_with("b_sst"))
    
    forestry_sqrt_std = (forestry-mean(forestry))/sd(forestry)
    
  }else if(effect == "npgo"){
    forestry <- seq(-3,3,length.out=100)
    b <- posterior %>% select(ends_with("b_npgo"))
    
    forestry_sqrt_std = (forestry-mean(forestry))/sd(forestry)
    
  }
  
  
  
  no_forestry <- min(forestry_sqrt_std)
  
  # ch_chm_eca_sst=read.csv(here('stan models','outs','posterior',ch_chm_eca_sst),check.names=F)
  
  
  
  global_prediction <- apply(exp(as.matrix(b[,1])%*%
                                   (forestry_sqrt_std-no_forestry))*100 - 100,2,median)
  
  global_025 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.025), 
                      row.names = c("q025"))
  
  global_975 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.975),
                      row.names = c("q975"))
  global_750 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.75),
                      row.names = c("q750"))
  
  global_250 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.25),
                      row.names = c("q250"))
  
  global_900 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.9),
                      row.names = c("q900"))
  
  global_100 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.1),
                      row.names = c("q100"))
  
  
  
  
  
  
  global_df <- data.frame(forestry = forestry,
                          productivity_median = global_prediction,
                          q025 = global_025,
                          q975 = global_975,
                          q750 = global_750,
                          q250 = global_250,
                          q900 = global_900,
                          q100 = global_100)
  
  full_productivity <- NULL
  
  if(by_river){
    
    # no_forestry <- min(eca_std_sqrt)
    for (i in 1:length(unique(df$River_n))){
      river <- unique(df$River_n)[i]
      
      river_data <- df %>% filter(River_n == river)
      
      if(effect == "eca"){
        b_rv <- posterior %>% select(starts_with("b_for_rv")) %>%
          select(ends_with(paste0("[",river,"]")))
        
        forestry_sqrt_std_river <- seq(min(forestry_sqrt_std),max(river_data$sqrt.ECA.std),length.out=100)
        #minimum forestry possible -0
        forestry_river <- seq(0,max(river_data$ECA_age_proxy_forested_only),length.out=100)
        
        
      } else if(effect == "cpd"){
        b_rv <- posterior %>% select(starts_with("b_for_rv")) %>%
          select(ends_with(paste0("[",river,"]")))
        
        forestry_sqrt_std_river <- seq(min(forestry_sqrt_std),max(river_data$sqrt.CPD.std),length.out=100)
        #minimum forestry possible - 0
        forestry_river <- seq(0,max(river_data$disturbedarea_prct_cs),length.out=100)
        
        
      } else if(effect == "sst"){
        b_rv <- posterior %>% select(starts_with("b_sst_rv")) %>%
          select(ends_with(paste0("[",river,"]")))
        
        # use only spring_ersst
        
        forestry_sqrt_std_river <- seq(min(river_data$sst.std),max(river_data$sst.std),length.out=100)
        
        forestry_river <- seq(min(river_data$spring_ersst),max(river_data$spring_ersst),length.out=100)
        
      }
      
      
      
      
      
      
      
      productivity <- (exp(as.matrix(b_rv[,1])%*%
                             (forestry_sqrt_std-no_forestry)))*100 - 100
      
      productivity_median <- apply(productivity,2,median)
      
      productivity_median_df <- data.frame(River = unique(river_data$River),
                                           productivity_median = productivity_median,
                                           forestry = forestry_river) %>% 
        filter(forestry >= min(forestry_river), forestry <= max(forestry_river))
      
      full_productivity <- rbind(full_productivity, productivity_median_df)
      
    }
    
    # c <- posterior %>% select(ends_with("b_for"))
    
    median_prediction <- apply(exp(as.matrix(b[,1])%*%
                                     (forestry_sqrt_std-no_forestry))*100 - 100,2,median)
    
    median_df <- data.frame(forestry = forestry,
                            productivity_median = median_prediction)
    
    if(effect == "eca" || effect == "cpd"){
      p1 <- ggplot(full_productivity) +
        geom_line(aes(x = forestry, y = productivity_median, group = River,
                      color = "watershed level\nforestry effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = productivity_median, color = "global forestry effect"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = q025, ymax = q975),
                    alpha = 0.5, fill = "#ADcCA5") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("watershed level\nforestry effect" = "darkgray", 
                                         "global forestry effect" = "black")) +
        ylim(-100,100) +
        scale_x_continuous(n.breaks = 5) +
        labs(title = model,
             x = ifelse(effect == "eca", "Equivalent clearcut area", "Cumulative disturbance (%)"),
             y = "Median change\n in productivity (%)") +
        theme_classic() +
        theme(legend.position = c(0.8,0.8),
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    } else if(effect == "sst"){
      p1 <- ggplot(full_productivity) +
        geom_line(aes(x = forestry, y = productivity_median, group = River,
                      color = "watershed level\nSST effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = productivity_median, color = "global SST effect"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = q025, ymax = q975),
                    alpha = 0.5, fill = "gray") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("watershed level\nSST effect" = "darkgray", 
                                         "global SST effect" = "black")) +
        ylim(-100,100) +
        scale_x_continuous(n.breaks = 5, limits = c(10,15)) +
        labs(title = model,
             x = "Spring SST (°C)",
             y = "Median change\n in productivity (%)") +
        theme_classic() +
        theme(legend.position = c(0.8,0.8),
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    }
  } else {
    
    if(effect == "eca" || effect == "cpd"){
      p1 <- ggplot(full_productivity) +
        # geom_line(aes(x = forestry, y = productivity_median, group = River,
        #               color = "watershed level\nforestry effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = productivity_median, color = "Median coastwide\nchange"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = q025, ymax = q975,
                                          alpha = "95% credible\ninterval",
                                          fill = "95% credible\ninterval"))+
        geom_ribbon(data = global_df, aes(x = forestry, ymin = q250, ymax = q750,
                                          alpha = "50% credible\ninterval",
                                          fill = "50% credible\ninterval"))+
        geom_ribbon(data = global_df, aes(x = forestry, ymin = q100, ymax = q900,
                                          alpha = "80% credible\ninterval",
                                          fill = "80% credible\ninterval"))+
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("Median coastwide\nchange" = "#6F7c67")) +
        scale_fill_manual("",values  = c("95% credible\ninterval" =  "#ADcCA5",
                                         "50% credible\ninterval" = "#ADcCA5",
                                         "80% credible\ninterval" = "#ADcCA5"
        )) +
        scale_alpha_manual("",values  = c("95% credible\ninterval" = 0.25,
                                          "50% credible\ninterval" = 0.65,
                                          "80% credible\ninterval" = 0.45
        )) +
        ylim(-100,50) +
        scale_x_continuous(n.breaks = 5) +
        labs(title = model,
             x = ifelse(effect == "eca", "Equivalent clearcut area", "Cumulative disturbance (%)"),
             y = "Median change\n in productivity (%)") +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_blank(),
              legend.key.size = unit(0.2, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)),
               
        )
    } else if(effect == "sst"){
      p1 <- ggplot(full_productivity) +
        # geom_line(aes(x = forestry, y = productivity_median, group = River,
        #               color = "watershed level\nSST effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = productivity_median, color = "global SST effect"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = q025, ymax = q975),
                    alpha = 0.25, fill = "#C78c63") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = q250, ymax = q750),
                    alpha = 0.65, fill = "#C78c63") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = q100, ymax = q900),
                    alpha = 0.45, fill = "#C78c63") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("global SST effect" = "#C78c63")) +
        
        ylim(-50,100) +
        scale_x_continuous(n.breaks = 5, limits = c(10,15)) +
        labs(title = model,
             x = "Spring SST (°C)",
             y = "Median change\n in productivity (%)") +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    }else if(effect == "npgo"){
      p1 <- ggplot(full_productivity) +
        # geom_line(aes(x = forestry, y = productivity_median, group = River,
        #               color = "watershed level\nSST effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = productivity_median, color = "global NPGO effect"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = q025, ymax = q975),
                    alpha = 0.25, fill = "#829Dc6") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = q250, ymax = q750),
                    alpha = 0.65, fill = "#829Dc6") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = q100, ymax = q900),
                    alpha = 0.45, fill = "#829Dc6") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("global NPGO effect" = "#829Dc6")) +
        
        ylim(-50,100) +
        scale_x_continuous(n.breaks = 5, limits = c(-3,3)) +
        labs(title = model,
             x = "NPGO",
             y = "Median change\n in productivity (%)") +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    }
  }
  
  return(p1)
  
}




plot_productivity_decline_manuscript2 <- function(posterior = ric_chm_cpd_ocean_covariates_logR_long_chain, 
                                                 effect = "cpd", species = "chum", 
                                                 model = "Ricker", 
                                                 by_river = FALSE, hd = TRUE){
  
  
  if(species == "chum"){
    df <- ch20rsc 
    df$sst.std <- (ch20rsc$spring_ersst-mean(ch20rsc$spring_ersst))/sd(ch20rsc$spring_ersst)
    
  } else if(species == "pink"){
    df <- pk10r
    df$sst.std <- (pk10r$spring_ersst-mean(pk10r$spring_ersst))/sd(pk10r$spring_ersst)
  }
  
  if(effect == "eca"){
    forestry <- seq(0,1,length.out=100)
    b <- posterior %>% select(ends_with("b_for"))
    
    forestry_sqrt <- sqrt(forestry)
    
    forestry_sqrt_std = (forestry_sqrt-mean(forestry_sqrt))/sd(forestry_sqrt)
    
  }else if(effect == "cpd"){
    forestry <- seq(0,100,length.out=100)
    b <- posterior %>% select(ends_with("b_for"))
    forestry_sqrt <- sqrt(forestry)
    
    forestry_sqrt_std = (forestry_sqrt-mean(forestry_sqrt))/sd(forestry_sqrt)
    
  }else if(effect == "sst"){
    forestry <- seq(10,15,length.out=100)
    b <- posterior %>% select(ends_with("b_sst"))
    
    forestry_sqrt_std = (forestry-mean(forestry))/sd(forestry)
    
  }else if(effect == "npgo"){
    forestry <- seq(-3,3,length.out=100)
    b <- posterior %>% select(ends_with("b_npgo"))
    
    forestry_sqrt_std = (forestry-mean(forestry))/sd(forestry)
    
  }
  
  
  
  no_forestry <- min(forestry_sqrt_std)
  
  # ch_chm_eca_sst=read.csv(here('stan models','outs','posterior',ch_chm_eca_sst),check.names=F)
  
  
  
  global_prediction <- apply(exp(as.matrix(b[,1])%*%
                                   (forestry_sqrt_std-no_forestry))*100 - 100,2,median)
  
  global_025 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.025), 
                      row.names = c("q025"))
  
  global_975 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.975),
                      row.names = c("q975"))
  global_750 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.75),
                      row.names = c("q750"))
  
  global_250 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.25),
                      row.names = c("q250"))
  
  global_900 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.9),
                      row.names = c("q900"))
  
  global_100 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.1),
                      row.names = c("q100"))
  
  # calculate high density credible intervals
  
  hd_df <- apply(exp(as.matrix(b[,1])%*%
                       (forestry_sqrt_std-no_forestry))*100 - 100,2,bayestestR::hdi, 
                 ci = c(0.5,0.8,0.95))
  
  # extract hd intervals from hd_df
  
  for(i in 1:100){
    hd_50_lower <- hd_df[[i]]$CI_low[1]
    hd_50_upper <- hd_df[[i]]$CI_high[1]
    
    hd_80_lower <- hd_df[[i]]$CI_low[2]
    hd_80_upper <- hd_df[[i]]$CI_high[2]
    
    hd_95_lower <- hd_df[[i]]$CI_low[3]
    hd_95_upper <- hd_df[[i]]$CI_high[3]
    
    if(i == 1){
      hd_intervals <- data.frame(hd_50_lower = hd_50_lower,
                                 hd_50_upper = hd_50_upper,
                                 hd_80_lower = hd_80_lower,
                                 hd_80_upper = hd_80_upper,
                                 hd_95_lower = hd_95_lower,
                                 hd_95_upper = hd_95_upper)
    } else {
      hd_intervals <- rbind(hd_intervals, data.frame(hd_50_lower = hd_50_lower,
                                                     hd_50_upper = hd_50_upper,
                                                     hd_80_lower = hd_80_lower,
                                                     hd_80_upper = hd_80_upper,
                                                     hd_95_lower = hd_95_lower,
                                                     hd_95_upper = hd_95_upper))
    }
  }
  
  
  
  
  
  
  global_df <- data.frame(forestry = forestry,
                          productivity_median = global_prediction,
                          q025 = global_025,
                          q975 = global_975,
                          q750 = global_750,
                          q250 = global_250,
                          q900 = global_900,
                          q100 = global_100) %>% 
    cbind(hd_intervals)
  
  full_productivity <- NULL
  
  if(hd == TRUE){
    ymin_95 = "hd_95_lower"
    ymax_95 = "hd_95_upper"
    ymin_80 = "hd_80_lower"
    ymax_80 = "hd_80_upper"
    ymin_50 = "hd_50_lower"
    ymax_50 = "hd_50_upper"
    
  } else{
    ymin_95 = "q025"
    ymax_95 = "q975"
    ymin_80 = "q100"
    ymax_80 = "q900"
    ymin_50 = "q250"
    ymax_50 = "q750"
    
  }
  
  if(by_river){
    
    # no_forestry <- min(eca_std_sqrt)
    for (i in 1:length(unique(df$River_n))){
      river <- unique(df$River_n)[i]
      
      river_data <- df %>% filter(River_n == river)
      
      if(effect == "eca"){
        b_rv <- posterior %>% select(starts_with("b_for_rv")) %>%
          select(ends_with(paste0("[",river,"]")))
        
        forestry_sqrt_std_river <- seq(min(forestry_sqrt_std),max(river_data$sqrt.ECA.std),length.out=100)
        #minimum forestry possible -0
        forestry_river <- seq(0,max(river_data$ECA_age_proxy_forested_only),length.out=100)
        
        
      } else if(effect == "cpd"){
        b_rv <- posterior %>% select(starts_with("b_for_rv")) %>%
          select(ends_with(paste0("[",river,"]")))
        
        forestry_sqrt_std_river <- seq(min(forestry_sqrt_std),max(river_data$sqrt.CPD.std),length.out=100)
        #minimum forestry possible - 0
        forestry_river <- seq(0,max(river_data$disturbedarea_prct_cs),length.out=100)
        
        
      } else if(effect == "sst"){
        b_rv <- posterior %>% select(starts_with("b_sst_rv")) %>%
          select(ends_with(paste0("[",river,"]")))
        
        # use only spring_ersst
        
        forestry_sqrt_std_river <- seq(min(river_data$sst.std),max(river_data$sst.std),length.out=100)
        
        forestry_river <- seq(min(river_data$spring_ersst),max(river_data$spring_ersst),length.out=100)
        
      }
      
      
      
      
      
      
      
      productivity <- (exp(as.matrix(b_rv[,1])%*%
                             (forestry_sqrt_std-no_forestry)))*100 - 100
      
      productivity_median <- apply(productivity,2,median)
      
      productivity_median_df <- data.frame(River = unique(river_data$River),
                                           productivity_median = productivity_median,
                                           forestry = forestry_river) %>% 
        filter(forestry >= min(forestry_river), forestry <= max(forestry_river))
      
      full_productivity <- rbind(full_productivity, productivity_median_df)
      
    }
    
    # c <- posterior %>% select(ends_with("b_for"))
    
    median_prediction <- apply(exp(as.matrix(b[,1])%*%
                                     (forestry_sqrt_std-no_forestry))*100 - 100,2,median)
    
    median_df <- data.frame(forestry = forestry,
                            productivity_median = median_prediction)
    
    if(effect == "eca" || effect == "cpd"){
      p1 <- ggplot(full_productivity) +
        geom_line(aes(x = forestry, y = productivity_median, group = River,
                      color = "watershed level\nforestry effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = productivity_median, color = "global forestry effect"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = q025, ymax = q975),
                    alpha = 0.5, fill = "#ADcCA5") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("watershed level\nforestry effect" = "darkgray", 
                                         "global forestry effect" = "black")) +
        ylim(-100,100) +
        scale_x_continuous(n.breaks = 5) +
        labs(title = model,
             x = ifelse(effect == "eca", "Equivalent clearcut area", "Cumulative disturbance (%)"),
             y = "Median change\n in recruitment (%)") +
        theme_classic() +
        theme(legend.position = c(0.8,0.8),
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    } else if(effect == "sst"){
      p1 <- ggplot(full_productivity) +
        geom_line(aes(x = forestry, y = productivity_median, group = River,
                      color = "watershed level\nSST effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = productivity_median, color = "global SST effect"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = q025, ymax = q975),
                    alpha = 0.5, fill = "gray") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("watershed level\nSST effect" = "darkgray", 
                                         "global SST effect" = "black")) +
        ylim(-100,100) +
        scale_x_continuous(n.breaks = 5, limits = c(10,15)) +
        labs(title = model,
             x = "Spring SST (°C)",
             y = "Median change\n in recruitment (%)") +
        theme_classic() +
        theme(legend.position = c(0.8,0.8),
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    }
  } else {
    
    if(effect == "eca" || effect == "cpd"){
      p1 <- ggplot(full_productivity) +
        # geom_line(aes(x = forestry, y = productivity_median, group = River,
        #               color = "watershed level\nforestry effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = productivity_median, color = "Median coastwide\nchange"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_95), ymax = !!sym(ymax_95),
                                          alpha = "95% credible\ninterval",
                                          fill = "95% credible\ninterval"))+
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_50), ymax = !!sym(ymax_50),
                                          alpha = "50% credible\ninterval",
                                          fill = "50% credible\ninterval"))+
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_80), ymax = !!sym(ymax_80),
                                          alpha = "80% credible\ninterval",
                                          fill = "80% credible\ninterval"))+
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("Median coastwide\nchange" = "#6F7c67")) +
        scale_fill_manual("",values  = c("95% credible\ninterval" =  "#ADcCA5",
                                         "50% credible\ninterval" = "#ADcCA5",
                                         "80% credible\ninterval" = "#ADcCA5"
        )) +
        scale_alpha_manual("",values  = c("95% credible\ninterval" = 0.25,
                                          "50% credible\ninterval" = 0.65,
                                          "80% credible\ninterval" = 0.45
        )) +
        ylim(-100,50) +
        scale_x_continuous(n.breaks = 5) +
        labs(title = model,
             x = ifelse(effect == "eca", "Equivalent clearcut area", "Cumulative disturbance (%)"),
             y = "Median change\n in recruitment (%)") +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_blank(),
              legend.key.size = unit(0.2, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)),
               
        )
    } else if(effect == "sst"){
      p1 <- ggplot(full_productivity) +
        # geom_line(aes(x = forestry, y = productivity_median, group = River,
        #               color = "watershed level\nSST effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = productivity_median, color = "global SST effect"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_95), ymax = !!sym(ymax_95)),
                    alpha = 0.25, fill = "#C78c63") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_50), ymax = !!sym(ymax_50)),
                    alpha = 0.65, fill = "#C78c63") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_80), ymax = !!sym(ymax_80)),
                    alpha = 0.45, fill = "#C78c63") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("global SST effect" = "#C78c63")) +
        
        ylim(-50,100) +
        scale_x_continuous(n.breaks = 5, limits = c(10,15)) +
        labs(title = model,
             x = "Spring SST (°C)",
             y = "Median change\n in recruitment (%)") +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    }else if(effect == "npgo"){
      p1 <- ggplot(full_productivity) +
        # geom_line(aes(x = forestry, y = productivity_median, group = River,
        #               color = "watershed level\nSST effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = productivity_median, color = "global NPGO effect"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin =!!sym(ymin_95), ymax = !!sym(ymax_95)),
                    alpha = 0.25, fill = "#829Dc6") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_50), ymax = !!sym(ymax_50)),
                    alpha = 0.65, fill = "#829Dc6") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_80), ymax = !!sym(ymax_80)),
                    alpha = 0.45, fill = "#829Dc6") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("global NPGO effect" = "#829Dc6")) +
        
        ylim(-50,100) +
        scale_x_continuous(n.breaks = 5, limits = c(-3,3)) +
        labs(title = model,
             x = "NPGO",
             y = "Median change\n in recruitment (%)") +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    }
    
    
    
    
  }
  
  
  
  return(p1)
  
}


plot_productivity_decline_manuscript3 <- function(posterior = ric_chm_cpd_ocean_covariates_logR_long_chain, 
                                                  effect = "cpd", species = "chum", 
                                                  model = "Ricker", 
                                                  by_river = FALSE, hd = TRUE){
  
  
  if(species == "chum"){
    df <- ch20rsc 
    df$sst.std <- (ch20rsc$spring_ersst-mean(ch20rsc$spring_ersst))/sd(ch20rsc$spring_ersst)
    
  } else if(species == "pink"){
    df <- pk10r
    df$sst.std <- (pk10r$spring_ersst-mean(pk10r$spring_ersst))/sd(pk10r$spring_ersst)
  }
  
  if(effect == "eca"){
    forestry <- seq(0,1,length.out=100)
    b <- posterior %>% select(ends_with("b_for"))
    
    
    
    forestry_sqrt <- sqrt(forestry)
    
    forestry_sqrt_std = (forestry_sqrt-mean(forestry_sqrt))/sd(forestry_sqrt)
    
  }else if(effect == "cpd"){
    forestry <- seq(0,100,length.out=100)
    b <- posterior %>% select(ends_with("b_for"))
    
    
    forestry_sqrt <- sqrt(forestry)
    
    forestry_sqrt_std = (forestry_sqrt-mean(forestry_sqrt))/sd(forestry_sqrt)
    
  }else if(effect == "sst"){
    forestry <- seq(10,15,length.out=100)
    b <- posterior %>% select(ends_with("b_sst"))
    
    forestry_sqrt_std = (forestry-mean(forestry))/sd(forestry)
    
  }else if(effect == "npgo"){
    forestry <- seq(-3,3,length.out=100)
    b <- posterior %>% select(ends_with("b_npgo"))
    
    forestry_sqrt_std = (forestry-mean(forestry))/sd(forestry)
    
  }
  
  b_95 <- HDInterval::hdi(b, credMass = 0.95)
  b_80 <- HDInterval::hdi(b, credMass = 0.8)
  b_50 <- HDInterval::hdi(b, credMass = 0.5)
  
  
  no_forestry <- min(forestry_sqrt_std)
  
  
  
  # ch_chm_eca_sst=read.csv(here('stan models','outs','posterior',ch_chm_eca_sst),check.names=F)
  
  
  
  global_prediction <- apply(exp(as.matrix(b[,1])%*%
                                   (forestry_sqrt_std-no_forestry))*100 - 100,2,median)
  
  global_025 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.025), 
                      row.names = c("q025"))
  
  global_975 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.975),
                      row.names = c("q975"))
  global_750 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.75),
                      row.names = c("q750"))
  
  global_250 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.25),
                      row.names = c("q250"))
  
  global_900 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.9),
                      row.names = c("q900"))
  
  global_100 <- apply(exp(as.matrix(b[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.1),
                      row.names = c("q100"))
  
  # calculate high density credible intervals
  
  # hd_df <- apply(exp(as.matrix(b[,1])%*%
  #                      (forestry_sqrt_std-no_forestry))*100 - 100,2,HDInterval::hdi, 
  #                credMass = c(0.5,0.8,0.95))
  
  hd_50_lower <- (exp(as.matrix(b_50[,1])%*%
                        (forestry_sqrt_std-no_forestry))*100 - 100)[1,]
  
  hd_80_lower <- (exp(as.matrix(b_80[,1])%*%
                        (forestry_sqrt_std-no_forestry))*100 - 100)[1,]
  
  hd_95_lower <- (exp(as.matrix(b_95[,1])%*%
                        (forestry_sqrt_std-no_forestry))*100 - 100)[1,]
  
  hd_50_upper <- (exp(as.matrix(b_50[,1])%*%
                        (forestry_sqrt_std-no_forestry))*100 - 100)[2,]
  
  hd_80_upper <- (exp(as.matrix(b_80[,1])%*%
                        (forestry_sqrt_std-no_forestry))*100 - 100)[2,]
  
  hd_95_upper <- (exp(as.matrix(b_95[,1])%*%
                        (forestry_sqrt_std-no_forestry))*100 - 100)[2,]
                       
  hd_intervals <- data.frame(hd_50_lower = hd_50_lower,
                             hd_50_upper = hd_50_upper,
                             hd_80_lower = hd_80_lower,
                             hd_80_upper = hd_80_upper,
                             hd_95_lower = hd_95_lower,
                             hd_95_upper = hd_95_upper)
    
  global_df <- data.frame(forestry = forestry,
                          productivity_median = global_prediction,
                          q025 = global_025,
                          q975 = global_975,
                          q750 = global_750,
                          q250 = global_250,
                          q900 = global_900,
                          q100 = global_100) %>% 
    cbind(hd_intervals)
  
  full_productivity <- NULL
  
  if(hd == TRUE){
    ymin_95 = "hd_95_lower"
    ymax_95 = "hd_95_upper"
    ymin_80 = "hd_80_lower"
    ymax_80 = "hd_80_upper"
    ymin_50 = "hd_50_lower"
    ymax_50 = "hd_50_upper"
    
  } else{
    ymin_95 = "q025"
    ymax_95 = "q975"
    ymin_80 = "q100"
    ymax_80 = "q900"
    ymin_50 = "q250"
    ymax_50 = "q750"
    
  }
  
  if(by_river){
    
    # no_forestry <- min(eca_std_sqrt)
    for (i in 1:length(unique(df$River_n))){
      river <- unique(df$River_n)[i]
      
      river_data <- df %>% filter(River_n == river)
      
      if(effect == "eca"){
        b_rv <- posterior %>% select(starts_with("b_for_rv")) %>%
          select(ends_with(paste0("[",river,"]")))
        
        forestry_sqrt_std_river <- seq(min(forestry_sqrt_std),max(river_data$sqrt.ECA.std),length.out=100)
        #minimum forestry possible -0
        forestry_river <- seq(0,max(river_data$ECA_age_proxy_forested_only),length.out=100)
        
        
      } else if(effect == "cpd"){
        b_rv <- posterior %>% select(starts_with("b_for_rv")) %>%
          select(ends_with(paste0("[",river,"]")))
        
        forestry_sqrt_std_river <- seq(min(forestry_sqrt_std),max(river_data$sqrt.CPD.std),length.out=100)
        #minimum forestry possible - 0
        forestry_river <- seq(0,max(river_data$disturbedarea_prct_cs),length.out=100)
        
        
      } else if(effect == "sst"){
        b_rv <- posterior %>% select(starts_with("b_sst_rv")) %>%
          select(ends_with(paste0("[",river,"]")))
        
        # use only spring_ersst
        
        forestry_sqrt_std_river <- seq(min(river_data$sst.std),max(river_data$sst.std),length.out=100)
        
        forestry_river <- seq(min(river_data$spring_ersst),max(river_data$spring_ersst),length.out=100)
        
      }
      
      
      
      
      
      
      
      productivity <- (exp(as.matrix(b_rv[,1])%*%
                             (forestry_sqrt_std-no_forestry)))*100 - 100
      
      productivity_median <- apply(productivity,2,median)
      
      productivity_median_df <- data.frame(River = unique(river_data$River),
                                           productivity_median = productivity_median,
                                           forestry = forestry_river) %>% 
        filter(forestry >= min(forestry_river), forestry <= max(forestry_river))
      
      full_productivity <- rbind(full_productivity, productivity_median_df)
      
    }
    
    # c <- posterior %>% select(ends_with("b_for"))
    
    median_prediction <- apply(exp(as.matrix(b[,1])%*%
                                     (forestry_sqrt_std-no_forestry))*100 - 100,2,median)
    
    median_df <- data.frame(forestry = forestry,
                            productivity_median = median_prediction)
    
    if(effect == "eca" || effect == "cpd"){
      p1 <- ggplot(full_productivity) +
        geom_line(aes(x = forestry, y = productivity_median, group = River,
                      color = "watershed level\nforestry effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = productivity_median, color = "global forestry effect"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = q025, ymax = q975),
                    alpha = 0.5, fill = "#ADcCA5") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("watershed level\nforestry effect" = "darkgray", 
                                         "global forestry effect" = "black")) +
        ylim(-100,100) +
        scale_x_continuous(n.breaks = 5) +
        labs(title = model,
             x = ifelse(effect == "eca", "Equivalent clearcut area", "Cumulative disturbance (%)"),
             y = "Median change\n in recruitment (%)") +
        theme_classic() +
        theme(legend.position = c(0.8,0.8),
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    } else if(effect == "sst"){
      p1 <- ggplot(full_productivity) +
        geom_line(aes(x = forestry, y = productivity_median, group = River,
                      color = "watershed level\nSST effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = productivity_median, color = "global SST effect"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = q025, ymax = q975),
                    alpha = 0.5, fill = "gray") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("watershed level\nSST effect" = "darkgray", 
                                         "global SST effect" = "black")) +
        ylim(-100,100) +
        scale_x_continuous(n.breaks = 5, limits = c(10,15)) +
        labs(title = model,
             x = "Spring SST (°C)",
             y = "Median change\n in recruitment (%)") +
        theme_classic() +
        theme(legend.position = c(0.8,0.8),
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    }
  } else {
    
    if(effect == "eca" || effect == "cpd"){
      p1 <- ggplot(full_productivity) +
        # geom_line(aes(x = forestry, y = productivity_median, group = River,
        #               color = "watershed level\nforestry effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = productivity_median, color = "Median coastwide\nchange"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_95), ymax = !!sym(ymax_95),
                                          alpha = "95% credible\ninterval",
                                          fill = "95% credible\ninterval"))+
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_50), ymax = !!sym(ymax_50),
                                          alpha = "50% credible\ninterval",
                                          fill = "50% credible\ninterval"))+
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_80), ymax = !!sym(ymax_80),
                                          alpha = "80% credible\ninterval",
                                          fill = "80% credible\ninterval"))+
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("Median coastwide\nchange" = "#6F7c67")) +
        scale_fill_manual("",values  = c("95% credible\ninterval" =  "#ADcCA5",
                                         "50% credible\ninterval" = "#ADcCA5",
                                         "80% credible\ninterval" = "#ADcCA5"
        )) +
        scale_alpha_manual("",values  = c("95% credible\ninterval" = 0.25,
                                          "50% credible\ninterval" = 0.65,
                                          "80% credible\ninterval" = 0.45
        )) +
        ylim(-100,50) +
        scale_x_continuous(n.breaks = 5) +
        labs(title = model,
             x = ifelse(effect == "eca", "Equivalent clearcut area", "Cumulative disturbance (%)"),
             y = "Median change\n in recruitment (%)") +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_blank(),
              legend.key.size = unit(0.2, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)),
               
        )
    } else if(effect == "sst"){
      p1 <- ggplot(full_productivity) +
        # geom_line(aes(x = forestry, y = productivity_median, group = River,
        #               color = "watershed level\nSST effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = productivity_median, color = "global SST effect"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_95), ymax = !!sym(ymax_95)),
                    alpha = 0.25, fill = "#C78c63") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_50), ymax = !!sym(ymax_50)),
                    alpha = 0.65, fill = "#C78c63") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_80), ymax = !!sym(ymax_80)),
                    alpha = 0.45, fill = "#C78c63") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("global SST effect" = "#C78c63")) +
        
        ylim(-50,100) +
        scale_x_continuous(n.breaks = 5, limits = c(10,15)) +
        labs(title = model,
             x = "Spring SST (°C)",
             y = "Median change\n in recruitment (%)") +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    }else if(effect == "npgo"){
      p1 <- ggplot(full_productivity) +
        # geom_line(aes(x = forestry, y = productivity_median, group = River,
        #               color = "watershed level\nSST effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = productivity_median, color = "global NPGO effect"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin =!!sym(ymin_95), ymax = !!sym(ymax_95)),
                    alpha = 0.25, fill = "#829Dc6") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_50), ymax = !!sym(ymax_50)),
                    alpha = 0.65, fill = "#829Dc6") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_80), ymax = !!sym(ymax_80)),
                    alpha = 0.45, fill = "#829Dc6") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("global NPGO effect" = "#829Dc6")) +
        
        ylim(-50,100) +
        scale_x_continuous(n.breaks = 5, limits = c(-3,3)) +
        labs(title = model,
             x = "NPGO",
             y = "Median change\n in recruitment (%)") +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    }
    
    
    
    
  }
  
  
  
  return(p1)
  
}

plot_productivity_decline_manuscript4 <- function(posterior = ric_chm_cpd_ocean_covariates_logR_long_chain, 
                                                  effect = "cpd", species = "chum", 
                                                  model = "Ricker", 
                                                  by_river = FALSE, hd = TRUE){
  
  
  if(species == "chum"){
    df <- ch20rsc 
    df$sst.std <- (ch20rsc$spring_ersst-mean(ch20rsc$spring_ersst))/sd(ch20rsc$spring_ersst)
    
  } else if(species == "pink"){
    df <- pk10r
    df$sst.std <- (pk10r$spring_ersst-mean(pk10r$spring_ersst))/sd(pk10r$spring_ersst)
  }
  
  if(effect == "eca"){
    forestry <- seq(0,1,length.out=100)
    b <- posterior %>% select(ends_with("b_for"))
    
    b_posterior <- subset_draws(as_draws_df(posterior),"b_for") %>% 
      rename(b = "b_for")
    
    forestry_sqrt <- sqrt(forestry)
    
    forestry_sqrt_std = (forestry_sqrt-mean(forestry_sqrt))/sd(forestry_sqrt)
    
  }else if(effect == "cpd"){
    forestry <- seq(0,100,length.out=100)
    b <- posterior %>% select(ends_with("b_for"))
    
    b_posterior <- subset_draws(as_draws_df(posterior),"b_for") %>% 
      rename(b = "b_for")
    
    forestry_sqrt <- sqrt(forestry)
    
    forestry_sqrt_std = (forestry_sqrt-mean(forestry_sqrt))/sd(forestry_sqrt)
    
  }else if(effect == "sst"){
    forestry <- seq(10,15,length.out=100)
    b <- posterior %>% select(ends_with("b_sst"))
    
    b_posterior <- subset_draws(as_draws_df(posterior),"b_sst") %>% 
      rename(b = "b_sst")
    
    forestry_sqrt_std = (forestry-mean(forestry))/sd(forestry)
    
  }else if(effect == "npgo"){
    forestry <- seq(-3,3,length.out=100)
    b <- posterior %>% select(ends_with("b_npgo"))
    
    b_posterior <- subset_draws(as_draws_df(posterior),"b_npgo") %>% 
      rename(b = "b_npgo")
    
    forestry_sqrt_std = (forestry-mean(forestry))/sd(forestry)
    
  }
  
  # b_for_cu_posterior <- subset_draws(as_draws_df(posterior),paste0("b_for_cu[",cu,"]")) %>% 
  #   rename(b_for_cu = paste0("b_for_cu[",cu,"]"))
  recruitment_change_full_posterior = NULL
  for(forestry_values in forestry_sqrt_std){
    
    recruitment_change_posterior <- mutate_variables(b_posterior, 
                                                     recruitment_change = exp(b*(forestry_values - no_forestry))*100-100)
    
    recruitment_change_full_posterior <- cbind(recruitment_change_full_posterior, recruitment_change_posterior$recruitment_change)
    
    
  }
  recruitment_change_posterior <- mutate_variables(b_posterior, 
                                                   recruitment_change = exp(b*(forestry_sqrt_std - no_forestry))*100-100)
  
  recruitment_change_full_95 <- HDInterval::hdi(recruitment_change_full_posterior, credMass = 0.95)
  
  recruitment_change_full_80 <- HDInterval::hdi(recruitment_change_full_posterior, credMass = 0.80)
  
  recruitment_change_full_50 <- HDInterval::hdi(recruitment_change_full_posterior, credMass = 0.5)
  
  
  hd_intervals <- data.frame(ymin_50 = recruitment_change_full_50[1,],
                             ymax_50 = recruitment_change_full_50[2,],
                             ymin_80 = recruitment_change_full_80[1,],
                             ymax_80 = recruitment_change_full_80[2,],
                             ymin_95 = recruitment_change_full_95[1,],
                             ymax_95 = recruitment_change_full_95[2,])
  
  global_df <- data.frame(forestry = forestry,
                          recruitment_median = apply(recruitment_change_full_posterior,2, median) ) %>% 
    cbind(hd_intervals)
  
  
  if(effect == "eca" || effect == "cpd"){
      p1 <- ggplot() +
        # geom_line(aes(x = forestry, y = productivity_median, group = River,
        #               color = "watershed level\nforestry effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = recruitment_median, color = "Median coastwide\nchange"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = ymin_95, ymax = ymax_95,
                                          alpha = "95% credible\ninterval",
                                          fill = "95% credible\ninterval"))+
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_50), ymax = !!sym(ymax_50),
                                          alpha = "50% credible\ninterval",
                                          fill = "50% credible\ninterval"))+
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_80), ymax = !!sym(ymax_80),
                                          alpha = "80% credible\ninterval",
                                          fill = "80% credible\ninterval"))+
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("Median coastwide\nchange" = "#6F7c67")) +
        scale_fill_manual("",values  = c("95% credible\ninterval" =  "#ADcCA5",
                                         "50% credible\ninterval" = "#ADcCA5",
                                         "80% credible\ninterval" = "#ADcCA5"
        )) +
        scale_alpha_manual("",values  = c("95% credible\ninterval" = 0.25,
                                          "50% credible\ninterval" = 0.65,
                                          "80% credible\ninterval" = 0.45
        )) +
        ylim(-100,50) +
        scale_x_continuous(n.breaks = 5) +
        labs(title = model,
             x = ifelse(effect == "eca", "Equivalent clearcut area", "Cumulative disturbance (%)"),
             y = "Median change\n in recruitment (%)") +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_blank(),
              legend.key.size = unit(0.2, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)),
               
        )
    } else if(effect == "sst"){
      p1 <- ggplot(full_productivity) +
        # geom_line(aes(x = forestry, y = productivity_median, group = River,
        #               color = "watershed level\nSST effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = recruitment_median, color = "global SST effect"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_95), ymax = !!sym(ymax_95)),
                    alpha = 0.25, fill = "#C78c63") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_50), ymax = !!sym(ymax_50)),
                    alpha = 0.65, fill = "#C78c63") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_80), ymax = !!sym(ymax_80)),
                    alpha = 0.45, fill = "#C78c63") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("global SST effect" = "#C78c63")) +
        
        ylim(-50,100) +
        scale_x_continuous(n.breaks = 5, limits = c(10,15)) +
        labs(title = model,
             x = "Spring SST (°C)",
             y = "Median change\n in recruitment (%)") +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    }else if(effect == "npgo"){
      p1 <- ggplot(full_productivity) +
        # geom_line(aes(x = forestry, y = productivity_median, group = River,
        #               color = "watershed level\nSST effect"),alpha=0.5) +
        geom_line(data = global_df, 
                  aes(x = forestry, y = recruitment_median, color = "global NPGO effect"), linewidth = 1) +
        geom_ribbon(data = global_df, aes(x = forestry, ymin =!!sym(ymin_95), ymax = !!sym(ymax_95)),
                    alpha = 0.25, fill = "#829Dc6") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_50), ymax = !!sym(ymax_50)),
                    alpha = 0.65, fill = "#829Dc6") +
        geom_ribbon(data = global_df, aes(x = forestry, ymin = !!sym(ymin_80), ymax = !!sym(ymax_80)),
                    alpha = 0.45, fill = "#829Dc6") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("global NPGO effect" = "#829Dc6")) +
        
        ylim(-50,100) +
        scale_x_continuous(n.breaks = 5, limits = c(-3,3)) +
        labs(title = model,
             x = "NPGO",
             y = "Median change\n in recruitment (%)") +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(1, "cm"),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 10, hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5)))
    }
    
    
    
    
  
  
  
  
  return(p1)
  
}


plot_productivity_decline_manuscript1(model = "ET credible intervals")

plot_productivity_decline_manuscript2()

plot_productivity_decline_manuscript3(hd = TRUE, model = "HD credible intervals")
plot_productivity_decline_manuscript3(hd = FALSE, model = "ET credible intervals")


# plot full posterior of b and of recruitment change

ggplot() +
  stat_density(data= posterior, aes(x= b_for),
               geom = 'line', position = 'identity', 
               alpha = 0.05, linewidth = 1.5) +
  theme_classic()
  

posterior1 <- ric_chm_cpd_ocean_covariates_logR
posterior2 <- ric_chm_cpd_ocean_covariates_logR_long_chain
  
ggplot() +
  stat_density(data= posterior1, aes(x= exp(as.matrix(b_for)%*%
      (forestry_sqrt_std[100]-no_forestry))*100 - 100),
      geom = 'line', position = 'identity', 
      alpha = 0.5, linewidth = 1.5) +
  stat_density(data= posterior2, aes(x= exp(as.matrix(b_for)%*%
                                              (forestry_sqrt_std[100]-no_forestry))*100 - 100),
               geom = 'line', position = 'identity', color = "red",
               alpha = 0.5, linewidth = 1.5) +
  theme_classic()

# both posteriors are the same - long tail towards zero  
  



productivity_decline_cu_df <- function(posterior, effect, species){
  
  if(species == "chum"){
    df <- ch20rsc 
    
  } else if(species == "pink"){
    df <- pk10r
  }
  
  full_productivity <- NULL
  
  for (i in 1:length(unique(df$CU_n))){
    
    cu <- unique(df$CU_n)[i]
    
    cu_data <- df %>% filter(CU_n == cu)
    
    b_cu <- posterior %>% select(starts_with("b_for_cu")) %>%
      select(ends_with(paste0("[",cu,"]")))
    
    b_cu_l_95 <- HDInterval::hdi(b_cu, credMass = 0.95)[1]
    
    b_cu_u_95 <- HDInterval::hdi(b_cu, credMass = 0.95)[2]
    
    b_cu_l_50 <- HDInterval::hdi(b_cu, credMass = 0.5)[1]
    
    b_cu_u_50 <- HDInterval::hdi(b_cu, credMass = 0.5)[2]
    
    eca_sqrt_std_cu <- max(cu_data$sqrt.ECA.std)
    #minimum forestry possible -0
    eca_cu <- max(cu_data$ECA_age_proxy_forested_only)
    
    forestry_eca <- seq(0,1, length.out = 100)
    
    # cpd_sqrt_std_cu <- max(cu_data$sqrt.CPD.std) # should not be using max from CU
    #minimum forestry possible - 0
    cpd_cu <- cu_data %>% group_by(River) %>% 
      filter(disturbedarea_prct_cs == max(disturbedarea_prct_cs)) %>% 
      distinct(disturbedarea_prct_cs) %>% 
      ungroup %>% 
      summarize(mean = mean(disturbedarea_prct_cs))
    
    cpd_sqrt_std_cu <- cu_data %>% group_by(River) %>% 
      filter(sqrt.CPD.std == max(sqrt.CPD.std)) %>% 
      distinct(sqrt.CPD.std) %>% 
      ungroup %>% 
      summarize(mean = mean(sqrt.CPD.std))
    
    forestry_cpd <- seq(0,100, length.out = 100)
    
    #to calculate no forestry in the standardized scale
    forestry_sqrt <- sqrt(forestry_cpd)
    
    forestry_sqrt_std = (forestry_sqrt-mean(forestry_sqrt))/sd(forestry_sqrt)
    
    no_forestry <- min(forestry_sqrt_std)
    
    productivity <- (exp(as.matrix(b_cu[,1])%*%
                           (cpd_sqrt_std_cu$mean-no_forestry)))*100 - 100
    
    b_for_cu_posterior <- subset_draws(as_draws_df(posterior),paste0("b_for_cu[",cu,"]")) %>% 
      rename(b_for_cu = paste0("b_for_cu[",cu,"]"))
    
    recruitment_change_posterior <- mutate_variables(b_for_cu_posterior, 
                                              recruitment_change = exp(b_for_cu*(cpd_sqrt_std_cu$mean - no_forestry))*100-100)
    
    productivity_median <- apply(productivity,2,median)
    
    productivity_median_df <- data.frame(CU = unique(cu_data$CU_name),
                                         productivity_50 = apply(productivity,2,median),
                                         productivity_25 = apply(productivity,2,quantile, probs = 0.25),
                                         productivity_75 = apply(productivity,2,quantile, probs = 0.75),
                                         productivity_025 = apply(productivity,2,quantile, probs = 0.025),
                                         productivity_975 = apply(productivity,2,quantile, probs = 0.975),
                                         productivity_025_hdi = (exp(as.matrix(b_cu_l_95)%*%
                                                                       (cpd_sqrt_std_cu$mean-no_forestry)))*100 - 100,
                                         productivity_975_hdi = (exp(as.matrix(b_cu_u_95)%*%
                                                                       (cpd_sqrt_std_cu$mean-no_forestry)))*100 - 100,
                                         productivity_25_hdi = (exp(as.matrix(b_cu_l_50)%*%
                                                                     (cpd_sqrt_std_cu$mean-no_forestry)))*100 - 100,
                                         productivity_75_hdi = (exp(as.matrix(b_cu_u_50)%*%
                                                                     (cpd_sqrt_std_cu$mean-no_forestry)))*100 - 100,
                                         
                                         productivity_025_hdi_new = HDInterval::hdi(recruitment_change_posterior$recruitment_change, credMass = 0.95)[1],
                                         productivity_975_hdi_new = HDInterval::hdi(recruitment_change_posterior$recruitment_change, credMass = 0.95)[2],
                                         productivity_25_hdi_new = HDInterval::hdi(recruitment_change_posterior$recruitment_change, credMass = 0.5)[1],
                                         productivity_75_hdi_new = HDInterval::hdi(recruitment_change_posterior$recruitment_change, credMass = 0.5)[2],
                                         productivity_50_new = median(recruitment_change_posterior$recruitment_change),
                                         productivity_25_new = quantile(recruitment_change_posterior$recruitment_change, probs = 0.25),
                                         productivity_75_new = quantile(recruitment_change_posterior$recruitment_change, probs = 0.75),
                                         productivity_025_new = quantile(recruitment_change_posterior$recruitment_change, probs = 0.025),
                                         productivity_975_new = quantile(recruitment_change_posterior$recruitment_change, probs = 0.975),
                                         
                                         forestry = cpd_cu$mean,
                                         CU_n = unique(cu_data$CU_n))
    
    full_productivity <- rbind(full_productivity, productivity_median_df)
    
    
  }
  
  
  return(full_productivity)
  
  
}


ric_chm_cpd_productivity_decline_cu <- productivity_decline_cu_df(ric_chm_cpd_ocean_covariates_logR, 
                                                                  effect = "cpd", species = "chum")    

cu_forest_plot_compare <- ric_chm_cpd_productivity_decline_cu %>% 
  arrange(desc(productivity_50)) %>%
  mutate(CU2 = factor(CU, levels = CU)) %>% 
  ggplot(aes(x = CU, y = productivity_50)) +
  geom_point(aes(y = productivity_50, x = CU2), color = '#516479',fill = "white", size = 3, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_025, ymax = productivity_975, color = "95% ETI CI"), width = 0, alpha = 0.5, size = 1) +
  geom_errorbar(aes(ymin = productivity_025_hdi, ymax = productivity_975_hdi, color = "95% HDI CI"), width = 0, alpha = 0.5, size = 1) +
  geom_errorbar(aes(ymin = productivity_25, ymax = productivity_75 ), color = '#516479', width = 0, alpha = 0.7, size = 2) +
  geom_errorbar(aes(ymin = productivity_25_hdi, ymax = productivity_75_hdi ), color = 'hotpink', width = 0, alpha = 0.7, size = 2) +
  scale_color_manual(name = 'Interval type', values = c('95% ETI CI' = '#516479', 
                                                        '50% ETI CI' = '#516479', 
                                                        '95% HDI CI' = 'hotpink'))+
  #add estimated median decline to the right of each error bar
  geom_text(aes(label = paste(round(productivity_50,1),"%")), 
            hjust = -0.25, 
            vjust = -0.35,
            size = 3, color = "gray20") +
  #add dashed v line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  #label the 95%credible interval using annotate
  
  coord_flip() +
  # scale_color_manual(name = 'Model type', values = c('independent alpha' = 'cadetblue', 'hierarchical alpha' = 'coral', 'hierarchical alpha - ricker' = 'darkgoldenrod')) +
  labs(title = 'CI',
    x = 'Conservation Unit',
    y = 'Change in recruitment (%)') +
  theme_classic() +
  theme(legend.position = c(0.8,0.9),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent', color = "transparent"),
        legend.box.background = element_rect(fill='transparent', color = 'transparent'))+
  #increase xlim
  scale_x_discrete(expand = expansion(mult = 0.05))

cu_forest_plot_compare


cpd <- 1
no_forestry

recruitment_change_df <- mutate_variables(subset_draws(as_draws_df(posterior),"b_for"), 
                                          recruitment_change = exp(b_for*(cpd - no_forestry))*100-100)

# calculate 95% equal tailed and 95% hd credible intervals

eti_ci <- recruitment_change_df %>%
  summarise(eti_lower = quantile(recruitment_change, probs = 0.025),
            eti_upper = quantile(recruitment_change, probs = 0.975),
            eti_median = median(recruitment_change))

hd_ci <- HDInterval::hdi(recruitment_change_df$recruitment_change, credMass = 0.95)

ric_chm_cpd_productivity_decline_cu_new <- productivity_decline_cu_df(ric_chm_cpd_ocean_covariates_logR, 
                                                                  effect = "cpd", species = "chum")    


cu_forest_plot_compare_new <- ric_chm_cpd_productivity_decline_cu_new %>% 
  arrange(desc(productivity_50)) %>%
  mutate(CU2 = factor(CU, levels = CU)) %>% 
  ggplot(aes(x = CU, y = productivity_50_new)) +
  geom_point(aes(y = productivity_50_new, x = CU2), color = '#516479',fill = "white", size = 3, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_025_new, ymax = productivity_975_new, color = "95% ETI CI"), width = 0, alpha = 0.5, size = 1) +
  geom_errorbar(aes(ymin = productivity_025_hdi_new, ymax = productivity_975_hdi_new, color = "95% HDI CI"), width = 0, alpha = 0.5, size = 1) +
  geom_errorbar(aes(ymin = productivity_25_new, ymax = productivity_75_new ), color = '#516479', width = 0, alpha = 0.7, size = 2) +
  geom_errorbar(aes(ymin = productivity_25_hdi_new, ymax = productivity_75_hdi_new ), color = 'hotpink', width = 0, alpha = 0.7, size = 2) +
  scale_color_manual(name = 'Interval type', values = c('95% ETI CI' = '#516479', 
                                                        '50% ETI CI' = '#516479', 
                                                        '95% HDI CI' = 'hotpink'))+
  #add estimated median decline to the right of each error bar
  geom_text(aes(label = paste(round(productivity_50_new,1),"%")), 
            hjust = -0.25, 
            vjust = -0.35,
            size = 3, color = "gray20") +
  #add dashed v line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  #label the 95%credible interval using annotate
  
  coord_flip() +
  # scale_color_manual(name = 'Model type', values = c('independent alpha' = 'cadetblue', 'hierarchical alpha' = 'coral', 'hierarchical alpha - ricker' = 'darkgoldenrod')) +
  labs(title = 'CI from posterior package',
    x = 'Conservation Unit',
    y = 'Change in recruitment (%)') +
  theme_classic() +
  theme(legend.position = c(0.8,0.9),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent', color = "transparent"),
        legend.box.background = element_rect(fill='transparent', color = 'transparent'))+
  #increase xlim
  scale_x_discrete(expand = expansion(mult = 0.05))

cu_forest_plot_compare_new
  
# conclusion - use equal tailed distributions 
# find out why the estimates of median are different with new posterior in the supplement.






