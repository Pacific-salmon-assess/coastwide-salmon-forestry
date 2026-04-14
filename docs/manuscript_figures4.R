# plot figures with new pink model with sigma vector


library(here);library(dplyr); library(stringr)
library(ggplot2)
library(tidyverse)
library(bayesplot)
library(patchwork)
library(hues)
library(GGally)
library(latex2exp)
library(ggrepel)
library(latex2exp)
library(ggpubr)
library(bayestestR)
library(ggrepel)


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

# change "East Hg" to"East Haida Gwaii"
ch20rsc$CU_name <- ifelse(ch20rsc$CU_name == "East Hg", "East Haida Gwaii", ch20rsc$CU_name)


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



plot_forestry_effect_manuscript <- function(posterior = bh_chm_eca, 
                                            effect = "eca", 
                                            species = "chum", 
                                            model = "BH model with NPGO", 
                                            xlim = c(-2, 2), 
                                            by_cu = FALSE, 
                                            global_only = FALSE){
  x_name = "Standardized coefficients"
  
  if(effect == "eca"){
    
    y_name = "Posterior density of\nECA effect"
    # color_var = "eca_level"
    # color_name = "Max ECA level\nin River (%)"
    # scale_color = scale_color_manual(name = "Max ECA level\nin River (%)",
    #                                  values = c('low' = '#35978f', 'medium' = 'gray', 'high' = '#cf812d'))
    # 
    
    
  } else if(effect == "cpd"){
    # x_name = "Standardized coefficients of CPD"
    y_name = "Posterior density of\ncumulative disturbance effect"
    # color_var = "cpd_max"
    # color_name = "Max CPD\nin River (%)"
    # scale_color = scale_color_gradient2(name = 'Max CPD\nin River (%)',
    #                                     low = '#35978f', mid = 'gray', high = '#cf812d', midpoint = 50)
    
  }
  if(by_cu){
    posterior_df <- posterior %>%
      select(starts_with('b_for_cu')) %>%
      pivot_longer(cols = everything(), 
                   names_to = 'CU', 
                   names_prefix = 'b_for_cu',
                   values_to = effect) %>%
      mutate(CU_n = as.numeric(str_extract(CU, '\\d+'))) %>% 
      select(-CU) %>%
      left_join(ch20rsc %>% select(CU, CU_n, CU_name) %>% distinct(), by = 'CU_n') 
    
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
    
    #color by CU and lacel the CU which has minimum and maximum median effect
    plot1 <- ggplot() +
      stat_density(data= posterior_df, aes(!!sym(effect), color = !!sym(color_var), group = CU_n),
                   geom = 'line', position = 'identity', 
                   alpha = 0.05, linewidth = 1.5) +
      geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
      labs(title = model, x = x_name, y = y_name) +
      xlim(xlim[1], xlim[2]) +
      scale_color +
      geom_density(aes(posterior$b_for), color = 'black', linewidth = 1.2, alpha = 0.2)+
      # geom_text(data = min_cu, aes(x = median_effect, y = 1, lacel = CU_name, color = !!sym(color_var)), size = 4)+
      # geom_text(data = max_cu, aes(x = median_effect, y = 1, lacel = CU_name, color = !!sym(color_var)), size = 4)+
      # #vline at the median value of the posterior
      geom_vline(xintercept = median(posterior$b_for), color = 'black', linetype = 'dashed')+
      theme_classic()+
      theme(legend.position = "right",
            legend.key.width = unit(0.5, "cm"),
            legend.key.height = unit(0.8, "lines"),
            #1 column for legend
            legend.cox = "vertical",
            legend.text = element_text(size = 7),
            legend.spacing.y = unit(0.001, "cm"),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            # plot.title = element_text(size = 10, hjust = 0.5)
            plot.title = element_text(size = 10, hjust = 0, vjust = -0.1)
      )+
      guides(color = guide_legend(override.aes = list(alpha = 1), ncol =1))
    
  } else {
    
    
    posterior_df <- posterior %>%
      select(starts_with('b_for_rv')) %>%
      pivot_longer(cols = everything(), 
                   names_to = 'River', 
                   names_prefix = 'b_for_rv',
                   values_to = effect) %>%
      mutate(River_n = as.numeric(str_extract(River, '\\d+')))# %>% 
    # select(-River) 
    
    #color by river
    plot1 <- ggplot() +
      stat_density(data= posterior_df, aes(!!sym(effect), 
                                           group = River), 
                   color = "#ADcCA5",
                   geom = 'line', position = 'identity', 
                   alpha = 0.05, linewidth = 1) + 
      geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
      labs(title = model, x = x_name, y = y_name) +
      xlim(xlim[1], xlim[2]) +
      # scale_color +
      geom_density(aes(posterior$b_for), color = 'black', linewidth = 1.2, alpha = 0.2)+
      #vline at the median value of the posterior
      geom_vline(xintercept = median(posterior$b_for), color = 'black', linetype = 'dashed')+
      theme_classic()+
      theme(legend.position = 'right',
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            legend.text = element_text(size = 7),
            legend.spacing.y = unit(0.001, "cm"),
            legend.key.width = unit(0.5, "cm"),
            legend.key.height = unit(0.8, "lines"),
            # plot.title = element_text(size = 10, hjust = 0.5)
            plot.title = element_text(size = 10, hjust = 0, vjust = -0.1)
      )+
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.6), ncol = 1))
    
    if(global_only){
      
      
      plot1 <- ggplot() +
        geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
        labs(title = model, x = x_name, y = y_name) +
        xlim(xlim[1], xlim[2]) +
        # scale_color +
        geom_density(aes(posterior$b_for), color = 'black', 
                     linewidth = 1.2, alpha = 0.1)+
        geom_text(x = median(posterior$b_for), y = 1, 
                  lacel = paste("Median:", round(median(posterior$b_for), 2)),
                  color = 'black', size = 4, hjust = 0.5)+
        #vline at the median value of the posterior
        geom_vline(xintercept = median(posterior$b_for), color = 'black', linetype = 'dashed')+
        theme_classic()+
        theme(legend.position = c(0.8,0.5),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              legend.text = element_text(size = 7),
              legend.spacing.y = unit(0.001, "cm"),
              legend.key.width = unit(0.5, "cm"),
              legend.key.height = unit(0.8, "lines"),
              # plot.title = element_text(size = 10, hjust = 0.5)
              plot.title = element_text(size = 10, hjust = 0, vjust = -0.1)
        )+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.6), ncol = 1))
      
    } else{
      
      
      posterior_df <- posterior %>%
        select(starts_with('b_for_rv')) %>%
        pivot_longer(cols = everything(), 
                     names_to = 'River', 
                     names_prefix = 'b_for_rv',
                     values_to = effect) %>%
        mutate(River_n = as.numeric(str_extract(River, '\\d+'))) #%>% 
      # select(-River) 
      
      #color by river
      plot1 <- ggplot() +
        stat_density(data= posterior_df, aes(!!sym(effect), 
                                             group = River, 
                                             color = 'River'),
                     geom = 'line', position = 'identity', 
                     alpha = 0.05, linewidth = 1) +
        geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
        labs(title = model, x = x_name, y = y_name) +
        xlim(xlim[1], xlim[2]) +
        # scale_color +
        stat_density(aes(posterior$b_for, color = 'Coastwide'), geom = 'line', position = 'identity', 
                     linewidth = 1.2, alpha = 1)+
        #vline at the median value of the posterior
        geom_vline(xintercept = median(posterior$b_for), color = 'black', linetype = 'dashed')+
        scale_color_manual(values = c('River' = '#ADcCA5', 'Coastwide' = 'black'))+
        theme_classic()+
        theme(legend.position = c(0.8,0.5),
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              legend.text = element_text(size = 7),
              legend.spacing.y = unit(0.001, "cm"),
              legend.key.width = unit(0.5, "cm"),
              legend.key.height = unit(0.8, "lines"),
              # plot.title = element_text(size = 10, hjust = 0.5)
              plot.title = element_text(size = 10, hjust = 0, vjust = -0.1)
        )+
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.6), ncol = 1))
      
    }
    
  }
  
  return(plot1)
}


plot_productivity_decline_manuscript <- function(posterior = bh_chm_eca_sst, 
                                                 effect = "eca", species = "chum", 
                                                 model = "BH model with NPGO, LH SST", 
                                                 by_river = FALSE, hd = FALSE){
  
  
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
                       (forestry_sqrt_std-no_forestry))*100 - 100,2,hdi, 
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

plot_sst_effect_manuscript <- function(posterior = bh_chm_eca_sst, effect = "sst", species = "chum", 
                                       model = "BH model with NPGO", 
                                       color_by = "none", xlim = c(-1, 1)) {
  
  if(species == "chum"){
    df <- ch20rsc
    posterior_rv <- posterior %>% 
      dplyr::select(starts_with('b_sst_rv')) %>%
      pivot_longer(cols = everything(), 
                   names_to = 'River', 
                   names_prefix = 'b_sst_rv',
                   values_to = 'sst_effect') %>%
      mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
      dplyr::select(-River) %>%
      left_join(df %>% dplyr::select(River_n, CU, CU_name, X_LONG, Y_LAT) %>% 
                  distinct(), by = 'River_n')
    
  } else if(species == "pink"){
    df <- pk10r
    posterior_rv <- posterior %>% 
      dplyr::select(starts_with('b_sst_rv')) %>%
      pivot_longer(cols = everything(), 
                   names_to = 'River', 
                   names_prefix = 'b_sst_rv',
                   values_to = 'sst_effect') %>%
      mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
      dplyr::select(-River) %>%
      left_join(df %>% dplyr::select(River_n, X_LONG, Y_LAT) %>% 
                  distinct(), by = 'River_n')
  }
  
  
  
  # if(color_by == "CU_name"){
  #   color_var = "CU_name"
  #   color_name = "CU"
  #   scale_color = scale_color_iwanthue(name = 'Conservation Unit', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
  #                                      lmin = 10, lmax = 95)
  #   
  #   
  # } else if(color_by == "Y_LAT"){
  #   color_var = "Y_LAT"
  #   color_name = "Latitude"
  #   scale_color = scale_color_viridis_c(name = 'Latitude', direction = -1)
  #   
  #   
  # } else if(color_by == "CU"){
  #   
  #   color_var = "CU"
  #   color_name = "CU"
  #   scale_color = scale_color_iwanthue(name = 'Conservation Unit', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
  #                                      lmin = 10, lmax = 95)
  #   
  # } else if(color_by == "none"){
  #   scale_color = scale_color_manual()
  # }
  
  if(species == "chum"){
    p1 <- ggplot()+
      stat_density(data= posterior_rv, aes(sst_effect, 
                                           # color = !!sym(color_by), 
                                           group = River_n),
                   color = "#C78c63",
                   geom = 'line', position = 'identity', 
                   alpha = 0.03, linewidth = 1) +
      geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
      labs(title = model, 
           x = 'Standardized coefficients', y = 'Posterior density of\nSST effect') +
      xlim(xlim[1], xlim[2]) +
      geom_density(aes(posterior$b_sst), color = 'black', linewidth = 1.2, alpha = 0.2)+
      # scale_color_iwanthue(name = 'CU', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
      #                      lmin = 10, lmax = 95) +
      # scale_color +
      geom_vline(xintercept = median(posterior$b_sst), color = 'black', linetype = 'dashed') + 
      theme_classic()+
      theme(legend.position = "right",
            legend.key.width = unit(0.5, "cm"),
            legend.key.height = unit(0.5, "lines"),
            legend.text = element_text(size = 7),
            legend.spacing.y = unit(0.002, "cm"),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            # plot.title = element_text(size = 10, hjust = 0.5)
            plot.title = element_text(size = 10, hjust = 0, vjust = -0.1)
      )+
      
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 2), ncol = 1))
  } else if(species == "pink"){
    p1 <- ggplot()+
      stat_density(data= posterior_rv, aes(sst_effect, group = River_n),
                   geom = 'line', position = 'identity', color = "#C78c63",
                   alpha = 0.03, linewidth = 1) +
      geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
      labs(title = model, 
           x = 'Standardized coefficients', y = 'Posterior density of\nSST effect') +
      xlim(xlim[1], xlim[2]) +
      geom_density(aes(posterior$b_sst), color = 'black', linewidth = 1.2, alpha = 0.2)+
      # scale_color_iwanthue(name = 'CU', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
      #                      lmin = 10, lmax = 95) +
      geom_vline(xintercept = median(posterior$b_sst), color = 'black', linetype = 'dashed') + 
      theme_classic()+
      theme(legend.position = "right",
            legend.key.width = unit(0.5, "cm"),
            legend.key.height = unit(0.5, "lines"),
            legend.text = element_text(size = 7),
            legend.spacing.y = unit(0.002, "cm"),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            # plot.title = element_text(size = 10, hjust = 0.5)
            plot.title = element_text(size = 10, hjust = 0, vjust = -0.1)
      )+
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 2), ncol = 1))
  }
  
  return(p1)
  
}


plot_npgo_effect_manuscript <- function(posterior = bh_chm_eca_npgo, 
                                        effect = "npgo", species = "chum", 
                                        model = "BH model with NPGO", 
                                        xlim = c(-0.5,0.5)) {
  
  
  
  if(species == "chum"){
    df <- ch20rsc
    posterior_rv <- posterior %>% 
      select(starts_with('b_npgo_rv')) %>%
      pivot_longer(cols = everything(), 
                   names_to = 'River', 
                   names_prefix = 'b_npgo_rv',
                   values_to = 'npgo_effect') %>%
      mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
      select(-River) %>%
      left_join(df %>% select(River_n, CU_name, X_LONG, Y_LAT) %>% distinct(), by = 'River_n')
  } else if(species == "pink"){
    df <- pk10r
    posterior_rv <- posterior %>% 
      select(starts_with('b_npgo_rv')) %>%
      pivot_longer(cols = everything(), 
                   names_to = 'River', 
                   names_prefix = 'b_npgo_rv',
                   values_to = 'npgo_effect') %>%
      mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
      select(-River) %>%
      left_join(df %>% select(River_n, X_LONG, Y_LAT) %>% distinct(), by = 'River_n')
  }
  
  
  
  if(species == "chum"){
    p1 <- ggplot()+
      stat_density(data= posterior_rv, aes(npgo_effect, 
                                           # color = CU_name, 
                                           group = River_n),
                   color = "#829Dc6",
                   geom = 'line', position = 'identity', 
                   alpha = 0.03, linewidth = 1) +
      geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
      labs(title = model, 
           x = 'Standardized coefficients', y = 'Posterior density of\nNPGO effect') +
      xlim(xlim[1], xlim[2]) +
      geom_density(aes(posterior$b_npgo), color = 'black', linewidth = 1.2, alpha = 0.2)+
      scale_color_iwanthue(name = 'Conservation Unit', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
                           lmin = 10, lmax = 95) +
      geom_vline(xintercept = median(posterior$b_npgo), color = 'black', linetype = 'dashed') + 
      theme_classic()+
      theme(legend.position = "right",
            legend.key.width = unit(0.5, "cm"),
            legend.key.height = unit(0.5, "lines"),
            legend.text = element_text(size = 7),
            legend.spacing.y = unit(0.002, "cm"),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            # plot.title = element_text(size = 10, hjust = 0.5)
            plot.title = element_text(size = 10, hjust = 0, vjust = -0.1)
      )+
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 2), ncol = 1))
    
    
  } else if(species == "pink"){
    p1 <- ggplot()+
      stat_density(data= posterior_rv, aes(npgo_effect, group = River_n),
                   geom = 'line', position = 'identity', color ="#829Dc6",
                   alpha = 0.03, linewidth = 1) +
      geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
      labs(title = model, 
           x = 'Standardized coefficients', y = 'Posterior density of\nNPGO effect') +
      xlim(xlim[1], xlim[2]) +
      geom_density(aes(posterior$b_npgo), color = 'black', linewidth = 1.2, alpha = 0.2)+
      # scale_color_iwanthue(name = 'CU', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
      #                      lmin = 10, lmax = 95) +
      geom_vline(xintercept = median(posterior$b_npgo), color = 'black', linetype = 'dashed') + 
      theme_classic()+
      theme(legend.position = "right",
            legend.key.width = unit(0.5, "cm"),
            legend.key.height = unit(0.5, "lines"),
            legend.text = element_text(size = 7),
            legend.spacing.y = unit(0.002, "cm"),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            # plot.title = element_text(size = 10, hjust = 0.5)
            plot.title = element_text(size = 10, hjust = 0, vjust = -0.1)
      )+
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 2), ncol = 1))
  }
  
  
  
  return(p1)
  
}



ric_chm_eca_ocean_covariates_logR=read.csv(here('stan models','outs','posterior',
                                                'ric_chm_eca_ocean_covariates_logR_long_chain.csv'),check.names=F)
ric_chm_cpd_ocean_covariates_logR=read.csv(here('stan models','outs','posterior',
                                                'ric_chm_cpd_ocean_covariates_logR_long_chain.csv'),check.names=F)



ric_pk_eca_ersst = read.csv(here('stan models',
                                 'outs',
                                 'posterior',
                                 'ric_pk_eca_st_noac_ocean_covariates_logR_long_chain_sigma_vector.csv'),check.names=F)

ric_pk_cpd_ersst = read.csv(here('stan models',
                                 'outs',
                                 'posterior',
                                 'ric_pk_cpd_st_noac_ocean_covariates_logR_long_chain_sigma_vector.csv'),check.names=F)





# Pink 

plot_pink_G <- plot_forestry_effect_manuscript(posterior = ric_pk_cpd_ersst, 
                                               effect = "cpd", species = "pink", 
                                               # model = "Ricker model with CPD", xlim = c(-1, 1))
                                               model = "Cumulative disturbance", xlim = c(-0.5, 0.5))+
  theme(legend.position = "none")+
  labs(y = "Posterior density")+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))

plot_pink_H <- plot_forestry_effect_manuscript(posterior = ric_pk_eca_ersst,
                                               effect = "eca", species = "pink", 
                                               # model = "Ricker model with ECA", xlim = c(-1, 1))
                                               model = "Equivalent clearcut area", xlim = c(-0.5, 0.5))+
  theme(legend.position = "none")+
  labs(y = "Posterior density")+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))

plot_pink_I <- plot_sst_effect_manuscript(posterior = ric_pk_cpd_ersst,
                                          effect = "sst", species = "pink", 
                                          # model = "Ricker model with CPD", 
                                          model = "Sea-surface temperature",
                                          # xlim = c(-1, 1))
                                          xlim = c(-0.5, 0.5))+
  labs(y = "Posterior density")+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))


plot_pink_J <- plot_npgo_effect_manuscript(posterior = ric_pk_cpd_ersst,
                                           effect = "npgo", species = "pink", 
                                           # model = "Ricker model with CPD",
                                           model = "North Pacific Gyre Oscillation",
                                           # xlim = c(-1, 1))
                                           xlim = c(-0.5, 0.5))+
  labs(y = "Posterior density")+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))


plot_pink_K <- plot_productivity_decline_manuscript(posterior = ric_pk_cpd_ersst,  
                                                    species = "pink",
                                                    effect = "cpd", 
                                                    # model = "Ricker model with CPD")
                                                    model = "", hd = TRUE)+
  ylim(c(-75,75))+
  labs(y = "Change in recruitment (%)")+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))

plot_pink_L <- plot_productivity_decline_manuscript(posterior = ric_pk_eca_ersst,
                                                    species = "pink",
                                                    effect = "eca", 
                                                    # model = "Ricker model with ECA")
                                                    model = "", hd = TRUE)+
  ylim(c(-75,75))+
  labs(y = "Change in recruitment (%)")+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))



plot_pink2 <- ((plot_pink_G / plot_pink_H / plot_pink_I / plot_pink_J)+ plot_layout(axes = 'collect') |
                 ((plot_pink_K  / plot_pink_L ) + plot_layout(axes = 'collect'))) + plot_layout(widths = c(1.7,1))

plot_pink2


