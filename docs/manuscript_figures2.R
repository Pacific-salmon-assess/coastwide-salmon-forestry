# changes need to be made to the current mansucript figures (as of nov 2025)

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

# make multiple versions of figure results for the manuscript
# make one version with no colours of cpd or eca
# make one version with only cu level forestry districutions



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

L_o=pk10r_o%>%group_by(River)%>%summarize(l=n(),by=min(BroodYear),tmin=(min(BroodYear)-min(pk10r_o$croodYear))/2+1,tmax=(max(BroodYear)-min(pk10r_o$BroodYear))/2)
L_e=pk10r_e%>%group_by(River)%>%summarize(l=n(),by=min(BroodYear),tmin=(min(BroodYear)-min(pk10r_e$croodYear))/2+1,tmax=(max(BroodYear)-min(pk10r_e$BroodYear))/2)
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



plot_forestry_effect_manuscript <- function(posterior = ch_chm_eca, 
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

plot_sst_effect_manuscript <- function(posterior = ch_chm_eca_sst, effect = "sst", species = "chum", 
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


plot_npgo_effect_manuscript <- function(posterior = ch_chm_eca_npgo, 
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
                                                'ric_chm_eca_ocean_covariates_logR.csv'),check.names=F)
ric_chm_cpd_ocean_covariates_logR=read.csv(here('stan models','outs','posterior',
                                                'ric_chm_cpd_ocean_covariates_logR.csv'),check.names=F)



ric_pk_eca_ersst = read.csv(here('stan models',
                                 'outs',
                                 'posterior',
                                 'ric_pk_eca_st_noac_ocean_covariates_logR.csv'),check.names=F)

ric_pk_cpd_ersst = read.csv(here('stan models',
                                 'outs',
                                 'posterior',
                                 'ric_pk_cpd_st_noac_ocean_covariates_logR.csv'),check.names=F)






#new plot for manuscript without plot numbers, with decline in recruitment, with same y axis for % change



plot_chum_A <- plot_forestry_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                               effect = "cpd", species = "chum", 
                                               model = "Cumulative disturbance", xlim = c(-0.5, 0.5))+
  theme(legend.position = c(0.8,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 8))+
  labs(y = "Posterior density")+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))

#add legend for credible intervals for one figure
plot_chum_B <- plot_forestry_effect_manuscript(posterior = ric_chm_eca_ocean_covariates_logR, 
                                               effect = "eca", species = "chum", 
                                               model = "Equivalent clearcut area", xlim = c(-0.5, 0.5))+
  theme(legend.position = "none")+
  labs(y = "Posterior density")+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))

plot_chum_C <- plot_sst_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                          effect = "sst", species = "chum", 
                                          model = "Sea-surface temperature", 
                                          xlim = c(-0.5, 0.5))+
  theme(legend.position = "none")+
  labs(y = "Posterior density")+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))

plot_chum_D <- plot_npgo_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                           effect = "npgo", species = "chum", 
                                           model = "North Pacific Gyre Oscillation", 
                                           xlim = c(-0.5, 0.5))+
  theme(legend.position = "none")+
  labs(y = "Posterior density")+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))


plot_chum_E <- plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                                    effect = "cpd", model = "", hd = TRUE)+
  ylim(c(-75,75))+
  labs(y = "Change in recruitment (%)")+
  theme(legend.position = c(0.95,0.99),
        legend.direction = "vertical",
        legend.key.size = unit(0.001, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.spacing = unit(0.00, "cm"),
        legend.justification = c("right", "top"),
        legend.box.background = element_rect(color = "transparent", fill = "transparent"),
        legend.box.margin = margin(0,0,0,0),
        #make it horizontal
        legend.box = "horizontal",
        legend.text = element_text(size = 7))+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))
# guides(fill=guide_legend(nrow=1, byrow=TRUE))

plot_chum_F <- plot_productivity_decline_manuscript(posterior = ric_chm_eca_ocean_covariates_logR,
                                                    effect = "eca", model = "", hd = TRUE)+
  ylim(c(-75,75))+
  labs(y = "Change in recruitment (%)")+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))
# 
# 
# plot_chum_G <- plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR,
#                                                     effect = "sst", model = "Sea-surface temperature")+
#   ylim(c(-75,75))+
#   labs(y = "Median change in recruitment (%)")
# 
# plot_chum_H <- plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR,
#                                                     effect = "npgo", model = "North Pacific Gyre Oscillation")+
#   ylim(c(-75,75))+
#   labs(y = "Median change in recruitment (%)")



# plot_chum <- (plot_chum_a1 / plot_chum_a2 / plot_chum_a3 / plot_chum_a4)+ plot_layout(axes = 'collect') | 
#   ((plot_chum_c1  + plot_chum_c2 ) + plot_layout(axes = 'collect'))/(
#     (plot_chum_c3  + plot_chum_c4 ) + plot_layout(axes = 'collect')) 




plot_chum2 <- ((plot_chum_A / plot_chum_B / plot_chum_C / plot_chum_D)+ plot_layout(axes = 'collect') | 
  ((plot_chum_E  / plot_chum_F ) + plot_layout(axes = 'collect')) ) + plot_layout(widths = c(1.7,1))
# plot_chum[[1]] <- plot_chum[[1]] + plot_layout(tag_level = 'new') 
# 

plot_chum2


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

# plot_pink_O <- plot_productivity_decline_manuscript(posterior = ric_pk_cpd_ersst,
#                                                     species = "pink",
#                                                     effect = "sst", 
#                                                     # model = "Ricker model with CPD")
#                                                     model = "Sea-surface temperature")+
#   ylim(c(-75,75))+
#   labs(y = "Median change in recruitment (%)")
# 
# plot_pink_P <- plot_productivity_decline_manuscript(posterior = ric_pk_cpd_ersst,
#                                                     species = "pink",
#                                                     effect = "npgo", 
#                                                     # model = "Ricker model with CPD")
#                                                     model = "North Pacific Gyre Oscillation")+
#   ylim(c(-75,75))+
#   labs(y = "Median change in recruitment (%)")






# plot_pink <- (plot_pink_a1 / plot_pink_a2 / plot_pink_a3 / plot_pink_a4)+ plot_layout(axes = 'collect') |
#   ((plot_pink_c1  + plot_pink_c2 ) + plot_layout(axes = 'collect'))/(
#     (plot_pink_c3  + plot_pink_c4 ) + plot_layout(axes = 'collect'))


plot_pink2 <- ((plot_pink_G / plot_pink_H / plot_pink_I / plot_pink_J)+ plot_layout(axes = 'collect') |
  ((plot_pink_K  / plot_pink_L ) + plot_layout(axes = 'collect'))) + plot_layout(widths = c(1.7,1))

plot_pink_chum_w_title2 <- ((wrap_elements(panel = plot_chum2 + plot_annotation(tag_levels = list(c('A','B','C','D','I','J')), 
                                                                               title = "Chum salmon")& 
                                             theme(plot.tag.position = c(0.05, 1),
                                                   plot.tag = element_text(size = 10, hjust = 0, vjust = 0, face = "bold")))) / (wrap_elements(panel = plot_pink2+  
                                                                                                           plot_annotation(tag_levels = list(c('E','F','G','H','K','L')),
                                                                                                                           # tag_suffix = ')',
                                                                                                                           title = "Pink salmon")& 
                                                                                                           theme(plot.tag.position = c(0.05, 1),
                                                                                                                 plot.tag = element_text(size = 10, hjust = 0, vjust = 0, face = "bold"))))) 
plot_pink_chum_w_title2

##save

ggsave(here('figures','manuscript_dec2025_chum_pink_ricker_only_w_hd_CI.png'),
       plot = plot_pink_chum_w_title2,
       width = 8,
       height = 12,
       units = 'in',
       dpi = 500)



#make same figure as above but with cu level estimates, instead of river level estimates - using b_for_cu, instead of b_for_rv


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
    
    productivity_median <- apply(productivity,2,median)
    
    productivity_median_df <- data.frame(CU = unique(cu_data$CU_name),
                                         productivity_50 = apply(productivity,2,median),
                                         productivity_25 = apply(productivity,2,quantile, probs = 0.25),
                                         productivity_75 = apply(productivity,2,quantile, probs = 0.75),
                                         productivity_025 = apply(productivity,2,quantile, probs = 0.025),
                                         productivity_975 = apply(productivity,2,quantile, probs = 0.975),
                                         productivity_025_hdi = apply(productivity,2, hdi, ci = 0.95)[[1]]$CI_low,
                                         productivity_975_hdi = apply(productivity,2, hdi, ci = 0.95)[[1]]$CI_high,
                                         productivity_25_hdi = apply(productivity,2, hdi, ci = 0.5)[[1]]$CI_low,
                                         productivity_75_hdi = apply(productivity,2, hdi, ci = 0.5)[[1]]$CI_high,
                                         forestry = cpd_cu$mean,
                                         CU_n = unique(cu_data$CU_n))
    
    full_productivity <- rbind(full_productivity, productivity_median_df)
    
    
  }
  
  
  return(full_productivity)
  
  
}

ric_chm_cpd_productivity_decline_cu <- productivity_decline_cu_df(ric_chm_cpd_ocean_covariates_logR, 
                                                                  effect = "cpd", species = "chum")    

# make forest plot of estimates of decline from highest to lowest CU

cu_forest_plot_compare <- ric_chm_cpd_productivity_decline_cu %>% 
  arrange(desc(productivity_50)) %>%
  mutate(CU2 = factor(CU, levels = CU)) %>% 
  ggplot(aes(x = CU, y = productivity_50)) +
  geom_point(aes(y = productivity_50, x = CU2), color = '#516479',fill = "white", size = 3, alpha = 0.5) +
  # geom_errorbar(aes(ymin = productivity_025, ymax = productivity_975, color = "95% ET\ncredible interval"), width = 0, alpha = 0.5, size = 1) +
  geom_errorbar(aes(ymin = productivity_025_hdi, ymax = productivity_975_hdi, color = "95% credible\ninterval"), width = 0, alpha = 0.5, size = 1) +
  geom_errorbar(aes(ymin = productivity_25_hdi, ymax = productivity_75_hdi, color = "50% credible\ninterval"), width = 0, alpha = 0.7, size = 2) +
  # geom_errorbar(aes(ymin = productivity_25, ymax = productivity_75, color = '50% ET\ncredible interval'), width = 0, alpha = 0.7, size = 2) +
  scale_color_manual(name = '', values = c('95% credible\ninterval' = '#516479', 
                                           '50% credible\ninterval' = '#516479'))+
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
  labs(#title = 'Estimated percent change in CU-level productivity',
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

ggsave(here("figures","manuscript_dec2025_chum_ricker_cpd_productivity_decline_by_cu_forest_plot_hd_intervals.png"),
       cu_forest_plot_compare, width = 5, height = 6, bg = "white")


#make table of cu level effect sizes of 1) cumulative disturbance, 2) SST, and 3) NPGO

effect_sizes_cu_df <- function(posterior, effect, species){
  
  if(species == "chum"){
    df <- ch20rsc 
    
  } else if(species == "pink"){
    df <- pk10r
  }
  
  effect_df <- NULL
  
  for (i in 1:length(unique(df$CU_n))){
    
    cu <- unique(df$CU_n)[i]
    
    cu_data <- df %>% filter(CU_n == cu)
    
    if(effect == "cpd"){
      b_cu <- posterior %>% select(starts_with("b_for_cu")) %>%
        select(ends_with(paste0("[",cu,"]")))
      
    } else if(effect == "sst"){
      b_cu <- posterior %>% select(starts_with("b_sst_cu")) %>%
        select(ends_with(paste0("[",cu,"]")))
    } else if(effect == "npgo"){
      b_cu <- posterior %>% select(starts_with("b_npgo_cu")) %>%
        select(ends_with(paste0("[",cu,"]")))
    }
    
    effect_df_cu <- data.frame(CU = unique(cu_data$CU_name),
                               effect_median = round(apply(as.matrix(b_cu[,1]),2,median),2))
                               # sym(effect)_25 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.25),
                               # sym(effect)_75 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.75),
                               # sym(effect)_025 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.025),
                               # sym(effect)_975 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.975))
    
    effect_df <- rbind(effect_df, effect_df_cu)
  }
  
  return(effect_df)
}

effect_chum_cpd_sst_npgo <- effect_sizes_cu_df(ric_chm_cpd_ocean_covariates_logR, effect = "cpd", species = "chum") %>% 
  rename(cpd_effect_median = effect_median) %>%
  left_join(effect_sizes_cu_df(ric_chm_cpd_ocean_covariates_logR, effect = "sst", species = "chum") %>%
              rename(sst_effect_median = effect_median), by = "CU") %>%
  left_join(effect_sizes_cu_df(ric_chm_cpd_ocean_covariates_logR, effect = "npgo", species = "chum") %>%
              rename(npgo_effect_median = effect_median), by = "CU") 

#save table

write.csv(effect_chum_cpd_sst_npgo,
          here("tables","manuscript_nov2025_chum_ricker_cpd_sst_npgo_effect_sizes_by_cu_table.csv"),
          row.names = FALSE)

#do same for river level effects



effect_sizes_river_df <- function(posterior, effect, species){
  
  if(species == "chum"){
    df <- ch20rsc 
    
  } else if(species == "pink"){
    df <- pk10r
  }
  
  effect_df <- NULL
  
  for (i in 1:length(unique(df$River_n))){
    
    river <- unique(df$River_n)[i]
    
    river_data <- df %>% filter(River_n == river)
    
    if(effect == "cpd" || effect == "eca"){
      b_rv <- posterior %>% select(starts_with("b_for_rv")) %>%
        select(ends_with(paste0("[",river,"]")))
      
    } else if(effect == "sst"){
      b_rv <- posterior %>% select(starts_with("b_sst_rv")) %>%
        select(ends_with(paste0("[",river,"]")))
    } else if(effect == "npgo"){
      b_rv <- posterior %>% select(starts_with("b_npgo_rv")) %>%
        select(ends_with(paste0("[",river,"]")))
    }
    
    effect_df_rv <- data.frame(River = unique(river_data$River),
                               effect_median = round(apply(as.matrix(b_rv[,1]),2,median),2))
    # sym(effect)_25 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.25),
    # sym(effect)_75 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.75),
    # sym(effect)_025 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.025),
    # sym(effect)_975 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.975))
    
    effect_df <- rbind(effect_df, effect_df_rv)
  }
  
  return(effect_df)
}

effect_chum_cpd_sst_npgo_rv <- effect_sizes_river_df(ric_chm_cpd_ocean_covariates_logR, effect = "cpd", species = "chum") %>% 
  rename(cpd_effect_median = effect_median) %>%
  left_join(effect_sizes_river_df(ric_chm_cpd_ocean_covariates_logR, effect = "sst", species = "chum") %>%
              rename(sst_effect_median = effect_median), by = "River") %>%
  left_join(effect_sizes_river_df(ric_chm_cpd_ocean_covariates_logR, effect = "npgo", species = "chum") %>%
              rename(npgo_effect_median = effect_median), by = "River") 

#calculate the proportion of rivers in which the magnitude of effect of cpd is greater 
#than the magnitude of effect of sst and the magnitude of effect of npgo

sum((abs(effect_chum_cpd_sst_npgo_rv$cpd_effect_median)>abs(effect_chum_cpd_sst_npgo_rv$sst_effect_median))&
      (abs(effect_chum_cpd_sst_npgo_rv$cpd_effect_median)>abs(effect_chum_cpd_sst_npgo_rv$npgo_effect_median)))


effect_chum_eca_sst_npgo_rv <- effect_sizes_river_df(ric_chm_eca_ocean_covariates_logR, effect = "eca", species = "chum") %>% 
  rename(eca_effect_median = effect_median) %>%
  left_join(effect_sizes_river_df(ric_chm_eca_ocean_covariates_logR, effect = "sst", species = "chum") %>%
              rename(sst_effect_median = effect_median), by = "River") %>%
  left_join(effect_sizes_river_df(ric_chm_eca_ocean_covariates_logR, effect = "npgo", species = "chum") %>%
              rename(npgo_effect_median = effect_median), by = "River") 

#calculate the proportion of rivers in which the magnitude of effect of eca is greater 
#than the magnitude of effect of sst and the magnitude of effect of npgo

sum((abs(effect_chum_eca_sst_npgo_rv$eca_effect_median)>abs(effect_chum_eca_sst_npgo_rv$sst_effect_median))&
      (abs(effect_chum_eca_sst_npgo_rv$eca_effect_median)>abs(effect_chum_eca_sst_npgo_rv$npgo_effect_median)))/length(unique(ch20rsc$River))


range(effect_chum_cpd_sst_npgo_rv$cpd_effect_median/effect_chum_cpd_sst_npgo_rv$sst_effect_median)

# calculate the coastwide effect of forestry and ocean conditions.

median(ric_chm_cpd_ocean_covariates_logR$b_for)/median(ric_chm_cpd_ocean_covariates_logR$b_sst)
median(ric_chm_cpd_ocean_covariates_logR$b_for)/median(ric_chm_cpd_ocean_covariates_logR$b_npgo)
median(ric_chm_eca_ocean_covariates_logR$b_for)/median(ric_chm_eca_ocean_covariates_logR$b_sst)
median(ric_chm_eca_ocean_covariates_logR$b_for)/median(ric_chm_eca_ocean_covariates_logR$b_npgo)


#theoretical figures of decline

ricker_alpha <- function(St, alpha, beta, Smax, Forestry){
  log_Rt_St = alpha - St/Smax + beta*Forestry
  return(cbind(exp(log_Rt_St)*St, log_Rt_St))
}


St <- seq(0,1000,10)

alpha <- 1.2

beta <- -0.5

Smax <- 500

Forestry <- c(-1,0,1)

df_ricker <- data.frame(St = St)

for (i in 1:length(Forestry)){
  df_ricker[[paste0("Rt_",Forestry[i])]] <- ricker_alpha(St, alpha, beta, Smax, Forestry[i])[,1]
  df_ricker[[paste0("ln_Rt_St_",Forestry[i])]] <- ricker_alpha(St, alpha, beta, Smax, Forestry[i])[,2]
}

df_long_ricker <- df_ricker %>% 
  pivot_longer(cols = starts_with("Rt_"), names_to = "Forestry", values_to = "Rt", names_prefix = "Rt_") %>% 
  mutate(Forestry = as.numeric(Forestry)) #%>% 
  # pivot_longer(cols = starts_with("ln_"), names_to = "Forestry2", values_to = "ln_Rt_St", names_prefix = "ln_Rt_St_") %>% 
  # mutate(Forestry2 = as.numeric(Forestry2))

#calculate (Rt(forestry = 1)/Rt(forestry = -1) -1)*100

df_ricker_percent_change <- df_ricker %>% 
  rename("Rt_low_forestry" = "Rt_-1",
         "Rt_mid_forestry" = "Rt_0",
         "Rt_high_forestry" = "Rt_1",
         "ln_Rt_low_forestry" = "ln_Rt_St_-1",
         "ln_Rt_mid_forestry" = "ln_Rt_St_0",
         "ln_Rt_high_forestry" = "ln_Rt_St_1") %>%
  mutate(percent_change_recruit = ((Rt_high_forestry/Rt_low_forestry)-1)*100,
         percent_productivity = (exp(ln_Rt_high_forestry)/exp(ln_Rt_low_forestry)-1)*100)


# plot decline for different values of St
ggplot(df_ricker_percent_change, aes(x = St)) +
  geom_line(aes(y = percent_change_recruit), color = 'darkgreen', size = 1) +
  labs(title = "Theoretical Ricker model decline with forestry",
       x = "Spawning stock (St)",
       y = "Percent change in recruitment\n(high forestry vs low forestry)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))


ggplot(df_ricker_percent_change, aes(x = St)) +
  geom_line(aes(y = percent_productivity), color = 'darkgreen', size = 1) +
  labs(title = "Theoretical Ricker model decline with forestry",
       x = "Spawning stock (St)",
       y = "Percent change in productivty\n(high forestry vs low forestry)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))




plot_productivity_decline_manuscript(ric_chm_cpd_ocean_covariates_logR, 
                                     effect = "cpd", species = "chum",
                                     model = "Ricker",
                                     by_river = TRUE)

recruitment_decline_river_df <- function(posterior, effect, species){
  
  if(species == "chum"){
    df <- ch20rsc 
    
  } else if(species == "pink"){
    df <- pk10r
  }
  
  full_productivity <- NULL
  
  for (i in 1:length(unique(df$River_n))){
    
    river <- unique(df$River_n)[i]
    
    river_data <- df %>% filter(River_n == river)
    
    b_rv <- posterior %>% select(starts_with("b_for_rv")) %>%
      select(ends_with(paste0("[",river,"]")))
    
    eca_sqrt_std_river <- max(river_data$sqrt.ECA.std)
    #minimum forestry possible -0
    eca_river <- max(river_data$ECA_age_proxy_forested_only)
    
    forestry_eca <- seq(0,1, length.out = 100)
    
    cpd_sqrt_std_river <- max(river_data$sqrt.CPD.std)
    #minimum forestry possible - 0
    cpd_river <- max(river_data$disturbedarea_prct_cs)
    
    forestry_cpd <- seq(0,100, length.out = 100)
    
    #to calculate no forestry in the standardized scale
    forestry_sqrt <- sqrt(forestry_cpd)
    
    forestry_sqrt_std = (forestry_sqrt-mean(forestry_sqrt))/sd(forestry_sqrt)
    
    no_forestry <- min(forestry_sqrt_std)
    
    productivity <- (exp(as.matrix(b_rv[,1])%*%
                           (cpd_sqrt_std_river-no_forestry)))*100 - 100
    
    productivity_median <- apply(productivity,2,median)
    
    productivity_median_df <- data.frame(River = unique(river_data$River),
                                         productivity_50 = apply(productivity,2,median),
                                         productivity_25 = apply(productivity,2,quantile, probs = 0.25),
                                         productivity_75 = apply(productivity,2,quantile, probs = 0.75),
                                         productivity_025 = apply(productivity,2,quantile, probs = 0.025),
                                         productivity_975 = apply(productivity,2,quantile, probs = 0.975),
                                         productivity_025_hdi = apply(productivity,2, HDInterval::hdi, credMass = 0.95)[1,],
                                         productivity_975_hdi = apply(productivity,2, HDInterval::hdi, credMass = 0.95)[2,],
                                         productivity_25_hdi = apply(productivity,2, HDInterval::hdi, credMass = 0.5)[1,],
                                         productivity_75_hdi = apply(productivity,2, HDInterval::hdi, credMass = 0.5)[2,],
                          
                                         forestry = cpd_river,
                                         CU = unique(river_data$CU_name)
    )
    
    full_productivity <- rbind(full_productivity, productivity_median_df)
    
    
  }
  
  
  return(full_productivity)
  
  
  
  
  
  
  
}

ric_chm_cpd_recruitment_decline_river <- recruitment_decline_river_df(ric_chm_cpd_ocean_covariates_logR_long_chain, 
                                                                        effect = "cpd", species = "chum")


important_rivers <- c("Nimpkish River", "Skeena River", "Fraser River", "Capilano River",
                      "Squamish River", "Cheakamus River", "Pitt River", "Alouette River",
                      "Chilliwack River", "Vedder River", "Cowichan River", "Koksilah River",
                      "Goldstream River", "Campbell River", "Qualicum River", "Kingcome River", "Shoal Harbour Creek")

casestudy_watersheds <- c("Carnation Creek", "Viner Sound Creek", 
                          "Neekas Creek", "Deena Creek", "Phillips River")



# make forest plot of estimates of decline from highest to lowest rivers
ric_chm_cpd_recruitment_decline_river %>% 
  arrange(desc(productivity_50)) %>%
  mutate(River = str_to_title(River)) %>% 
  #change River to title case and then compare to important rivers
  mutate(important = ifelse((River %in% important_rivers | River %in% casestudy_watersheds), "yes", "no")) %>%
  mutate(River2 = factor(River, levels = River)) %>% 
  ggplot(aes(x = River, y = productivity_50)) +
  geom_point(aes(y = productivity_50, x = River2), color = '#516479',fill = "white", size = 1, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_025_hdi, ymax = productivity_975_hdi ), color = '#516479', width = 0, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_25_hdi, ymax = productivity_75_hdi ), color = '#516479', width = 0, alpha = 0.7) +
  geom_text_repel(color = "gray20", aes(label = ifelse(important == "yes", paste(River,CU, sep = ", "), NA)), 
                  size = 3, max.overlaps = 20,
                  direction    = "y", 
                  box.padding = 0.3, hjust = -1.5) + 
  coord_flip() +
  # scale_color_manual(name = 'Model type', values = c('independent alpha' = 'cadetblue', 'hierarchical alpha' = 'coral', 'hierarchical alpha - ricker' = 'darkgoldenrod')) +
  labs(#title = 'Estimated percent change in river-level productivity',
    x = 'River',
    y = 'Change in recruitment (%)') +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size = 4),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

#save 
ggsave(here("figures","manuscript_dec2025_chum_ricker_cpd_recruitment_decline_by_river_forest_plot.png"), width = 8, height = 10)


