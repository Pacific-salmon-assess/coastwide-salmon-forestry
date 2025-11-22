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


bh_chm_eca_ocean_covariates=read.csv(here('stan models','outs','posterior',
                                          'bh_chm_eca_ocean_covariates.csv'),check.names=F)
bh_chm_cpd_ocean_covariates=read.csv(here('stan models','outs','posterior',
                                          'bh_chm_cpd_ocean_covariates.csv'),check.names=F)


ric_chm_eca_ocean_covariates_logR=read.csv(here('stan models','outs','posterior',
                                           'ric_chm_eca_ocean_covariates_logR.csv'),check.names=F)
ric_chm_cpd_ocean_covariates_logR=read.csv(here('stan models','outs','posterior',
                                           'ric_chm_cpd_ocean_covariates_logR.csv'),check.names=F)

(((plot_forestry_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                 effect = "cpd", species = "chum", 
                                 model = "Ricker model with CPD", xlim = c(-1.3, 1.3))/
    plot_forestry_effect_manuscript(posterior = ric_chm_eca_ocean_covariates_logR, 
                                    effect = "eca", species = "chum", 
                                    model = "Ricker model with ECA", xlim = c(-1.3, 1.3)))+plot_layout(axes = 'collect')|
    (plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                         effect = "cpd", model = "Ricker model with CPD")+
    plot_productivity_decline_manuscript(posterior = ric_chm_eca_ocean_covariates_logR, 
                                         effect = "eca", model = "Ricker model with ECA")) + plot_layout(axes = 'collect'))/
    (plot_sst_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                               effect = "sst", species = "chum", 
                               model = "Ricker model with CPD", 
                               xlim = c(-1.3, 1.3))/
    plot_npgo_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                effect = "npgo", species = "chum", 
                                model = "Ricker model with CPD", 
                                xlim = c(-1.3, 1.3))+plot_layout(axes = 'collect')|
      (plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                            effect = "sst", model = "Ricker model with CPD")+
         plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                              effect = "npgo", model = "Ricker model with CPD")) + plot_layout(axes = 'collect'))/
    (plot_forestry_effect_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                    effect = "cpd", species = "chum", 
                                    model = "BH model with CPD", xlim = c(-1.3, 1.3))/
    plot_forestry_effect_manuscript(posterior = ch_chm_eca_ocean_covariates, 
                                    effect = "eca", species = "chum", 
                                    model = "BH model with ECA", xlim = c(-1.3, 1.3))+plot_layout(axes = 'collect')|
      (plot_productivity_decline_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                            effect = "cpd", model = "BH model with CPD")+
         plot_productivity_decline_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                              effect = "eca", model = "BH model with CPD")) + plot_layout(axes = 'collect'))/
    (plot_sst_effect_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                               effect = "sst", species = "chum", 
                               model = "BH model with CPD", 
                               xlim = c(-1.3, 1.3))/
    plot_npgo_effect_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                effect = "npgo", species = "chum", 
                                model = "BH model with CPD", 
                                xlim = c(-1.3, 1.3)) + plot_layout(axes = 'collect')|
      (plot_productivity_decline_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                            effect = "sst", model = "BH model with CPD")+
         plot_productivity_decline_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                              effect = "npgo", model = "BH model with CPD")) + plot_layout(axes = 'collect'))+
    plot_layout(axes = 'collect')
  
    ) 
      
#save figure

ggsave(here("figures","manuscript_aug2025_chum_results.png"), width = 10, height = 12)


(plot_forestry_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                 effect = "cpd", species = "chum", 
                                 model = "Ricker model with CPD", xlim = c(-1.3, 1.3))/
    plot_forestry_effect_manuscript(posterior = ric_chm_eca_ocean_covariates_logR, 
                                    effect = "eca", species = "chum", 
                                    model = "Ricker model with ECA", xlim = c(-1.3, 1.3))/
    plot_sst_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                               effect = "sst", species = "chum", 
                               model = "Ricker model with CPD", 
                               xlim = c(-1.3, 1.3))/
    plot_npgo_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                effect = "npgo", species = "chum", 
                                model = "Ricker model with CPD", 
                                xlim = c(-1.3, 1.3))/
    plot_forestry_effect_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                    effect = "cpd", species = "chum", 
                                    model = "BH model with CPD", xlim = c(-1.3, 1.3))/
    plot_forestry_effect_manuscript(posterior = ch_chm_eca_ocean_covariates, 
                                    effect = "eca", species = "chum", 
                                    model = "BH model with ECA", xlim = c(-1.3, 1.3))/
    plot_sst_effect_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                               effect = "sst", species = "chum", 
                               model = "BH model with CPD", 
                               xlim = c(-1.3, 1.3))/
    plot_npgo_effect_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                effect = "npgo", species = "chum", 
                                model = "BH model with CPD", 
                              xlim = c(-1.3, 1.3))+ plot_layout(axes = 'collect') |
  (((plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                          species = "chum",odel = "Ricker model with CPD")+
     plot_productivity_decline_manuscript(posterior = ric_chm_eca_ocean_covariates_logR, 
                                          effect = "eca", model = "Ricker model with ECA")) + plot_layout(axes = 'collect'))/
     ((plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                           effect = "sst", model = "Ricker model with CPD")+
        plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                             effect = "npgo", model = "Ricker model with CPD")) + plot_layout(axes = 'collect'))/
     ((plot_productivity_decline_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                            effect = "cpd", model = "BH model with CPD")+
         plot_productivity_decline_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                              effect = "eca", model = "BH model with CPD")) + plot_layout(axes = 'collect'))/
     ((plot_productivity_decline_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                           effect = "sst", model = "BH model with CPD")+
        plot_productivity_decline_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                             effect = "npgo", model = "BH model with CPD")) + plot_layout(axes = 'collect'))
     
     
  )
)

#save figure

ggsave(here("figures","manuscript_aug2025_chum_results.png"), width = 10, height = 12)

ric_pk_eca_ersst = read.csv(here('stan models',
                                 'outs',
                                 'posterior',
                                 'ric_pk_eca_st_noac_ocean_covariates_logR.csv'),check.names=F)

ric_pk_cpd_ersst = read.csv(here('stan models',
                                 'outs',
                                 'posterior',
                                 'ric_pk_cpd_st_noac_ocean_covariates_logR.csv'),check.names=F)

bh_pk_cpd_ersst = read.csv(here('stan models',
                                'outs',
                                'posterior',
                                'bh_pk_cpd_noac_w_npgo_ersst.csv'),check.names=F)

bh_pk_eca_ersst = read.csv(here('stan models',
                                'outs',
                                'posterior',
                                'bh_pk_eca_noac_w_npgo_ersst.csv'),check.names=F)

# make same plot for pink salmon

(plot_forestry_effect_manuscript(posterior = ric_pk_cpd_ersst, 
                                 effect = "cpd", species = "pink", 
                                 model = "Ricker model with CPD", xlim = c(-1.3, 1.3))/
    plot_forestry_effect_manuscript(posterior = ric_pk_eca_ersst, 
                                    effect = "eca", species = "pink", 
                                    model = "Ricker model with ECA", xlim = c(-1.3, 1.3))/
    plot_sst_effect_manuscript(posterior = ric_pk_cpd_ersst,
                               effect = "sst", species = "pink", 
                               model = "Ricker model with CPD", 
                               xlim = c(-1.3, 1.3))/
    plot_npgo_effect_manuscript(posterior = ric_pk_cpd_ersst, 
                                effect = "npgo", species = "pink", 
                                model = "Ricker model with CPD", 
                                xlim = c(-1.3, 1.3))/
    plot_forestry_effect_manuscript(posterior = ch_pk_cpd_ersst, 
                                    effect = "cpd", species = "pink", 
                                    model = "BH model with CPD", xlim = c(-1.3, 1.3))/
    plot_forestry_effect_manuscript(posterior = ch_pk_eca_ersst,
                                    effect = "eca", species = "pink", 
                                    model = "BH model with ECA", xlim = c(-1.3, 1.3))/
    plot_sst_effect_manuscript(posterior = ch_pk_cpd_ersst,
                               effect = "sst", species = "pink", 
                               model = "BH model with CPD", 
                               xlim = c(-1.3, 1.3))/
    plot_npgo_effect_manuscript(posterior = ch_pk_cpd_ersst, 
                                effect = "npgo", species = "pink", 
                                model = "BH model with CPD", 
                                xlim = c(-1.3, 1.3))+ plot_layout(axes = 'collect'))|
  (((plot_productivity_decline_manuscript(posterior = ric_pk_cpd_ersst,  
                                          species = "pink",
                                        effect = "cpd", model = "Ricker model with CPD")+
     plot_productivity_decline_manuscript(posterior = ric_pk_eca_ersst,  
                                          species = "pink",
                                          effect = "eca", model = "Ricker model with ECA")) + plot_layout(axes = 'collect'))/
     ((plot_productivity_decline_manuscript(posterior = ric_pk_cpd_ersst,  
                                            species = "pink",
                                            effect = "sst", model = "Ricker model with CPD")+
         plot_productivity_decline_manuscript(posterior = ric_pk_cpd_ersst,  
                                              species = "pink",
                                              effect = "npgo", model = "Ricker model with CPD")) + plot_layout(axes = 'collect'))/
     ((plot_productivity_decline_manuscript(posterior = ch_pk_cpd_ersst,  
                                            species = "pink",
                                            effect = "cpd", model = "BH model with CPD")+
         plot_productivity_decline_manuscript(posterior = ch_pk_eca_ersst, 
                                              species = "pink", 
                                              effect = "eca", model = "BH model with ECA")) + plot_layout(axes = 'collect'))/
     ((plot_productivity_decline_manuscript(posterior = ch_pk_cpd_ersst,  
                                            species = "pink",
                                            effect = "sst", model = "BH model with CPD")+
         plot_productivity_decline_manuscript(posterior = ch_pk_cpd_ersst, 
                                              species = "pink", 
                                              effect = "npgo", model = "BH model with CPD")) + plot_layout(axes = 'collect'))
     
     
  

)


ggsave(here("figures","manuscript_aug2025_pink_results.png"), width = 10, height = 12)(here("figures","manuscript_aug2025_pink_results.png"), width = 10, height = 12)




# Figure with just Ricker models - Chum at the cottom and Pinks at the top


plot_chum_a1 <- plot_forestry_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                                effect = "cpd", species = "chum", 
                                                model = "Cumulative disturbance", xlim = c(-0.5, 0.5))

plot_chum_a2 <- plot_forestry_effect_manuscript(posterior = ric_chm_eca_ocean_covariates_logR, 
                                                effect = "eca", species = "chum", 
                                                model = "Equivalent clearcut area", xlim = c(-0.5, 0.5))

plot_chum_a3 <- plot_sst_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                                effect = "sst", species = "chum", 
                                                model = "Sea-surface temperature", 
                                                xlim = c(-0.5, 0.5))

plot_chum_a4 <- plot_npgo_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                                effect = "npgo", species = "chum", 
                                                model = "North Pacific Gyre Oscillation", 
                                                xlim = c(-0.5, 0.5))


plot_chum_c1 <- plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                                 effect = "cpd", model = "Cumulative disturbance")

plot_chum_c2 <- plot_productivity_decline_manuscript(posterior = ric_chm_eca_ocean_covariates_logR,
                                                 effect = "eca", model = "Equivalent clearcut area")

plot_chum_c3 <- plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR,
                                                 effect = "sst", model = "Sea-surface temperature")

plot_chum_c4 <- plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR,
                                                 effect = "npgo", model = "North Pacific Gyre Oscillation")



plot_chum <- (plot_chum_a1 / plot_chum_a2 / plot_chum_a3 / plot_chum_a4)+ plot_layout(axes = 'collect') | 
  ((plot_chum_c1  + plot_chum_c2 ) + plot_layout(axes = 'collect'))/(
  (plot_chum_c3  + plot_chum_c4 ) + plot_layout(axes = 'collect')) 

# plot_chum[[1]] <- plot_chum[[1]] + plot_layout(tag_level = 'new') 
# 
# 
# plot_chum + plot_annotation(tag_levels = c('a', '1'), title = "Chum salmon")& 
#   theme(plot.tag = element_text(size = 8))


plot_chum_w_title <- (plot_chum + plot_annotation(tag_levels = list(c('a1','a2', 'a3','a4','b5', 'b6','b7','b8')), title = "Chum salmon")& 
  theme(plot.tag = element_text(size = 8)))


# Pink 

plot_pink_a1 <- plot_forestry_effect_manuscript(posterior = ric_pk_cpd_ersst, 
                                                effect = "cpd", species = "pink", 
                                                # model = "Ricker model with CPD", xlim = c(-1, 1))
                                                model = "Cumulative disturbance", xlim = c(-0.5, 0.5))

plot_pink_a2 <- plot_forestry_effect_manuscript(posterior = ric_pk_eca_ersst,
                                                effect = "eca", species = "pink", 
                                                # model = "Ricker model with ECA", xlim = c(-1, 1))
                                                model = "Equivalent clearcut area", xlim = c(-0.5, 0.5))

plot_pink_a3 <- plot_sst_effect_manuscript(posterior = ric_pk_cpd_ersst,
                                                effect = "sst", species = "pink", 
                                                # model = "Ricker model with CPD", 
                                                model = "Sea-surface temperature",
                                                # xlim = c(-1, 1))
                                                xlim = c(-0.5, 0.5))

plot_pink_a4 <- plot_npgo_effect_manuscript(posterior = ric_pk_cpd_ersst,
                                                effect = "npgo", species = "pink", 
                                                # model = "Ricker model with CPD",
                                                model = "North Pacific Gyre Oscillation",
                                                # xlim = c(-1, 1))
                                                xlim = c(-0.5, 0.5))


plot_pink_c1 <- plot_productivity_decline_manuscript(posterior = ric_pk_cpd_ersst,  
                                                 species = "pink",
                                                 effect = "cpd", 
                                                 # model = "Ricker model with CPD")
                                                 model = "Cumulative disturbance")

plot_pink_c2 <- plot_productivity_decline_manuscript(posterior = ric_pk_eca_ersst,
                                                 species = "pink",
                                                 effect = "eca", 
                                                 # model = "Ricker model with ECA")
                                                 model = "Equivalent clearcut area")

plot_pink_c3 <- plot_productivity_decline_manuscript(posterior = ric_pk_cpd_ersst,
                                                 species = "pink",
                                                 effect = "sst", 
                                                 # model = "Ricker model with CPD")
                                                 model = "Sea-surface temperature")

plot_pink_c4 <- plot_productivity_decline_manuscript(posterior = ric_pk_cpd_ersst,
                                                 species = "pink",
                                                 effect = "npgo", 
                                                 # model = "Ricker model with CPD")
                                                 model = "North Pacific Gyre Oscillation")






plot_pink <- (plot_pink_a1 / plot_pink_a2 / plot_pink_a3 / plot_pink_a4)+ plot_layout(axes = 'collect') |
  ((plot_pink_c1  + plot_pink_c2 ) + plot_layout(axes = 'collect'))/(
    (plot_pink_c3  + plot_pink_c4 ) + plot_layout(axes = 'collect'))
plot_pink + plot_annotation(tag_levels = c('a', '1'), title = "Pink salmon")& 
  theme(plot.tag = element_text(size = 8))
# plot_pink[[1]]  <- plot_pink[[1]]  + plot_layout(tag_level = 'new')
# 
# plot_pink[[2]]  <- plot_pink[[2]]  + plot_layout(tag_level = 'new')
# 
# 
# plot_pink + plot_annotation(tag_levels = c('a', '1'), title = "Pink salmon")& 
#   theme(plot.tag = element_text(size = 8))

plot_pink_w_title <- (plot_pink+ plot_annotation(tag_levels = list(c('c1','c2', 'c3','c4','d5', 'd6','d7','d8')), title = "Pink salmon")& 
  theme(plot.tag = element_text(size = 8)))


plot_pink_chum <- ((plot_chum + plot_annotation(title = "Chum salmon")) / (plot_pink + plot_annotation(title = "Pink salmon")) + 
                     plot_annotation(tag_levels = list(c('a1','a2', 'a3','a4','b5', 'b6','b7','b8','c1','c2', 'c3','c4','d5', 'd6','d7','d8')))& 
                     theme(plot.tag = element_text(size = 8)))

plot_pink_chum

# add title to chum plots and chum plots using wrap_elements(panel = p1 + ggtitle('Look at me shrink')

plot_pink_chum_w_title <- ((wrap_elements(panel = plot_chum + plot_annotation(tag_levels = list(c('A1','A2', 'A3','A4','B1', 'B2','B3','B4')), title = "Chum salmon")& 
                                            theme(plot.tag = element_text(size = 8)))) / (wrap_elements(panel = plot_pink+  plot_annotation(tag_levels = list(c('C1','C2', 'C3','C4','D1', 'D2','D3','D4')), title = "Pink salmon")& 
                                                                                                          theme(plot.tag = element_text(size = 8))))) 
plot_pink_chum_w_title

ggsave(here("figures","manuscript_oct2025_chum_pink_ricker_only.png"), plot = plot_pink_chum_w_title, width = 10, height = 12)

# forest plot of the predicted decline by CU and watershed

#new plot for manuscript without plot numbers, with decline in recruitment, with same y axis for % change






plot_productivity_decline_manuscript(ric_chm_cpd_ocean_covariates_logR, 
                          effect = "cpd", species = "chum",
                          model = "Ricker",
                          by_river = TRUE)

productivity_decline_river_df <- function(posterior, effect, species){
  
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
                                         forestry = cpd_river,
                                         CU = unique(river_data$CU_name)
                                         )
    
    full_productivity <- rbind(full_productivity, productivity_median_df)
    
  
  }


  return(full_productivity)
  
  
  
  
  
  
  
}

ric_chm_cpd_productivity_decline_river <- productivity_decline_river_df(ric_chm_cpd_ocean_covariates_logR, 
                                                                          effect = "cpd", species = "chum")

# look up river that has "SKEENA" in the name

ch20rsc %>% 
  filter(str_detect(CU_NAME, "SKEENA")) %>%
  select(River, CU_name) %>%
  distinct()


#make a list of important rivers in BC



important_rivers <- c("Nimpkish River", "Skeena River", "Fraser River", "Capilano River",
                      "Squamish River", "Cheakamus River", "Pitt River", "Alouette River",
                      "Chilliwack River", "Vedder River", "Cowichan River", "Koksilah River",
                      "Goldstream River", "Campbell River", "Quinsam River", "Qualicum River",
                      "Englishman River", "Big Qualicum River", "Nanaimo River", "Chemainus River",
                      "Crofton Creek", "Sproat River", "Somass River", "Nitinat River",
                      "Jordan River", "Cox Creek", "Carnation Creek")

#make a list of important rivers for first nations - include Kingcome, Viner, Shoal, Mackenzie, 
# Wakeman, Ahta, Kakweikan
#and 
# River       CU_name
# 1  ANDESITE CREEK  Lower Skeena
# 2      CLAY CREEK  Lower Skeena
# 3   DOG TAG CREEK  Lower Skeena
# 4   ECSTALL RIVER  Lower Skeena
# 5 GITNADOIX RIVER  Lower Skeena
# 6    KASIKS RIVER  Lower Skeena
# 7  KITWANGA RIVER Middle Skeena

important_rivers_first_nations <- c("Kingcome River", "Viner Sound Creek", "Shoal Harbour Creek", 
                                    "Mackenzie River","Wakeman River", "Ahta River", 
                                    "Kakweiken River", "Andesite Creek", "Clay Creek",
                                    "Dog Tag Creek", "Ecstall River", "Gitnadoix River",
                                    "Kasiks River", "Kitwanga River")


# make forest plot of estimates of decline from highest to lowest rivers
ric_chm_cpd_productivity_decline_river %>% 
  arrange(desc(productivity_50)) %>%
  mutate(River = str_to_title(River)) %>% 
  #change River to title case and then compare to important rivers
  mutate(important = ifelse((River %in% important_rivers | River %in% important_rivers_first_nations), "yes", "no")) %>%
  mutate(River2 = factor(River, levels = River)) %>% 
  ggplot(aes(x = River, y = productivity_50)) +
  geom_point(aes(y = productivity_50, x = River2), color = 'cadetblue',fill = "white", size = 1, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_025, ymax = productivity_975 ), color = 'cadetblue', width = 0.2, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_25, ymax = productivity_75 ), color = 'cadetblue', width = 0.4, alpha = 0.7) +
  geom_text_repel(color = "gray20", aes(label = ifelse(important == "yes", paste(River,CU, sep = ", "), NA)), 
                  size = 3, max.overlaps = 20,
                  direction    = "y", 
                  box.padding = 0.3, hjust = -1.5) + 
  coord_flip() +
  # scale_color_manual(name = 'Model type', values = c('independent alpha' = 'cadetblue', 'hierarchical alpha' = 'coral', 'hierarchical alpha - ricker' = 'darkgoldenrod')) +
  labs(#title = 'Estimated percent change in river-level productivity',
       x = 'River',
       y = 'Percent change in productivity') +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size = 4),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

#save 
ggsave(here("figures","manuscript_sep2025_chum_ricker_cpd_productivity_decline_by_river_forest_plot.png"), width = 8, height = 10)


#make same figure as above but with cu level estimates, instead of river level estimates - using b_for_cu, instead of b_for_rv
library(bayestestR)

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
                                         forestry = cpd_cu$mean,
                                         CU_n = unique(cu_data$CU_n))
    
    full_productivity <- rbind(full_productivity, productivity_median_df)
    
    
  }
  
  
  return(full_productivity)
  
  
}


ric_chm_cpd_productivity_decline_cu <- productivity_decline_cu_df(ric_chm_cpd_ocean_covariates_logR, 
                                                                    effect = "cpd", species = "chum")    

# make forest plot of estimates of decline from highest to lowest CU

ric_chm_cpd_productivity_decline_cu %>% 
  arrange(desc(productivity_50)) %>%
  mutate(CU2 = factor(CU, levels = CU)) %>% 
  ggplot(aes(x = CU, y = productivity_50)) +
  geom_point(aes(y = productivity_50, x = CU2), color = '#516479',fill = "white", size = 3, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_025, ymax = productivity_975 ), color = '#516479', width = 0, alpha = 0.5, size = 1) +
  geom_errorbar(aes(ymin = productivity_25, ymax = productivity_75 ), color = '#516479', width = 0, alpha = 0.7, size = 2) +
  #add estimated median decline to the right of each error bar
  geom_text(aes(label = paste(round(productivity_50,1),"%")), 
            hjust = -0.25, 
            vjust = -0.35,
            size = 3, color = "gray20") +
  coord_flip() +
  # scale_color_manual(name = 'Model type', values = c('independent alpha' = 'cadetblue', 'hierarchical alpha' = 'coral', 'hierarchical alpha - ricker' = 'darkgoldenrod')) +
  labs(#title = 'Estimated percent change in CU-level productivity',
       x = 'Conservation Unit',
       y = 'Percent change in productivity') +
  theme_classic() +
  theme(legend.position = "none",
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

#save

ggsave(here("figures","manuscript_oct2025_chum_ricker_cpd_productivity_decline_by_cu_forest_plot.png"), width = 6, height = 8)


#make a list of CUs ordered in ascending order of effect of b_for_cu

if(species == "chum"){
  df <- ch20rsc 
  
} else if(species == "pink"){
  df <- pk10r
}

full_b_cu <- NULL

for (i in 1:length(unique(df$CU_n))){
  
  cu <- unique(df$CU_n)[i]
  
  cu_data <- df %>% filter(CU_n == cu)
  
  b_cu <- posterior %>% select(starts_with("b_for_cu")) %>%
    select(ends_with(paste0("[",cu,"]"))) 
  
  b_cu_median_df <- data.frame(CU = unique(cu_data$CU_name),
                                       b_cu_50 = apply(b_cu,2,median),
                                       b_cu_25 = apply(b_cu,2,quantile, probs = 0.25),
                                       b_cu_75 = apply(b_cu,2,quantile, probs = 0.75),
                                       b_cu_025 = apply(b_cu,2,quantile, probs = 0.025),
                                       b_cu_975 = apply(b_cu,2,quantile, probs = 0.975),
                                       CU_n = unique(cu_data$CU_n))
  
  full_b_cu <- rbind(full_b_cu, b_cu_median_df)
  
}

full_b_cu_ascend <- full_b_cu %>%
  arrange(-b_cu_50)


#forest plot of effects of forestry at CU level

full_b_cu_ascend %>% 
  mutate(CU2 = factor(CU, levels = CU)) %>% 
  ggplot(aes(x = CU, y = b_cu_50)) +
  geom_point(aes(y = b_cu_50, x = CU2), color = 'cadetblue',fill = "white", size = 3, alpha = 0.5) +
  geom_errorbar(aes(ymin = b_cu_025, ymax = b_cu_975 ), color = 'cadetblue', width = 0, alpha = 0.5, size = 1) +
  geom_errorbar(aes(ymin = b_cu_25, ymax = b_cu_75 ), color = 'cadetblue', width = 0, alpha = 0.7, size = 2) +
  #add estimated median decline to the right of each error bar
  geom_text(aes(label = round(b_cu_50,3)), 
            hjust = -2,
            vjust = -0.25,
            size = 3, color = "gray20") +
  coord_flip() +
  # scale_color_manual(name = 'Model type', values = c('independent alpha' = 'cadetblue', 'hierarchical alpha' = 'coral', 'hierarchical alpha - ricker' = 'darkgoldenrod')) +
  labs(#title = 'Estimated effect of forestry on CU-level b_cu',
       x = 'Conservation Unit',
       y = 'Effect of forestry on productivity') +
  theme_classic() +
  theme(legend.position = "none",
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

#save

ggsave(here("figures","manuscript_sep2025_chum_ricker_cpd_forestry_effect_by_cu_forest_plot.png"), width = 6, height = 8)


# make forest plot of estimates of decline from highest to lowest rivers
# and keep all the river from a CU together
# arrange CU according to full_b_cu_ascend$CU

important_rivers <- c("Nimpkish River", "Skeena River", "Fraser River", "Capilano River",
                      "Squamish River", "Cheakamus River", "Pitt River", "Alouette River",
                      "Chilliwack River", "Vedder River", "Cowichan River", "Koksilah River",
                      "Goldstream River", "Campbell River", "Qualicum River", "Kingcome River", "Shoal Harbour Creek")

casestudy_watersheds <- c("Carnation Creek", "Viner Sound Creek", 
                         "Neekas Creek", "Deena Creek", "Phillips River")

ric_chm_cpd_productivity_decline_river_arrange <- ric_chm_cpd_productivity_decline_river %>% 
  arrange(desc(productivity_50)) %>%
  left_join(full_b_cu_ascend %>% select(CU, CU_n, b_cu_50) %>% distinct(), by = "CU") %>%
  arrange(desc(b_cu_50)) %>% 
  mutate(River = str_to_title(River)) %>% 
  #change River to title case and then compare to important rivers
  mutate(important = ifelse((River %in% important_rivers | River %in% casestudy_watersheds), "yes", "no")) %>%
  mutate(River2 = factor(River, levels = River))
  

ric_chm_cpd_productivity_decline_river_arrange %>%
  ggplot(aes(x = River, y = productivity_50)) +
  geom_point(aes(y = productivity_50, x = River2), color = 'cadetblue',fill = "white", size = 1, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_025, ymax = productivity_975 ), color = 'cadetblue', width = 0.2, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_25, ymax = productivity_75 ), color = 'cadetblue', width = 0.4, alpha = 0.7) +
  geom_text_repel(color = "gray20", aes(label = ifelse(important == "yes", paste(River,CU, sep = ", "), NA)), 
                  size = 3, max.overlaps = 20,
                  direction    = "y", 
                  box.padding = 0.3, hjust = -1.5) + 
  #add label of the CU with long bracket
  # geom_bracket(aes(xmin = 5.5, xmax = 8.5, y.position = 15, label = "Lower Skeena CU"), 
  #              label.size = 3, tip.length = 0.02, label.size.text = 3) +
  coord_flip() +
  # scale_color_manual(name = 'Model type', values = c('independent alpha' = 'cadetblue', 'hierarchical alpha' = 'coral', 'hierarchical alpha - ricker' = 'darkgoldenrod')) +
  labs(#title = 'Estimated percent change in river-level productivity',
    x = 'River',
    y = 'Percent change in productivity') +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size = 4),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

#save 
ggsave(here("figures","manuscript_sep2025_chum_ricker_cpd_productivity_decline_by_river_and_cu_forest_plot.png"), width = 8, height = 10)


# theoretical ricker curves for supplementary -----------------------------




#theoretical Ricker functions

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
  mutate(Forestry = as.numeric(Forestry)) %>% 
  pivot_longer(cols = starts_with("ln_"), names_to = "Forestry2", values_to = "ln_Rt_St", names_prefix = "ln_Rt_St_") %>% 
  mutate(Forestry2 = as.numeric(Forestry2))

p4_ric_1 <- df_long_ricker %>% 
  as.data.frame() %>% 
  ggplot(aes(x=St, y=Rt, group = Forestry, color = Forestry))+
  geom_line(size = 2, alpha = 0.5)+
  ylim(0, 1000)+
  scale_color_gradient2(name = 'Forestry',
                        low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 0) + 
  labs(#title = "Ricker Spawner Recruit Curve",
    x = "Spawners",
    y = "Recruits")+
  #annotate equation
  annotate("text", x =500, y = 200, label = TeX("$R = S e^{\\alpha - \\frac{S}{Smax} + \\beta_{forestry} Forestry$}"), parse = TRUE, size = 4)+
  
  theme_classic()


q4_ric_1 <- df_long_ricker %>%
  as.data.frame() %>% 
  select(ln_Rt_St, Forestry2, St) %>% 
  distinct() %>% 
  ggplot(aes(x=St, y=ln_Rt_St, group = Forestry2, color = Forestry2))+
  geom_line(size = 2, alpha = 0.5)+
  ylim(-2,2)+
  scale_color_gradient2(name = 'Forestry',
                        low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 0) + 
  labs(#title = "Ricker Spawner Recruit Curve",
    x = "Spawners",
    y = "log (Recruits/Spawners)")+
  #annotate equation
  annotate("text", x = 300, y = -1, label = TeX("$\\ln{\\frac{R}{S}} = \\alpha - \\frac{S}{Smax} + \\beta_{forestry} Forestry$"), parse = TRUE, size = 3)+
  theme_classic()

manuscript_ric_1 <- df_long_ricker %>% 
  # filter(Forestry == 0) %>%
  as.data.frame() %>%
  ggplot(aes(x=St, y=Rt, group = Forestry, color = Forestry))+
  geom_line(size = 2, alpha = 0.5)+
  ylim(0, 1100)+
  scale_color_gradient2(name = 'Forestry',
                        low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 0) + 
  labs(#title = "Ricker Spawner Recruit Curve",
    x = "Spawners",
    y = "Recruits")+
  
  # annotate("text", x =500, y = 200, label = TeX("$R = S e^{\\alpha - \\frac{S}{Smax} + \\beta^{forestry} Forestry$}"), parse = TRUE, size = 7)+
  geom_text(aes(label = "High Forestry", x = 500, y= 250), color = '#bf812d',
            size = 5) + 
  geom_text(aes(label = "Low Forestry", x = 500, y= 1100), color = '#35978f', size = 5) + 
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

manuscript_ric_1



ggsave(here("figures", "manuscript_oct2025_ricker_w_forestry.png"), 
       plot = manuscript_ric_1, 
       width = 8, height = 6, dpi = 300)



manuscript_linear_ric_1 <- df_long_ricker %>% 
  # filter(Forestry2 == 0) %>%
  as.data.frame() %>%
  ggplot(aes(x=St, y=ln_Rt_St, group = Forestry2, color = Forestry2))+
  geom_line(size = 2, alpha = 0.5)+
  ylim(-2,2)+
  scale_color_gradient2(name = 'Forestry',
                        low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 0) + 
  labs(#title = "Ricker Spawner Recruit Curve",
    x = "Spawners",
    y = TeX(r'($\ln{\frac{Recruits}{Spawners}}$)'))+
    # y = "ln (Recruits/Spawners)")+
  # annotate("text", x = 350, y = -1, label = TeX("$\\ln{\\frac{R}{S}} = \\alpha - \\frac{S}{Smax} + \\beta_{forestry} Forestry$"), parse = TRUE, size = 7)+
  geom_text(aes(label = "High Forestry", x = 200, y= -0.2), color = '#bf812d',
            size = 5) + 
  geom_text(aes(label = "Low Forestry", x = 200, y= 2), color = '#35978f', size = 5) + 
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

manuscript_linear_ric_1

ggsave(here("figures", "manuscript_oct2025_ricker_linear_w_forestry.png"),
       plot = manuscript_linear_ric_1, 
       width = 8, height = 6, dpi = 300)

manuscript_ric <- manuscript_ric_1 + manuscript_linear_ric_1

ggsave(here("figures", "manuscript_oct2025_ricker_both_w_forestry.png"),
       plot = manuscript_ric, 
       width = 12, height = 6, dpi = 300)


manuscript_logR_ric_1 <- df_ricker %>% 
  pivot_longer(cols = starts_with("Rt_"), names_to = "Forestry", values_to = "Rt", names_prefix = "Rt_") %>% 
  mutate(Forestry = as.numeric(Forestry)) %>% 
  # View()
  ggplot(aes(x=St, y=log(Rt), group = Forestry, color = Forestry))+
  geom_line(size = 2, alpha = 0.5)+
  # ylim(-2,2)+
  scale_color_gradient2(name = 'Forestry',
                        low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 0) + 
  labs(#title = "Ricker Spawner Recruit Curve",
    x = "Spawners",
    y = TeX(r'($\log {(Recruits)}$)'))+
  # y = "ln (Recruits/Spawners)")+
  # annotate("text", x = 350, y = -1, label = TeX("$\\ln{\\frac{R}{S}} = \\alpha - \\frac{S}{Smax} + \\beta_{forestry} Forestry$"), parse = TRUE, size = 7)+
  geom_text(aes(label = "High Forestry", x = 400, y= 5.25), color = '#bf812d',
            size = 5) + 
  geom_text(aes(label = "Low Forestry", x = 250, y= 7.25), color = '#35978f', size = 5) + 
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

manuscript_logR_ric_1

manuscript_ric3 <- manuscript_ric_1 + manuscript_linear_ric_1 + manuscript_logR_ric_1 + plot_annotation(tag_levels = 'A')

ggsave(here("figures", "manuscript_oct2025_ricker_all3_w_forestry.png"),
       plot = manuscript_ric3, 
       width = 12, height = 4, dpi = 300)

#make theoretical prediction of estimate change in recruitment (percent) for various level of forestry

# calculate

df_ricker_change_recruits <- df_ricker %>% 
  select("St", "Rt_1", "Rt_-1") %>%
  rename_with(~str_replace_all(., "-1", "low_forestry")) %>%
  rename_with(~str_replace_all(., "1", "high_forestry")) %>%
  mutate(recruits_change = (Rt_high_forestry - Rt_low_forestry)*100/Rt_low_forestry) %>% View()


St <- seq(0,1000,10)

alpha <- 1.2

beta <- -0.5

Smax <- 500

Forestry <- seq(-2, 2,100)#c(-1,0,1)

df_ricker <- data.frame(St = St)

for (i in 1:length(Forestry)){
  df_ricker[[paste0("Rt_",Forestry[i])]] <- ricker_alpha(St, alpha, beta, Smax, Forestry[i])[,1]
  df_ricker[[paste0("ln_Rt_St_",Forestry[i])]] <- ricker_alpha(St, alpha, beta, Smax, Forestry[i])[,2]
}




# make df for manuscript with productivity declines at different values --------

# also get median values of covaraite effect from posterior

ric_chm_eca_ocean_covariates_logR %>% 
  select("b_for") %>% 
  summarize(median = median(b_for),
            quantile_025 = quantile(b_for, probs = 0.025),
            quantile_975 = quantile(b_for, probs = 0.975))
# median quantile_025 quantile_975
# 1 -0.09862975   -0.2147934 -0.008607823

ric_pk_eca_ersst %>% 
  select("b_for") %>% 
  summarize(median = median(b_for),
            quantile_025 = quantile(b_for, probs = 0.025),
            quantile_975 = quantile(b_for, probs = 0.975))
# median quantile_025 quantile_975
# 1 0.007171355  -0.04269134   0.05549505

ric_chm_cpd_ocean_covariates_logR %>% 
  select("b_for") %>% 
  summarize(median = median(b_for),
            quantile_025 = quantile(b_for, probs = 0.025),
            quantile_975 = quantile(b_for, probs = 0.975))
# median quantile_025 quantile_975
# 1 -0.145166   -0.2423172  -0.06860073

ric_pk_cpd_ersst %>% 
  select("b_for") %>% 
  summarize(median = median(b_for),
            quantile_025 = quantile(b_for, probs = 0.025),
            quantile_975 = quantile(b_for, probs = 0.975))
# median quantile_025 quantile_975
# 1 -0.02796995  -0.07715119   0.02082885



ric_chm_cpd_ocean_covariates_logR %>% 
  select("b_sst") %>% 
  summarize(median = median(b_sst),
            quantile_025 = quantile(b_sst, probs = 0.025),
            quantile_975 = quantile(b_sst, probs = 0.975))
# median quantile_025 quantile_975
# 1 -0.02498645   -0.0833777   0.04199203

ric_chm_cpd_ocean_covariates_logR %>% 
  select("b_npgo") %>% 
  summarize(median = median(b_npgo),
            quantile_025 = quantile(b_npgo, probs = 0.025),
            quantile_975 = quantile(b_npgo, probs = 0.975))
# median quantile_025 quantile_975
# 1 0.0593236   0.01312529    0.1039678


# calculate the average (and max) covariate value and then calculate productivity
# decline for that value

productivity_decline_global_df <- function(posterior = ric_chm_cpd_ocean_covariates_logR, 
                              effect = "cpd", species = "chum"){
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
  
  return(global_df)
  
  
}

ric_chm_cpd_global_productivity_decline <- productivity_decline_global_df(ric_chm_cpd_ocean_covariates_logR, 
                                                                            effect = "cpd", species = "chum") %>% 
  # find the value closest to the mean and max disturbance
  filter((forestry - mean(ch20rsc$disturbedarea_prct_cs)) == min(abs(forestry - mean(ch20rsc$disturbedarea_prct_cs))) |
           (forestry - max(ch20rsc$disturbedarea_prct_cs)) == min(abs(forestry - max(ch20rsc$disturbedarea_prct_cs))))

ric_pk_cpd_global_productivity_decline <- productivity_decline_global_df(ric_pk_cpd_ersst, 
                                                                          effect = "cpd", species = "pink") %>% 
  # find the value closest to the mean and max disturbance
  filter(abs(forestry - mean(pk10r$disturbedarea_prct_cs)) == min(abs(forestry - mean(pk10r$disturbedarea_prct_cs))) |
           abs(forestry - max(pk10r$disturbedarea_prct_cs)) == min(abs(forestry - max(pk10r$disturbedarea_prct_cs))))

ric_chm_eca_global_productivity_decline <- productivity_decline_global_df(ric_chm_eca_ocean_covariates_logR, 
                                                                            effect = "eca", species = "chum") %>% 
  # find the value closest to the mean and max ECA
  filter(abs(forestry - mean(ch20rsc$ECA_age_proxy_forested_only)) == min(abs(forestry - mean(ch20rsc$ECA_age_proxy_forested_only))) |
           abs(forestry - max(ch20rsc$ECA_age_proxy_forested_only)) == min(abs(forestry - max(ch20rsc$ECA_age_proxy_forested_only))))

ric_pk_eca_global_productivity_decline <- productivity_decline_global_df(ric_pk_eca_ersst, 
                                                                          effect = "eca", species = "pink") %>% 
  # find the value closest to the mean and max ECA
  filter(abs(forestry - mean(pk10r$ECA_age_proxy_forested_only)) == min(abs(forestry - mean(pk10r$ECA_age_proxy_forested_only))) |
           abs(forestry - max(pk10r$ECA_age_proxy_forested_only)) == min(abs(forestry - max(pk10r$ECA_age_proxy_forested_only))))


ric_chm_sst_global_productivity_decline <- productivity_decline_global_df(ric_chm_cpd_ocean_covariates_logR, 
                                                                            effect = "sst", species = "chum") %>% 
  # find the value closest to the mean and max sst
  filter(abs(forestry - mean(ch20rsc$spring_ersst)) == min(abs(forestry - mean(ch20rsc$spring_ersst))) |
           abs(forestry - max(ch20rsc$spring_ersst)) == min(abs(forestry - max(ch20rsc$spring_ersst))))
ric_chm_sst_global_productivity_decline


ric_pk_sst_global_productivity_decline <- productivity_decline_global_df(ric_pk_cpd_ersst, 
                                                                          effect = "sst", species = "pink") %>% 
  # find the value closest to the mean and max sst
  filter(abs(forestry - mean(pk10r$spring_ersst)) == min(abs(forestry - mean(pk10r$spring_ersst))) |
           abs(forestry - max(pk10r$spring_ersst)) == min(abs(forestry - max(pk10r$spring_ersst))))
ric_pk_sst_global_productivity_decline


ric_chm_npgo_global_productivity_decline <- productivity_decline_global_df(ric_chm_cpd_ocean_covariates_logR, 
                                                                            effect = "npgo", species = "chum") %>% 
  # find the value closest to the mean and max npgo
  filter(abs(forestry - mean(ch20rsc$npgo)) == min(abs(forestry - mean(ch20rsc$npgo))) |
           abs(forestry - max(ch20rsc$npgo)) == min(abs(forestry - max(ch20rsc$npgo))))
ric_chm_npgo_global_productivity_decline


ric_pk_npgo_global_productivity_decline <- productivity_decline_global_df(ric_pk_cpd_ersst, 
                                                                          effect = "npgo", species = "pink") %>% 
  # find the value closest to the mean and max npgo
  filter(abs(forestry - mean(pk10r$npgo)) == min(abs(forestry - mean(pk10r$npgo))) |
           abs(forestry - max(pk10r$npgo)) == min(abs(forestry - max(pk10r$npgo))))
ric_pk_npgo_global_productivity_decline





# what is the mean ECA of dataset

mean(ch20rsc$ECA_age_proxy_forested_only)
# 0.085

max(ch20rsc$ECA_age_proxy_forested_only)
# 0.9216

# most recent mean ECA

ch20rsc %>% select(ECA_age_proxy_forested_only, River, BroodYear) %>% group_by(River) %>% filter(BroodYear == max(BroodYear)) %>% ungroup() %>% summarize(mean_max_ECA = mean(ECA_age_proxy_forested_only))
# 0.10


# mean CPD of dataset

mean(ch20rsc$disturbedarea_prct_cs)
# 19.77

# most recent mean CPD

ch20rsc %>% select(disturbedarea_prct_cs, River, BroodYear) %>% group_by(River) %>% filter(BroodYear == max(BroodYear)) %>% ungroup() %>% summarize(mean_max_CPD = mean(disturbedarea_prct_cs))
# 28.7

# mean SST

mean(ch20rsc$spring_ersst)
#10.94


# most recent mean sst

ch20rsc %>% select(spring_ersst, River, BroodYear) %>% group_by(River) %>% filter(BroodYear == max(BroodYear)) %>% ungroup() %>% summarize(mean_max_sst = mean(spring_ersst))
# 11.1


# mean NPGO

mean(ch20rsc$npgo)
#-0.012

# most recent mean npgo

ch20rsc %>% select(npgo, River, BroodYear) %>% group_by(River) %>% filter(BroodYear == max(BroodYear)) %>% ungroup() %>% summarize(mean_max_npgo = mean(npgo))


#0.356



# map of productivity decline ---------------------------------------------



# Required packages 
list.of.packages <- c("tidyverse", "extrafont", "giscoR", "ggplot2", "sf", "rmapshaper", "rnaturalearth", "patchwork", "here",
                      "bcmaps", "ggsflabel") 

# What packages need to be installed?
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])] 

## Install missing packages
if (length(new.packages)) {
  install.packages(new.packages, dependencies = TRUE, repos = c("http://cran.rstudio.com/", "http://R-Forge.R-project.org"))
}

## Loading libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))

# Load population sheds ----
# pop_sheds <- st_read(dsn = "data/processed/Max_ECA_sheds_Nov_2024.gpkg")  # Note that there are 1746 polygons, but only 1745 have VRI information in 2022.

pop_sheds <- st_read(dsn = here("origional-ecofish-data-models","Data","Spatial",
                                "bc_plots","Max_ECA_sheds_Nov_2024.gpkg"))  # Note that there are 1746 polygons, but only 1745 have VRI information in 2022.

# Get Canada 
world <- ne_countries(scale='medium',returnclass = 'sf')
Canada <- subset(world, admin == "Canada")
USA <- subset(world, admin == "United States of America")

# Watershed by region plot for forest  appendix ----
# Custom extent
b <- st_bbox(pop_sheds) 
bbox <- st_as_sfc(b)

# # Custom colors based on development districts
# colrs <- c("#e64b35", "#3c5388", "#2ca02c", "#8a4198", "#eea236", "#8f4c2d")
# names(colrs) <- c("North Coast - Skeena", 
#                   "South Coast", 
#                   "Haida Gwaii",
#                   "North Island - Central Coast",
#                   "Campbell River", 
#                   "South Island")


# Custom colors based on productivity change

# read posterior
posterior <- ric_chm_cpd_ocean_covariates_logR
lookup <- read.csv(here("origional-ecofish-data-models","Data","Spatial",
                        "bc_plots", "lookup.csv"))


# merge lookup and ch20rsc based on GFE_ID

lookup <- lookup %>% 
  left_join(ch20rsc %>% select(GFE_ID, River_n, BroodYear, sqrt.CPD.std, Species, disturbedarea_prct_cs) %>% 
              group_by(GFE_ID, River_n, Species) %>% #filter only max year
              filter(BroodYear == max(BroodYear)) %>% unique(), by = c("GFE_ID","Species"))

rivers <- lookup %>% select(River_n) %>%  filter(!is.na(River_n)) %>% pull(River_n) %>% unique()

no_forestry <- -2.745


productivity_df <- data.frame()

for(river in rivers){
  posterior_b_rv <- posterior %>%
    select(starts_with("b_for_rv")) %>% 
    # select(ends_with(paste0(".",river,".")))
    select(ends_with(paste0("[",river,"]"))) 
  
  sqrt.cpd.std <- lookup %>% filter(River_n == river) %>% pull(sqrt.CPD.std)
  cpd <- lookup %>% filter(River_n == river) %>% pull(disturbedarea_prct_cs)
  
  # posterior_alpha_rv <- posterior %>%
  #   select(starts_with("alpha_j")) %>% 
  #   select(ends_with(paste0(".",river,".")))
  
  productivity <- (exp(as.matrix(posterior_b_rv[,1])%*%
                         (sqrt.cpd.std-no_forestry)))*100 - 100
  
  productivity_median <- apply(productivity,2,median)
  
  productivity_df <- productivity_df %>% rbind(data.frame(River_n = river, sqrt.cpd.std = sqrt.cpd.std, cpd = cpd,
                                                          productivity_median = productivity_median))
  
}

# merge productivity_df with lookup and then pop_sheds

lookup_productivity_df <- lookup %>% 
  left_join(productivity_df, by = "River_n") %>% 
  mutate(LINEAR_FEATURE_ID = as.character(LINEAR_FEATURE_ID))

pop_sheds <- pop_sheds %>% 
  left_join(lookup_productivity_df %>% select(LINEAR_FEATURE_ID, River,productivity_median, sqrt.cpd.std, cpd), by = c("outlet_lfid"="LINEAR_FEATURE_ID"))



# Globe Inset
wld <- gisco_get_countries(resolution = "20")
ortho_crs <-'+proj=ortho +lat_0=45 +lon_0=-126 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs'

ocean <- st_point(x = c(0,0)) %>%
  st_buffer(dist = 6371000) %>% # radio Tierra
  st_sfc(crs = ortho_crs)

world <-   st_intersection(wld, st_transform(ocean, 4326)) %>%
  st_transform(crs = ortho_crs)

world_line <- ms_innerlines(world)

wld_map <- ggplot(world) +
  geom_sf(data = ocean, fill = "#deebf7", linewidth = .2) +
  geom_sf(fill = "grey50", 
          colour = NA,
          show.legend = F) +
  geom_sf(data = world_line, linewidth = .05, colour = "white") +
  geom_sf(data = bbox, col = "red", fill = NA, linewidth = 1) +
  scale_fill_manual(values = c("grey50", "red")) + 
  theme_void()


bc_boundary <- bc_bound() %>% st_transform(4326)


bc_region_productivity <- ggplot() +
  geom_sf(data = bc_boundary, fill = "gray90", color = alpha("slategray",0.1), linewidth = 0.5) +
  geom_sf(data = pop_sheds %>% 
            filter(!is.na(productivity_median)), 
          # aes(fill = productivity_median, color = cpd), linewidth = 0.5) +
          aes(fill = productivity_median, color = cpd), linewidth = 0.5, alpha = 0.8) +
  # geom_sf_label_repel(data = pop_sheds %>% 
  #                       filter(productivity_median < -60), aes(label = River), fill = alpha("gray90",0.2),
  #                     size = 2, color = "black", fontface = "bold", alpha = 0.5) +
  # Plot Window
  coord_sf(crs = st_crs(pop_sheds), xlim = c(b["xmin"], b["xmax"] -70000) , ylim = c(b["ymin"] - 4000, b["ymax"]-70000)) +
  # Theme
  # scale_fill_viridis_c() +
  scale_color_gradient2(name = 'Cumulative\ndisturbance (%)',
                        low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 50) +
  scale_fill_distiller(name = "Productivity\nchange (%)", 
                       type = "seq",
                       direction = -1,
                       palette = "Greys",
                       limits = c(-70, 0))+
  # scale_color_distiller(name = "Productivity\nchange (%)", 
  #                       type = "seq",
  #                       direction = -1,
  #                       palette = "Greys",
  #                       limits = c(-70, 0))+
  theme(panel.grid.major = element_line(colour = "aliceblue", linetype = "dashed", 
                                        size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"), 
        panel.border = element_rect(fill = NA, color = NA),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        base_family = "ArcherPro Book",
        legend.position = c(0.9, 0.5),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.font = "ArcherPro Book",
        legend.key.size = unit(0.8, "cm"),
        legend.key.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.background = element_rect(fill = "transparent", size = 0.5)
  ) +
  labs(x = "", y = "", color = "", fill = "") +
  # Spatial annotation
  ggspatial::annotation_scale(
    location = "bl",
    bar_cols = c("black", "white"),
    text_family = "ArcherPro Book"
  ) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    # pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "black",
      text_family = "ArcherPro Book"
    )
  )  + 
  inset_element(wld_map, left = 0.05, bottom = 0.1, right = 0.3, top = 0.35, align_to = "plot")

bc_region_productivity


ggsave(here("figures", "coastwide_productivity_change_map_cpd_border_hres_ricker_oct_2025.png"), bc_region_productivity, width = 10, height = 12, dpi = 1000)


# forest plot of ricers with CU label


ric_chm_cpd_productivity_decline_river_arrange %>%
  ggplot(aes(x = River, y = productivity_50)) +
  geom_point(aes(y = productivity_50, x = River2), color = 'cadetblue',fill = "white", size = 1, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_025, ymax = productivity_975 ), color = 'cadetblue', width = 0.2, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_25, ymax = productivity_75 ), color = 'cadetblue', width = 0.4, alpha = 0.7) +
  # facet_wrap(~factor(CU, levels = unique(ric_chm_cpd_productivity_decline_river_arrange$CU)),
  #            strip.position = "bottom", 
  #            nrow = 1) +
  geom_text_repel(color = "gray20", aes(label = ifelse(important == "yes", paste(River,CU, sep = ", "), NA)), 
                  size = 3, max.overlaps = 20,
                  direction    = "y", 
                  nudge_y      = 150,
                  box.padding = 0.3, hjust = 0) + 
  #add label of the CU
  coord_flip() +
  # scale_color_manual(name = 'Model type', values = c('independent alpha' = 'cadetblue', 'hierarchical alpha' = 'coral', 'hierarchical alpha - ricker' = 'darkgoldenrod')) +
  labs(#title = 'Estimated percent change in river-level productivity',
    x = 'River',
    y = 'Percent change in productivity') +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size = 4),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.placement = 'outside')

#save 
ggsave(here("figures","manuscript_sep2025_chum_ricker_cpd_productivity_decline_by_river_and_cu_forest_plot.png"), width = 8, height = 10)



ric_chm_cpd_productivity_decline_river_arrange %>%
  ggplot(aes(y = River, x = productivity_50)) +
  geom_point(aes(x = productivity_50, y = River2), color = 'cadetblue',fill = "white", size = 1, alpha = 0.5) +
  geom_errorbar(aes(xmin = productivity_025, xmax = productivity_975 ), color = 'cadetblue', width = 0.2, alpha = 0.5) +
  geom_errorbar(aes(xmin = productivity_25, xmax = productivity_75 ), color = 'cadetblue', width = 0.4, alpha = 0.7) + 
  facet_wrap(~factor(CU, levels = unique(ric_chm_cpd_productivity_decline_river_arrange$CU)),
             ncol = 1,
             strip.position = "right",
             scales = "free_y"
             )+
  labs(#title = 'Estimated percent change in river-level productivity',
    y = 'River',
    x = 'Percent change in productivity') +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size = 4),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.placement = 'outside',
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        strip.text.y.right = element_text(angle = 0)
        )




# make forest plot of estimates of decline from highest to lowest rivers
ric_chm_cpd_productivity_decline_river %>% 
  arrange(desc(productivity_50)) %>%
  mutate(River = str_to_title(River)) %>% 
  #change River to title case and then compare to important rivers
  mutate(important = ifelse((River %in% important_rivers | River %in% casestudy_watersheds), "yes", "no")) %>%
  mutate(River2 = factor(River, levels = River)) %>% 
  ggplot(aes(x = River, y = productivity_50)) +
  geom_point(aes(y = productivity_50, x = River2), color = '#516479',fill = "white", size = 1, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_025, ymax = productivity_975 ), color = '#516479', width = 0, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_25, ymax = productivity_75 ), color = '#516479', width = 0, alpha = 0.7) +
  geom_text_repel(color = "gray20", aes(label = ifelse(important == "yes", paste(River,CU, sep = ", "), NA)), 
                  size = 3, max.overlaps = 20,
                  direction    = "y", 
                  box.padding = 0.3, hjust = -1.5) + 
  coord_flip() +
  # scale_color_manual(name = 'Model type', values = c('independent alpha' = 'cadetblue', 'hierarchical alpha' = 'coral', 'hierarchical alpha - ricker' = 'darkgoldenrod')) +
  labs(#title = 'Estimated percent change in river-level productivity',
    x = 'River',
    y = 'Percent change in productivity') +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size = 4),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

#save 
ggsave(here("figures","manuscript_oct2025_chum_ricker_cpd_productivity_decline_by_river_forest_plot.png"), width = 8, height = 10)


# have map figure bc_region_productivity and forest plot with CU as inset

cu_forest_plot_compare <- ric_chm_cpd_productivity_decline_cu %>% 
  arrange(desc(productivity_50)) %>%
  mutate(CU2 = factor(CU, levels = CU)) %>% 
  ggplot(aes(x = CU, y = productivity_50)) +
  geom_point(aes(y = productivity_50, x = CU2), color = '#516479',fill = "white", size = 3, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_025, ymax = productivity_975, color = "95% ETI CI"), width = 0, alpha = 0.5, size = 1) +
  geom_errorbar(aes(ymin = productivity_025_hdi, ymax = productivity_975_hdi, color = "95% HDI CI"), width = 0, alpha = 0.5, size = 1) +
  geom_errorbar(aes(ymin = productivity_25, ymax = productivity_75 ), color = '#516479', width = 0, alpha = 0.7, size = 2) +
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

ggsave(here("figures","manuscript_oct2025_chum_ricker_cpd_productivity_decline_by_cu_forest_plot_compare_intervals.png"),
       cu_forest_plot_compare, width = 4, height = 6, bg = "white")


cu_forest_plot <- ric_chm_cpd_productivity_decline_cu %>% 
  arrange(desc(productivity_50)) %>%
  mutate(CU2 = factor(CU, levels = CU)) %>% 
  ggplot(aes(x = CU, y = productivity_50)) +
  geom_point(aes(y = productivity_50, x = CU2), color = '#516479',fill = "white", size = 3, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_025, ymax = productivity_975, color = "95% confidence\ninterval"), width = 0, alpha = 0.5, size = 1) +
  # geom_errorbar(aes(ymin = productivity_025_hdi, ymax = productivity_975_hdi, color = "95% HDI CI"), width = 0, alpha = 0.5, size = 1) +
  geom_errorbar(aes(ymin = productivity_25, ymax = productivity_75,color = '50% confidence\ninterval' ), width = 0, alpha = 0.7, size = 2) +
  scale_color_manual(name = '', values = c('95% confidence\ninterval' = '#516479', 
                                                        '50% confidence\ninterval' = '#516479', 
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

cu_forest_plot

ggsave(here("figures","manuscript_nov2025_chum_ricker_cpd_productivity_decline_by_cu_forest_plot.png"),
       cu_forest_plot, width = 5, height = 6, bg = "white")



bc_region_productivity2 <- ggplot() +
  geom_sf(data = bc_boundary, fill = "gray90", color = alpha("slategray",0.1), linewidth = 0.5) +
  geom_sf(data = pop_sheds %>% 
            filter(!is.na(productivity_median)), 
          # aes(fill = productivity_median, color = cpd), linewidth = 0.5) +
          aes(fill = productivity_median, color = cpd), linewidth = 0.5, alpha = 0.8) +
  # geom_sf_label_repel(data = pop_sheds %>% 
  #                       filter(productivity_median < -60), aes(label = River), fill = alpha("gray90",0.2),
  #                     size = 2, color = "black", fontface = "bold", alpha = 0.5) +
  # Plot Window
  coord_sf(crs = st_crs(pop_sheds), xlim = c(b["xmin"], b["xmax"] -70000) , ylim = c(b["ymin"] - 4000, b["ymax"]-70000)) +
  # Theme
  # scale_fill_viridis_c() +
  scale_color_gradient2(name = 'Cumulative\ndisturbance (%)',
                        low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 50) +
  scale_fill_distiller(name = "Productivity\nchange (%)", 
                       type = "seq",
                       direction = -1,
                       palette = "Greys",
                       limits = c(-70, 0))+
  # scale_color_distiller(name = "Productivity\nchange (%)", 
  #                       type = "seq",
  #                       direction = -1,
  #                       palette = "Greys",
  #                       limits = c(-70, 0))+
  theme(panel.grid.major = element_line(colour = "aliceblue", linetype = "dashed", 
                                        size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"), 
        panel.border = element_rect(fill = NA, color = NA),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        base_family = "ArcherPro Book",
        legend.position = c(0.9, 0.4),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.font = "ArcherPro Book",
        legend.key.size = unit(0.8, "cm"),
        legend.key.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.background = element_rect(fill = "transparent", size = 0.5)
  ) +
  labs(x = "", y = "", color = "", fill = "") +
  # Spatial annotation
  ggspatial::annotation_scale(
    location = "bl",
    bar_cols = c("black", "white"),
    text_family = "ArcherPro Book"
  ) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    # pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "black",
      text_family = "ArcherPro Book"
    )
  )  + 
  inset_element(wld_map, left = 0.05, bottom = 0.1, right = 0.3, top = 0.35, align_to = "plot")

bc_region_productivity2

bc_region_productivity2_w_inset <- bc_region_productivity2 + 
  inset_element(cu_forest_plot, left = 0.55, bottom = 0.6, right = 0.99, top = 0.99, align_to = "plot")


ggsave(here("figures", "coastwide_productivity_change_map_cpd_border_hres_ricker_oct_2025_w_inset.png"), bc_region_productivity2_w_inset, width = 10, height = 12, dpi = 1000)


# make table for manuscript
## table should have CU name as column 1, productivity change by Cumulative disturbance, 
## productivity change by SST, productivity change by NPGO

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
    
    if(effect == "cpd"){
      b_cu <- posterior %>% select(starts_with("b_for_cu")) %>%
        select(ends_with(paste0("[",cu,"]")))
      
      # cpd_sqrt_std_cu <- max(cu_data$sqrt.CPD.std) # should not be using max from CU
      #minimum forestry possible - 0
      
      covariate_cu <- cu_data %>% group_by(River) %>% 
        filter(disturbedarea_prct_cs == max(disturbedarea_prct_cs)) %>% 
        distinct(disturbedarea_prct_cs) %>% 
        ungroup %>% 
        summarize(mean = mean(disturbedarea_prct_cs))
      
      covariate_sqrt_std_cu <- cu_data %>% group_by(River) %>% 
        filter(sqrt.CPD.std == max(sqrt.CPD.std)) %>% 
        distinct(sqrt.CPD.std) %>% 
        ungroup %>% 
        summarize(mean = mean(sqrt.CPD.std))
      
      predicted_covariate <- seq(0,100, length.out = 100)
      
      #to calculate no forestry in the standardized scale
      predicted_covariate_sqrt <- sqrt(predicted_covariate)
      
      predicted_covariate_std = (predicted_covariate_sqrt-mean(predicted_covariate_sqrt))/sd(predicted_covariate_sqrt)
      
      min_covariate <- min(predicted_covariate_std)
      
      
    } else if(effect == "sst"){
      b_cu <- posterior %>% select(starts_with("b_sst_cu")) %>%
        select(ends_with(paste0("[",cu,"]")))
      
      covariate_cu <- cu_data %>% group_by(River) %>% 
        filter(spring_ersst == max(spring_ersst)) %>% 
        distinct(spring_ersst) %>% 
        ungroup %>% 
        summarize(mean = mean(spring_ersst))
      
      covariate_sqrt_std_cu <- cu_data %>% group_by(River) %>% 
        filter(sst.std  == max(sst.std)) %>% 
        distinct(sst.std) %>% 
        ungroup %>% 
        summarize(mean = mean(sst.std))
      
      
      predicted_covariate <- seq(min(cu_data$spring_ersst), max(cu_data$spring_ersst), length.out = 100)
      
      
      
      
      predicted_covariate_std = (predicted_covariate-mean(predicted_covariate))/sd(predicted_covariate)
      
      min_covariate <- min(predicted_covariate_std)
      
    } else if(effect == "npgo"){
      b_cu <- posterior %>% select(starts_with("b_npgo_cu")) %>%
        select(ends_with(paste0("[",cu,"]")))
      
      covariate_cu <- cu_data %>% group_by(River) %>% 
        filter(npgo == max(npgo)) %>% 
        distinct(npgo) %>% 
        ungroup %>% 
        summarize(mean = mean(npgo))
      
      covariate_sqrt_std_cu <- cu_data %>% group_by(River) %>% 
        filter(npgo.std == max(npgo.std)) %>% 
        distinct(npgo.std) %>% 
        ungroup %>% 
        summarize(mean = mean(npgo.std))
      
      predicted_covariate <- seq(min(cu_data$npgo), max(cu_data$npgo), length.out = 100)
      
      predicted_covariate_std <- (predicted_covariate-mean(predicted_covariate))/sd(predicted_covariate)
      
      min_covariate <- min(predicted_covariate_std)
      
    
      
    } else if(effect == "eca"){
      b_cu <- posterior %>% select(starts_with("b_for_cu")) %>%
        select(ends_with(paste0("[",cu,"]")))
      
      covariate_cu <- cu_data %>% group_by(River) %>%
        filter(ECA_age_proxy_forested_only == max(ECA_age_proxy_forested_only)) %>% 
        distinct(ECA_age_proxy_forested_only) %>% 
        ungroup %>% 
        summarize(mean = mean(ECA_age_proxy_forested_only))
      
      covariate_sqrt_std_cu <- cu_data %>% group_by(River) %>%
        filter(sqrt.ECA.std  == max(sqrt.ECA.std )) %>% 
        distinct(sqrt.ECA.std ) %>% 
        ungroup %>% 
        summarize(mean = mean(sqrt.ECA.std ))
      
      predicted_covariate <- seq(0,1, length.out = 100)
      
      predicted_covariate_sqrt = sqrt(predicted_covariate)
      
      predicted_covariate_std = (predicted_covariate_sqrt-mean(predicted_covariate_sqrt))/sd(predicted_covariate_sqrt)
      
      min_covariate <- min(forestry_sqrt_std)
      
      
    }
    
    
    productivity <- (exp(as.matrix(b_cu[,1])%*%
                           (covariate_sqrt_std_cu$mean-min_covariate)))*100 - 100
    
    productivity_median <- apply(productivity,2,median)
    
    productivity_median_df <- data.frame(CU = unique(cu_data$CU_name),
                                         productivity_50 = apply(productivity,2,median),
                                         productivity_25 = apply(productivity,2,quantile, probs = 0.25),
                                         productivity_75 = apply(productivity,2,quantile, probs = 0.75),
                                         productivity_025 = apply(productivity,2,quantile, probs = 0.025),
                                         productivity_975 = apply(productivity,2,quantile, probs = 0.975),
                                         covariate = covariate_cu$mean,
                                         CU_n = unique(cu_data$CU_n))
    
    full_productivity <- rbind(full_productivity, productivity_median_df)
    
    
  }
  
  
  return(full_productivity)
  
  
}


ric_chm_sst_productivity_decline_cu <- productivity_decline_cu_df(ric_chm_cpd_ocean_covariates_logR, 
                                                              effect = "sst", species = "chum")

ric_chm_npgo_productivity_decline_cu <- productivity_decline_cu_df(ric_chm_cpd_ocean_covariates_logR, 
                                                              effect = "npgo", species = "chum")

manuscript_productivity_change_table <- ric_chm_cpd_productivity_decline_cu %>% 
  select(CU, productivity_50, productivity_025, productivity_975) %>% 
  rename(productivity_change_cpd_50 = productivity_50,
         productivity_change_cpd_025 = productivity_025,
         productivity_change_cpd_975 = productivity_975) %>% 
  left_join(ric_chm_sst_productivity_decline_cu %>% 
              select(CU, productivity_50, productivity_025, productivity_975) %>% 
              rename(productivity_change_sst_50 = productivity_50,
                     productivity_change_sst_025 = productivity_025,
                     productivity_change_sst_975 = productivity_975),
            by = "CU") %>% 
  left_join(ric_chm_npgo_productivity_decline_cu %>% 
              select(CU, productivity_50, productivity_025, productivity_975) %>% 
              rename(productivity_change_npgo_50 = productivity_50,
                     productivity_change_npgo_025 = productivity_025,
                     productivity_change_npgo_975 = productivity_975),
            by = "CU")

manuscript_productivity_change_table %>% View()



#save csv

write.csv(manuscript_productivity_change_table, here("tables","manuscript_productivity_change_table_by_cu_oct_2025.csv"), row.names = FALSE)



#round all values and only select median columns

manuscript_productivity_change_table_rounded <- manuscript_productivity_change_table %>% 
  select(CU, productivity_change_cpd_50, productivity_change_sst_50, productivity_change_npgo_50) %>% 
  mutate(across(where(is.numeric), ~ round(., 1)))

write.csv(manuscript_productivity_change_table_rounded, here("tables","manuscript_productivity_change_table_by_cu_oct_2025_subset.csv"), row.names = FALSE)





#new plot for manuscript without plot numbers, with decline in recruitment, with same y axis for % change



plot_chum_A <- plot_forestry_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                                effect = "cpd", species = "chum", 
                                                model = "Cumulative disturbance", xlim = c(-0.5, 0.5))+
  theme(legend.position = c(0.9,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 6))
#add legend for credible intervals for one figure
plot_chum_B <- plot_forestry_effect_manuscript(posterior = ric_chm_eca_ocean_covariates_logR, 
                                                effect = "eca", species = "chum", 
                                                model = "Equivalent clearcut area", xlim = c(-0.5, 0.5))+
  theme(legend.position = "none")

plot_chum_C <- plot_sst_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                           effect = "sst", species = "chum", 
                                           model = "Sea-surface temperature", 
                                           xlim = c(-0.5, 0.5))+
  theme(legend.position = "none")

plot_chum_D <- plot_npgo_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                            effect = "npgo", species = "chum", 
                                            model = "North Pacific Gyre Oscillation", 
                                            xlim = c(-0.5, 0.5))+
  theme(legend.position = "none")


plot_chum_E <- plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                                     effect = "cpd", model = "Cumulative disturbance")+
  ylim(c(-75,75))+
  labs(y = "Median change in recruitment (%)")+
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
        legend.text = element_text(size = 7))#+
# guides(fill=guide_legend(nrow=1, byrow=TRUE))

plot_chum_F <- plot_productivity_decline_manuscript(posterior = ric_chm_eca_ocean_covariates_logR,
                                                     effect = "eca", model = "Equivalent clearcut area")+
  ylim(c(-75,75))+
  labs(y = "Median change in recruitment (%)")


plot_chum_G <- plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR,
                                                     effect = "sst", model = "Sea-surface temperature")+
  ylim(c(-75,75))+
  labs(y = "Median change in recruitment (%)")

plot_chum_H <- plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR,
                                                     effect = "npgo", model = "North Pacific Gyre Oscillation")+
  ylim(c(-75,75))+
  labs(y = "Median change in recruitment (%)")



# plot_chum <- (plot_chum_a1 / plot_chum_a2 / plot_chum_a3 / plot_chum_a4)+ plot_layout(axes = 'collect') | 
#   ((plot_chum_c1  + plot_chum_c2 ) + plot_layout(axes = 'collect'))/(
#     (plot_chum_c3  + plot_chum_c4 ) + plot_layout(axes = 'collect')) 

plot_chum <- (plot_chum_A / plot_chum_B / plot_chum_C / plot_chum_D)+ plot_layout(axes = 'collect') | 
  ((plot_chum_E  + plot_chum_F ) + plot_layout(axes = 'collect'))/(
    (plot_chum_G  + plot_chum_H ) + plot_layout(axes = 'collect'))

plot_chum


plot_chum2 <- (plot_chum_A / plot_chum_B / plot_chum_C / plot_chum_D)+ plot_layout(axes = 'collect') | 
  ((plot_chum_E  / plot_chum_F ) + plot_layout(axes = 'collect'))
# plot_chum[[1]] <- plot_chum[[1]] + plot_layout(tag_level = 'new') 
# 
# 
# plot_chum + plot_annotation(tag_levels = c('a', '1'), title = "Chum salmon")& 
#   theme(plot.tag = element_text(size = 8))


plot_chum_w_title <- (plot_chum + plot_annotation(tag_levels = "A",
                                                  title = "Chum salmon")& 
                        theme(plot.tag = element_text(size = 8)))


# Pink 

plot_pink_I <- plot_forestry_effect_manuscript(posterior = ric_pk_cpd_ersst, 
                                                effect = "cpd", species = "pink", 
                                                # model = "Ricker model with CPD", xlim = c(-1, 1))
                                                model = "Cumulative disturbance", xlim = c(-0.5, 0.5))+
  theme(legend.position = "none")

plot_pink_J <- plot_forestry_effect_manuscript(posterior = ric_pk_eca_ersst,
                                                effect = "eca", species = "pink", 
                                                # model = "Ricker model with ECA", xlim = c(-1, 1))
                                                model = "Equivalent clearcut area", xlim = c(-0.5, 0.5))+
  theme(legend.position = "none")

plot_pink_K <- plot_sst_effect_manuscript(posterior = ric_pk_cpd_ersst,
                                           effect = "sst", species = "pink", 
                                           # model = "Ricker model with CPD", 
                                           model = "Sea-surface temperature",
                                           # xlim = c(-1, 1))
                                           xlim = c(-0.5, 0.5))

plot_pink_L <- plot_npgo_effect_manuscript(posterior = ric_pk_cpd_ersst,
                                            effect = "npgo", species = "pink", 
                                            # model = "Ricker model with CPD",
                                            model = "North Pacific Gyre Oscillation",
                                            # xlim = c(-1, 1))
                                            xlim = c(-0.5, 0.5))


plot_pink_M <- plot_productivity_decline_manuscript(posterior = ric_pk_cpd_ersst,  
                                                     species = "pink",
                                                     effect = "cpd", 
                                                     # model = "Ricker model with CPD")
                                                     model = "Cumulative disturbance")+
  ylim(c(-75,75))+
  labs(y = "Median change in recruitment (%)")

plot_pink_N <- plot_productivity_decline_manuscript(posterior = ric_pk_eca_ersst,
                                                     species = "pink",
                                                     effect = "eca", 
                                                     # model = "Ricker model with ECA")
                                                     model = "Equivalent clearcut area")+
  ylim(c(-75,75))+
  labs(y = "Median change in recruitment (%)")

plot_pink_O <- plot_productivity_decline_manuscript(posterior = ric_pk_cpd_ersst,
                                                     species = "pink",
                                                     effect = "sst", 
                                                     # model = "Ricker model with CPD")
                                                     model = "Sea-surface temperature")+
  ylim(c(-75,75))+
  labs(y = "Median change in recruitment (%)")

plot_pink_P <- plot_productivity_decline_manuscript(posterior = ric_pk_cpd_ersst,
                                                     species = "pink",
                                                     effect = "npgo", 
                                                     # model = "Ricker model with CPD")
                                                     model = "North Pacific Gyre Oscillation")+
  ylim(c(-75,75))+
  labs(y = "Median change in recruitment (%)")






# plot_pink <- (plot_pink_a1 / plot_pink_a2 / plot_pink_a3 / plot_pink_a4)+ plot_layout(axes = 'collect') |
#   ((plot_pink_c1  + plot_pink_c2 ) + plot_layout(axes = 'collect'))/(
#     (plot_pink_c3  + plot_pink_c4 ) + plot_layout(axes = 'collect'))

plot_pink <- (plot_pink_I / plot_pink_J / plot_pink_K / plot_pink_L)+ plot_layout(axes = 'collect') |
  ((plot_pink_M  + plot_pink_N ) + plot_layout(axes = 'collect'))/(
    (plot_pink_O  + plot_pink_P ) + plot_layout(axes = 'collect'))


plot_pink2 <- (plot_pink_I / plot_pink_J / plot_pink_K / plot_pink_L)+ plot_layout(axes = 'collect') |
  ((plot_pink_M  / plot_pink_N ) + plot_layout(axes = 'collect'))

plot_pink_chum_w_title2 <- ((wrap_elements(panel = plot_chum + plot_annotation(tag_levels = "A", 
                                                                              title = "Chum salmon")& 
                                            theme(plot.tag = element_text(size = 8)))) / (wrap_elements(panel = plot_pink+  
                                                                                                          plot_annotation(tag_levels = list(c('I','J','K','L','M','N','O','P')), 
                                                                                                                          title = "Pink salmon")& 
                                                                                                          theme(plot.tag = element_text(size = 8))))) 
plot_pink_chum_w_title2



ggsave(here("figures","manuscript_nov2025_chum_pink_ricker_only_w_CI.png"), plot = plot_pink_chum_w_title2, width = 10, height = 12)

plot_pink_chum_w_title3 <- ((wrap_elements(panel = plot_chum2
                                           
                                           + plot_annotation(tag_levels = "A", 
                                                                               title = "Chum salmon")& 
                                             theme(plot.tag = element_text(size = 8)))) / (wrap_elements(panel = plot_pink2+  
                                                                                                           plot_annotation(tag_levels = list(c('I','J','K','L','M','N','O','P')), 
                                                                                                                           title = "Pink salmon")& 
                                                                                                           theme(plot.tag = element_text(size = 8))))) 
plot_pink_chum_w_title3

#figure of posterior distribution of ratio of effect sizes

posterior1 <- ric_chm_cpd_ocean_covariates_logR

ratio_forestry_sst_npgo <-  posterior1 %>% 
  select("b_for","b_sst","b_npgo") %>% 
  mutate(ratio_for_sst = b_for/b_sst,
         ratio_for_npgo = b_for/b_npgo)

#plot distribution of ratio
ggplot(ratio_forestry_sst_npgo, aes(x = ratio_for_sst)) +
  geom_density(fill = "cadetblue", alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  xlim(-20,20)+
  labs(title = "Posterior distribution of ratio of forestry to SST effect sizes",
       x = "Ratio of forestry to SST effect sizes",
       y = "Density") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))


#plot distribution of ratio
ggplot(ratio_forestry_sst_npgo, aes(x = ratio_for_npgo)) +
  geom_density(fill = "cadetblue", alpha = 0.5) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "gray40") +
  xlim(-20,20)+
  labs(title = "Posterior distribution of ratio of forestry to NPGO effect sizes",
       x = "Ratio of forestry to NPGO effect sizes",
       y = "Density") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))


# calculate ratios of effects for each CU and plot distributions for all CUS

ratio_forestry_sst_npgo_cu <-  posterior1 %>% 
  select(starts_with(c("b_for_cu","b_sst_cu","b_npgo_cu"))) %>% 
  mutate(mcmc_number = as.numeric(rownames(.))) %>%
  pivot_longer(-mcmc_number, names_to = "parameter", values_to = "value") %>%
  # pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(effect = case_when(str_detect(parameter, "b_for_cu") ~ "forestry",
                            str_detect(parameter, "b_sst_cu") ~ "sst",
                            str_detect(parameter, "b_npgo_cu") ~ "npgo"),
         CU_n = str_extract(parameter, "\\[(.*?)\\]") %>% str_remove_all("\\[|\\]")) %>%
  select(-c("parameter")) %>% 
  pivot_wider(names_from = effect, values_from = value) %>%
  mutate(ratio_for_sst = forestry/sst,
         ratio_for_npgo = forestry/npgo)


#plot density plots of ratio for each CU
sst_ratio <- ggplot() +
  stat_density(data = ratio_forestry_sst_npgo_cu, aes(x = ratio_for_sst, group = CU_n),
               color = "#C9AE9F", geom = 'line', position = 'identity', 
               alpha = 0.2, linewidth = 1.5)+
  # geom_vline(xintercept = -1, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 1, color = 'slategray', linewidth = 0.8) +
  
  geom_vline(xintercept = median(ratio_forestry_sst_npgo$ratio_for_sst), color = 'black', linetype = 'dashed') + 
  #also include overall ratio as black bold line
  geom_density(data = ratio_forestry_sst_npgo, aes(x = ratio_for_sst), 
               color = "black", size = 1.5)+
  xlim(-10,10)+
  labs(title = "",
       x = "Ratio of cumulative disturbance to SST effect sizes",
       y = "Posterior density") +
  theme_classic() +
  theme(legend.position = "right",
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.8, "lines"),
        #1 column for legend
        legend.cox = "vertical",
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.001, "cm"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        # plot.title = element_text(size = 10, hjust = 0.5)
        plot.title = element_text(size = 12, hjust = 0, vjust = -0.1))

#ggsave
ggsave(here("figures","manuscript_oct2025_chum_ricker_ratio_forestry_sst_effect_sizes_by_cu.png"),
       width = 6, height = 4, bg = "white")



#plot density plots of ratio for each CU
npgo_ratio <- ggplot() +
  stat_density(data = ratio_forestry_sst_npgo_cu, aes(x = ratio_for_npgo, group = CU_n),
               color = "#A9B7CC", geom = 'line', position = 'identity', 
               alpha = 0.2, linewidth = 1.5)+
  # geom_vline(xintercept = -1, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = -1, color = 'slategray', linewidth = 0.8) +
  
  geom_vline(xintercept = median(ratio_forestry_sst_npgo$ratio_for_npgo), color = 'black', linetype = 'dashed') + 
  #also include overall ratio as black bold line
  geom_density(data = ratio_forestry_sst_npgo, aes(x = ratio_for_npgo), 
               color = "black", size = 1.5)+
  xlim(-10,10)+
  labs(title = "",
       x = "Ratio of cumulative disturbance to NPGO effect sizes",
       y = "Posterior density") +
  theme_classic() +
  theme(legend.position = "right",
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.8, "lines"),
        #1 column for legend
        legend.cox = "vertical",
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.001, "cm"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        # plot.title = element_text(size = 10, hjust = 0.5)
        plot.title = element_text(size = 12, hjust = 0, vjust = -0.1))

#ggsave

ggsave(here("figures","manuscript_oct2025_chum_ricker_ratio_forestry_npgo_effect_sizes_by_cu.png"),
       width = 6, height = 4, bg = "white")




sst_ratio + npgo_ratio

ggsave(here("figures","manuscript_oct2025_chum_ricker_ratio_forestry_sst_npgo_effect_sizes_by_cu.png"),
       width = 8, height = 4, bg = "white")


#do same for pink

posterior_pink <- ric_pk_cpd_ersst

ratio_forestry_sst_npgo_pk <- posterior_pink %>% 
  select("b_for","b_sst","b_npgo") %>% 
  mutate(ratio_for_sst = b_for/b_sst,
         ratio_for_npgo = b_for/b_npgo)

# calculate ratios of effects for each CU and plot distributions for all CUS

ratio_forestry_sst_npgo_cu_pk <-  posterior_pink %>% 
  select(starts_with(c("b_for_cu","b_sst_cu","b_npgo_cu"))) %>% 
  mutate(mcmc_number = as.numeric(rownames(.))) %>%
  pivot_longer(-mcmc_number, names_to = "parameter", values_to = "value") %>%
  # pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(effect = case_when(str_detect(parameter, "b_for_cu") ~ "forestry",
                            str_detect(parameter, "b_sst_cu") ~ "sst",
                            str_detect(parameter, "b_npgo_cu") ~ "npgo"),
         CU_n = str_extract(parameter, "\\[(.*?)\\]") %>% str_remove_all("\\[|\\]")) %>%
  select(-c("parameter")) %>% 
  pivot_wider(names_from = effect, values_from = value) %>%
  mutate(ratio_for_sst = forestry/sst,
         ratio_for_npgo = forestry/npgo)


#plot density plots of ratio for each CU
sst_ratio_pk <- ggplot() +
  stat_density(data = ratio_forestry_sst_npgo_cu_pk, aes(x = ratio_for_sst, group = CU_n),
               color = "#C9AE9F", geom = 'line', position = 'identity', 
               alpha = 0.2, linewidth = 1.5)+
  # geom_vline(xintercept = -1, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 1, color = 'slategray', linewidth = 0.8) +
  
  geom_vline(xintercept = median(ratio_forestry_sst_npgo_pk$ratio_for_sst), color = 'black', linetype = 'dashed') + 
  #also include overall ratio as black bold line
  geom_density(data = ratio_forestry_sst_npgo_pk, aes(x = ratio_for_sst), 
               color = "black", size = 1.5)+
  xlim(-10,10)+
  labs(title = "",
       x = "Ratio of cumulative disturbance to SST effect sizes",
       y = "Posterior density") +
  theme_classic() +
  theme(legend.position = "right",
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.8, "lines"),
        #1 column for legend
        legend.cox = "vertical",
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.001, "cm"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        # plot.title = element_text(size = 10, hjust = 0.5)
        plot.title = element_text(size = 12, hjust = 0, vjust = -0.1))



#plot density plots of ratio for each CU
npgo_ratio_pk <- ggplot() +
  stat_density(data = ratio_forestry_sst_npgo_cu_pk, aes(x = ratio_for_npgo, group = CU_n),
               color = "#A9B7CC", geom = 'line', position = 'identity', 
               alpha = 0.2, linewidth = 1.5)+
  # geom_vline(xintercept = -1, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = -1, color = 'slategray', linewidth = 0.8) +
  
  geom_vline(xintercept = median(ratio_forestry_sst_npgo_pk$ratio_for_npgo), color = 'black', linetype = 'dashed') + 
  #also include overall ratio as black bold line
  geom_density(data = ratio_forestry_sst_npgo_pk, aes(x = ratio_for_npgo), 
               color = "black", size = 1.5)+
  xlim(-10,10)+
  labs(title = "",
       x = "Ratio of cumulative disturbance to NPGO effect sizes",
       y = "Posterior density") +
  theme_classic() +
  theme(legend.position = "right",
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.8, "lines"),
        #1 column for legend
        legend.cox = "vertical",
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0.001, "cm"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        # plot.title = element_text(size = 10, hjust = 0.5)
        plot.title = element_text(size = 12, hjust = 0, vjust = -0.1))

sst_ratio_pk+npgo_ratio_pk

#ggsave
ggsave(here("figures","manuscript_oct2025_pink_ricker_ratio_forestry_sst_npgo_effect_sizes_by_cu.png"),
       width = 8, height = 4, bg = "white")

# make a table of medians of  ratios of forestry effect to sst effect 
# and forestry effect to sst effect by CU

table_ratio_chum <- ratio_forestry_sst_npgo_cu %>% 
  select(-mcmc_number) %>% 
  mutate(CU_n = as.numeric(CU_n)) %>% 
  group_by(CU_n) %>% 
  summarise(median_ratio_for_sst = median(ratio_for_sst),
            median_ratio_for_npgo = median(ratio_for_npgo)
            # median_forestry = median(forestry),
            # median_sst = median(sst),
            # median_npgo = median(npgo)
            ) %>% 
  left_join(ch20rsc %>% select(CU_n, CU_name) %>% distinct(), by = "CU_n") %>% 
  #round all values to two decimal places
  mutate(across(where(is.numeric), ~ round(., 2))) %>% 
  select(-CU_n, CU_name, median_ratio_for_sst, median_ratio_for_npgo)
  

#save the table
write.csv(table_ratio_chum, here("tables",
                                 "manuscript_nov2025_chum_ricker_ratio_forestry_sst_npgo_effect_sizes_by_cu.csv"), row.names = FALSE)


# calculate ratios of effects for each CU and plot distributions for all CUS

ratio_forestry_sst_npgo_river <-  posterior1 %>% 
  select(starts_with(c("b_for_rv","b_sst_rv","b_npgo_rv"))) %>% 
  mutate(mcmc_number = as.numeric(rownames(.))) %>%
  pivot_longer(-mcmc_number, names_to = "parameter", values_to = "value") %>%
  # pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(effect = case_when(str_detect(parameter, "b_for_rv") ~ "forestry",
                            str_detect(parameter, "b_sst_rv") ~ "sst",
                            str_detect(parameter, "b_npgo_rv") ~ "npgo"),
         River_n = str_extract(parameter, "\\[(.*?)\\]") %>% str_remove_all("\\[|\\]")) %>%
  select(-c("parameter")) %>% 
  pivot_wider(names_from = effect, values_from = value) %>%
  mutate(ratio_for_sst = forestry/sst,
         ratio_for_npgo = forestry/npgo)


#plot density plots of ratio for each River
sst_ratio_river <- ggplot() +
  stat_density(data = ratio_forestry_sst_npgo_river, aes(x = ratio_for_sst, group = River_n,
               color = "River"), geom = 'line', position = 'identity', 
               alpha = 0.03, linewidth = 1)+
  # geom_vline(xintercept = -1, linetype = "dashed", color = "gray40") +
  # geom_vline(xintercept = 1, color = 'slategray', linewidth = 0.8) +
  
  geom_vline(xintercept = median(ratio_forestry_sst_npgo$ratio_for_sst), color = "black", linetype = 'dashed') + 
  #also include overall ratio as black bold line
  stat_density(data = ratio_forestry_sst_npgo, aes(x = ratio_for_sst, 
               color = "Coastwide"), geom = 'line', position = 'identity', 
               linewidth = 1.2, alpha =1)+
  scale_color_manual(name = "",values = c("River" = "#C9AE9F", "Coastwide" = "black"))+
  xlim(-10,10)+
  labs(title = "",
       x = "Ratio of effect of cumulative disturbance to SST",
       y = "Posterior density") +
  theme_classic()+
  theme(legend.position = c(0.9,0.9),
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
  
  guides(color = guide_legend(override.aes = list(alpha = 0.9, linewidth = 1), ncol = 1))


npgo_ratio_river <- ggplot() +
  stat_density(data = ratio_forestry_sst_npgo_river, aes(x = ratio_for_npgo, group = River_n,
               color = "River"), geom = 'line', position = 'identity', 
               alpha = 0.03, linewidth = 1)+
  # geom_vline(xintercept = -1, linetype = "dashed", color = "gray40") +
  # geom_vline(xintercept = -1, color = 'slategray', linewidth = 0.8) +
  
  geom_vline(xintercept = median(ratio_forestry_sst_npgo$ratio_for_npgo), color = "black", linetype = 'dashed') + 
  #also include overall ratio as black bold line
  stat_density(data = ratio_forestry_sst_npgo, aes(x = ratio_for_npgo, 
               color = "Coastwide"), geom = 'line', position = 'identity', 
               linewidth = 1.2, alpha =1)+
  scale_color_manual(name = "",values = c("River" = "#A9B7CC", "Coastwide" = "black"))+
  xlim(-10,10)+
  labs(title = "",
       x = "Ratio of effect of cumulative disturbance to NPGO",
       y = "Posterior density") +
  theme_classic()+
  theme(legend.position = c(0.9,0.9),
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
  
  guides(color = guide_legend(override.aes = list(alpha = 0.9, linewidth = 1), ncol = 1))

sst_ratio_river + npgo_ratio_river

#save

ggsave(here("figures","manuscript_nov2025_chum_ricker_ratio_forestry_sst_npgo_effect_sizes_by_river.png"),
       width = 8, height = 4, bg = "white")

ratio_forestry_sst_npgo_river_pink <-  posterior_pink %>% 
  select(starts_with(c("b_for_rv","b_sst_rv","b_npgo_rv"))) %>% 
  mutate(mcmc_number = as.numeric(rownames(.))) %>%
  pivot_longer(-mcmc_number, names_to = "parameter", values_to = "value") %>%
  # pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(effect = case_when(str_detect(parameter, "b_for_rv") ~ "forestry",
                            str_detect(parameter, "b_sst_rv") ~ "sst",
                            str_detect(parameter, "b_npgo_rv") ~ "npgo"),
         River_n = str_extract(parameter, "\\[(.*?)\\]") %>% str_remove_all("\\[|\\]")) %>%
  select(-c("parameter")) %>% 
  pivot_wider(names_from = effect, values_from = value) %>%
  mutate(ratio_for_sst = forestry/sst,
         ratio_for_npgo = forestry/npgo)

#plot density plots of ratio for each River

sst_ratio_river_pk <- ggplot() +
  stat_density(data = ratio_forestry_sst_npgo_river_pink, aes(x = ratio_for_sst, group = River_n,
               color = "River"), geom = 'line', position = 'identity', 
               alpha = 0.03, linewidth = 1)+
  # geom_vline(xintercept = -1, linetype = "dashed", color = "gray40") +
  # geom_vline(xintercept = 1, color = 'slategray', linewidth = 0.8) +
  
  geom_vline(xintercept = median(ratio_forestry_sst_npgo_pk$ratio_for_sst, color = "Coastwide"), linetype = 'dashed') + 
  #also include overall ratio as black bold line
  stat_density(data = ratio_forestry_sst_npgo_pk, aes(x = ratio_for_sst, 
               color = "Coastwide"), geom = 'line', position = 'identity', 
               linewidth = 1.2, alpha =1)+
  scale_color_manual(name = "",values = c("River" = "#C9AE9F", "Coastwide" = "black"))+
  xlim(-10,10)+
  labs(title = "",
       x = "Ratio of effect of cumulative disturbance to SST",
       y = "Posterior density") +
  theme_classic()+
  theme(legend.position = c(0.9,0.9),
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
  
  guides(color = guide_legend(override.aes = list(alpha = 0.9, linewidth = 1), ncol = 1))

npgo_ratio_river_pk <- ggplot() +
  stat_density(data = ratio_forestry_sst_npgo_river_pink, aes(x = ratio_for_npgo, group = River_n,
               color = "River"), geom = 'line', position = 'identity', 
               alpha = 0.03, linewidth = 1)+
  # geom_vline(xintercept = -1, linetype = "dashed", color = "gray40") +
  # geom_vline(xintercept = -1, color = 'slategray', linewidth = 0.8) +
  
  geom_vline(xintercept = median(ratio_forestry_sst_npgo_pk$ratio_for_npgo, color = "Coastwide"), linetype = 'dashed') + 
  #also include overall ratio as black bold line
  stat_density(data = ratio_forestry_sst_npgo_pk, aes(x = ratio_for_npgo, 
               color = "Coastwide"), geom = 'line', position = 'identity', 
               linewidth = 1.2, alpha =1)+
  scale_color_manual(name = "",values = c("River" = "#A9B7CC", "Coastwide" = "black"))+
  xlim(-10,10)+
  labs(title = "",
       x = "Ratio of effect of cumulative disturbance to NPGO",
       y = "Posterior density") +
  theme_classic()+
  theme(legend.position = c(0.9,0.9),
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
  
  guides(color = guide_legend(override.aes = list(alpha = 0.9, linewidth = 1), ncol = 1))

sst_ratio_river_pk + npgo_ratio_river_pk

#save
ggsave(here("figures","manuscript_nov2025_pink_ricker_ratio_forestry_sst_npgo_effect_sizes_by_river.png"),
       width = 8, height = 4, bg = "white")


