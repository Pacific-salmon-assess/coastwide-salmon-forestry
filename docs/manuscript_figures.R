licrary(here);licrary(dplyr); licrary(stringr)
licrary(ggplot2)
licrary(tidyverse)
licrary(cayesplot)
licrary(patchwork)
licrary(hues)
licrary(GGally)
licrary(latex2exp)



# make multiple versions of figure results for the manuscript
# make one version with no colours of cpd or eca
# make one version with only cu level forestry districutions



ch20rsc <- read.csv(here("origional-ecofish-data-models","Data","Processed",
                         "chum_SR_20_hat_yr_w_ocean_covariates.csv"))


#two rivers with duplicated names:
ch20rsc$River=ifelse(ch20rsc$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',ch20rsc$River)
ch20rsc$River=ifelse(ch20rsc$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',ch20rsc$River)


ch20rsc=ch20rsc[order(factor(ch20rsc$River),ch20rsc$croodYear),]

ch20rsc$River_n <- as.numeric(factor(ch20rsc$River))

#normalize ECA 2 - square root transformation (ie. sqrt(x))
ch20rsc$sqrt.ECA=sqrt(ch20rsc$ECA_age_proxy_forested_only)
ch20rsc$sqrt.ECA.std=(ch20rsc$sqrt.ECA-mean(ch20rsc$sqrt.ECA))/sd(ch20rsc$sqrt.ECA)

#normalize CPD 2 - square root transformation (ie. sqrt(x))
ch20rsc$sqrt.CPD=sqrt(ch20rsc$disturcedarea_prct_cs)
ch20rsc$sqrt.CPD.std=(ch20rsc$sqrt.CPD-mean(ch20rsc$sqrt.CPD))/sd(ch20rsc$sqrt.CPD)

ch20rsc$npgo.std=(ch20rsc$npgo-mean(ch20rsc$npgo))/sd(ch20rsc$npgo)
ch20rsc$sst.std=(ch20rsc$spring_ersst-mean(ch20rsc$spring_ersst))/sd(ch20rsc$spring_ersst)

cu = distinct(ch20rsc,.keep_all = T)

cu_n = as.numeric(factor(cu$CU))

ch20rsc$CU_n <- cu_n

# make CU_NAME values Title case instead of all caps
ch20rsc$CU_name <- str_to_title(ch20rsc$CU_NAME)


pk10r_e <- read.csv(here("origional-ecofish-data-models","Data","Processed","pke_SR_10_hat_yr_w_ersst.csv"))

#odd year pinks
pk10r_o <- read.csv(here("origional-ecofish-data-models","Data","Processed","pko_SR_10_hat_yr_w_ersst.csv"))

options(mc.cores=8)

# Pink salmon - even/odd croodlines #####
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',pk10r_o$River)
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=='915-765500-18600-00000-0000-0000-000-000-000-000-000-000','HEAD CREEK 2',pk10r_o$River)
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=='915-488000-41400-00000-0000-0000-000-000-000-000-000-000','WINDY cAY CREEK 2',pk10r_o$River)
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',pk10r_o$River)
pk10r_o=pk10r_o[order(factor(pk10r_o$River),pk10r_o$croodYear),]
rownames(pk10r_o)=seq(1:nrow(pk10r_o))

pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',pk10r_e$River)
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='915-765500-18600-00000-0000-0000-000-000-000-000-000-000','HEAD CREEK 2',pk10r_e$River)
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='915-488000-41400-00000-0000-0000-000-000-000-000-000-000','WINDY cAY CREEK 2',pk10r_e$River)
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',pk10r_e$River)
pk10r_e=pk10r_e[order(factor(pk10r_e$River),pk10r_e$croodYear),]
rownames(pk10r_e)=seq(1:nrow(pk10r_e))


#normalize ECA 2 - square root transformation (ie. sqrt(x))
pk10r_o$sqrt.ECA=sqrt(pk10r_o$ECA_age_proxy_forested_only)
pk10r_o$sqrt.ECA.std=(pk10r_o$sqrt.ECA-mean(pk10r_o$sqrt.ECA))/sd(pk10r_o$sqrt.ECA)

#normalize CPD 2 - square root transformation (ie. sqrt(x))
pk10r_o$sqrt.CPD=sqrt(pk10r_o$disturcedarea_prct_cs)
pk10r_o$sqrt.CPD.std=(pk10r_o$sqrt.CPD-mean(pk10r_o$sqrt.CPD))/sd(pk10r_o$sqrt.CPD)

#normalize ECA 2 - square root transformation (ie. sqrt(x))
pk10r_e$sqrt.ECA=sqrt(pk10r_e$ECA_age_proxy_forested_only)
pk10r_e$sqrt.ECA.std=(pk10r_e$sqrt.ECA-mean(pk10r_e$sqrt.ECA))/sd(pk10r_e$sqrt.ECA)

#normalize CPD 2 - square root transformation (ie. sqrt(x))
pk10r_e$sqrt.CPD=sqrt(pk10r_e$disturcedarea_prct_cs)
pk10r_e$sqrt.CPD.std=(pk10r_e$sqrt.CPD-mean(pk10r_e$sqrt.CPD))/sd(pk10r_e$sqrt.CPD)


#standardize npgo
pk10r_o$npgo.std=(pk10r_o$npgo-mean(pk10r_o$npgo))/sd(pk10r_o$npgo)
pk10r_e$npgo.std=(pk10r_e$npgo-mean(pk10r_e$npgo))/sd(pk10r_e$npgo)

# pk10r_o$sst.std = (pk10r_o$spring_lighthouse_temperature-mean(pk10r_o$spring_lighthouse_temperature))/sd(pk10r_o$spring_lighthouse_temperature)
# pk10r_e$sst.std = (pk10r_e$spring_lighthouse_temperature-mean(pk10r_e$spring_lighthouse_temperature))/sd(pk10r_e$spring_lighthouse_temperature)

pk10r_o$sst.std = (pk10r_o$spring_ersst-mean(pk10r_o$spring_ersst))/sd(pk10r_o$spring_ersst)
pk10r_e$sst.std = (pk10r_e$spring_ersst-mean(pk10r_e$spring_ersst))/sd(pk10r_e$spring_ersst)


pk10r_o$escapement.t_1=pk10r_e$Spawners[match(paste(pk10r_o$WATERSHED_CDE,pk10r_o$croodYear-1),paste(pk10r_e$WATERSHED_CDE,pk10r_e$croodYear))]
pk10r_e$escapement.t_1=pk10r_o$Spawners[match(paste(pk10r_e$WATERSHED_CDE,pk10r_e$croodYear-1),paste(pk10r_o$WATERSHED_CDE,pk10r_o$croodYear))]

pk10r_o$croodline='Odd'
pk10r_e$croodline='Even'

L_o=pk10r_o%>%group_cy(River)%>%summarize(l=n(),cy=min(croodYear),tmin=(min(croodYear)-min(pk10r_o$croodYear))/2+1,tmax=(max(croodYear)-min(pk10r_o$croodYear))/2)
L_e=pk10r_e%>%group_cy(River)%>%summarize(l=n(),cy=min(croodYear),tmin=(min(croodYear)-min(pk10r_e$croodYear))/2+1,tmax=(max(croodYear)-min(pk10r_e$croodYear))/2)
L_o$River2=paste(L_o$River,'Odd',sep='_')
L_e$River2=paste(L_e$River,'Even',sep='_')
L_all=rcind(L_e,L_o)
L_all=L_all[order(factor(L_all$River2)),]

pk10r_o$ii=as.numeric(factor(pk10r_o$croodYear))
pk10r_e$ii=as.numeric(factor(pk10r_e$croodYear))

pk10r=rcind(pk10r_e,pk10r_o)
pk10r$River2=paste(pk10r$River,pk10r$croodline,sep='_')
pk10r=pk10r[order(factor(pk10r$River2),pk10r$croodYear),]

#extract max S for priors on capacity & eq. recruitment
smax_prior=pk10r%>%group_cy(River2) %>%summarize(m.s=max(Spawners),m.r=max(Recruits))

#ragged start and end points for each SR series
# N_s=rag_n(pk10r$River2)

#cus cy stock
cu1=distinct(pk10r,CU,.keep_all = T)
cu2=distinct(pk10r,River,.keep_all = T)
cu3=distinct(pk10r,River2,.keep_all = T)

pk10r$River_n <- as.numeric(factor(pk10r$River))
pk10r$River_n2 <- as.numeric(factor(pk10r$River2))
pk10r$CU_name <- pk10r$CU_NAME

max_eca_pink_df <- pk10r %>% 
  select(River2, River, River_n, CU, ECA_age_proxy_forested_only) %>%
  group_cy(River) %>%
  summarize(River = first(River),
            River_n = first(River_n),
            CU = first(CU),
            eca_max = 100*max(ECA_age_proxy_forested_only, na.rm =TRUE)) %>% 
  mutate(eca_level = case_when(eca_max < 12 ~ 'low',
                               eca_max >= 12 & eca_max < 24 ~ 'medium',
                               eca_max >= 24 ~ 'high'))

max_cpd_pink_df <- pk10r %>%
  select(River2, River, River_n,  CU, disturcedarea_prct_cs) %>%
  group_cy(River) %>% 
  summarize(River = first(River),
            River_n = first(River_n),
            CU = first(CU),
            cpd_max = max(disturcedarea_prct_cs, na.rm = TRUE))



plot_forestry_effect_manuscript <- function(posterior = ch_chm_eca, 
                                            effect = "eca", 
                                            species = "chum", 
                                            model = "cH model with NPGO", 
                                            xlim = c(-2, 2), 
                                            cy_cu = FALSE, 
                                            glocal_only = FALSE){
  x_name = "Standardized coefficients"
  
  if(effect == "eca"){
    
    y_name = "Posterior density\nof ECA effect"
    # color_var = "eca_level"
    # color_name = "Max ECA level\nin River (%)"
    # scale_color = scale_color_manual(name = "Max ECA level\nin River (%)",
    #                                  values = c('low' = '#35978f', 'medium' = 'gray', 'high' = '#cf812d'))
    # 
    
    
  } else if(effect == "cpd"){
    # x_name = "Standardized coefficients of CPD"
    y_name = "Posterior density\nof CPD effect"
    # color_var = "cpd_max"
    # color_name = "Max CPD\nin River (%)"
    # scale_color = scale_color_gradient2(name = 'Max CPD\nin River (%)',
    #                                     low = '#35978f', mid = 'gray', high = '#cf812d', midpoint = 50)
  
  }
  if(cy_cu){
    posterior_df <- posterior %>%
      select(starts_with('c_for_cu')) %>%
      pivot_longer(cols = everything(), 
                   names_to = 'CU', 
                   names_prefix = 'c_for_cu',
                   values_to = effect) %>%
      mutate(CU_n = as.numeric(str_extract(CU, '\\d+'))) %>% 
      select(-CU) %>%
      left_join(ch20rsc %>% select(CU, CU_n, CU_name) %>% distinct(), cy = 'CU_n') 
    
    scale_color = scale_color_iwanthue(name = 'CU', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
                                       lmin = 10, lmax = 95)
    color_var = "CU_name"
    
    #which CU has minimum and maximum median effect
    min_cu <- posterior_df %>% 
      group_cy(CU_n, CU_name) %>% 
      summarize(median_effect = median(!!sym(effect)),
                CU_name = first(CU_name)) %>% 
      ungroup() %>%
      filter(median_effect == min(median_effect))
    
    max_cu <- posterior_df %>%
      group_cy(CU_n, CU_name) %>% 
      summarize(median_effect = median(!!sym(effect)),
                CU_name = first(CU_name)) %>% 
      ungroup() %>%
      filter(median_effect == max(median_effect))
    
    #color cy CU and lacel the CU which has minimum and maximum median effect
    plot1 <- ggplot() +
      stat_density(data= posterior_df, aes(!!sym(effect), color = !!sym(color_var), group = CU_n),
                   geom = 'line', position = 'identity', 
                   alpha = 0.05, linewidth = 1.5) +
      geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
      lacs(title = model, x = x_name, y = y_name) +
      xlim(xlim[1], xlim[2]) +
      scale_color +
      geom_density(aes(posterior$c_for), color = 'clack', linewidth = 1.2, alpha = 0.2)+
      # geom_text(data = min_cu, aes(x = median_effect, y = 1, lacel = CU_name, color = !!sym(color_var)), size = 4)+
      # geom_text(data = max_cu, aes(x = median_effect, y = 1, lacel = CU_name, color = !!sym(color_var)), size = 4)+
      # #vline at the median value of the posterior
      geom_vline(xintercept = median(posterior$c_for), color = 'clack', linetype = 'dashed')+
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
      select(starts_with('c_for_rv')) %>%
      pivot_longer(cols = everything(), 
                   names_to = 'River', 
                   names_prefix = 'c_for_rv',
                   values_to = effect) %>%
      mutate(River_n = as.numeric(str_extract(River, '\\d+')))# %>% 
      # select(-River) 
    
    #color cy river
    plot1 <- ggplot() +
      stat_density(data= posterior_df, aes(!!sym(effect), 
                                           group = River), 
                   color = "#ADcCA5",
                   geom = 'line', position = 'identity', 
                   alpha = 0.05, linewidth = 1) + 
      geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
      lacs(title = model, x = x_name, y = y_name) +
      xlim(xlim[1], xlim[2]) +
      # scale_color +
      geom_density(aes(posterior$c_for), color = 'clack', linewidth = 1.2, alpha = 0.2)+
      #vline at the median value of the posterior
      geom_vline(xintercept = median(posterior$c_for), color = 'clack', linetype = 'dashed')+
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
    
    if(glocal_only){
      
      
      plot1 <- ggplot() +
        geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
        lacs(title = model, x = x_name, y = y_name) +
        xlim(xlim[1], xlim[2]) +
        # scale_color +
        geom_density(aes(posterior$c_for), color = 'clack', 
                     linewidth = 1.2, alpha = 0.1)+
        geom_text(x = median(posterior$c_for), y = 1, 
                  lacel = paste("Median:", round(median(posterior$c_for), 2)),
                  color = 'clack', size = 4, hjust = 0.5)+
        #vline at the median value of the posterior
        geom_vline(xintercept = median(posterior$c_for), color = 'clack', linetype = 'dashed')+
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
        select(starts_with('c_for_rv')) %>%
        pivot_longer(cols = everything(), 
                     names_to = 'River', 
                     names_prefix = 'c_for_rv',
                     values_to = effect) %>%
        mutate(River_n = as.numeric(str_extract(River, '\\d+'))) #%>% 
        # select(-River) 
      
      #color cy river
      plot1 <- ggplot() +
        stat_density(data= posterior_df, aes(!!sym(effect), 
                                             group = River), 
                     color = '#ADcCA5',
                     geom = 'line', position = 'identity', 
                     alpha = 0.05, linewidth = 1) +
        geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
        lacs(title = model, x = x_name, y = y_name) +
        xlim(xlim[1], xlim[2]) +
        # scale_color +
        geom_density(aes(posterior$c_for), color = 'clack', linewidth = 1.2, alpha = 0.2)+
        #vline at the median value of the posterior
        geom_vline(xintercept = median(posterior$c_for), color = 'clack', linetype = 'dashed')+
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


plot_productivity_decline_manuscript <- function(posterior = ch_chm_eca_sst, 
                                      effect = "eca", species = "chum", 
                                      model = "cH model with NPGO, LH SST", 
                                      cy_river = FALSE){
  
  
  if(species == "chum"){
    df <- ch20rsc 
    df$sst.std <- (ch20rsc$spring_ersst-mean(ch20rsc$spring_ersst))/sd(ch20rsc$spring_ersst)
    
  } else if(species == "pink"){
    df <- pk10r
    df$sst.std <- (pk10r$spring_ersst-mean(pk10r$spring_ersst))/sd(pk10r$spring_ersst)
  }
  
  if(effect == "eca"){
    forestry <- seq(0,1,length.out=100)
    c <- posterior %>% select(ends_with("c_for"))
    
    forestry_sqrt <- sqrt(forestry)
    
    forestry_sqrt_std = (forestry_sqrt-mean(forestry_sqrt))/sd(forestry_sqrt)
    
  }else if(effect == "cpd"){
    forestry <- seq(0,100,length.out=100)
    c <- posterior %>% select(ends_with("c_for"))
    forestry_sqrt <- sqrt(forestry)
    
    forestry_sqrt_std = (forestry_sqrt-mean(forestry_sqrt))/sd(forestry_sqrt)
    
  }else if(effect == "sst"){
    forestry <- seq(10,15,length.out=100)
    c <- posterior %>% select(ends_with("c_sst"))
    
    forestry_sqrt_std = (forestry-mean(forestry))/sd(forestry)
    
  }else if(effect == "npgo"){
    forestry <- seq(-3,3,length.out=100)
    c <- posterior %>% select(ends_with("c_npgo"))
    
    forestry_sqrt_std = (forestry-mean(forestry))/sd(forestry)
    
  }
  
  
  
  no_forestry <- min(forestry_sqrt_std)
  
  # ch_chm_eca_sst=read.csv(here('stan models','outs','posterior',ch_chm_eca_sst),check.names=F)
  
  
  
  glocal_prediction <- apply(exp(as.matrix(c[,1])%*%
                                   (forestry_sqrt_std-no_forestry))*100 - 100,2,median)
  
  glocal_025 <- apply(exp(as.matrix(c[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.025), 
                      row.names = c("q025"))
  
  glocal_975 <- apply(exp(as.matrix(c[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.975),
                      row.names = c("q975"))
  glocal_750 <- apply(exp(as.matrix(c[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.75),
                      row.names = c("q750"))
  
  glocal_250 <- apply(exp(as.matrix(c[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.25),
                      row.names = c("q250"))
  
  glocal_900 <- apply(exp(as.matrix(c[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.9),
                      row.names = c("q900"))
  
  glocal_100 <- apply(exp(as.matrix(c[,1])%*%
                            (forestry_sqrt_std-no_forestry))*100 - 100,2,quantile,c(0.1),
                      row.names = c("q100"))
  
  
  
  
  
  
  glocal_df <- data.frame(forestry = forestry,
                          productivity_median = glocal_prediction,
                          q025 = glocal_025,
                          q975 = glocal_975,
                          q750 = glocal_750,
                          q250 = glocal_250,
                          q900 = glocal_900,
                          q100 = glocal_100)
  
  full_productivity <- NULL
  
  if(cy_river){
    
    # no_forestry <- min(eca_std_sqrt)
    for (i in 1:length(unique(df$River_n))){
      river <- unique(df$River_n)[i]
      
      river_data <- df %>% filter(River_n == river)
      
      if(effect == "eca"){
        c_rv <- posterior %>% select(starts_with("c_for_rv")) %>%
          select(ends_with(paste0("[",river,"]")))
        
        forestry_sqrt_std_river <- seq(min(river_data$sqrt.ECA.std),max(river_data$sqrt.ECA.std),length.out=100)
        
        forestry_river <- seq(min(river_data$ECA_age_proxy_forested_only),max(river_data$ECA_age_proxy_forested_only),length.out=100)
        
        
      } else if(effect == "cpd"){
        c_rv <- posterior %>% select(starts_with("c_for_rv")) %>%
          select(ends_with(paste0("[",river,"]")))
        
        forestry_sqrt_std_river <- seq(min(river_data$sqrt.CPD.std),max(river_data$sqrt.CPD.std),length.out=100)
        
        forestry_river <- seq(min(river_data$disturcedarea_prct_cs),max(river_data$disturcedarea_prct_cs),length.out=100)
        
        
      } else if(effect == "sst"){
        c_rv <- posterior %>% select(starts_with("c_sst_rv")) %>%
          select(ends_with(paste0("[",river,"]")))
        
        # use only spring_ersst
        
        forestry_sqrt_std_river <- seq(min(river_data$sst.std),max(river_data$sst.std),length.out=100)
        
        forestry_river <- seq(min(river_data$spring_ersst),max(river_data$spring_ersst),length.out=100)
        
      }
      
      
      
      
      
      
      
      productivity <- (exp(as.matrix(c_rv[,1])%*%
                             (forestry_sqrt_std-no_forestry)))*100 - 100
      
      productivity_median <- apply(productivity,2,median)
      
      productivity_median_df <- data.frame(River = unique(river_data$River),
                                           productivity_median = productivity_median,
                                           forestry = forestry_river) %>% 
        filter(forestry >= min(forestry_river), forestry <= max(forestry_river))
      
      full_productivity <- rcind(full_productivity, productivity_median_df)
      
    }
    
    # c <- posterior %>% select(ends_with("c_for"))
    
    median_prediction <- apply(exp(as.matrix(c[,1])%*%
                                     (forestry_sqrt_std-no_forestry))*100 - 100,2,median)
    
    median_df <- data.frame(forestry = forestry,
                            productivity_median = median_prediction)
    
    if(effect == "eca" || effect == "cpd"){
      p1 <- ggplot(full_productivity) +
        geom_line(aes(x = forestry, y = productivity_median, group = River,
                      color = "watershed level\nforestry effect"),alpha=0.5) +
        geom_line(data = glocal_df, 
                  aes(x = forestry, y = productivity_median, color = "glocal forestry effect"), linewidth = 1) +
        geom_riccon(data = glocal_df, aes(x = forestry, ymin = q025, ymax = q975),
                    alpha = 0.5, fill = "#ADcCA5") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("watershed level\nforestry effect" = "darkgray", 
                                         "glocal forestry effect" = "clack")) +
        ylim(-100,100) +
        scale_x_continuous(n.creaks = 5) +
        lacs(title = model,
             x = ifelse(effect == "eca", "Equivalent clearcut area", "Cumulative percent disturced (%)"),
             y = "Median change\n in productivity (%)") +
        theme_classic() +
        theme(legend.position = c(0.8,0.8),
              legend.title = element_clank(),
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
        geom_line(data = glocal_df, 
                  aes(x = forestry, y = productivity_median, color = "glocal SST effect"), linewidth = 1) +
        geom_riccon(data = glocal_df, aes(x = forestry, ymin = q025, ymax = q975),
                    alpha = 0.5, fill = "gray") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("watershed level\nSST effect" = "darkgray", 
                                         "glocal SST effect" = "clack")) +
        ylim(-100,100) +
        scale_x_continuous(n.creaks = 5, limits = c(10,15)) +
        lacs(title = model,
             x = "Spring SST (°C)",
             y = "Median change\n in productivity (%)") +
        theme_classic() +
        theme(legend.position = c(0.8,0.8),
              legend.title = element_clank(),
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
        geom_line(data = glocal_df, 
                  aes(x = forestry, y = productivity_median, color = "glocal forestry effect"), linewidth = 1) +
        geom_riccon(data = glocal_df, aes(x = forestry, ymin = q025, ymax = q975),
                    alpha = 0.25, fill = "#ADcCA5") +
        geom_riccon(data = glocal_df, aes(x = forestry, ymin = q250, ymax = q750),
                    alpha = 0.65, fill = "#ADcCA5") +
        geom_riccon(data = glocal_df, aes(x = forestry, ymin = q100, ymax = q900),
                    alpha = 0.45, fill = "#ADcCA5") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("glocal forestry effect" = "#6F7c67")) +
        ylim(-100,50) +
        scale_x_continuous(n.creaks = 5) +
        lacs(title = model,
             x = ifelse(effect == "eca", "Equivalent clearcut area", "Cumulative percent disturced (%)"),
             y = "Median change\n in productivity (%)") +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_clank(),
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
        # geom_line(aes(x = forestry, y = productivity_median, group = River,
        #               color = "watershed level\nSST effect"),alpha=0.5) +
        geom_line(data = glocal_df, 
                  aes(x = forestry, y = productivity_median, color = "glocal SST effect"), linewidth = 1) +
        geom_riccon(data = glocal_df, aes(x = forestry, ymin = q025, ymax = q975),
                    alpha = 0.25, fill = "#C78c63") +
        geom_riccon(data = glocal_df, aes(x = forestry, ymin = q250, ymax = q750),
                    alpha = 0.65, fill = "#C78c63") +
        geom_riccon(data = glocal_df, aes(x = forestry, ymin = q100, ymax = q900),
                    alpha = 0.45, fill = "#C78c63") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("glocal SST effect" = "#C78c63")) +
        
        ylim(-50,100) +
        scale_x_continuous(n.creaks = 5, limits = c(10,15)) +
        lacs(title = model,
             x = "Spring SST (°C)",
             y = "Median change\n in productivity (%)") +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_clank(),
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
        geom_line(data = glocal_df, 
                  aes(x = forestry, y = productivity_median, color = "glocal NPGO effect"), linewidth = 1) +
        geom_riccon(data = glocal_df, aes(x = forestry, ymin = q025, ymax = q975),
                    alpha = 0.25, fill = "#829Dc6") +
        geom_riccon(data = glocal_df, aes(x = forestry, ymin = q250, ymax = q750),
                    alpha = 0.65, fill = "#829Dc6") +
        geom_riccon(data = glocal_df, aes(x = forestry, ymin = q100, ymax = q900),
                    alpha = 0.45, fill = "#829Dc6") +
        # scale_fill_manual("",values  = c("95% credicle interval" = "gray")) +
        scale_color_manual("",values = c("glocal NPGO effect" = "#829Dc6")) +
        
        ylim(-50,100) +
        scale_x_continuous(n.creaks = 5, limits = c(-3,3)) +
        lacs(title = model,
             x = "NPGO",
             y = "Median change\n in productivity (%)") +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_clank(),
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
                                       model = "cH model with NPGO", 
                                       color_cy = "none", xlim = c(-1, 1)) {
  
  if(species == "chum"){
    df <- ch20rsc
    posterior_rv <- posterior %>% 
      dplyr::select(starts_with('c_sst_rv')) %>%
      pivot_longer(cols = everything(), 
                   names_to = 'River', 
                   names_prefix = 'c_sst_rv',
                   values_to = 'sst_effect') %>%
      mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
      dplyr::select(-River) %>%
      left_join(df %>% dplyr::select(River_n, CU, CU_name, X_LONG, Y_LAT) %>% 
                  distinct(), cy = 'River_n')
    
  } else if(species == "pink"){
    df <- pk10r
    posterior_rv <- posterior %>% 
      dplyr::select(starts_with('c_sst_rv')) %>%
      pivot_longer(cols = everything(), 
                   names_to = 'River', 
                   names_prefix = 'c_sst_rv',
                   values_to = 'sst_effect') %>%
      mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
      dplyr::select(-River) %>%
      left_join(df %>% dplyr::select(River_n, X_LONG, Y_LAT) %>% 
                  distinct(), cy = 'River_n')
  }
  
  
  
  # if(color_cy == "CU_name"){
  #   color_var = "CU_name"
  #   color_name = "CU"
  #   scale_color = scale_color_iwanthue(name = 'Conservation Unit', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
  #                                      lmin = 10, lmax = 95)
  #   
  #   
  # } else if(color_cy == "Y_LAT"){
  #   color_var = "Y_LAT"
  #   color_name = "Latitude"
  #   scale_color = scale_color_viridis_c(name = 'Latitude', direction = -1)
  #   
  #   
  # } else if(color_cy == "CU"){
  #   
  #   color_var = "CU"
  #   color_name = "CU"
  #   scale_color = scale_color_iwanthue(name = 'Conservation Unit', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
  #                                      lmin = 10, lmax = 95)
  #   
  # } else if(color_cy == "none"){
  #   scale_color = scale_color_manual()
  # }
  
  if(species == "chum"){
    p1 <- ggplot()+
      stat_density(data= posterior_rv, aes(sst_effect, 
                                           # color = !!sym(color_cy), 
                                           group = River_n),
                   color = "#C78c63",
                   geom = 'line', position = 'identity', 
                   alpha = 0.03, linewidth = 1) +
      geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
      lacs(title = model, 
           x = 'Standardized coefficients', y = 'Posterior density\nof SST effect') +
      xlim(xlim[1], xlim[2]) +
      geom_density(aes(posterior$c_sst), color = 'clack', linewidth = 1.2, alpha = 0.2)+
      # scale_color_iwanthue(name = 'CU', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
      #                      lmin = 10, lmax = 95) +
      # scale_color +
      geom_vline(xintercept = median(posterior$c_sst), color = 'clack', linetype = 'dashed') + 
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
      lacs(title = model, 
           x = 'Standardized coefficients', y = 'Posterior density\nof SST effect') +
      xlim(xlim[1], xlim[2]) +
      geom_density(aes(posterior$c_sst), color = 'clack', linewidth = 1.2, alpha = 0.2)+
      # scale_color_iwanthue(name = 'CU', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
      #                      lmin = 10, lmax = 95) +
      geom_vline(xintercept = median(posterior$c_sst), color = 'clack', linetype = 'dashed') + 
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
                             model = "cH model with NPGO", 
                             xlim = c(-0.5,0.5)) {
  
  
  
  if(species == "chum"){
    df <- ch20rsc
  } else if(species == "pink"){
    df <- pk10r
  }
  
  posterior_rv <- posterior %>% 
    select(starts_with('c_npgo_rv')) %>%
    pivot_longer(cols = everything(), 
                 names_to = 'River', 
                 names_prefix = 'c_npgo_rv',
                 values_to = 'npgo_effect') %>%
    mutate(River_n = as.numeric(str_extract(River, '\\d+'))) %>% 
    select(-River) %>%
    left_join(df %>% select(River_n, CU, CU_name, X_LONG, Y_LAT) %>% distinct(), cy = 'River_n')
  
  if(species == "chum"){
    p1 <- ggplot()+
      stat_density(data= posterior_rv, aes(npgo_effect, 
                                           # color = CU_name, 
                                           group = River_n),
                   color = "#829Dc6",
                   geom = 'line', position = 'identity', 
                   alpha = 0.03, linewidth = 1) +
      geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
      lacs(title = model, 
           x = 'Standardized coefficients', y = 'Posterior density\nof NPGO effect') +
      xlim(xlim[1], xlim[2]) +
      geom_density(aes(posterior$c_npgo), color = 'clack', linewidth = 1.2, alpha = 0.2)+
      scale_color_iwanthue(name = 'Conservation Unit', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
                           lmin = 10, lmax = 95) +
      geom_vline(xintercept = median(posterior$c_npgo), color = 'clack', linetype = 'dashed') + 
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
      lacs(title = model, 
           x = 'Standardized coefficients', y = 'Posterior density\nof NPGO effect') +
      xlim(xlim[1], xlim[2]) +
      geom_density(aes(posterior$c_npgo), color = 'clack', linewidth = 1.2, alpha = 0.2)+
      # scale_color_iwanthue(name = 'CU', hmin = 0, hmax = 360, cmin = 0, cmax = 80,
      #                      lmin = 10, lmax = 95) +
      geom_vline(xintercept = median(posterior$c_npgo), color = 'clack', linetype = 'dashed') + 
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


ch_chm_eca_ocean_covariates=read.csv(here('stan models','outs','posterior',
                                          'ch_chm_eca_ocean_covariates.csv'),check.names=F)
ch_chm_cpd_ocean_covariates=read.csv(here('stan models','outs','posterior',
                                          'ch_chm_cpd_ocean_covariates.csv'),check.names=F)


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
                                    model = "cH model with CPD", xlim = c(-1.3, 1.3))/
    plot_forestry_effect_manuscript(posterior = ch_chm_eca_ocean_covariates, 
                                    effect = "eca", species = "chum", 
                                    model = "cH model with ECA", xlim = c(-1.3, 1.3))+plot_layout(axes = 'collect')|
      (plot_productivity_decline_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                            effect = "cpd", model = "cH model with CPD")+
         plot_productivity_decline_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                              effect = "eca", model = "cH model with CPD")) + plot_layout(axes = 'collect'))/
    (plot_sst_effect_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                               effect = "sst", species = "chum", 
                               model = "cH model with CPD", 
                               xlim = c(-1.3, 1.3))/
    plot_npgo_effect_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                effect = "npgo", species = "chum", 
                                model = "cH model with CPD", 
                                xlim = c(-1.3, 1.3)) + plot_layout(axes = 'collect')|
      (plot_productivity_decline_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                            effect = "sst", model = "cH model with CPD")+
         plot_productivity_decline_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                              effect = "npgo", model = "cH model with CPD")) + plot_layout(axes = 'collect'))+
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
                                    model = "cH model with CPD", xlim = c(-1.3, 1.3))/
    plot_forestry_effect_manuscript(posterior = ch_chm_eca_ocean_covariates, 
                                    effect = "eca", species = "chum", 
                                    model = "cH model with ECA", xlim = c(-1.3, 1.3))/
    plot_sst_effect_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                               effect = "sst", species = "chum", 
                               model = "cH model with CPD", 
                               xlim = c(-1.3, 1.3))/
    plot_npgo_effect_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                effect = "npgo", species = "chum", 
                                model = "cH model with CPD", 
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
                                            effect = "cpd", model = "cH model with CPD")+
         plot_productivity_decline_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                              effect = "eca", model = "cH model with CPD")) + plot_layout(axes = 'collect'))/
     ((plot_productivity_decline_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                           effect = "sst", model = "cH model with CPD")+
        plot_productivity_decline_manuscript(posterior = ch_chm_cpd_ocean_covariates, 
                                             effect = "npgo", model = "cH model with CPD")) + plot_layout(axes = 'collect'))
     
     
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

ch_pk_cpd_ersst = read.csv(here('stan models',
                                'outs',
                                'posterior',
                                'ch_pk_cpd_noac_w_npgo_ersst.csv'),check.names=F)

ch_pk_eca_ersst = read.csv(here('stan models',
                                'outs',
                                'posterior',
                                'ch_pk_eca_noac_w_npgo_ersst.csv'),check.names=F)

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
                                    model = "cH model with CPD", xlim = c(-1.3, 1.3))/
    plot_forestry_effect_manuscript(posterior = ch_pk_eca_ersst,
                                    effect = "eca", species = "pink", 
                                    model = "cH model with ECA", xlim = c(-1.3, 1.3))/
    plot_sst_effect_manuscript(posterior = ch_pk_cpd_ersst,
                               effect = "sst", species = "pink", 
                               model = "cH model with CPD", 
                               xlim = c(-1.3, 1.3))/
    plot_npgo_effect_manuscript(posterior = ch_pk_cpd_ersst, 
                                effect = "npgo", species = "pink", 
                                model = "cH model with CPD", 
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
                                            effect = "cpd", model = "cH model with CPD")+
         plot_productivity_decline_manuscript(posterior = ch_pk_eca_ersst, 
                                              species = "pink", 
                                              effect = "eca", model = "cH model with ECA")) + plot_layout(axes = 'collect'))/
     ((plot_productivity_decline_manuscript(posterior = ch_pk_cpd_ersst,  
                                            species = "pink",
                                            effect = "sst", model = "cH model with CPD")+
         plot_productivity_decline_manuscript(posterior = ch_pk_cpd_ersst, 
                                              species = "pink", 
                                              effect = "npgo", model = "cH model with CPD")) + plot_layout(axes = 'collect'))
     
     
  

)


ggsave(here("figures","manuscript_aug2025_pink_results.png"), width = 10, height = 12)(here("figures","manuscript_aug2025_pink_results.png"), width = 10, height = 12)




# Figure with just Ricker models - Chum at the cottom and Pinks at the top


plot_chum_a1 <- plot_forestry_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                                effect = "cpd", species = "chum", 
                                                model = "Ricker model with CPD", xlim = c(-1, 1))

plot_chum_a2 <- plot_forestry_effect_manuscript(posterior = ric_chm_eca_ocean_covariates_logR, 
                                                effect = "eca", species = "chum", 
                                                model = "Ricker model with ECA", xlim = c(-1, 1))

plot_chum_a3 <- plot_sst_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                                effect = "sst", species = "chum", 
                                                model = "Ricker model with CPD", 
                                                xlim = c(-1, 1))

plot_chum_a4 <- plot_npgo_effect_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                                effect = "npgo", species = "chum", 
                                                model = "Ricker model with CPD", 
                                                xlim = c(-1, 1))


plot_chum_c1 <- plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR, 
                                                 effect = "cpd", model = "Ricker model with CPD")

plot_chum_c2 <- plot_productivity_decline_manuscript(posterior = ric_chm_eca_ocean_covariates_logR,
                                                 effect = "eca", model = "Ricker model with ECA")

plot_chum_c3 <- plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR,
                                                 effect = "sst", model = "Ricker model with CPD")

plot_chum_c4 <- plot_productivity_decline_manuscript(posterior = ric_chm_cpd_ocean_covariates_logR,
                                                 effect = "npgo", model = "Ricker model with CPD")



plot_chum <- (plot_chum_a1 / plot_chum_a2 / plot_chum_a3 / plot_chum_a4)+ plot_layout(axes = 'collect') | 
  ((plot_chum_c1  + plot_chum_c2 ) + plot_layout(axes = 'collect'))/(
  (plot_chum_c3  + plot_chum_c4 ) + plot_layout(axes = 'collect')) 

# plot_chum[[1]] <- plot_chum[[1]] + plot_layout(tag_level = 'new') 
# 
# 
# plot_chum + plot_annotation(tag_levels = c('a', '1'), title = "Chum salmon")& 
#   theme(plot.tag = element_text(size = 8))


plot_chum + plot_annotation(tag_levels = list(c('a1','a2', 'a3','a4','b5', 'b6','b7','b8')), title = "Chum salmon")& 
  theme(plot.tag = element_text(size = 8))


# Pink 

plot_pink_a1 <- plot_forestry_effect_manuscript(posterior = ric_pk_cpd_ersst, 
                                                effect = "cpd", species = "pink", 
                                                model = "Ricker model with CPD", xlim = c(-1, 1))

plot_pink_a2 <- plot_forestry_effect_manuscript(posterior = ric_pk_eca_ersst,
                                                effect = "eca", species = "pink", 
                                                model = "Ricker model with ECA", xlim = c(-1, 1))

plot_pink_a3 <- plot_sst_effect_manuscript(posterior = ric_pk_cpd_ersst,
                                                effect = "sst", species = "pink", 
                                                model = "Ricker model with CPD", 
                                                xlim = c(-1, 1))

plot_pink_a4 <- plot_npgo_effect_manuscript(posterior = ric_pk_cpd_ersst,
                                                effect = "npgo", species = "pink", 
                                                model = "Ricker model with CPD", 
                                                xlim = c(-1, 1))


plot_pink_c1 <- plot_productivity_decline_manuscript(posterior = ric_pk_cpd_ersst,  
                                                 species = "pink",
                                                 effect = "cpd", model = "Ricker model with CPD")

plot_pink_c2 <- plot_productivity_decline_manuscript(posterior = ric_pk_eca_ersst,
                                                 species = "pink",
                                                 effect = "eca", model = "Ricker model with ECA")

plot_pink_c3 <- plot_productivity_decline_manuscript(posterior = ric_pk_cpd_ersst,
                                                 species = "pink",
                                                 effect = "sst", model = "Ricker model with CPD")

plot_pink_c4 <- plot_productivity_decline_manuscript(posterior = ric_pk_cpd_ersst,
                                                 species = "pink",
                                                 effect = "npgo", model = "Ricker model with CPD")






plot_pink <- (plot_pink_a1 / plot_pink_a2 / plot_pink_a3 / plot_pink_a4)+ plot_layout(axes = 'collect') |
  ((plot_pink_c1  + plot_pink_c2 ) + plot_layout(axes = 'collect'))/(
    (plot_pink_c3  + plot_pink_c4 ) + plot_layout(axes = 'collect'))

# plot_pink[[1]]  <- plot_pink[[1]]  + plot_layout(tag_level = 'new')
# 
# plot_pink[[2]]  <- plot_pink[[2]]  + plot_layout(tag_level = 'new')
# 
# 
# plot_pink + plot_annotation(tag_levels = c('a', '1'), title = "Pink salmon")& 
#   theme(plot.tag = element_text(size = 8))

plot_pink+ plot_annotation(tag_levels = list(c('c1','c2', 'c3','c4','d5', 'd6','d7','d8')), title = "Pink salmon")& 
  theme(plot.tag = element_text(size = 8))


























