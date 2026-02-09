# Goal - change figures to have equal tailed credible intervals and CDA instead of cumulative disturbance
# And use long chain posteriors


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
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=='915-488000-41400-00000-0000-0000-000-000-000-000-000-000','WINDY BAY CREEK 2',pk10r_o$River)
pk10r_o$River=ifelse(pk10r_o$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',pk10r_o$River)
pk10r_o=pk10r_o[order(factor(pk10r_o$River),pk10r_o$BroodYear),]
rownames(pk10r_o)=seq(1:nrow(pk10r_o))

pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',pk10r_e$River)
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='915-765500-18600-00000-0000-0000-000-000-000-000-000-000','HEAD CREEK 2',pk10r_e$River)
pk10r_e$River=ifelse(pk10r_e$WATERSHED_CDE=='915-488000-41400-00000-0000-0000-000-000-000-000-000-000','WINDY BAY CREEK 2',pk10r_e$River)
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

ric_chm_eca_ocean_covariates_logR_long_chain=read.csv(here('stan models','outs','posterior',
                                                           'ric_chm_eca_ocean_covariates_logR_long_chain.csv'),check.names=F)
ric_chm_cpd_ocean_covariates_logR_long_chain=read.csv(here('stan models','outs','posterior',
                                                           'ric_chm_cpd_ocean_covariates_logR_long_chain.csv'),check.names=F)



ric_pk_eca_ersst_long_chain = read.csv(here('stan models',
                                            'outs',
                                            'posterior',
                                            'ric_pk_eca_st_noac_ocean_covariates_logR_long_chain.csv'),check.names=F)

ric_pk_cpd_ersst_long_chain = read.csv(here('stan models',
                                            'outs',
                                            'posterior',
                                            'ric_pk_cpd_st_noac_ocean_covariates_logR_long_chain.csv'),check.names=F)

#starting with supplementary figures


recruitment_decline_river_df_new <- function(posterior, effect, species){
  
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
    
    forestry_cpd_df = data.frame(forestry_cpd, forestry_sqrt,forestry_sqrt_std)
    
    high_forestry = forestry_cpd_df$forestry_sqrt_std[which.min(abs(forestry_cpd_df$forestry_cpd - cpd_river))]
    
    no_forestry <- min(forestry_sqrt_std)
    
    productivity <- (exp(as.matrix(b_rv[,1])%*%
                           (high_forestry-no_forestry)))*100 - 100
    
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

ric_chm_cpd_recruitment_decline_river <- recruitment_decline_river_df_new(ric_chm_cpd_ocean_covariates_logR_long_chain, 
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
  geom_errorbar(aes(ymin = productivity_025, ymax = productivity_975 ), color = '#516479', width = 0, alpha = 0.5) +
  geom_errorbar(aes(ymin = productivity_25, ymax = productivity_75 ), color = '#516479', width = 0, alpha = 0.7) +
  geom_text_repel(color = "gray20", aes(label = ifelse(important == "yes", paste(River,CU,
                                                                                 paste0(round(productivity_50,1),"%"), sep = ", "), NA)), 
                  size = 3, max.overlaps = 20,
                  direction    = "y", 
                  box.padding = 0.3, hjust = -1.5) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
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
ggsave(here("figures","manuscript_feb2026_chum_ricker_cpd_recruitment_decline_by_river_forest_plot.png"), width = 8, height = 10)


posterior_mu2 <- read.csv(here('stan models','outs','posterior',
                           'ric_chm_cpd_ocean_covariates_logR_long_chain_mu2.csv'),check.names=F)

glimpse(posterior_mu2)


n_rows <- nrow(ch20rsc)

#make df

residual_df <- data.frame(observed = ch20rsc$ln_RS, residual = NA)



residual_df <- data.frame(observed = log(ch20rsc$Recruits), forestry = ch20rsc$disturbedarea_prct_cs) %>% 
  mutate(predicted = posterior_mu2 %>% 
           select(starts_with("mu2")) %>%
           apply(., 2, median),
         residual = observed - predicted)

#plot residuals as a function of fitted

chum_residuals <- ggplot(residual_df)+
  geom_point(aes(x = predicted, y = residual, color = forestry), alpha = 0.2, size = 2) +
  geom_hline(yintercept = 0, color = 'black', linetype = 'dashed') +
  labs(title = "Chum", x = TeX(r"(Predicted $\log (Recruits)$)"), y = "Residuals") +
  ylim(-6, 6) +
  theme_classic() +
  scale_color_gradient2(name = 'CDA (%)',
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
        plot.title = element_text(size = 16, hjust = 0)
  )

# save

ggsave(here("figures", "residuals_w_autocorrelation_logR_chum_long_chain.png"),
       width = 6, height = 4, dpi = 300, units = "in")

# residuals for pink

posterior_mu2_pk <- read.csv(here('stan models','outs','posterior',
                               'ric_pk_cpd_st_noac_ocean_covariates_logR_long_chain_mu2.csv'),check.names=F)

#glimpse(posterior_mu2_pk)

n_rows_pk <- nrow(pk10r)

#make df

residual_df_pk <- data.frame(observed = pk10r$ln_RS, residual = NA)

mu_cols <- grep("^mu2\\[", names(posterior_mu2_pk), value = TRUE)
medians <- vapply(posterior_mu2_pk[mu_cols], median, numeric(1))




residual_df_pk <- data.frame(observed = log(pk10r$Recruits), forestry = pk10r$disturbedarea_prct_cs) %>% 
  mutate(predicted = posterior_mu2_pk %>% 
           select(starts_with("mu2")) %>%
           apply(., 2, median),
         residual = observed - predicted)

#plot residuals as a function of fitted

pink_residuals <- ggplot(residual_df_pk)+
  geom_point(aes(x = predicted, y = residual, color = forestry), alpha = 0.2, size = 2) +
  geom_hline(yintercept = 0, color = 'black', linetype = 'dashed') +
  labs(title = "Pink", x = TeX(r"(Predicted $\log (Recruits)$)"), y = "Residuals") +
  ylim(-6, 6) +
  theme_classic() +
  scale_color_gradient2(name = 'CDA (%)',
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
        plot.title = element_text(size = 16, hjust = 0)
  )

# save

ggsave(here("figures", "residuals_w_autocorrelation_logR_pink_long_chain.png"),
       width = 6, height = 4, dpi = 300, units = "in")


#put them together with titles and save

combined_residuals <- (chum_residuals)/(pink_residuals) + plot_layout(guides = 'collect', axis_titles = 'collect_x')

combined_residuals

ggsave(here("figures", "residuals_w_autocorrelation_logR_chum_pink_long_chain.png"),
       width = 6, height = 6, dpi = 300, units = "in")



# redo CU-level estimated recruitment decline figure  ---------------------

productivity_decline_cu_df_new <- function(posterior, effect, species){
  
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
    
    
    # cpd_sqrt_std_cu <- max(cu_data$sqrt.CPD.std) # should not be using max from CU
    #minimum forestry possible - 0
    real_cpd_cu <- cu_data %>% group_by(River) %>% 
      filter(disturbedarea_prct_cs == max(disturbedarea_prct_cs)) %>% 
      distinct(disturbedarea_prct_cs) %>% 
      ungroup %>% 
      summarize(mean = mean(disturbedarea_prct_cs))
    
    real_cpd_sqrt_std_cu <- cu_data %>% group_by(River) %>% 
      filter(sqrt.CPD.std == max(sqrt.CPD.std)) %>% 
      distinct(sqrt.CPD.std) %>% 
      ungroup %>% 
      summarize(mean = mean(sqrt.CPD.std))
    
    theoretical_cpd <- seq(0,100, length.out = 100)
    
    #to calculate no forestry in the standardized scale
    theoretical_cpd_sqrt <- sqrt(theoretical_cpd)
    
    theoretical_cpd_sqrt_std = (theoretical_cpd_sqrt-mean(theoretical_cpd_sqrt))/sd(theoretical_cpd_sqrt)
    
    #make df for theoretical values
    theoretical_df <- data.frame(theoretical_cpd, theoretical_cpd_sqrt_std)
    
    current_forestry <- theoretical_df$theoretical_cpd_sqrt_std[which.min(abs(theoretical_df$theoretical_cpd - real_cpd_cu$mean))]
    
    no_forestry <- min(theoretical_df$theoretical_cpd_sqrt_std)
    
    productivity <- (exp(as.matrix(b_cu[,1])%*%
                           (current_forestry-no_forestry)))*100 - 100
    
    productivity_median <- apply(productivity,2,median)
    
    productivity_median_df <- data.frame(CU = unique(cu_data$CU_name),
                                         productivity_50 = apply(productivity,2,median),
                                         productivity_25 = apply(productivity,2,quantile, probs = 0.25),
                                         productivity_75 = apply(productivity,2,quantile, probs = 0.75),
                                         productivity_025 = apply(productivity,2,quantile, probs = 0.025),
                                         productivity_975 = apply(productivity,2,quantile, probs = 0.975),
                                         # productivity_025_hdi = apply(productivity,2, hdi, ci = 0.95)[[1]]$CI_low,
                                         # productivity_975_hdi = apply(productivity,2, hdi, ci = 0.95)[[1]]$CI_high,
                                         forestry = real_cpd_cu$mean,
                                         CU_n = unique(cu_data$CU_n))
    
    full_productivity <- rbind(full_productivity, productivity_median_df)
    
    
  }
  
  
  return(full_productivity)
  
  
}


ric_chm_cpd_productivity_decline_cu_new <- productivity_decline_cu_df_new(ric_chm_cpd_ocean_covariates_logR_long_chain, 
                                                                  effect = "cpd", species = "chum")    

cu_forest_plot_new <- ric_chm_cpd_productivity_decline_cu_new %>% 
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
  #add dashed v line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  coord_flip() +
  # scale_color_manual(name = 'Model type', values = c('independent alpha' = 'cadetblue', 'hierarchical alpha' = 'coral', 'hierarchical alpha - ricker' = 'darkgoldenrod')) +
  labs(#title = 'Estimated percent change in CU-level productivity',
    x = 'Conservation Unit',
    y = 'Change in recruitment (%)') +
  theme_classic() +
  theme(legend.position = "none",
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave(here("figures","manuscript_feb2026_chum_ricker_cpd_recruitment_decline_by_cu_forest_plot.png"),
       cu_forest_plot_new, width = 5, height = 6, bg = "white")











