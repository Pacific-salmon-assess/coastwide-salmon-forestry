# add results to supplementary materials of the manuscript -

## coefficient for CDA, ECA, SST, NPGO
## decline predicted at current average CDA for chum
## decline predicted at current average ECA for chum
## decline predicted at current average CDA for pink  
## decline predicted at current average ECA for pink
## decline predicted at current average SST for pink  
## decline predicted at current average NPGO for pink


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



effect_sizes_df <- function(posterior, species){
  
  if(species == "chum"){
    df <- ch20rsc 
    
  } else if(species == "pink"){
    df <- pk10r
  }
  
  effect_df <- NULL
  b_effect <- posterior %>% select("b_for", "b_sst", "b_npgo")
  
  
  effect_df <- data.frame(effect = names(b_effect),
                          effect_median = round(apply(as.matrix(b_effect),2,median),2),
                          effect_ci_lower = round(apply(as.matrix(b_effect),2,HDInterval::hdi),2)[1,],
                          effect_ci_upper = round(apply(as.matrix(b_effect),2,HDInterval::hdi),2)[2,]
                          )
  return(effect_df)
  
}


chum_effect_sizes_eca <- effect_sizes_df(ric_chm_eca_ocean_covariates_logR_long_chain, species = "chum" )
chum_effect_sizes_cpd <- effect_sizes_df(ric_chm_cpd_ocean_covariates_logR_long_chain, species = "chum" )
pink_effect_sizes_eca <- effect_sizes_df(ric_pk_eca_ersst_long_chain, species = "pink" )
pink_effect_sizes_cpd <- effect_sizes_df(ric_pk_cpd_ersst_long_chain, species = "pink" )

all_effect_sizes <- rbind(
  chum_effect_sizes_eca %>% mutate(species = "Chum", forestry_metric = "ECA"),
  chum_effect_sizes_cpd %>% mutate(species = "Chum", forestry_metric = "CDA"),
  pink_effect_sizes_eca %>% mutate(species = "Pink", forestry_metric = "ECA"),
  pink_effect_sizes_cpd %>% mutate(species = "Pink", forestry_metric = "CDA")
) %>% 
  mutate(effect = case_when(effect == "b_for" ~ "Effect of forestry",
                            effect == "b_sst" ~ "Effect of SST",
                            effect == "b_npgo" ~ "Effect of NPGO")) %>% 
  mutate(effect_size = paste0(as.character(effect_median), " [",as.character(effect_ci_lower), ", ",as.character(effect_ci_upper),"]")) %>%
  select(-effect_median, -effect_ci_lower, -effect_ci_upper) %>%
  pivot_wider(names_from = effect,
              values_from = effect_size)

all_effect_sizes
#save

write.csv(all_effect_sizes, here("tables",
                                 "forestry_ocean_effect_sizes_by_species_and_metric.csv"),
          row.names = FALSE)


#calculate recruitment decline at most recent average CDA for chum and for pink


recruitment_decline_df <- function(posterior, effect, species, covariate_value){
  
  if(species == "chum"){
    df <- ch20rsc 
    
  } else if(species == "pink"){
    df <- pk10r
  }
  
  recruitment_df <- NULL
  
  b <- posterior %>% select(starts_with("b_for"))
  
  if(effect == "eca"){
    covariate_std = (sqrt(covariate_value)-mean(df$sqrt.ECA))/sd(df$sqrt.ECA)
    
    forestry_eca <- seq(0,1, length.out = 100)
    
    forestry_sqrt <- sqrt(forestry_eca)
    
    forestry_sqrt_std = (forestry_sqrt-mean(forestry_sqrt))/sd(forestry_sqrt)
    
    no_forestry <- min(forestry_sqrt_std)
    
    
    recruitment = (exp(as.matrix(b[,1])%*%
                           (covariate_std - no_forestry)))*100 - 100
    
    
  } else if(effect == "cpd"){
    covariate_std = (sqrt(covariate_value)-mean(df$sqrt.CPD))/sd(df$sqrt.CPD)
    
    forestry_cpd <- seq(0,100, length.out = 100)
    
    #to calculate no forestry in the standardized scale
    forestry_sqrt <- sqrt(forestry_cpd)
    
    forestry_sqrt_std = (forestry_sqrt-mean(forestry_sqrt))/sd(forestry_sqrt)
    
    no_forestry <- min(forestry_sqrt_std)
    
    recruitment_decline = (exp(as.matrix(b[,1])%*%
                           (covariate_std - no_forestry)))*100 - 100
  }
  
  
  recruitment_df = data.frame(species = species,
                               effect = effect,
                               covariate_value = covariate_value,
                               recruitment_median = round(median(recruitment),2),
                               recruitment_025 = HDInterval::hdi(recruitment)[1],
                               recruitment_975 = HDInterval::hdi(recruitment)[2]
  )
  
  return(recruitment_df)
}
  

# most recent average CPD

current_average_cpd <- df %>% 
  group_by(River_n) %>% 
  summarize(max_cpd = max(disturbedarea_prct_cs)) %>% 
  summarize(current_average_cpd = mean(max_cpd))

max_average_cpd <- df %>% 
  group_by(River_n) %>% 
  summarize(max_eca = max(ECA_age_proxy_forested_only)) %>% 
  summarize(max_average_eca = mean(max_eca))


#calculate recruitment decline for current_average_cpd and max_average_eca

chum_recruitment_decline_cpd <- recruitment_decline_df(ric_chm_cpd_ocean_covariates_logR_long_chain,
                                                         effect = "cpd",
                                                         species = "chum",
                                                         covariate_value = current_average_cpd$current_average_cpd)

#something is wrong - inconsistent with previously calculated values


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
    
    if(effect == "cpd" || effect == "eca"){
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



