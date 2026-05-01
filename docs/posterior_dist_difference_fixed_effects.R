# goal - to plot the posterior distribution of the difference between |effect of forestry| and |effect of NPGO or SST|

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



posterior1 <- ric_chm_cpd_ocean_covariates_logR_long_chain

diff_forestry_sst_npgo <-  posterior1 %>% 
  select("b_for","b_sst","b_npgo") %>% 
  mutate(diff_for_sst = abs(b_for) - abs(b_sst),
         diff_for_npgo = abs(b_for) - abs(b_npgo))

#print median, and 95% CI for diff_for_sst and diff_for_npgo
diff_for_sst_npgo_summary <- diff_forestry_sst_npgo %>% 
  summarise(median_diff_for_sst = median(diff_for_sst),
            ci_lower_diff_for_sst = quantile(diff_for_sst, 0.025),
            ci_upper_diff_for_sst = quantile(diff_for_sst, 0.975),
            median_diff_for_npgo = median(diff_for_npgo),
            ci_lower_diff_for_npgo = quantile(diff_for_npgo, 0.025),
            ci_upper_diff_for_npgo = quantile(diff_for_npgo, 0.975))



#plot distribution of difference
ggplot(diff_forestry_sst_npgo) +
  geom_density(aes(x = diff_for_sst), fill = "orange", alpha = 0.5) +
  geom_density(aes(x = diff_for_npgo), fill = "cadetblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  xlim(-1,1)+
  labs(title = "Posterior distribution of difference of forestry and SST effect sizes",
       x = "Difference of forestry to SST effect sizes",
       y = "Density") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

# calculate ratios of effects for each CU and plot distributions for all CUS

diff_forestry_sst_npgo_cu <-  posterior1 %>% 
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
  mutate(diff_for_sst = abs(forestry) - abs(sst),
         diff_for_npgo = abs(forestry) - abs(npgo))

# make a table with CU, median diff_for_sst, 95% CI, median diff_for_npgo, 95% CI
diff_forestry_sst_npgo_cu_summary <- diff_forestry_sst_npgo_cu %>% 
  group_by(CU_n) %>% 
  summarise(median_diff_for_sst = median(diff_for_sst),
            ci_lower_diff_for_sst = quantile(diff_for_sst, 0.025),
            ci_upper_diff_for_sst = quantile(diff_for_sst, 0.975),
            median_diff_for_npgo = median(diff_for_npgo),
            ci_lower_diff_for_npgo = quantile(diff_for_npgo, 0.025),
            ci_upper_diff_for_npgo = quantile(diff_for_npgo, 0.975)) %>% 
  left_join(cu %>% select(CU_n, CU_name), by = "CU_n")





#plot density plots of diff for each CU
sst_diff <- ggplot() +
  stat_density(data = diff_forestry_sst_npgo_cu, aes(x = diff_for_sst, group = CU_n, color = "CU"),
               geom = 'line', position = 'identity', 
               alpha = 0.2, linewidth = 1.5)+
  # geom_vline(xintercept = -1, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
  
  geom_vline(xintercept = median(diff_forestry_sst_npgo$diff_for_sst), color = 'black', linetype = 'dashed') + 
  #also include overall diff as black bold line
  stat_density(data = diff_forestry_sst_npgo, aes(x = diff_for_sst, color = "coastwide"), 
               geom = 'line', position = 'identity', 
               linewidth = 1.5)+
  xlim(-0.5,0.5)+
  scale_color_manual(values = c("CU" = "#C9AE9F", "coastwide" = "black"))+
  labs(title = "",
       x = expression(abs( phantom(" ") * beta^CDA* phantom(" ") )  - abs( phantom(" ") * beta^SST* phantom(" ") )),
       y = "Posterior density") +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(alpha = 1) ) ) +
  theme(legend.position = c(0.8,0.8),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.8, "lines"),
        #1 column for legend
        legend.cox = "vertical",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.001, "cm"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        # plot.title = element_text(size = 10, hjust = 0.5)
        plot.title = element_text(size = 12, hjust = 0, vjust = -0.1))


sst_diff


# do same for npgo

npgo_diff <- ggplot() +
  stat_density(data = diff_forestry_sst_npgo_cu, aes(x = diff_for_npgo, group = CU_n, color = 'CU'),
               geom = 'line', position = 'identity', 
               alpha = 0.2, linewidth = 1.5)+
  # geom_vline(xintercept = -1, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
  
  geom_vline(xintercept = median(diff_forestry_sst_npgo$diff_for_npgo), color = 'black', linetype = 'dashed') + 
  #also include overall diff as black bold line
  stat_density(data = diff_forestry_sst_npgo, aes(x = diff_for_npgo, color = 'coastwide'), 
               geom = 'line', position = 'identity', 
               linewidth = 1.5)+
  xlim(-0.5,0.5)+
  labs(title = "",
       x = expression(abs( phantom(" ") * beta^CDA* phantom(" ") )  - abs( phantom(" ") * beta^NPGO* phantom(" ") )),
       y = "Posterior density") +
  scale_color_manual(values = c("CU" = "#A9B7CC", "coastwide" = "black"))+
  theme_classic() +
  guides(color = guide_legend(override.aes = list(alpha = 1) ) ) +
  theme(legend.position = c(0.9,0.9),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.8, "lines"),
        #1 column for legend
        legend.cox = "vertical",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.001, "cm"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        # plot.title = element_text(size = 10, hjust = 0.5)
        plot.title = element_text(size = 12, hjust = 0, vjust = -0.1))

npgo_diff

#save sst_dff + npgo_diff
both_plots_cpd <- sst_diff + npgo_diff + plot_layout(axes = 'collect')+
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 12))

ggsave(here('figures','difference_effect_sizes_sst.png'), plot = both_plots_cpd, width = 6, height = 3, dpi = 300)

# do the same with ECA

diff_eca_sst_npgo_cu <-  ric_chm_eca_ocean_covariates_logR_long_chain %>% 
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
  mutate(diff_for_sst = abs(forestry) - abs(sst),
         diff_for_npgo = abs(forestry) - abs(npgo))

diff_eca_sst_npgo <-  ric_chm_eca_ocean_covariates_logR_long_chain %>%  
  select("b_for","b_sst","b_npgo") %>% 
  mutate(diff_for_sst = abs(b_for) - abs(b_sst),
         diff_for_npgo = abs(b_for) - abs(b_npgo))



#plot density plots of diff for each CU
sst_diff_eca <- ggplot() +
  stat_density(data = diff_eca_sst_npgo_cu, aes(x = diff_for_sst, group = CU_n, color = "CU"),
               geom = 'line', position = 'identity', 
               alpha = 0.2, linewidth = 1.5)+
  # geom_vline(xintercept = -1, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
  
  geom_vline(xintercept = median(diff_eca_sst_npgo$diff_for_sst), color = 'black', linetype = 'dashed') + 
  #also include overall diff as black bold line
  stat_density(data = diff_eca_sst_npgo, aes(x = diff_for_sst, color = "coastwide"), 
               geom = 'line', position = 'identity', 
               linewidth = 1.5)+
  xlim(-0.5,0.5)+
  scale_color_manual(values = c("CU" = "#C9AE9F", "coastwide" = "black"))+
  labs(title = "",
       x = expression(abs( phantom(" ") * beta^ECA* phantom(" ") )  - abs( phantom(" ") * beta^SST* phantom(" ") )),
       y = "Posterior density") +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(alpha = 1) ) ) +
  theme(legend.position = "none",
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.8, "lines"),
        #1 column for legend
        legend.cox = "vertical",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.001, "cm"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        # plot.title = element_text(size = 10, hjust = 0.5)
        plot.title = element_text(size = 12, hjust = 0, vjust = -0.1))


sst_diff_eca


#npgo

npgo_diff_eca <- ggplot() +
  stat_density(data = diff_eca_sst_npgo_cu, aes(x = diff_for_npgo, group = CU_n, color = 'CU'),
               geom = 'line', position = 'identity', 
               alpha = 0.2, linewidth = 1.5)+
  # geom_vline(xintercept = -1, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, color = 'slategray', linewidth = 0.8) +
  
  geom_vline(xintercept = median(diff_eca_sst_npgo$diff_for_npgo), color = 'black', linetype = 'dashed') + 
  #also include overall diff as black bold line
  stat_density(data = diff_eca_sst_npgo, aes(x = diff_for_npgo, color = 'coastwide'), 
               geom = 'line', position = 'identity', 
               linewidth = 1.5)+
  xlim(-0.5,0.5)+
  labs(title = "",
       x = expression(abs( phantom(" ") * beta^ECA* phantom(" ") )  - abs( phantom(" ") * beta^NPGO* phantom(" ") )),
       y = "Posterior density") +
  scale_color_manual(values = c("CU" = "#A9B7CC", "coastwide" = "black"))+
  theme_classic() +
  guides(color = guide_legend(override.aes = list(alpha = 1) ) ) +
  theme(legend.position = "none",
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.8, "lines"),
        #1 column for legend
        legend.cox = "vertical",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.001, "cm"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        # plot.title = element_text(size = 10, hjust = 0.5
        plot.title = element_text(size = 12, hjust = 0, vjust = -0.1))

npgo_diff_eca

# put all the figures together

both_plots_eca_cda <- (sst_diff + npgo_diff ) / (sst_diff_eca + npgo_diff_eca ) + plot_layout(axes = 'collect')+
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 12))


both_plots_eca_cda 


#save
ggsave(here('figures','difference_effect_sizes_eca_cda.png'), plot = both_plots_eca_cda, width = 8, height = 6, dpi = 300)









# get median and 95% CI for diff_for_npgo
diff_for_npgo_summary <- diff_forestry_sst_npgo %>% 
  summarise(median_diff_for_npgo = median(diff_for_npgo),
            ci_lower_diff_for_npgo = quantile(diff_for_npgo, 0.025),
            ci_upper_diff_for_npgo = quantile(diff_for_npgo, 0.975))


# for sst

diff_for_sst_summary <- diff_forestry_sst_npgo %>% 
  summarise(median_diff_for_sst = median(diff_for_sst),
            ci_lower_diff_for_sst = quantile(diff_for_sst, 0.025),
            ci_upper_diff_for_sst = quantile(diff_for_sst, 0.975))

# get CU level medians and 95% CI

diff_forestry_sst_npgo_cu_summary <- diff_forestry_sst_npgo_cu %>% 
  group_by(CU_n) %>% 
  summarise(median_diff_for_sst = median(diff_for_sst),
            ci_lower_diff_for_sst = quantile(diff_for_sst, 0.025),
            ci_upper_diff_for_sst = quantile(diff_for_sst, 0.975),
            median_diff_for_npgo = median(diff_for_npgo),
            ci_lower_diff_for_npgo = quantile(diff_for_npgo, 0.025),
            ci_upper_diff_for_npgo = quantile(diff_for_npgo, 0.975))


diff_forestry_sst_npgo_cu_summary


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

effect_sizes_difference_cu_df <- function(posterior,  species){
  
  if(species == "chum"){
    df <- ch20rsc 
    
  } else if(species == "pink"){
    df <- pk10r
  }
  
  effect_df <- NULL
  
  for (i in 1:length(unique(df$CU_n))){
    
    cu <- unique(df$CU_n)[i]
    
    cu_data <- df %>% filter(CU_n == cu)
    
    b_cu_for <- posterior %>% select(starts_with("b_for_cu")) %>%
        select(ends_with(paste0("[",cu,"]")))
      
    b_cu_sst <- posterior %>% select(starts_with("b_sst_cu")) %>%
        select(ends_with(paste0("[",cu,"]")))
    
    b_cu_npgo <- posterior %>% select(starts_with("b_npgo_cu")) %>%
        select(ends_with(paste0("[",cu,"]")))
    
    
    effect_df_cu <- data.frame(CU = unique(cu_data$CU_name),
                               diff_for_sst_median = round(apply(as.matrix(abs(b_cu_for[,1]) - abs(b_cu_sst[,1])),2,median),2),
                               diff_for_sst_025 = round(apply(as.matrix(abs(b_cu_for[,1]) - abs(b_cu_sst[,1])),2,quantile, probs = 0.025),2),
                               diff_for_sst_975 = round(apply(as.matrix(abs(b_cu_for[,1]) - abs(b_cu_sst[,1])),2,quantile, probs = 0.975),2),
                               diff_for_npgo_median = round(apply(as.matrix(abs(b_cu_for[,1]) - abs(b_cu_npgo[,1])),2,median),2),
                               diff_for_npgo_025 = round(apply(as.matrix(abs(b_cu_for[,1]) - abs(b_cu_npgo[,1])),2,quantile, probs = 0.025),2),
                               diff_for_npgo_975 = round(apply(as.matrix(abs(b_cu_for[,1]) - abs(b_cu_npgo[,1])),2,quantile, probs = 0.975),2))
                               
    # sym(effect)_25 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.25),
    # sym(effect)_75 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.75),
    # sym(effect)_025 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.025),
    # sym(effect)_975 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.975))
    
    effect_df <- rbind(effect_df, effect_df_cu)
  }
  
  return(effect_df)
}

# format dataframe "median (lower CI, upper CI)"

effect_sizes_difference_cu_df(ric_chm_cpd_ocean_covariates_logR_long_chain, species = "chum") %>% 
  mutate(diff_for_sst = paste0(diff_for_sst_median, " (", diff_for_sst_025, ", ", diff_for_sst_975, ")"),
         diff_for_npgo = paste0(diff_for_npgo_median, " (", diff_for_npgo_025, ", ", diff_for_npgo_975, ")")) %>% 
  select(CU, diff_for_sst, diff_for_npgo)

#do same for eca

effect_sizes_difference_cu_df(ric_chm_eca_ocean_covariates_logR_long_chain, species = "chum") %>% 
  mutate(diff_for_sst = paste0(diff_for_sst_median, " (", diff_for_sst_025, ", ", diff_for_sst_975, ")"),
         diff_for_npgo = paste0(diff_for_npgo_median, " (", diff_for_npgo_025, ", ", diff_for_npgo_975, ")")) %>% 
  select(CU, diff_for_sst, diff_for_npgo)

#save csv

effect_sizes_difference_cu_df(ric_chm_cpd_ocean_covariates_logR_long_chain, species = "chum") %>% 
  mutate(diff_cda_sst = paste0(diff_for_sst_median, " (", diff_for_sst_025, ", ", diff_for_sst_975, ")"),
         diff_cda_npgo = paste0(diff_for_npgo_median, " (", diff_for_npgo_025, ", ", diff_for_npgo_975, ")")) %>% 
  select(CU, diff_cda_sst, diff_cda_npgo) %>% 
  left_join(
    effect_sizes_difference_cu_df(ric_chm_eca_ocean_covariates_logR_long_chain, species = "chum") %>% 
      mutate(diff_eca_sst = paste0(diff_for_sst_median, " (", diff_for_sst_025, ", ", diff_for_sst_975, ")"),
             diff_eca_npgo = paste0(diff_for_npgo_median, " (", diff_for_npgo_025, ", ", diff_for_npgo_975, ")")) %>% 
      select(CU, diff_eca_sst, diff_eca_npgo),
    by = "CU"
  ) %>% 
  write.csv(here('tables','difference_effect_sizes_for_cu_chum.csv'), row.names = FALSE)


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
                               effect_median = round(apply(as.matrix(b_rv[,1]),2,median),2),
                               CU = unique(river_data$CU_name))
    # sym(effect)_25 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.25),
    # sym(effect)_75 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.75),
    # sym(effect)_025 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.025),
    # sym(effect)_975 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.975))
    
    effect_df <- rbind(effect_df, effect_df_rv)
  }
  
  return(effect_df)
}


effect_sizes_difference_river_df <- function(posterior,  species){
  
  if(species == "chum"){
    df <- ch20rsc 
    
  } else if(species == "pink"){
    df <- pk10r
  }
  
  effect_df <- NULL
  
  for (i in 1:length(unique(df$River_n))){
    
    river <- unique(df$River_n)[i]
    
    river_data <- df %>% filter(River_n == river)
    
    b_rv_for <- posterior %>% select(starts_with("b_for_rv")) %>%
        select(ends_with(paste0("[",river,"]")))
      
    b_rv_sst <- posterior %>% select(starts_with("b_sst_rv")) %>%
        select(ends_with(paste0("[",river,"]")))
    
    b_rv_npgo <- posterior %>% select(starts_with("b_npgo_rv")) %>%
        select(ends_with(paste0("[",river,"]")))
    
    
    effect_df_rv <- data.frame(River = unique(river_data$River),
                               diff_for_sst_median = round(apply(as.matrix(abs(b_rv_for[,1]) - abs(b_rv_sst[,1])),2,median),2),
                               diff_for_sst_025 = round(apply(as.matrix(abs(b_rv_for[,1]) - abs(b_rv_sst[,1])),2,quantile, probs = 0.025),2),
                               diff_for_sst_975 = round(apply(as.matrix(abs(b_rv_for[,1]) - abs(b_rv_sst[,1])),2,quantile, probs = 0.975),2),
                               diff_for_npgo_median = round(apply(as.matrix(abs(b_rv_for[,1]) - abs(b_rv_npgo[,1])),2,median),2),
                               diff_for_npgo_025 = round(apply(as.matrix(abs(b_rv_for[,1]) - abs(b_rv_npgo[,1])),2,quantile, probs = 0.025),2),
                               diff_for_npgo_975 = round(apply(as.matrix(abs(b_rv_for[,1]) - abs(b_rv_npgo[,1])),2,quantile, probs = 0.975),2),
                               CU = unique(river_data$CU_name))
                               
    # sym(effect)_25
    # sym(effect)_75 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.75),
    # sym(effect)_025 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.025),
    # sym(effect)_975 = apply(as.matrix(b_cu[,1]),2,quantile, probs = 0.975))
    
    effect_df <- rbind(effect_df, effect_df_rv)
  }
  
  return(effect_df)
  
}


effect_sizes_difference_river_df(ric_chm_cpd_ocean_covariates_logR_long_chain, species = "chum") %>% 
  mutate(diff_for_sst = paste0(diff_for_sst_median, " (", diff_for_sst_025, ", ", diff_for_sst_975, ")"),
         diff_for_npgo = paste0(diff_for_npgo_median, " (", diff_for_npgo_025, ", ", diff_for_npgo_975, ")")) %>% 
  select(CU, diff_for_sst, diff_for_npgo)




# 
# prop_rivers <- effect_chum_cpd_sst_npgo_rv %>% 
#   mutate(flag = (abs(cpd_effect_median) >= abs(sst_effect_median) & abs(cpd_effect_median) >= abs(npgo_effect_median))) %>%
#   # filter(CU == "Southwest Vancouver Island") %>% 
#   # View()
#   group_by(CU) %>% 
#   summarize(n_rivers_forestry = sum(flag),
#             n_rivers = n()) %>% 
#   mutate(proportion_forestry_greater = round(n_rivers_forestry/n_rivers,2)) %>%
#   select(CU, proportion_forestry_greater)
# 
# 

#ratios


ratio_cda_sst_npgo <-  ric_chm_cpd_ocean_covariates_logR_long_chain %>% 
  select("b_for","b_sst","b_npgo") %>% 
  #check if there are any zeros
  # filter(b_sst != 0, b_npgo != 0)  #none
  mutate(ratio_cda_sst = b_for/b_sst,
         ratio_cda_npgo = b_for/b_npgo)
  
ratio_cda_sst_npgo_cu <- ric_chm_cpd_ocean_covariates_logR_long_chain %>% 
    select(starts_with(c("b_for_cu","b_sst_cu","b_npgo_cu"))) %>% 
    mutate(mcmc_number = as.numeric(rownames(.))) %>%
    pivot_longer(-mcmc_number, names_to = "parameter", values_to = "value") %>%
    # pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(effect = case_when(str_detect(parameter, "b_for_cu") ~ "cda",
                              str_detect(parameter, "b_sst_cu") ~ "sst",
                              str_detect(parameter, "b_npgo_cu") ~ "npgo"),
           CU_n = str_extract(parameter, "\\[(.*?)\\]") %>% str_remove_all("\\[|\\]")) %>%
    select(-c("parameter")) %>% 
    pivot_wider(names_from = effect, values_from = value) %>%
    mutate(ratio_cda_sst = cda/sst,
           ratio_cda_npgo = cda/npgo)
  


  
  
  #plot density plots of ratio for each CU
sst_ratio <- ggplot() +
  stat_density(data = ratio_cda_sst_npgo_cu, aes(x = ratio_cda_sst, group = CU_n, color = "CU"),
               geom = 'line', position = 'identity', 
               alpha = 0.2, linewidth = 1.5)+
  # geom_vline(xintercept = -1, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 1, color = 'slategray', linewidth = 0.8) +
  
  geom_vline(xintercept = median(ratio_cda_sst_npgo$ratio_cda_sst), color = 'black', linetype = 'dashed') + 
  #also include overall ratio as black bold line
  stat_density(data = ratio_cda_sst_npgo, aes(x = ratio_cda_sst, color = "coastwide"), 
               geom = 'line', position = 'identity', 
               linewidth = 1.5)+
  xlim(-20,20)+
  scale_color_manual(values = c("CU" = "#C9AE9F", "coastwide" = "black"))+
  labs(title = "",
       x = expression(abs( phantom(" ") * beta^CDA* phantom(" ") )  - abs( phantom(" ") * beta^SST* phantom(" ") )),
       y = "Posterior density") +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(alpha = 1) ) ) +
  theme(legend.position = c(0.8,0.8),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.8, "lines"),
        #1 column cda legend
        legend.cox = "vertical",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.001, "cm"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        # plot.title = element_text(size = 10, hjust = 0.5)
        plot.title = element_text(size = 12, hjust = 0, vjust = -0.1))


sst_ratio


# do same cda npgo

npgo_ratio <- ggplot() +
  stat_density(data = ratio_cda_sst_npgo_cu, aes(x = ratio_cda_npgo, group = CU_n, color = 'CU'),
               geom = 'line', position = 'identity', 
               alpha = 0.2, linewidth = 1.5)+
  # geom_vline(xintercept = -1, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 1, color = 'slategray', linewidth = 0.8) +
  
  geom_vline(xintercept = median(ratio_cda_sst_npgo$ratio_cda_npgo), color = 'black', linetype = 'dashed') + 
  #also include overall ratio as black bold line
  stat_density(data = ratio_cda_sst_npgo, aes(x = ratio_cda_npgo, color = 'coastwide'), 
               geom = 'line', position = 'identity', 
               linewidth = 1.5)+
  xlim(-20,20)+
  labs(title = "",
       x = expression(abs( phantom(" ") * beta^CDA* phantom(" ") )  - abs( phantom(" ") * beta^NPGO* phantom(" ") )),
       y = "Posterior density") +
  scale_color_manual(values = c("CU" = "#A9B7CC", "coastwide" = "black"))+
  theme_classic() +
  guides(color = guide_legend(override.aes = list(alpha = 1) ) ) +
  theme(legend.position = c(0.9,0.9),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.8, "lines"),
        #1 column for legend
        legend.cox = "vertical",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.001, "cm"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        # plot.title = element_text(size = 10, hjust = 0.5)
        plot.title = element_text(size = 12, hjust = 0, vjust = -0.1))

npgo_ratio
  
















