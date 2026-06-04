# goal make correlation plots between forestry disturbance within watersheds in each CU

#libraries

library(ggplot2)
library(tidyverse)
library(here)
library(GGally)
library(hues)
library(sf)
library(patchwork)



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







#make list of unique CU
cu_list <- ch20rsc %>% 
  group_by(CU_NAME) %>% 
  summarise(n = n_distinct(River)) %>% 
  arrange(n) %>% 
  filter(n>1,n<15) %>%
  # select(CU) %>%
  unique()


pop_sheds <- st_read(dsn = here("origional-ecofish-data-models","Data","Spatial","bc_plots",
                                "Max_ECA_sheds_Nov_2024.gpkg"))  # Note that there are 1746 polygons, but only 1745 have VRI information in 2022.

lookup <- read.csv(here("origional-ecofish-data-models","Data","Spatial","bc_plots","lookup.csv"))

ch20rsc_region <- ch20rsc %>% 
  left_join(lookup %>% select(GFE_ID,LINEAR_FEATURE_ID, Species), by = c("GFE_ID" = "GFE_ID", "Species" = "Species")) %>%
  left_join(pop_sheds %>% select(Region,outlet_lfid)%>% mutate(outlet_lfid = as.integer(outlet_lfid)), by = c("LINEAR_FEATURE_ID" = "outlet_lfid"))


panel.ts <- function(data, mapping){
  # print(data$BroodYear)
  # print(data[,mapping$x])
  x <- unique(unlist(lapply(mapping, all.vars)))
  # print(x)
  # print(data[,x])
  new_data <- cbind(values = data[,x], time=data$BroodYear)
  # print(x)
  ggplot(new_data, aes(x = time, y = values))+
    geom_line(size = 1.2, alpha = 0.5, color = "#C9C5BA")
}

panel_two_ts <- function(data, mapping){
  x <- unique(unlist(lapply(mapping, all.vars)))[1]
  y <- unique(unlist(lapply(mapping, all.vars)))[2]
  # print(x)
  # print(y)
  new_data <- cbind(x = data[,x], y = data[,y], time=data$BroodYear)
  ggplot(new_data)+
    geom_line(aes(x = time, y = x), color = "#B89898", size = 1, alpha = 0.8)+
    geom_line(aes(x = time, y = y), color = "#97B1A6", size = 1, alpha =0.8)
}

# plot time series of all watersheds forestry disturbance in the same plot

ggplot(data = ch20rsc)+
  geom_line(aes( x = BroodYear, y = disturbedarea_prct_cs, group = River, color = CU), size = 1, alpha = 0.2)+
  scale_color_iwanthue()+
  labs(x = "Year", y = "CDA")+
  theme_classic()+
  theme(legend.position = "rught")

#save
ggsave(here("figures","timeseries_cda_all_rivers.png"), width = 6, height = 4)

ggplot(data = ch20rsc)+
  geom_line(aes( x = BroodYear, y = ECA_age_proxy_forested_only, group = River, color = CU), size = 1, alpha = 0.2)+
  scale_color_iwanthue()+
  labs(x = "Year", y = "ECA")+
  theme_classic()+
  theme(legend.position = "none")

#save

ggsave(here("figures","timeseries_eca_all_rivers.png"), width = 6, height = 4)

colrs <- c("#e64b35", "#3c5388", "#2ca02c", "#8a4198", "#eea236", "#8f4c2d")



names(colrs) <- c("North Coast - Skeena",
                  "South Coast",
                  "Haida Gwaii",
                  "North Island - Central Coast",
                  "Campbell River",
                  "South Island")

cda_region <- ggplot(data = ch20rsc_region)+
  geom_line(aes( x = BroodYear, y = disturbedarea_prct_cs, group = River, color = Region), size = 0.1, alpha = 0.6)+
  # color = "gray")+
  scale_color_manual(values = colrs)+
  labs(x = "Year", y = "CDA")+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.title = element_blank(),
        legend.direction = "horizontal"
  )+
  guides(colour = guide_legend(override.aes = list(linewidth = 2, alpha = 1), ncol =3))

#save

ggsave(here("figures","timeseries_cda_all_rivers_region.png"), width = 6, height = 4)

cda_gray <- ggplot(data = ch20rsc_region)+
  geom_line(aes( x = BroodYear, y = disturbedarea_prct_cs, group = River, 
                 # color = Region
                 ), 
            size = 0.1, alpha = 0.7,
            color = "gray"
            )+
  scale_color_manual(values = colrs)+
  labs(x = "Year", y = "CDA")+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.title = element_blank(),
        legend.direction = "horizontal"
  )+
  guides(colour = guide_legend(override.aes = list(linewidth = 2, alpha = 1), ncol =3))

#save

ggsave(here("figures","timeseries_cda_all_rivers_gray.png"), width = 6, height = 4)


eca_region <- ggplot(data = ch20rsc_region)+
  geom_line(aes( x = BroodYear, y = ECA_age_proxy_forested_only*100, group = River, color = Region), size = 0.1, alpha = 0.6)+
            # color = "gray")+
  scale_color_manual(values = colrs)+
  labs(x = "Year", y = "ECA")+
  ylim(0,100)+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.title = element_blank(),
        legend.direction = "horizontal"
        )+
  guides(colour = guide_legend(override.aes = list(linewidth = 2, alpha = 1), ncol =3))

#save

ggsave(here("figures","timeseries_eca_all_rivers_region.png"), width = 6, height = 4)

cda_region + eca_region + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(here("figures","timeseries_cda_eca_all_rivers_region.png"), width = 8, height = 4)

eca_gray <- ggplot(data = ch20rsc_region)+
  geom_line(aes( x = BroodYear, y = ECA_age_proxy_forested_only, group = River, 
                 color = Region), 
            size = 0.1, alpha = 0.7,
            color = "gray"
  )+
  scale_color_manual(values = colrs)+
  labs(x = "Year", y = "ECA")+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.title = element_blank(),
        legend.direction = "horizontal"
  )+
  guides(colour = guide_legend(override.aes = list(linewidth = 2, alpha = 1), ncol =6))

#save

ggsave(here("figures","timeseries_eca_all_rivers_gray.png"), width = 6, height = 4)


cda_gray + eca_gray

ggsave(here("figures","timeseries_cda_eca_all_rivers_gray.png"), width = 8, height = 4)



#make a matrix of the all the watershed disturbances. The number of rows are the maximum number of years with data
# the number of columns are the number of ricers

matrix_cda <- ch20rsc %>% 
  select(BroodYear, disturbedarea_prct_cs, River) %>% 
  # group_by(River) %>% 
  # group_modify(~add_row(.x,BroodYear = 1700, disturbedarea_prct_cs = 0)) %>% 
  pivot_wider(names_from = River, values_from = disturbedarea_prct_cs) %>% 
  # add a row with BroodYear 1800 and value 0 for all rivers
  # add_row(c(1700,rep(0,374)))
  mutate(BroodYear = as.integer(BroodYear)) %>% 
  arrange(BroodYear) 
  

matrix_eca <- ch20rsc %>% 
  select(BroodYear, River, ECA_age_proxy_forested_only) %>% 
  # group_by(River) %>% 
  # group_modify(~add_row(.x,BroodYear = 1700, disturbedarea_prct_cs = 0)) %>% 
  pivot_wider(names_from = River, values_from = ECA_age_proxy_forested_only) %>% 
  # add a row with BroodYear 1800 and value 0 for all rivers
  # add_row(c(1700,rep(0,374)))
  mutate(BroodYear = as.integer(BroodYear)) %>% 
  arrange(BroodYear) 




correlation_matrix = cor(matrix_cda[,2:ncol(matrix_cda)], use = "pairwise.complete.obs")

correlation_matrix_eca = cor(matrix_eca[,2:ncol(matrix_eca)], use = "pairwise.complete.obs")


# make histogram of all the pairwise correlations - lower triangle of the correlation matrix

hist(correlation_matrix[lower.tri(correlation_matrix, diag = FALSE)])

cda_matrix <- data.frame(x=correlation_matrix[lower.tri(correlation_matrix, diag = FALSE)])
  # glimpse()
ggplot(cda_matrix) +
  geom_histogram(aes(x), bins = 50, 
                 fill = "#B89898", color = "black") +
  geom_vline(xintercept = mean(cda_matrix$x, na.rm = TRUE))+
  labs(x = "Pairwise correlation between CDA time series across watersheds")+
  theme_classic()

#save plot
ggsave(here("figures","histogram_correlation_cda.png"), width = 6, height = 4)



eca_matrix <- data.frame(x=correlation_matrix_eca[lower.tri(correlation_matrix_eca, diag = FALSE)])
# glimpse()
ggplot(eca_matrix) +
  geom_histogram(aes(x), bins = 50, 
                 fill = "#B89898", color = "black") +
  geom_vline(xintercept = mean(eca_matrix$x, na.rm = TRUE))+
  labs(x = "Pairwise correlation between ECA time series across watersheds")+
  theme_classic()


#save plot
ggsave(here("figures","histogram_correlation_eca.png"), width = 6, height = 4)




# do same for CU


matrix_cda_CU <- ch20rsc %>% 
  select(BroodYear, disturbedarea_prct_cs, River, CU_NAME) %>% 
  group_by(CU_NAME, BroodYear) %>%
  summarize(mean_CDA = mean(disturbedarea_prct_cs, na.rm = TRUE)) %>% 
  # group_modify(~add_row(.x,BroodYear = 1700, disturbedarea_prct_cs = 0)) %>% 
  pivot_wider(names_from = CU_NAME, values_from = mean_CDA) %>% 
  # add a row with BroodYear 1800 and value 0 for all rivers
  # add_row(c(1700,rep(0,374)))
  mutate(BroodYear = as.integer(BroodYear)) %>% 
  arrange(BroodYear) 



correlation_matrix_cda_CU = cor(matrix_cda_CU[,2:ncol(matrix_cda_CU)], use = "pairwise.complete.obs")

cda_matrix_CU <- data.frame(x=correlation_matrix_cda_CU[lower.tri(correlation_matrix_cda_CU, diag = FALSE)])
# glimpse()
ggplot(cda_matrix_CU) +
  geom_histogram(aes(x), bins = 20, 
                 fill = "#B89898", color = "black") +
  geom_vline(xintercept = mean(cda_matrix$x, na.rm = TRUE))+
  labs(x = "Pairwise correlation between mean CDA time series across CUs")+
  theme_classic()

#save plot
ggsave(here("figures","histogram_correlation_cda_CU.png"), width = 6, height = 4)


#do same for eca

matrix_eca_CU <- ch20rsc %>% 
  select(BroodYear, ECA_age_proxy_forested_only, River, CU_NAME) %>% 
  group_by(CU_NAME, BroodYear) %>%
  summarize(mean_ECA = mean(ECA_age_proxy_forested_only, na.rm = TRUE)) %>% 
  # group_modify(~add_row(.x,BroodYear = 1700, disturbedarea_prct_cs = 0)) %>% 
  pivot_wider(names_from = CU_NAME, values_from = mean_ECA) %>% 
  # add a row with BroodYear 1800 and value 0 for all rivers
  # add_row(c(1700,rep(0,374)))
  mutate(BroodYear = as.integer(BroodYear)) %>% 
  arrange(BroodYear)

correlation_matrix_eca_CU = cor(matrix_eca_CU[,2:ncol(matrix_eca_CU)], use = "pairwise.complete.obs")

eca_matrix_CU <- data.frame(x=correlation_matrix_eca_CU[lower.tri(correlation_matrix_eca_CU, diag = FALSE)])

ggplot(eca_matrix_CU) +
  geom_histogram(aes(x), bins = 20, 
                 fill = "#B89898", color = "black") +
  geom_vline(xintercept = mean(eca_matrix$x, na.rm = TRUE))+
  labs(x = "Pairwise correlation between mean ECA time series across CUs")+
  theme_classic()

#save plot
ggsave(here("figures","histogram_correlation_eca_CU.png"), width = 6, height = 4)


# Do same for pinks

eca_gray_pink <- ggplot(data = pk10r)+
  geom_line(aes( x = BroodYear, y = ECA_age_proxy_forested_only, group = River, 
                 # color = Region
                 ), 
            size = 0.1, alpha = 0.6,
            color = "gray"
            )+
  scale_color_manual(values = colrs)+
  labs(x = "Year", y = "ECA")+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.title = element_blank(),
        legend.direction = "horizontal"
  )+
  guides(colour = guide_legend(override.aes = list(linewidth = 2, alpha = 1), ncol =3))

cda_gray_pink <- ggplot(data = pk10r)+
  geom_line(aes( x = BroodYear, y = disturbedarea_prct_cs, group = River, 
                 # color = Region
  ), 
  size = 0.1, alpha = 0.6,
  color = "gray"
  )+
  scale_color_manual(values = colrs)+
  labs(x = "Year", y = "CDA")+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.title = element_blank(),
        legend.direction = "horizontal"
  )+
  guides(colour = guide_legend(override.aes = list(linewidth = 2, alpha = 1), ncol =3))


cda_gray_pink + eca_gray_pink

ggsave(here("figures","timeseries_cda_eca_all_rivers_gray_pink.png"), width = 8, height = 4)

# colour by region

pk10r_region <- pk10r %>% 
  left_join(lookup %>% select(GFE_ID,LINEAR_FEATURE_ID, Species), by = c("GFE_ID" = "GFE_ID", "Species" = "Species")) %>%
  left_join(pop_sheds %>% select(Region,outlet_lfid)%>% mutate(outlet_lfid = as.integer(outlet_lfid)), by = c("LINEAR_FEATURE_ID" = "outlet_lfid"))


cda_region_pink <- ggplot(data = pk10r_region)+
  geom_line(aes( x = BroodYear, y = disturbedarea_prct_cs, group = River, color = Region), size = 0.1, alpha = 0.6)+
  # color = "gray")+
  scale_color_manual(values = colrs)+
  labs(x = "Year", y = "CDA")+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.title = element_blank(),
        legend.direction = "horizontal"
  )+
  guides(colour = guide_legend(override.aes = list(linewidth = 2, alpha = 1), ncol =3))

cda_region_pink

eca_region_pink <- ggplot(data = pk10r_region)+
  geom_line(aes( x = BroodYear, y = ECA_age_proxy_forested_only*100, group = River, color = Region), size = 0.1, alpha = 0.6)+
  # color = "gray")+
  scale_color_manual(values = colrs)+
  labs(x = "Year", y = "ECA")+
  ylim(0,100)+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.title = element_blank(),
        legend.direction = "horizontal"
  )+
  guides(colour = guide_legend(override.aes = list(linewidth = 2, alpha = 1), ncol =3))

eca_region_pink

cda_region_pink + eca_region_pink + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(here("figures","timeseries_cda_eca_all_rivers_region_pink.png"), width = 8, height = 4)



((cda_region + eca_region) + plot_annotation(title = "Chum")) /((cda_region_pink + eca_region_pink) + plot_annotation(title = "Pink")) + 
  plot_layout(guides = "collect") & theme(legend.position = "bottom") & plot_annotation(tag_levels = "A") & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 8, hjust = 0, vjust = 0))


ggsave(here("figures","timeseries_cda_eca_all_rivers_region_chum_pink.png"), width = 8, height = 6)















