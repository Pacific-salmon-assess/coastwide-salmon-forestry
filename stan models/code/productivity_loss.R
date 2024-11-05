#goal - plot the loss in productivity at current cpd levels

#edits - make colors depend on cpd values and have vertical lines at the means


# Libraries 
library(here)
library(tidyverse)
library(ggplot2)

#data
ch20r <- read.csv(here("origional-ecofish-data-models","Data","Processed","chum_SR_20_hat_yr_w_npgo.csv"))
ch20r$cuid=as.numeric(factor(ch20r$CU))


#two rivers named salmon river... recode
ch20r$River=ifelse(ch20r$WATERSHED_CDE=='950-169400-00000-00000-0000-0000-000-000-000-000-000-000','SALMON RIVER 2',ch20r$River)
ch20r$River=ifelse(ch20r$WATERSHED_CDE=="915-486500-05300-00000-0000-0000-000-000-000-000-000-000",'LAGOON CREEK 2',ch20r$River)


ch20r$disturbedarea_prct_cs.std=scale(ch20r$disturbedarea_prct_cs)

#normalize ECA 2 - square root transformation (ie. sqrt(x))
ch20r$sqrt.ECA=sqrt(ch20r$ECA_age_proxy_forested_only)
ch20r$sqrt.ECA.std=(ch20r$sqrt.ECA-mean(ch20r$sqrt.ECA))/sd(ch20r$sqrt.ECA)

#normalize CPD 2 - square root transformation (ie. sqrt(x))
ch20r$sqrt.CPD=sqrt(ch20r$disturbedarea_prct_cs)
ch20r$sqrt.CPD.std=(ch20r$sqrt.CPD-mean(ch20r$sqrt.CPD))/sd(ch20r$sqrt.CPD)

#add numeric values for each river

ch20r <- ch20r %>% 
  mutate(River_n = as.numeric(as.factor(River)))

#for every river, add a column for maximum cpd

ch20r <- ch20r %>% 
  group_by(River) %>% 
  mutate(max_cpd = max(sqrt.CPD.std), max_year = max(BroodYear))

#plot cpd by year for every river

ch20r %>% 
  ggplot(aes(x=BroodYear,y=sqrt.CPD.std,group=River))+
  geom_line(aes(alpha = 0.2))+
  labs(title='CPD by Year for Every River',
       x='Year',
       y='CPD (standardized)')+
  theme_minimal()



posterior =read.csv(here('stan models',
                        'outs','posterior',
                        'bh_chm_cpd_ac_oct24.csv'),check.names=F)


#time varying alpha and river specific alpha

alpha <- posterior %>% 
  select(starts_with("alpha"))
#alpha_t[1], alpha_t[2] etc and alpha_j[1], alpha_j[2] etc

#calculate alpha_t + alpha_j for every combination of t (year) and j (river)

#from ch20r, select river and max_year, print unique values of river

alpha_combos <- ch20r %>% 
  select(River_n, max_year) %>% 
  unique() %>% 
  mutate(alpha_t = paste0("alpha_t[",unique(max_year) - 1953,"]") , alpha_j = paste0("alpha_j[",River_n,"]"),
         alpha_t_j = paste0("alpha_t[",unique(max_year) - 1953,"]_alpha_j[",River_n,"]"))



#create empty df

df <- data.frame()


# Create new columns for each combination of alpha_t and alpha_j
for (n_col in 1:length(alpha_combos$alpha_t_j)) {
  col <- alpha_combos$alpha_t_j[n_col]
  t_col <- alpha_combos$alpha_t[n_col]
  j_col <- alpha_combos$alpha_j[n_col]
  new_col_name <- paste0(t_col, "_", j_col)
  print(new_col_name)
  alpha <- alpha %>% mutate(!!col := !!sym(t_col) + !!sym(j_col))
  
}




#calculate productivity (R/S) at latest levels of CPD

#loop over all watersheds
median_productivity <- data.frame()
#make empty df with 4000 rows


for (i in 1:length(unique(ch20r$River_n))) {
  river <- unique(ch20r$River_n)[i]
  print(river)
  # print(river)
  #filter data for river
  river_data <- ch20r %>% filter(River_n == river)
  # #filter alpha for river
  river_alpha <- alpha %>% select(ends_with(paste0("_alpha_j[",river,"]")))
  #b effect of forestry
  b <- posterior %>% select(starts_with("b_for_rv")) %>% select(ends_with(paste0("[",river,"]")))
  #cpd
  max_cpd <- max(river_data$sqrt.CPD.std)
  no_forestry <- min(ch20r$sqrt.CPD.std)
  
  # #calculate productivity
  if(i==1){
  productivity <- (exp(as.matrix(river_alpha) + 
                        as.matrix(b)%*%
                        max_cpd)/exp(as.matrix(river_alpha) + 
                                       as.matrix(b)%*%
                                       no_forestry))*100 - 100}
  else{
    productivity <- cbind(productivity, (exp(as.matrix(river_alpha) + 
                        as.matrix(b)%*%
                        max_cpd)/exp(as.matrix(river_alpha) + 
                                       as.matrix(b)%*%
                                       no_forestry))*100 - 100)
  }
  median_productivity <- rbind(median_productivity, 
  data.frame(river = river, median_productivity = median(productivity[,i])))
  # river_data <- river_data %>% mutate(productivity = exp(river_alpha %*% t(as.matrix(river_data$sqrt.CPD.std))) / (1 + exp(river_alpha %*% t(as.matrix(river_data$sqrt.CPD.std))))
  # #add productivity to df
  # df <- rbind(df, river_data)
}

#change the colnames of productivity to river names, make df long

colnames(productivity) <- unique(ch20r$River_n)
# productivity <- data.frame(productivity) %>% 
#   mutate(sample = rownames(productivity))

productivity_long <- data.frame(productivity) %>% 
  pivot_longer(everything(), names_to = 'River_n', 
               values_to = 'productivity_loss') %>% 
  mutate(River_n = as.numeric(substr(River_n, 2, nchar(River_n)))) %>% 
  left_join(ch20r %>% select(River_n, River, CU, max_cpd) %>% 
              distinct(), 
            by = 'River_n', relationship = "many-to-one")



ggplot(median_productivity, 
       aes(x=median_productivity))+
  geom_histogram(aes(alpha = 0.5))+
  # xlim(-250,250)+
  labs(title='Productivity Loss at Current CPD Levels',
       x='Productivity Loss (%)',
       y='Count')+
  theme_minimal()

ggplot(productivity_long %>% filter(productivity_loss < 100), 
       aes(x=productivity_loss))+
  geom_density(aes(alpha = 0.5))+
  facet_wrap(~CU, ncol = 4)+
  # xlim(-250,250)+
  labs(title='Productivity Loss at Current CPD Levels',
       x='Productivity Loss (%)',
       y='Count')+
  theme_classic()+
  theme(legend.position="none",strip.background = element_blank())

  
cu_names <- data.frame(CU = c("CM-1","CM-2","CM-3","CM-4","CM-5","CM-6",
                              "CM-7","CM-8","CM-9","CM-10","CM-11","CM-12",
                              "CM-13","CM-14","CM-15","CM-16","CM-17","CM-18",
                              "CM-19","CM-20","CM-21","CM-22","CM-23","CM-24",
                              "CM-25","CM-26","CM-27","CM-28","CM-29","CM-30",
                              "CM-31","CM-32","CM-33","CM-34","CM-35","CM-36",
                              "CM-37","CM-38", "CM-39"),
                        CU_name = c("Fraser Canyon",
                                    "Lower Fraser",
                                    "Howe Sound-Burrard Inlet",
                                    "Georgia Strait",
                                    "East Vancouver Island",
                                    "Loughborough",
                                    "Bute Inlet",
                                    "Southern Coastal Streams",
                                    "Upper Knight",
                                    "Southwest Vancouver Island",
                                    "Northwest Vancouver Island",
                                    "Smith Inlet",
                                    "Rivers Inlet",
                                    "Wannock",
                                    "Spiller-Fitz Hugh-Burke",
                                    "Bella Coola - Dean Rivers",
                                    "Bella Coola River - Late",
                                    "Hecate Lowlands",
                                    "Mussel-Kynoch",
                                    "Douglas-Gardner",
                                    "East Haida Gwaii",
                                    "Skidegate",
                                    "West Haida Gwaii",
                                    "North Haida Gwaii",
                                    "North Haida Gwaii-Stanley Creek",
                                    "Skeena Estuary",
                                    "Lower Skeena",
                                    "Middle Skeena",
                                    "Upper Skeena",
                                    "Portland Inlet",
                                    "Lower Nass",
                                    "Portland Canal-Observatory",
                                    "Unuk",
                                    "Lower Stikine",
                                    "Whiting",
                                    "Taku",
                                    "Lynn Canal",
                                    "Teslin",
                                    "Lower Liard"))

productivity_long <- productivity_long %>%
  left_join(cu_names, by = 'CU')


ggplot(productivity_long %>% filter(productivity_loss < 100), 
       aes(x=productivity_loss, fill = River))+
  geom_density(linewidth = 0, alpha = 0.2, color = NA)+
  facet_wrap(~CU_name, ncol = 4, scales = "free")+
  scale_colour_hue(l = 70, c = 50)+
  scale_x_continuous(n.breaks = 3, breaks = waiver())+
  labs(x='Change in productivity (%)',
       y='Density')+
  theme_classic()+
  theme(legend.position="none",strip.background = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 13))

ggsave(here('figures','productivity_loss.png'), width = 12, height = 10, dpi = 300)

#make same figure with colors depending on cpd values
means <- productivity_long %>%
  filter(productivity_loss < 100) %>%
  group_by(CU_name, River) %>%
  summarize(mean_loss = mean(productivity_loss, na.rm = TRUE))
calculate_density_at_mean <- function(df) {
  dens <- density(df$productivity_loss)
  mean_value <- mean(df$productivity_loss, na.rm = TRUE)
  # Find the density height closest to the mean value
  closest_x <- which.min(abs(dens$x - mean_value))
  tibble(mean_loss = mean_value, density_height = dens$y[closest_x])
}

# Add density height to means
density_means <- productivity_long %>%
  filter(productivity_loss < 100) %>%
  group_by(CU_name, River) %>%
  group_modify(~ calculate_density_at_mean(.x))

ggplot(productivity_long %>% filter(productivity_loss < 100), 
       aes(x=productivity_loss, fill = max_cpd, group = River))+
  geom_density(linewidth = 0, alpha = 0.5, color = NA)+
  # Add vertical lines for the mean of each distribution
  # geom_vline(data = means, aes(xintercept = mean_loss, group = River), 
  #            color = "gray", linewidth = 0.5, alpha = 0.5) +
  geom_segment(data = density_means, aes(x = mean_loss, xend = mean_loss, 
                                         y = 0, yend = density_height), 
               color = "slategray", linewidth = 0.5, alpha = 0.7) +
  facet_wrap(~CU_name, ncol = 4, scales = "free")+
  scale_x_continuous(n.breaks = 3, breaks = waiver())+
  # scale_y_continuous(n.breaks = 3, breaks = waiver(), limits = c(0, 0.1))+
  scale_fill_gradient2(low = "#5ab4ac",
                       mid = "gray",
                       high = "#d8b365", name = "CPD",
                       guide = guide_colorbar(direction = "horizontal"))+#low = "darkgreen", high = "gray")+
  labs(x='Change in productivity (%)',
       y='Density')+
  theme_classic()+
  theme(legend.position="bottom",
        #legend.position = c(0.94, 0.01),  # Position at bottom right inside the plot
        legend.justification = c(1, 0),   # Align the legend
        legend.background = element_rect(fill = "white", color = "black"), # Optional: add a border and background for clarity
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 13))


ggsave(here('figures','productivity_loss2.png'), width = 12, height = 10, dpi = 300)


#trying without alpha because it should cancel out


#loop over all watersheds
median_productivity <- data.frame()
#make empty df with 4000 rows

for (i in 1:length(unique(ch20r$River_n))) {
  river <- unique(ch20r$River_n)[i]
  print(river)
  # print(river)
  #filter data for river
  river_data <- ch20r %>% filter(River_n == river)
  # #filter alpha for river
  # river_alpha <- alpha %>% select(ends_with(paste0("_alpha_j[",river,"]")))
  #b effect of forestry
  b <- posterior %>% select(starts_with("b_for_rv")) %>% select(ends_with(paste0("[",river,"]")))
  #cpd
  max_cpd <- max(river_data$sqrt.CPD.std)
  no_forestry <- min(ch20r$sqrt.CPD.std)
  
  # #calculate productivity
  if(i==1){
    productivity <- (exp(as.matrix(b)%*%
                           max_cpd)/exp(as.matrix(b)%*%
                                          no_forestry))*100 - 100}
  else{
    productivity <- cbind(productivity, (exp(as.matrix(b)%*%
                                               max_cpd)/exp(as.matrix(b)%*%
                                                              no_forestry))*100 - 100)
  }
  median_productivity <- rbind(median_productivity, 
                               data.frame(river = river, median_productivity = median(productivity[,i])))
  # river_data <- river_data %>% mutate(productivity = exp(river_alpha %*% t(as.matrix(river_data$sqrt.CPD.std))) / (1 + exp(river_alpha %*% t(as.matrix(river_data$sqrt.CPD.std))))
  # #add productivity to df
  # df <- rbind(df, river_data)
}

#change the colnames of productivity to river names, make df long

colnames(productivity) <- unique(ch20r$River_n)
productivity <- data.frame(productivity) %>% 
  mutate(sample = rownames(productivity))

productivity_long <- data.frame(productivity) %>% 
  pivot_longer(everything(), names_to = 'River_n', 
               values_to = 'productivity_loss') %>% 
  mutate(River_n = as.numeric(substr(River_n, 2, nchar(River_n)))) %>% 
  left_join(ch20r %>% select(River_n, River, CU) %>% 
              distinct(), 
            by = 'River_n', relationship = "many-to-one")

#looks like the dataset now has more rivers. need to update posterior 
#to proceed. Going to run the models on the server
