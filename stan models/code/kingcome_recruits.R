library(here)
library(tidyverse)
library(ggplot2)
library(zoo)

posterior <- read.csv(here(
  'stan models','outs','fits',
  'posterior','fit4bh_chm_cpd_sqrt_mk.csv'),check.names=F)

ch20r <- read.csv("./origional-ecofish-data-models/Data/Processed/chum_SR_20_hat_yr_reduced_VRI90.csv")
ch20r$sqrt.ECA=sqrt(ch20r$ECA_age_proxy_forested_only)
ch20r$sqrt.CPD=sqrt(ch20r$disturbedarea_prct_cs)
ch20r$sqrt.CPD.std=(ch20r$sqrt.CPD-mean(ch20r$sqrt.CPD))/sd(ch20r$sqrt.CPD)

#add numeric values for each river

ch20r <- ch20r %>% 
  mutate(River_n = as.numeric(as.factor(River)))

chum_kingcome <- ch20r %>% 
  filter(River == "KINGCOME RIVER") %>% 
  mutate(no_forestry = min(ch20r$sqrt.CPD.std))

#subset all columns in posterior whose name starts with "alpha"
#then taken median of all columns keeping all the rows

alpha <- posterior %>% 
  select(starts_with("alpha_t")) %>% 
  rowwise() %>%
  mutate(alpha_median = median(c_across(everything()),na.rm = T))

rkr <- posterior %>% 
  select(starts_with("Rk[135]"))


b_eca <- posterior %>% 
  select(b_ECA)

b_eca_cu <- posterior %>% 
  select("b_ECA_cu[21]")

predicted_recruits_kingcome <- exp(matrix(alpha$alpha_median, 
                                          nrow = nrow(alpha),
                                          ncol = nrow(chum_kingcome)) + 
                                     as.matrix(b_eca$b_ECA)%*%
                                     chum_kingcome$sqrt.CPD.std)%*%diag(chum_kingcome$Spawners)/(
                                       1 + as.matrix(exp(alpha$alpha_median)/rkr)%*%chum_kingcome$Spawners
                                    ) 

predicted_recruits_kingcome_wo_forestry <- exp(matrix(alpha$alpha_median, 
                                          nrow = nrow(alpha),
                                          ncol = nrow(chum_kingcome)) + 
                                     as.matrix(b_eca$b_ECA)%*%
                                     chum_kingcome$no_forestry)%*%diag(chum_kingcome$Spawners)/(
                                       1 + as.matrix(exp(alpha$alpha_median)/rkr)%*%chum_kingcome$Spawners
                                     ) 

#column names should be the BroodYear from chum_kingcome
predicted_recruits_kingcome <- as.data.frame(predicted_recruits_kingcome)
colnames(predicted_recruits_kingcome) <- chum_kingcome$BroodYear


predicted_recruits_kingcome_wo_forestry <- as.data.frame(predicted_recruits_kingcome_wo_forestry)
  
colnames(predicted_recruits_kingcome_wo_forestry) <- chum_kingcome$BroodYear

#plotting the predicted recruits for kingcome river, lines for each year, 3000 lines total

#make the data frame long for ggplot

predicted_recruits_kingcome_long <- predicted_recruits_kingcome %>% 
  mutate(row = row_number()) %>% 
  pivot_longer(cols = -row, names_to = "BroodYear", values_to = "Recruits_predicted") %>% 
  mutate(BroodYear = as.numeric(BroodYear))
  
predicted_recruits_kingcome_wo_forestry_long <- predicted_recruits_kingcome_wo_forestry %>%
  mutate(row = row_number()) %>% 
  pivot_longer(cols = -row, names_to = "BroodYear", values_to = "Recruits_predicted_no_forestry") %>% 
  mutate(BroodYear = as.numeric(BroodYear))

predicted_recruits_kingcome_combined <- predicted_recruits_kingcome_long %>% 
  left_join(predicted_recruits_kingcome_wo_forestry_long, by = c("row" = "row", "BroodYear" = "BroodYear"))


#add actual recruit numbers from chum_kingcome to the predicted_recruits_kingcome_long data frame, join by year

predicted_recruits_kingcome_long <- predicted_recruits_kingcome_long %>% 
  left_join(chum_kingcome %>% select(Recruits,BroodYear), by = c("BroodYear" = "BroodYear"))


ggplot(predicted_recruits_kingcome_combined, aes(x = BroodYear, y = Recruits_predicted, group = row)) +
  geom_line(color = "#aba0b1", alpha = 0.01) +
  geom_line(aes(y = Recruits_predicted_no_forestry),color = "#a0b1a1" , alpha = 0.01) +
  geom_line(aes(y = Recruits), color = "slategray", linewidth = 1.5, alpha = 0.5) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Predicted Recruits for Kingcome River",
       x = "Brood year",
       y = "Recruits")

ggsave(here('figures','predicted_recruits_kingcome.jpeg'), width = 10, height = 6)                              


#try with time varying alpha and river specific alpha

alpha <- posterior %>% 
  select(starts_with("alpha_t"),"alpha_j[135]") 


colnames(alpha)[-60] <- paste0("alpha_t_",as.character(as.numeric(ifelse(substr(colnames(alpha)[-60],10,10)=="]", 
                                     substr(colnames(alpha)[-60],9,9), 
                                     substr(colnames(alpha)[-60],9,10))) + 1953))
  
  
#calculate alpha_1953 = alpha_t_1953 - alpha_j[135], alpha_1954 = alpha_t_1954 - alpha_j[135], etc 

alpha_t <- alpha %>% 
  select(-starts_with("alpha_j"))

alpha_j <- alpha %>% 
  select(starts_with("alpha_j"))

alpha <- alpha_t %>% 
  mutate(across(starts_with("alpha_t"), ~ . + alpha_j$"alpha_j[135]"))
  

#keep only the columns of alpha with year that is in the BroodYear of chum_kingcome

alpha_kingcome <- alpha %>% 
  select(ends_with(as.character(chum_kingcome$BroodYear)))


predicted_recruits_kingcome_2 <- exp(as.matrix(alpha_kingcome) + 
                                       as.matrix(b_eca$b_ECA)%*%
                                     chum_kingcome$sqrt.CPD.std)%*%diag(chum_kingcome$Spawners)/(
                                       1 + as.matrix(exp(alpha_kingcome)/
                                                       matrix(rkr$`Rk[135]`,
                                                              nrow = nrow(rkr),
                                                              ncol = ncol(alpha_kingcome)))%*%diag(chum_kingcome$Spawners)
                                     )

predicted_recruits_kingcome_no_forestry_2 <-  exp(as.matrix(alpha_kingcome) + 
                                                                                    as.matrix(b_eca$b_ECA)%*%
                                                                                    chum_kingcome$no_forestry)%*%diag(chum_kingcome$Spawners)/(
                                                                                      1 + as.matrix(exp(alpha_kingcome)/
                                                                                                      matrix(rkr$`Rk[135]`,
                                                                                                             nrow = nrow(rkr),
                                                                                                             ncol = ncol(alpha_kingcome)))%*%diag(chum_kingcome$Spawners)
                                                                                    )


predicted_recruits_kingcome_2 <- as.data.frame(predicted_recruits_kingcome_2)
colnames(predicted_recruits_kingcome_2) <- chum_kingcome$BroodYear

predicted_recruits_kingcome_no_forestry_2 <- as.data.frame(predicted_recruits_kingcome_no_forestry_2)
colnames(predicted_recruits_kingcome_no_forestry_2) <- chum_kingcome$BroodYear


predicted_recruits_kingcome_long_2 <- predicted_recruits_kingcome_2 %>% 
  mutate(row = row_number()) %>% 
  pivot_longer(cols = -row, names_to = "BroodYear", values_to = "Recruits_predicted") %>% 
  mutate(BroodYear = as.numeric(BroodYear))


predicted_recruits_kingcome_no_forestry_long_2 <- predicted_recruits_kingcome_no_forestry_2 %>%
  mutate(row = row_number()) %>% 
  pivot_longer(cols = -row, names_to = "BroodYear", values_to = "Recruits_predicted_no_forestry") %>% 
  mutate(BroodYear = as.numeric(BroodYear))

predicted_recruits_kingcome_long_2 <- predicted_recruits_kingcome_long_2 %>%
  left_join(chum_kingcome %>% select(Recruits,BroodYear), by = c("BroodYear" = "BroodYear"))

predicted_recruits_kingcome_combined_2 <- predicted_recruits_kingcome_long_2 %>%
  left_join(predicted_recruits_kingcome_no_forestry_long_2, by = c("row" = "row", "BroodYear" = "BroodYear"))


ggplot(predicted_recruits_kingcome_combined_2, aes(x = BroodYear, y = Recruits_predicted, group = row)) +
  geom_line(color = "plum4", alpha = 0.01) +
  geom_line(aes(y = Recruits_predicted_no_forestry),color = "aquamarine4" , alpha = 0.01) +
  geom_line(aes(y = Recruits), color = "#6e706e", linewidth = 1.5, alpha = 0.5) +
  theme(legend.position = "right") +
  labs(title = "Predicted Recruits for Kingcome River",
       x = "Brood year",
       y = "Recruits")+
  theme_classic()

ggsave(here('figures','predicted_recruits_kingcome_time_varying_alpha.jpeg'), width = 10, height = 6)
png(filename=here('figures','predicted_recruits_kingcome_time_varying_alpha.png'))

ggplot(predicted_recruits_kingcome_combined_2, aes(x = BroodYear, y = Recruits_predicted, group = row)) +
  # geom_line(color = "plum4", alpha = 0.01) +
  geom_line(aes(y = Recruits_predicted_no_forestry),color = "aquamarine4" , alpha = 0.01) +
  geom_line(aes(y = Recruits), color = "#6e706e", linewidth = 1.5, alpha = 0.5) +
  theme(legend.position = "right") +
  labs(title = "Predicted Recruits for Kingcome River",
       x = "Brood year",
       y = "Recruits")+
  theme_classic()

ggsave(here('figures','predicted_recruits_kingcome_time_varying_alpha_wo_forestry.jpeg'), width = 10, height = 6)


ggplot(predicted_recruits_kingcome_combined_2, aes(x = BroodYear, y = Recruits_predicted, group = row)) +
  geom_line(color = "plum4", alpha = 0.01) +
  # geom_line(aes(y = Recruits_predicted_no_forestry),color = "aquamarine4" , alpha = 0.01) +
  geom_line(aes(y = Recruits), color = "#6e706e", linewidth = 1.5, alpha = 0.5) +
  theme(legend.position = "right") +
  labs(title = "Predicted Recruits for Kingcome River",
       x = "Brood year",
       y = "Recruits")+
  theme_classic()
ggsave(here('figures','predicted_recruits_kingcome_time_varying_alpha_w_forestry.jpeg'), width = 10, height = 6)


#check how many times BroodYear 2000 is in the BroodYear of predicted_recruits_kingcome_combined_2

predicted_recruits_kingcome_combined_2 %>% 
  filter(BroodYear == 1960) %>% 
  nrow()

#using b_eca[21]

predicted_recruits_kingcome_3 <- exp(as.matrix(alpha_kingcome) + 
                                       as.matrix(b_eca_cu$'b_ECA_cu[21]')%*%
                                       chum_kingcome$sqrt.CPD.std)%*%diag(chum_kingcome$Spawners)/(
                                         1 + as.matrix(exp(alpha_kingcome)/
                                                         matrix(rkr$`Rk[135]`,
                                                                nrow = nrow(rkr),
                                                                ncol = ncol(alpha_kingcome)))%*%diag(chum_kingcome$Spawners)
                                       )


predicted_recruits_kingcome_no_forestry_3 <-  exp(as.matrix(alpha_kingcome) + 
                                                    as.matrix(b_eca_cu$'b_ECA_cu[21]')%*%
                                                    chum_kingcome$no_forestry)%*%diag(chum_kingcome$Spawners)/(
                                                      1 + as.matrix(exp(alpha_kingcome)/
                                                                      matrix(rkr$`Rk[135]`,
                                                                             nrow = nrow(rkr),
                                                                             ncol = ncol(alpha_kingcome)))%*%diag(chum_kingcome$Spawners)
                                                    )

predicted_recruits_kingcome_3 <- as.data.frame(predicted_recruits_kingcome_3)
colnames(predicted_recruits_kingcome_3) <- chum_kingcome$BroodYear

predicted_recruits_kingcome_no_forestry_3 <- as.data.frame(predicted_recruits_kingcome_no_forestry_3)
colnames(predicted_recruits_kingcome_no_forestry_3) <- chum_kingcome$BroodYear


predicted_recruits_kingcome_long_3 <- predicted_recruits_kingcome_3 %>% 
  mutate(row = row_number()) %>% 
  pivot_longer(cols = -row, names_to = "BroodYear", values_to = "Recruits_predicted") %>% 
  mutate(BroodYear = as.numeric(BroodYear))


predicted_recruits_kingcome_no_forestry_long_3 <- predicted_recruits_kingcome_no_forestry_3 %>%
  mutate(row = row_number()) %>% 
  pivot_longer(cols = -row, names_to = "BroodYear", values_to = "Recruits_predicted_no_forestry") %>% 
  mutate(BroodYear = as.numeric(BroodYear))

predicted_recruits_kingcome_long_3 <- predicted_recruits_kingcome_long_3 %>%
  left_join(chum_kingcome %>% select(Recruits,BroodYear), by = c("BroodYear" = "BroodYear"))

predicted_recruits_kingcome_combined_3 <- predicted_recruits_kingcome_long_3 %>%
  left_join(predicted_recruits_kingcome_no_forestry_long_3, by = c("row" = "row", "BroodYear" = "BroodYear"))

plot1 <- ggplot(predicted_recruits_kingcome_combined_3, aes(x = BroodYear, y = Recruits_predicted, group = row)) +
  geom_line(color = "plum4", alpha = 0.01) +
  geom_line(aes(y = Recruits_predicted_no_forestry),color = "aquamarine4" , alpha = 0.01) +
  geom_line(aes(y = Recruits), color = "#6e706e", linewidth = 1.5, alpha = 0.5) +
  theme(legend.position = "right") +
  labs(title = "Predicted Recruits for Kingcome River",
       x = "Brood year",
       y = "Recruits")+
  theme_classic()


ggsave(here('figures','predicted_recruits_kingcome_time_varying_alpha_3.png'), 
       width = 10, height = 6)


#make new dataframe with 3000 rows for each BroodYear from 1955 to 2008 and add column for each replicate 1 to 3000

new_df <- data.frame(BroodYear = rep(1955:2008, each = 3000)) %>% 
  mutate(row = rep(1:3000, times = 54))

new_df_combined <- new_df %>%
  left_join(predicted_recruits_kingcome_combined_3, by = c("row" = "row", "BroodYear" = "BroodYear"))


ggplot(new_df_combined, aes(x = BroodYear, y = Recruits_predicted, group = row)) +
  geom_line(color = "plum4", alpha = 0.01) +
  geom_line(aes(y = Recruits_predicted_no_forestry),color = "aquamarine4" , alpha = 0.01) +
  geom_line(aes(y = Recruits), color = "#6e706e", linewidth = 1.5, alpha = 0.5) +
  theme(legend.position = "right") +
  labs(title = "Predicted Recruits for Kingcome River",
       x = "Brood year",
       y = "Recruits")+
  theme_classic()

#calculate difference between median Recruits_predicted and Recruits_predicted_no_forestry

diff <- predicted_recruits_kingcome_combined_3 %>% 
  select(BroodYear, Recruits_predicted, Recruits_predicted_no_forestry) %>%
  group_by(BroodYear) %>% 
  summarise(Recruits_predicted = median(Recruits_predicted), 
            Recruits_predicted_no_forestry = median(Recruits_predicted_no_forestry)) %>% 
  left_join(chum_kingcome %>% select(Recruits,BroodYear), by = c("BroodYear" = "BroodYear")) %>%
  mutate(diff_predictions = Recruits_predicted_no_forestry - Recruits_predicted,
         diff_observed = Recruits_predicted_no_forestry - Recruits)


diff_interpolated <- predicted_recruits_kingcome_combined_3 %>% 
  select(BroodYear, Recruits_predicted, Recruits_predicted_no_forestry) %>%
  group_by(BroodYear) %>% 
  summarise(Recruits_predicted = median(Recruits_predicted), 
            Recruits_predicted_no_forestry = median(Recruits_predicted_no_forestry)) %>% 
  left_join(chum_kingcome %>% select(Recruits,BroodYear), by = c("BroodYear" = "BroodYear")) %>%
  #add rows for missing Broodyears
  complete(BroodYear = 1955:2008) %>%
  #linear interpolation
  mutate(Recruits_predicted_inp = na.approx(Recruits_predicted),
         Recruits_predicted_no_forestry_inp = na.approx(Recruits_predicted_no_forestry),
         Recruits_inp = na.approx(Recruits)) %>%
  mutate(diff_predictions = Recruits_predicted_no_forestry_inp - Recruits_predicted_inp,
         diff_observed = Recruits_predicted_no_forestry_inp - Recruits_inp) %>% 
  mutate(cumulative_diff_predictions = cumsum(diff_predictions),
         cumulative_diff_observed = cumsum(diff_observed))

#plot the cumulative difference

ggplot(diff_interpolated, aes(x = BroodYear, y = cumulative_diff_predictions)) +
  geom_line(aes(color = "Predicted")) +
  geom_line(aes(y = cumulative_diff_observed, color ="Observed")) +
  scale_color_manual(values = c("Predicted" = "blue4", "Observed" = "green4"), name = "") +
  labs(x = "Brood year",
       y = "Cumulative difference in recruits",
       legend = "")+
  theme_classic()

ggsave(here('figures','cumulative_diff_kingcome_interpolated.png'), 
       width = 10, height = 6)


