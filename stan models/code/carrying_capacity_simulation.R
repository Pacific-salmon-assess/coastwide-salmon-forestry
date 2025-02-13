# code to simulate data from ricker model with varying carrying capacities
# plotting log R/S vs S to test for the effect of carrying capacity on the slope

# load libraries

library(ggplot2)
library(tidyverse)


# define the ricker model

ricker <- function(N, r, K) {
  N * exp(r * (1 - N / K))
}

# ricker(10, 0.5, 100)

#simulate data

set.seed(123)

# parameters

r <- 0.5
K <- c(10, 50, 100, 500, 1000)
N0 <- 100
time <- 100

for (i in 1:length(K)) {
  N <- numeric(time)
  N[1] <- N0
  for (t in 2:time) {
    N[t] <- ricker(N[t - 1], r, K[i])
  }
  data <- data.frame(time = 1:time, N = N, K = K[i])
  if (i == 1) {
    df <- data
  } else {
    df <- rbind(df, data)
  }
}

#define R as the population size at time t+1

df$R <- c(df$N[-1], NA)
df$S <- df$N

# plot log R/S vs S

ggplot(df, aes(x = S, y = log(R/S))) +
  geom_point(alpha = 0.2) +
  # geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~K) +
  theme_minimal() +
  labs(x = "S", y = "log(R/S)") +
  ggtitle("Log R/S vs S for different carrying capacities")



