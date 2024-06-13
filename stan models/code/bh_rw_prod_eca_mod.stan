data{
  int<lower=1> N;//number of observations
  int L; //length covered by all time-series
  int C; //number of CUs
  int J;//number of river-level stocks
  int C_i[J]; //CU index for each stock
  int C_j[N]; //index of CUs
  int C_ii[N]; //index of CU-year combinations
  int J_i[N]; //index of stocks
  int ii[N]; //index of brood years
  vector[N] R_S; //matrix of productivity among stocks - log(recruits per spawner)
  vector[N] S; //design matrix of spawners
  vector[N] ECA; //design matrix of stock-specific ECA through time
  int<lower=0> start_y[J];       // ragged start point for observations (N)
  int<lower=0> end_y[J];         // ragged end points for observations (N)
  vector[J] pRk_mean; //priors on equilibrium recruitment - based on obs R
  vector[J] pRk_sig; //priors on sigma for equilibrium recruitment - based on obs R
}
transformed data{
vector[J] logRk_pr;
vector[J] logRk_pr_sig;

for(t in 1:J){
logRk_pr_sig[t]=sqrt(log(1+((pRk_sig[t])*(pRk_sig[t]))/((pRk_mean[t])*(pRk_mean[t])))); //this converts sigma on the untransformed scale to a log scale
logRk_pr[t]=log(pRk_mean[t])-0.5*logRk_pr_sig[t]*logRk_pr_sig[t]; //convert smax prior to per capita slope - transform to log scale with bias correction
}

}
parameters{
  //initial productivity for each CU
  real alpha_0; //initial global productivity
  vector[C] alpha_cu; //persistent CU-level differences in productivity
  vector[J] alpha_j; //persistent river-level differences in productivity
  real<lower=0> sigma_a_cu; //variance in long-term productivity among CUs
  real<lower=0> sigma_a_rv; //variance in long-term productivity among rivers
  real<lower=0> sigma_a_t; //temporal variance in coastwide productivity term
  vector[L-1] a_dev; //stock-level deviations in year-to-year productivity

//BH eq. recruitment
  vector<lower=0>[J] Rk; //equilibrium recruitment for BH

//covariate effects
real b_ECA; //global (across stock) mean effect of ECA
vector[C] b_ECA_cu; //CU-specific ECA effect
real<lower=0> sd_ECA; //variance in CU-level ECA effect

 //variance components
 real<lower=0> mu_sigma; ///mean sigma among all stocks
 vector<lower=0>[C] cu_sigma; ///CU-level sigma
 real<lower=0> sd_sigma; ///variance in sigma within CUs (pooled among CUs)
 real<lower=0> sd_sigma_cu; ///variance in sigma among CUs  
 vector<lower = 0>[J] sigma; ///stock-level sigma
 vector<lower = -1,upper=1>[J] rho; //autocorrelation parameter
}
transformed parameters{
  vector[L] alpha_t; //stock productivity over time
 
 //initial productivities
  alpha_t[1] = alpha_0;
  
  for(t in 2:L){
   alpha_t[t] = alpha_t[t-1] + a_dev[t-1]*sigma_a_t; //coastwide time-varying productivity
  }	
 }  
model{
  //priors
  //productivity
  alpha_0 ~ normal(1.5,5); //global intrinsic productivity for all stocks at time t=1;
  alpha_cu ~ normal(0,sd_a_cu); //CU-specific deviations in intrinsic productivity
  alpha_j ~ normal(alpha_cu[C_i],sd_a_rv); //river-specific static productivity adjustment (time-invariant)
  
  sd_a_cu ~ normal(0, 0.5); //variance in CU time-invariant productivities
  sd_a_rv ~ normal(0, 0.5); //variance in river-level time-invariant productivities
  sigma_a_t ~ normal(0, 0.5); //temporal variance in time-varying global productivity
  a_dev ~ std_normal(); //z-scores for productivity changes in each time-step
  
   //recruitment capacity for each stock - fit individually with weakly informative priors based on maximum observed recruitment
  for(j in 1:J) Rk[j] ~ lognormal(logRk_pr[j],logRk_pr_sig[j]); //stock-specific recruitment capacity
 
  //covariate effects
  b_ECA ~ normal(0,1); //standard normal prior for the effect of ECA
  b_ECA_cu ~ normal(b_ECA,sd_ECA); //CU-specific ECA effects
  
  //hierarchical variances
  sd_ECA ~ normal(0,0.5); //variance in stock-level ECA effects
  
  //variance terms
  mu_sigma ~ normal(0.5,0.5);
  cu_sigma ~ normal(mu_sigma,sd_sigma_cu);
  sd_sigma_cu ~ normal(0,0.5);
  sd_sigma ~ normal(0,0.5);
  sigma ~ normal(cu_sigma[C_i], sd_sigma);

 
  //likelihood sampling:

  
  R_S ~ normal(to_vector(alpha_t)[ii]+alpha_j[J_i] - log(1+(exp(alpha_t[ii]+alpha_j[J_i])./Rk[J_i])).*S+ b_ECA_cu[C_j].*ECA, sigma[J_i]); 
}

