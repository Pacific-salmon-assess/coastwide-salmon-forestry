data{
  int<lower=1> N;//number of observations
  int L; //length covered by all time-series
  int C; //number of CUs
  int J;//number of river-level stocks
  array[J] int C_i; //CU index for each stock
  array[N] int ii; //index of brood years
  array[N] real R_S; //vector of productivity for all stocks - log(recruits per spawner)
  array[N] real S; //vector of spawners
  array[N] real forest_loss; //vector of watershed forest loss through time
  array[J] int start_t; //start point brood year for each stock
  array[J] int  start_y; // ragged start point for observations (N)
  array[J] int end_y;   // ragged end points for observations (N)
  // array[J] real pSmax_mean; //priors on smax - based on observed spawner abundance
  // array[J] real pSmax_sig; //priors on sigma for smax - based on observed spawner abundance
  array[J] real p_k_mean; //priors on k - based on long term mean of observed spawner abundance
  array[J] real p_k_sig; //priors on sigma for k - based on sig for observed spawner abundance
}
transformed data{
// vector[J] logsmax_pr;
// vector[J] logsmax_pr_sig;
vector[J] log_p_k;
vector[J] log_p_k_sig;

for(t in 1:J){
// logsmax_pr_sig[t]=sqrt(log(1+((pSmax_sig[t])^2/(pSmax_mean[t])^2))); //this converts sigma on the untransformed scale to a log scale
// logsmax_pr[t]=log(pSmax_mean[t])-0.5*logsmax_pr_sig[t]^2; //convert smax prior to per capita slope - transform to log scale with bias correction

log_p_k_sig[t]=sqrt(log(1+((p_k_sig[t])^2/(p_k_mean[t])^2))); //this converts sigma on the untransformed scale to a log scale
log_p_k[t]=log(p_k_mean[t])-0.5*log_p_k_sig[t]^2; //convert k prior to per capita slope - transform to log scale with bias correction


}

}
parameters{
  //initial productivity for each CU
  real alpha0; //initial global productivity
  vector[C] z_a_cu; //CU-specific z-score deviation
  vector[J] z_a_rv; //River-specific  z-score deviation
  real<lower=0> sigma_a_cu; //variance in long-term productivity among CUs
  vector<lower=0>[C] sigma_a_rv; //variance in long-term productivity among rivers

//BH eq. recruitment
  // vector<lower=0>[J] Smax; //spawners where recruits are maximized
  vector<lower=0>[J] k; //carrying capacity

//covariate effects
 real b_for; //global (across stock) mean effect of forestry metrics
 vector[C] z_for_cu; //CU-specific forestry z-score deviation
 vector[J] z_for_rv; //River-specific forestry z-score deviation
 real<lower=0> sigma_for_cu; //variance in CU-level ECA effect
 real<lower=0> sigma_for_rv; //variance in river-level ECA effect

 //variance components
 real<lower=0> mu_sigma; ///mean sigma among all stocks
 vector[C] z_sig_cu; //persistent CU-level differences in productivity
 vector[J] z_sig_rv; //persistent river-level differences in productivity
 real<lower=0> sd_sigma; ///variance in sigma within CUs (pooled among CUs)
 real<lower=0> sd_sigma_cu; ///variance in sigma among CUs  

 vector<lower = -1,upper=1>[J] rho; //autocorrelation parameter
}
transformed parameters{

  //parameters estimated from non-centered parameterization
vector[C] alpha_cu; //persistent CU-level differences in productivity
vector[J] alpha_j; //persistent river-level differences in productivity
vector[C] b_for_cu; //CU-specific forestry effect
vector[J] b_for_rv; //River-specific forestry effect
vector<lower=0>[C] cu_sigma; ///CU-level sigma
vector<lower = 0>[J] sigma; ///stock-level sigma
vector<lower = 0>[J] sigmaAR; ///stock-level sigma
vector<lower=0>[N] b; //per capita density dependence term
	
  //productivity residuals through time
vector[N] e_t; //stock residual productivity at time t
vector[N] mu1; //initial expectation at each time for each stock
vector[N] mu2; //autocorr. adjusted expectation at each time for each stock
vector<lower=0>[N] k_rv; //time varying carrying capcity
vector[N] log_k_rv; //log of time varying carrying capcity for each river and year
vector[J] log_k; //log of time varying carrying capcity for each river
 

alpha_cu = alpha0+sigma_a_cu*z_a_cu; //non-centered estimate for CU time-invariant productivity
alpha_j = alpha_cu[C_i] + sigma_a_rv[C_i].*z_a_rv; //non-centered estimate for River time-invariant productivity

b_for_cu = b_for + sigma_for_cu*z_for_cu; //non-centered CU-varying estimate for forestry effects 
b_for_rv = b_for_cu[C_i] + sigma_for_rv*z_for_rv; //non-centered CU-varying estimate for forestry effects

cu_sigma = mu_sigma + sd_sigma_cu*z_sig_cu; //non-centered CU-varying estimate for sigma
sigma = cu_sigma[C_i] + sd_sigma*z_sig_rv; //non-centered CU-varying estimate for sigma
 
//residual productivity deviations
   for(j in 1:J){ //for every stock
   
    log_k[j] = log(k[j]);
    log_k_rv[start_y[j]] = log_k[j] + b_for_rv[j]*forest_loss[start_y[j]];
    k_rv[start_y[j]] = exp(log_k_rv[start_y[j]]);
    b[start_y[j]]=1/k_rv[start_y[j]];
    mu1[start_y[j]]=alpha_j[j]*(1-b[j]*S[start_y[j]]); //adjust expectation based on previous deviate - rho is raised to the power of the number of time steps (in years) between observations
    e_t[start_y[j]] = R_S[start_y[j]] - mu1[start_y[j]]; //first deviate for stock j

    for(t in (start_y[j]+1):(end_y[j])){ //adjust expectation based on autocorrelation
      log_k_rv[t] = log_k[j] + b_for_rv[j]*forest_loss[t];
      k_rv[t] = exp(log_k_rv[t]);
      b[t]=1/k_rv[t];
      mu2[t]  = alpha_j[j]*(1-b[t]*S[t])+ rho[j]^(ii[t]-ii[t-1])*e_t[t-1]; //adjust expectation based on previous deviate - rho is raised to the power of the number of time steps (in years) between observations
      e_t[t] = R_S[t] - (mu2[t]-(rho[j]^(ii[t]-ii[t-1]))*e_t[t-1]);  //residual for stock j at time t
	}
  }
  sigmaAR = sigma.*sqrt(1-rho^2); //correct stock sigma for autocorrelation (rho)	
 }  
model{
  //priors
  //productivity
   alpha0 ~ normal(1.5,2); //global intrinsic productivity for all stocks at time t=1;
   z_a_cu ~ std_normal(); //std normal prior for CU-deviations
   z_a_rv ~ std_normal(); //std normal prior for River-deviations
    
  sigma_a_cu ~  normal(0,1); //variance in CU time-invariant productivities
  sigma_a_rv ~  normal(0,1); //variance in river-level time-invariant productivities
  
   //recruitment capacity for each stock - fit individually with weakly informative priors based on maximum observed recruitment
  // for(j in 1:J) Smax[j] ~ lognormal(logsmax_pr[j],logsmax_pr_sig[j]); //stock-specific capacity
  for(j in 1:J){
  k[j] ~ lognormal(log_p_k[j],log_p_k_sig[j]); //stock-specific capacity
  
  }
 
 //covariate effects
  b_for ~ normal(0,1); //standard normal prior for the effect of forest loss
  z_for_cu ~ std_normal(); //std normal prior for CU-deviations
  z_for_rv ~ std_normal(); //std normal prior for River-deviations
   
  //hierarchical variances
 sigma_for_cu ~ normal(0,1); //variance in stock-level forest loss effects
  sigma_for_rv ~ normal(0,1); //variance in stock-level forest loss effects
    
  //variance terms
  mu_sigma ~ normal(1,1);
  z_sig_cu ~ std_normal();
  z_sig_rv ~ std_normal();
  sd_sigma_cu ~normal(0,1);
  sd_sigma ~ normal(0,1);
  
  rho ~ uniform(-1,1); //prior for autocorrelation

  //likelihood sampling:
  for(j in 1:J){
  R_S[start_y[j]]~ normal(mu1[start_y[j]],sigma[j]); //initial fit for each stock
  
  R_S[(start_y[j]+1):end_y[j]] ~ normal(mu2[(start_y[j]+1):end_y[j]], sigmaAR[j]); //subsequent samples including autocorrelation
  }
}
generated quantities{
vector[N] log_lik; //pointwise log likelihoods

for(j in 1:J){
 log_lik[start_y[j]]=normal_lpdf(R_S[start_y[j]]|mu1[start_y[j]],sigma[j]); //pointwise log likelihood calculation - initial estimates
 
 for(t in (start_y[j]+1):(end_y[j])){
 log_lik[t]=normal_lpdf(R_S[t]|mu2[t], sigmaAR[j]); //pointwise log likelihood calculation
 }
}
}
