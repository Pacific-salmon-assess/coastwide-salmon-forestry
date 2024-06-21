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
  int<lower=0> N_i[J]; //number of observations in each time-series
  int<lower=0> start_t[J]; //observations per stock
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
 vector[J] sigmaAR; //sigma - adjusted for rho
	
  //productivity residuals through time
  vector[N] e_t; //stock residual productivity at time t
  vector[N] mu1; //initial expectation at each time for each stock
  vector[N] mu2; //autocorr. adjusted expectation at each time for each stock
  
 
 //initial productivities
  alpha_t[1] = alpha_0;
  
  for(t in 2:L){
   alpha_t[t] = alpha_t[t-1] + a_dev[t-1]*sigma_a_t; //coastwide time-varying productivity
  }

//residual productivity deviations
   for(j in 1:J){ //for every stock
     	mu1[start_y[j]]=alpha_t[start_t[j]]+alpha_j[j]-log(1+(exp(alpha_t[start_t[j]]+alpha_j[j])/Rk[j])*S[start_y[j]])+b_ECA_cu[C_i[j]]*ECA[start_y[j]]; //expectation (unadjusted for autocorrelated residuals)
	
	e_t[start_y[j]] = R_S[start_y[j]] - mu1[start_y[j]]; //first deviate for stock j
	
	for(t in (start_y[j]+1):(end_y[j])){ //adjust expectation based on autocorrelation
	  mu2[t] = alpha_t[ii[t]]+alpha_j[j]-log(1+(exp(alpha_t[ii[t]]+alpha_j[j])/Rk[j])*S[t])+b_ECA_cu[C_i[j]]*ECA[t] + rho[j]*e_t[t-1]; //adjust expectation based on previous deviate - rho is raised to the power of the number of time steps (in years) between observations
         e_t[t] = R_S[t] - mu2[t];  //residual for stock j at time t
	}
  }
  sigmaAR = sigma.*sqrt(1-rho^2); //correct stock sigma for autocorrelation (rho)	
 }  
model{
  //priors
  //productivity
  alpha_0 ~ normal(1.5,5); //global intrinsic productivity for all stocks at time t=1;
  alpha_cu ~ normal(0,sigma_a_cu); //CU-specific deviations in intrinsic productivity
  alpha_j ~ normal(alpha_cu[C_i],sigma_a_rv); //river-specific static productivity adjustment (time-invariant)
  
  sigma_a_cu ~ normal(0, 0.5); //variance in CU time-invariant productivities
  sigma_a_rv ~ normal(0, 0.5); //variance in river-level time-invariant productivities
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
