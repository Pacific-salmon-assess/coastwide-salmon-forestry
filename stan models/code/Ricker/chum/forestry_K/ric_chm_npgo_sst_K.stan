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
  array[N] real npgo;
  array[N] real sst;
  array[J] int start_t; //start point brood year for each stock
  array[J] int start_y; // ragged start point for observations (N)
  array[J] int end_y;   // ragged end points for observations (N)
  array[J] real p_K_mean; //priors on equilibrium recruitment - based on obs R
  array[J] real p_K_sig; //priors on sigma for equilibrium recruitment - based on obs R
}

transformed data{
  vector[J] logK_pr;
  vector[J] logK_pr_sig;
  
  for (j in 1:J){
    // converting the standard deviation and mean of the prior to a log scale - ref: Quan & Zhang, Statistics in Medicine, 2003
    logK_pr_sig[j] = sqrt(log(1 +(p_K_sig[j]*p_K_sig[j])/(p_K_mean[j]*p_K_mean[j])));
    logK_pr[j] = log(p_K_mean[j]) - 0.5*logK_pr_sig[j]*logK_pr_sig[j];
  }
}



parameters{
  real alpha0; //initial global productivity
  vector[C] z_a_cu; //CU-specific z-score deviation
  vector[J] z_a_rv; //River-specific  z-score deviation
  real<lower=0> sigma_a_cu; //variance in long-term productivity among CUs
  real<lower=0> sigma_a_rv; //variance in long-term productivity among rivers
  
  vector<lower=0>[J] K; //carrying capacity
  
  //covariate effects
  real b_for; //global (across stock) mean effect of forestry metrics
  real b_npgo; //global (across stock) mean effect of npgo
  real b_sst; //global (across stock) mean effect of sst
  
  vector[C] z_for_cu; //CU-specific forestry z-score deviation
  vector[J] z_for_rv; //River-specific forestry z-score deviation
  real<lower=0> sigma_for_cu; //variance in CU-level ECA effect
  real<lower=0> sigma_for_rv; //variance in river-level ECA effect
  
  vector[C] z_npgo_cu; //CU-specific npgo z-score deviation
  vector[J] z_npgo_rv; //River-specific npgo z-score deviation
  real<lower=0> sigma_npgo_cu; //variance in CU-level npgo effect
  real<lower=0> sigma_npgo_rv; //variance in river-level npgo effect
  
  vector[C] z_sst_cu; //CU-specific sst z-score deviation
  vector[J] z_sst_rv; //River-specific sst z-score deviation
  real<lower=0> sigma_sst_cu; //variance in CU-level sst effect
  real<lower=0> sigma_sst_rv; //variance in river-level sst effect

  //variance components
  real<lower=0> mu_sigma; ///mean sigma among all stocks
  vector[C] z_sig_cu; //persistent CU-level differences in productivity
  vector[J] z_sig_rv; //persistent river-level differences in productivity
  real<lower=0> sd_sigma; ///variance in sigma within CUs (pooled among CUs)
  real<lower=0> sd_sigma_cu; ///variance in sigma among CUs  
  vector<lower = -1,upper=1>[J] rho; //autocorrelation parameter
}


transformed parameters{
  //estimated parameters
  vector[C] alpha_cu; //persistent CU-level differences in productivity
  vector[J] alpha_j; //persistent river-level differences in productivity
  vector[C] b_for_cu; //CU-specific forestry effect
  vector[J] b_for_rv; //River-specific forestry effect
  vector[C] b_npgo_cu; //CU-specific npgo effect
  vector[J] b_npgo_rv; //River-specific npgo effect
  vector[C] b_sst_cu; //CU-specific sst effect
  vector[J] b_sst_rv; //River-specific sst effect
  vector<lower=0>[C] cu_sigma; ///CU-level sigma
  vector<lower=0>[J] sigma; ///stock-level sigma
  vector<lower=0>[J] sigmaAR; ///stock-level sigma
  //productivity residuals through time
  vector[N] e_t; //stock residual productivity at time t
  vector[N] mu1; //initial expectation at each time for each stock
  vector[N] mu2; //autocorr. adjusted expectation at each time for each stock
  vector<lower=0>[N] K_for; //time varying K with forestry
  // vector[N] alpha_npgo_sst; //time varying alpha with npgo and sst
  vector<lower=0>[N] b; 
  vector[N] log_K_for; //log transformed K with forestry
  vector[J] log_K; //log transformed K
  
  alpha_cu = alpha0 + sigma_a_cu*z_a_cu; //non-centered estimate for CU time-invariant productivity
  alpha_j = alpha_cu[C_i] + sigma_a_rv*z_a_rv; //non-centered estimate for river time-invariant productivity
  
  b_for_cu = b_for + sigma_for_cu*z_for_cu; //non-centered CU-varying estimate for forestry effects 
  b_for_rv = b_for_cu[C_i] + sigma_for_rv*z_for_rv; //non-centered CU-varying estimate for forestry effects
  
  b_npgo_cu = b_npgo + sigma_npgo_cu*z_npgo_cu; //non-centered CU-varying estimate for npgo effects
  b_npgo_rv = b_npgo_cu[C_i] + sigma_npgo_rv*z_npgo_rv; //non-centered CU-varying estimate for npgo effects
  
  b_sst_cu = b_sst + sigma_sst_cu*z_sst_cu; //non-centered CU-varying estimate for sst effects
  b_sst_rv = b_sst_cu[C_i] + sigma_sst_rv*z_sst_rv; //non-centered CU-varying estimate for sst effects

  cu_sigma = mu_sigma + sd_sigma_cu*z_sig_cu; //non-centered CU-varying estimate for sigma
  sigma = cu_sigma[C_i] + sd_sigma*z_sig_rv; //non-centered CU-varying estimate for sigma
  
  //residual deviations in productivity
  for(j in 1:J){//for every stock
    // alpha_npgo_sst[start_y[j]] = alpha_j[j] + b_npgo_rv[j]*npgo[start_y[j]] + b_sst_rv[j]*sst[start_y[j]]; //adjust productivity based on  npgo, sst
    log_K[j] = log(K[j]); //log transform the carrying capacity
    log_K_for[start_y[j]] = log_K[j] + b_for_rv[j]*forest_loss[start_y[j]]; //log transform the carrying capacity adjusted for forest loss
    K_for[start_y[j]]=exp(log_K_for[start_y[j]]); 
    b[start_y[j]] = 1/K_for[start_y[j]];
    mu1[start_y[j]]= alpha_j[j]*(1-b[start_y[j]]*S[start_y[j]])+ b_npgo_rv[j]*npgo[start_y[j]] + b_sst_rv[j]*sst[start_y[j]]; 
    e_t[start_y[j]] = R_S[start_y[j]] - mu1[start_y[j]]; //first deviate for stock j
    
    for(t in (start_y[j]+1):(end_y[j])){ //adjust expectation based on autocorrelation
      // alpha_npgo_sst[t] = alpha_j[j] + b_npgo_rv[j]*npgo[t] + b_sst_rv[j]*sst[t]; //adjust productivity based on forest loss
      log_K_for[t] = log_K[j] + b_for_rv[j]*forest_loss[t]; //log transform the carrying capacity adjusted for forest loss
      K_for[t]=exp(log_K_for[t]); //exponentiate to get the carrying capacity
      b[t] = 1/K_for[t];
      mu2[t]  = alpha_j[j]*(1 - b[t]*S[t]) + b_npgo_rv[j]*npgo[t] + b_sst_rv[j]*sst[t] + rho[j]^(ii[t]-ii[t-1])*e_t[t-1]; //adjust expectation based on previous deviate - rho is raised to the power of the number of time steps (in years) between observations
      e_t[t] = R_S[t] - (mu2[t]-(rho[j]^(ii[t]-ii[t-1]))*e_t[t-1]);  //residual for stock j at time t
    }
  }
  sigmaAR = sigma.*sqrt(1-rho^2); //correct stock sigma for autocorrelation (rho)	
}

model{
  //priors
  alpha0 ~ normal(1.5,2); //global intrinsic productivity for all stocks at time t=1;
  z_a_cu ~ std_normal(); //std normal prior for CU-deviations
  z_a_rv ~ std_normal(); //std normal prior for River-deviations
    
  sigma_a_cu ~  normal(0,1); //variance in CU time-invariant productivities
  sigma_a_rv ~  normal(0,1); //variance in river-level time-invariant productivities
 
   //recruitment capacity for each stock - fit individually with weakly informative priors based on maximum observed recruitment
  for(j in 1:J) K[j] ~ lognormal(logK_pr[j],logK_pr_sig[j]); //stock-specific recruitment capacity
 
 //covariate effects
  b_for ~ normal(0,1); //standard normal prior for the effect of forest loss
  b_npgo ~ normal(0,1); //standard normal prior for the effect of npgo
  b_sst ~ normal(0,1); //standard normal prior for the effect of sst
  
  z_for_cu ~ std_normal(); //std normal prior for CU-deviations
  z_for_rv ~ std_normal(); //std normal prior for River-deviations
   
  //hierarchical variances
  sigma_for_cu ~ normal(0,1); //variance in stock-level forest loss effects
  sigma_for_rv ~ normal(0,1); //variance in stock-level forest loss effects
    
  z_npgo_cu ~ std_normal(); //std normal prior for CU-deviations
  z_npgo_rv ~ std_normal(); //std normal prior for River-deviations
  
  //hierarchical variances
  sigma_npgo_cu ~ normal(0,1); //variance in stock-level npgo effects
  sigma_npgo_rv ~ normal(0,1); //variance in stock-level npgo effects
  
  z_sst_cu ~ std_normal(); //std normal prior for CU-deviations
  z_sst_rv ~ std_normal(); //std normal prior for River-deviations
  
  //hierarchical variances
  sigma_sst_cu ~ normal(0,1); //variance in stock-level sst effects
  sigma_sst_rv ~ normal(0,1); //variance in stock-level sst effects
  
  
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

