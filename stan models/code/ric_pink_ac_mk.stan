data{
  int N;//number of observations
  int L; //length covered by all time-series
  int C; //number of CUs
  int R; //number of rivers
  int J;//number of river-level stocks
  array[C] int B_i; //Broodline index for each CU
  array [J] int C_i; //CU index for each stock
  array [R] int C_r; //CU index for each stock
  array [J] int R_i; //River index for each stock
  array [J] int BL; //index of odd (1) or even (2) broodines
  array [N] int ii; //index of brood years
  vector[N] R_S; //matrix of productivity among stocks - log(recruits per spawner)
  vector[N] S; //vector of spawners
  vector[N] forest_loss; //vector of stock-specific forest loss through time
  array [J] int start_t; //observations per stock
  array [J] int start_y;       // ragged start point for observations (N)
  array [J] int end_y;         // ragged end points for observations (N)
  array[J] real pSmax_mean; //priors on smax - based on observed spawner abundance
  array[J] real pSmax_sig; //priors on sigma for smax - based on observed spawner abundance
}
transformed data{
vector[J] logsmax_pr;
vector[J] logsmax_pr_sig;

for(t in 1:J){
logsmax_pr_sig[t]=sqrt(log(1+((pSmax_sig[t])^2/(pSmax_mean[t])^2))); //this converts sigma on the untransformed scale to a log scale
logsmax_pr[t]=log(pSmax_mean[t])-0.5*logsmax_pr_sig[t]^2; //convert smax prior to per capita slope - transform to log scale with bias correction
}

}
parameters{
  //initial productivity for each CU
  real alpha_01; //initial global productivity - even lineage
  real alpha_02; //initial global productivity - odd lineage
  vector[C] z_a_cu; //CU-specific z-score deviation
  vector[J] z_a_j; //River-specific  z-score deviation
  real<lower=0> sigma_a_cu; //variance in long-term productivity among CUs
  real<lower=0> sigma_a_j; //variance in long-term productivity among stocks
  vector<lower=0>[2] sigma_a_t; //temporal variance in coastwide productivity term - separate for odd/even broodlines
  matrix[2,L-1] a_dev; //stock-level deviations in year-to-year productivity



//covariate effects
real b_for; //global (across stock) mean effect of forestry metrics
 vector[C] z_for_cu; //CU-specific forestry z-score deviation
 vector[R] z_for_rv; //River-specific forestry z-score deviation
 real<lower=0> sigma_for_cu; //variance in CU-level ECA effect
 real<lower=0> sigma_for_rv; //variance in river-level ECA effect

 //variance components
 real<lower=0> mu_sigma; ///mean sigma among all stocks
 vector[C] z_sig_cu; //persistent CU-level differences in productivity
 vector[J] z_sig_j; //persistent river-level differences in productivity
 real<lower=0> sd_sigma; ///variance in sigma within CUs (pooled among CUs)
 real<lower=0> sd_sigma_cu; ///variance in sigma among CUs
 
 vector<lower = -1,upper=1>[J] rho; //autocorrelation parameter
 
 vector<lower=0>[J] Smax;
}
transformed parameters{
  matrix[2,L] alpha_t; //stock productivity over time
	
  //productivity residuals through time
  vector[N] e_t; //stock residual productivity at time t
  vector[N] mu1; //initial expectation at each time for each stock
  vector[N] mu2; //autocorr. adjusted expectation at each time for each stock
//Cushing density-dependence
  vector<lower=0>[J] b; //Cushing density term 
 //Transform non-centered parameters into parameter estimates: 
 vector[C] alpha_cu; //persistent CU-level differences in productivity
vector[J] alpha_j; //persistent river-level differences in productivity
vector[C] b_for_cu; //CU-specific forestry effect
vector[R] b_for_rv; //River-specific forestry effect
vector<lower=0>[C] cu_sigma; ///CU-level sigma
vector<lower = 0>[J] sigma; ///stock-level sigma
vector<lower = 0>[J] sigmaAR; ///stock-level sigma

alpha_cu = sigma_a_cu*z_a_cu; //non-centered estimate for CU time-invariant productivity
alpha_j = alpha_cu[C_i]+ sigma_a_j*z_a_j; //non-centered estimate for River time-invariant productivity

b_for_cu = b_for + sigma_for_cu*z_for_cu; //non-centered CU-varying estimate for forestry effects 
b_for_rv = b_for_cu[C_r] + sigma_for_rv*z_for_rv; //non-centered CU-varying estimate for forestry effects

cu_sigma = mu_sigma + sd_sigma_cu*z_sig_cu; //non-centered CU-varying estimate for sigma
sigma = cu_sigma[C_i] + sd_sigma*z_sig_j; //non-centered CU-varying estimate for sigma


 //initial productivities
  alpha_t[1,1] = alpha_01;
  alpha_t[2,1] = alpha_02;
  
  for(t in 2:L){
   alpha_t[1,t] = alpha_t[1,t-1] + a_dev[1,t-1]*sigma_a_t[1]; //coastwide time-varying productivity
   alpha_t[2,t] = alpha_t[2,t-1] + a_dev[2,t-1]*sigma_a_t[2]; //coastwide time-varying productivity
  }

//residual productivity deviations
   for(j in 1:J){ //for every stock
   
   b[j]=1/Smax[j];
   mu1[start_y[j]]=alpha_t[BL[j],start_t[j]]+alpha_j[j]-b[j]*S[start_y[j]]+b_for_rv[R_i[j]]*forest_loss[start_y[j]];
   e_t[start_y[j]] = R_S[start_y[j]] - mu1[start_y[j]]; //first deviate for stock j
   
	for(t in (start_y[j]+1):(end_y[j])){ //adjust expectation based on autocorrelation
	  mu2[t]  = alpha_t[BL[j],ii[t]]+alpha_j[j]-b[j]*S[t]+b_for_rv[R_i[j]]*forest_loss[t];
	  e_t[t] = R_S[t] - (mu2[t]-(rho[j]^(ii[t]-ii[t-1]))*e_t[t-1]);   //residual for stock j at time t
	
	}
  }
  sigmaAR = sigma.*sqrt(1-rho^2); //correct stock sigma for autocorrelation (rho)	

 }  
model{
  //priors
  //productivity
  alpha_01 ~ normal(1,5); //global intrinsic productivity for all stocks at time t=1;
  alpha_02 ~ normal(1,5); //global intrinsic productivity for all stocks at time t=1;
 
  z_a_cu ~ std_normal(); //CU-specific deviations in intrinsic productivity
 z_a_j ~ std_normal(); //river-specific static productivity adjustment (time-invariant)
  
  sigma_a_cu ~ normal(0, 1); //variance in CU time-invariant productivities
  sigma_a_j ~ normal(0, 1); //variance in river-level time-invariant productivities
  sigma_a_t ~ normal(0, 1); //temporal variance in time-varying global productivity
  
  for(i in 1:2) a_dev[i,] ~ std_normal(); //z-scores for productivity changes in each time-step
  
   //recruitment capacity for each stock - fit individually with weakly informative priors based on maximum observed recruitment
  for(j in 1:J) Smax[j] ~ lognormal(logsmax_pr[j],logsmax_pr_sig[j]); //stock-specific capacity
  
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
  z_sig_j ~ std_normal();
  sd_sigma_cu ~ normal(0,1);
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
