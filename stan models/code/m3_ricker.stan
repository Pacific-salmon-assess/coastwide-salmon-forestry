//static, autocorrelated
//no pooling
data{
  int<lower=1> N;//number of observations
  int L; //length covered by all time-series
  int C; //number of CUs
  int J;//number of river-level stocks
  int N_i[J]; //number of observations in each time-series
  int C_i[J]; //CU index for each stock
  int ii[N]; //index of brood years
  vector[N] R_S; //matrix of productivity among stocks - log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners
  matrix[N,J] ECA; //design matrix of stock-specific ECA 
  matrix[N,J] Area; //design matrix of stock-specific watershed area
  matrix[N,J] ExA; //design matrix of stock-specific ECA x watershed area
  int<lower=0> start_y[J];       // ragged start point for observations (N)
  int<lower=0> end_y[J];         // ragged end points for observations (N)
  vector[J] pSmax_mean; //priors on smax - based on observed spawner abundance
  vector[J] pSmax_sig;
}
transformed data{
vector[J] logbeta_pr;
vector[J] logbeta_pr_sig;

for(t in 1:J){
logbeta_pr_sig[t]=sqrt(log(1+((1/pSmax_sig[t])*(1/pSmax_sig[t]))/((1/pSmax_mean[t])*(1/pSmax_mean[t])))); //this converts sigma on the untransformed scale to a log scale
logbeta_pr[t]=log(1/pSmax_mean[t])-0.5*logbeta_pr_sig[t]*logbeta_pr_sig[t]; //convert smax prior to per capita slope - transform to log scale with bias correction
}

}
parameters{
  //intrinsic productivity
  real<lower=0> alpha_0; ///mean alpha among all stocks
  vector<lower=0>[C] alpha_cu; ///mean alpha among all stocks within each CU
  real<lower=0> sd_alpha_cu; ///among CU variance in alpha
  vector<lower=0>[J] alpha_j;//stock-specific alpha   
  real<lower=0> sd_alpha; //within CU variance

 //capacity
 vector<upper = 0>[J] log_b; // log transformed per capita density dependence by stock

//covariate effects
real b_ECA; //global (across stock) mean effect of ECA
vector[C] b_ECA_cu; //CU-specific ECA effect
real<lower=0> sd_ECA; //variance in CU-level ECA effect
real b_ECA_Area; //global interaction term for ECA x Area

 //variance components
 real<lower=0> mu_sigma; ///mean sigma among all stocks
 vector<lower=0>[C] cu_sigma; ///CU-level sigma
 real<lower=0> sd_sigma; ///variance in sigma within CUs (pooled among CUs)
 real<lower=0> sd_sigma_cu; ///variance in sigma among CUs  
 vector<lower = 0>[J] sigma; ///stock-level sigma
 vector<lower = -1,upper=1>[J] rho; //autocorrelation parameter
}
transformed parameters{
  vector[J] b; //capacity rate (untransformed)
  vector[J] sigmaAR; //sigma - adjusted for rho
	
  //productivity residuals through time
  vector[N] e_t; //stock residual productivity at time t
  vector[N] mu1; //initial expectation at each time for each stock
  vector[N] mu2; //autocorr. adjusted expectation at each time for each stock
  
  b=exp(log_b);
  
  //residual productivity deviations
   for(j in 1:J){ //for every stock
	mu1[start_y[j]:end_y[j]]=alpha_j[j] - b[j]*S[start_y[j]:end_y[j],j]+b_ECA_cu[C_i[j]]*ECA[start_y[j]:end_y[j],j]+b_ECA_Area*ExA[start_y[j]:end_y[j],j]; //expectation (unadjusted for autocorrelated residuals)

	e_t[start_y[j]] = R_S[start_y[j]] - mu1[start_y[j]]; //first deviate for stock j
	
	for(t in 1:(N_i[j]-1)){ //adjust expectation based on autocorrelation
		 mu2[start_y[j]+t] = mu1[start_y[j]+t] + (rho[j]^(ii[start_y[j]+t]-ii[start_y[j]+t-1])*e_t[start_y[j]+t-1]); //adjust expectation based on previous deviate - rho is raised to the power of the number of time steps (in years) between observations
         e_t[start_y[j]+t] = R_S[start_y[j]+t] - mu2[start_y[j]+t];  //residual for stock j at time t
	}
  
  }
  sigmaAR = sigma.*sqrt(1-rho^2); //correct stock sigma for autocorrelation (rho)
}  
model{
  //priors
  //prod
  alpha_0 ~ normal(1.5,1); //global intrinsic productivity for all stocks
  alpha_cu ~ normal(alpha_0,sd_alpha); //CU-specific deviations in intrinsic productivity
  alpha_j ~ normal(alpha_cu[C_i],sd_alpha_cu); //within CU-deviations among stocks
  
//hierarchical variances
  sd_alpha ~ normal(0,0.5); //among CU variance in productivity
  sd_alpha_cu ~ normal(0,0.5); //within CU variance in productivity
    
  //capacity for each stock - fit individually with weakly informative priors based on maximum observed spawners
  for(j in 1:J) log_b[j] ~ normal(logbeta_pr[j],logbeta_pr_sig[j]);
 
  //covariate effects
  b_ECA ~ normal(0,1); //standard normal prior
  b_ECA_cu ~ normal(b_ECA,sd_ECA); //CU-specific ECA effect
  
  b_Area ~ normal(0,1); //standard normal prior
  b_Area_cu ~ normal(b_Area,sd_Area); //CU-specific ECA effect
  
  b_ECA_Area ~ normal(0,1); //standard normal prior
  b_ECA_Area_cu ~ normal(b_ECA_Area,sd_ExA); //CU-specific ECA effect
 
  //hierarchical variances
  sd_ECA ~ normal(0,0.5); //variance in stock-level ECA effects
  sd_Area ~ normal(0,0.5); //variance in stock-level ECA effects
  sd_ExA ~ normal(0,0.5); //variance in stock-level ECA effects
  
  //variance terms
  mu_sigma ~ normal(0.5,1);
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
vector[J] Smax; //Smax - spawners where peak recruitment is achieved
vector[J] Smsy; //Smsy - spawners for max. sustainable yield
vector[J] Umsy; //Umsy - harvest rate corresponding to Smsy

for(j in 1:J){
 log_lik[start_y[j]]=normal_lpdf(R_S[start_y[j]]|mu1[start_y[j]],sigma[j]); //pointwise log likelihood calculation - initial estimates
 
 for(t in 1:(N_i[j]-1)){
 log_lik[start_y[j]+t]=normal_lpdf(R_S[start_y[j]+t]|mu2[start_y[j]+t], sigmaAR[j]); //pointwise log likelihood calculation
 }
 
 Smax[j] = 1/b[j];
 Umsy[j] = 1-lambert_w0(exp(1-alpha_j[j]));
 Smsy[j] = (1-lambert_w0(exp(1-alpha_j[j])))/b[j];
}
}


