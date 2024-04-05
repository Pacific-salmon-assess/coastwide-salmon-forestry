functions {
  vector gp_pred_rng(real[] x2,
                     vector y1, real[] x1,
                     real alpha, real rho, real sigma, real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] K =   cov_exp_quad(x1, alpha, rho)
                         + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);

      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho);
      vector[N2] f2_mu = (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      matrix[N2, N2] cov_f2 =   cov_exp_quad(x2, alpha, rho) - v_pred' * v_pred
                              + diag_matrix(rep_vector(delta, N2));
      f2 = multi_normal_rng(f2_mu, cov_f2);
    }
    return f2;
  }
}
data{
  int<lower=1> N;//number of observations
  int L; //length covered by all time-series
  int C; //number of CUs
  int J;//number of river-level stocks
  int N_i[J]; //number of observations in each time-series
  int C_i[J]; //CU index for each stock
  int ii[N]; //index of brood years
  vector[N] R_S; //matrix of productivity among stocks - log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners in time T
  vector[N] ECA_vec; //vector of stock-specific ECA through time
  int<lower=0> start_y[J];       // ragged start point for observations (N)
  int<lower=0> end_y[J];         // ragged end points for observations (N)
  vector[J] pRk_mean; //priors on equilibrium recruitment - based on obs R
  vector[J] pRk_sig; //priors on sigma for equilibrium recruitment - based on obs R
  int<lower=1> N_predict; //length of predicted function for GP
  vector[N_predict] x_pred; //vector to predict over for GP
}
transformed data{
vector[J] logRk_pr;
vector[J] logRk_pr_sig;

array[N] real ECA_vec_arr = to_array_1d(ECA_vec);
array[N_predict] x_predict = to_array_1d(x_pred);

for(t in 1:J){
logRk_pr_sig[t]=sqrt(log(1+((pRk_sig[t])*(pRk_sig[t]))/((pRk_mean[t])*(pRk_mean[t])))); //this converts sigma on the untransformed scale to a log scale
logRk_pr[t]=log(pRk_mean[t])-0.5*logRk_pr_sig[t]*logRk_pr_sig[t]; //convert smax prior to per capita slope - transform to log scale with bias correction
}

}
parameters{
  //intrinsic productivity
  real<lower=0> alpha_0; ///mean alpha among all stocks
  vector<lower=0>[C] alpha_cu; ///mean alpha among all stocks within each CU
  real<lower=0> sd_alpha_cu; ///among CU variance in alpha
  vector<lower=0>[J] alpha_j;//stock-specific alpha   
  real<lower=0> sd_alpha; //within CU variance

 //BH eq. recruitment
  vector<lower=0>[J] Rk; //equilibrium recruitment for BH

//covariate effects
vector[N] z_eca; //effect sizes
real<lower=0> gp_tau; //gaussian process scale parameter - controls variance in ECA response (vertical distance of non linear function)
real<lower=0> gp_rho;  //gaussian process (aka lengthscale) rho parameter - controls smoothness in ECA GP response
  
 //variance components
 real<lower=0> mu_sigma; ///mean sigma among all stocks
 vector<lower=0>[C] cu_sigma; ///CU-level sigma
 real<lower=0> sd_sigma; ///variance in sigma within CUs (pooled among CUs)
 real<lower=0> sd_sigma_cu; ///variance in sigma among CUs  
 vector<lower = 0>[J] sigma; ///stock-level sigma
 vector<lower = -1,upper=1>[J] rho; //autocorrelation parameter
}
transformed parameters{
  vector[J] sigmaAR; //sigma - adjusted for rho
  vector[N] B_eca; //realized ECA effect
	
  //productivity residuals through time
  vector[N] e_t; //stock residual productivity at time t
  vector[N] mu1; //initial expectation at each time for each stock
  vector[N] mu2; //autocorr. adjusted expectation at each time for each stock

	//Gaussian process for the effect of ECA
    matrix[N,N] L_Tmat;
	matrix[N,N] Tmat;
	Tmat = gp_exp_quad_cov(ECA_vec_arr, gp_tau, gp_rho);
	for(n in 1:N) Tmat[n, n] += 1e-12; //ensures positive for Cholesky decomp.
	L_Tmat = cholesky_decompose(Tmat);
	B_eca =L_Tmat*z_eca;
	
  //residual productivity deviations
   for(j in 1:J){ //for every stock
   
   //initial expectation (uncorrected for autocorrelation)
	mu1[start_y[j]:end_y[j]]=alpha_j[j]-log(1+(exp(alpha_j[j])/Rk[j])*S[start_y[j]:end_y[j],j])+B_eca[start_y[j]:end_y[j]]; //expectation (unadjusted for autocorrelated residuals)

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
  for(j in 1:J) Rk[j] ~ lognormal(logRk_pr[j],logRk_pr_sig[j]);
 
  //Gaussian process for ECA - nonlinear function generated from a MVN
  z_eca ~ normal(0,1);
  gp_tau ~ normal(0,2);
  gp_rho ~ inv_gamma(5, 5);
  
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
vector[J] Smsy; //Smsy - spawners for max. sustainable yield
vector[J] Umsy; //Umsy - harvest rate corresponding to Smsy
vector[N_predict] f_predict = gp_pred_rng(x_predict, y, x, alpha, rho, sigma, 1e-10);
vector[N_predict] y_predict;

for (n in 1:N_predict)
    y_predict[n] = normal_rng(f_predict[n], sigma);
}

for(j in 1:J){
 log_lik[start_y[j]]=normal_lpdf(R_S[start_y[j]]|mu1[start_y[j]],sigma[j]); //pointwise log likelihood calculation - initial estimates
 
 for(t in 1:(N_i[j]-1)){
 log_lik[start_y[j]+t]=normal_lpdf(R_S[start_y[j]+t]|mu2[start_y[j]+t], sigmaAR[j]); //pointwise log likelihood calculation
 }
 
 Smsy[j] = Rk[j]*sqrt(1/exp(alpha_j[j]))-(Rk[j]/exp(alpha_j[j]));
 Umsy[j] = 1-lambert_w0(exp(1-alpha_j[j]));
}
}


