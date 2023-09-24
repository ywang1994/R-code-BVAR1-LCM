 data {
  int<lower=0> N; // number of stocks
  int<lower=0> T; // time length
  int Y[2*N*T];
  cov_matrix[2] wishart_prior_mat; // for the wishart prior
  int wishart_nu; //for the df of wishart prior
  real<lower=0.0> ar_rho_prior_sigma;// sigma parameters ar_rho prior
  real<lower=0.0> beta0_prior_sigma; // sigma paramerers for beta0 prior
}
transformed data{
  vector[2] zeros;
  zeros = rep_vector(0,2);
}
parameters{
  matrix[2,N] beta0;
  matrix[2,T] Bi_AR_innov;
  matrix[2,N*T] LCM_eff;
  cov_matrix[2] Bi_AR_Prec;
  cov_matrix[2] LCM_Prec;
  vector<lower=-1.0,upper=1.0>[2] ar_rho;
}
transformed parameters{
  matrix[2,T] Bi_AR_eff;
  matrix[2,2] ar_rho_diagmat;
  cov_matrix[2] Bi_AR_Sigma; // Covariance matrix of the two AR1
  cov_matrix[2] LCM_Sigma; // Covariance matrix of the LCM
  real<lower=0.0> Bi_AR_prec1; // precision for the first ar1
  real<lower=0.0> Bi_AR_prec2; // precision for the second ar1
  real<lower=-1.0,upper=1.0> Bi_AR_rho; // correlation between two ar1'
  real<lower=0.0> LCM_prec1; // precision for the first LCM
  real<lower=0.0> LCM_prec2; // precision for the second LCM
  real<lower=-1.0,upper=1.0> LCM_rho; // correlation between two LCM's
  
  Bi_AR_Sigma = inverse_spd(Bi_AR_Prec);
  Bi_AR_prec1 = 1/Bi_AR_Sigma[1,1];
  Bi_AR_prec2 = 1/Bi_AR_Sigma[2,2];
  Bi_AR_rho = Bi_AR_Sigma[1,2]*Bi_AR_prec1^.5*Bi_AR_prec2^.5;
  
  LCM_Sigma = inverse_spd(LCM_Prec);
  LCM_prec1 = 1/LCM_Sigma[1,1];
  LCM_prec2 = 1/LCM_Sigma[2,2];
  LCM_rho = LCM_Sigma[1,2]*LCM_prec1^.5*LCM_prec2^.5;
  
  ar_rho_diagmat = diag_matrix(ar_rho);
  
  Bi_AR_eff[1:2,1] = Bi_AR_innov[1:2,1];
  for(t in 2:T){
    Bi_AR_eff[1:2,t] = ar_rho_diagmat*Bi_AR_eff[1:2,t-1] + Bi_AR_innov[1:2,t];
  }
  
}
model {
  for(i in 1:N){
      beta0[1:2,i] ~ normal(0,beta0_prior_sigma);
  }
  ar_rho ~ normal(0,ar_rho_prior_sigma);
  Bi_AR_Prec ~ wishart(wishart_nu,wishart_prior_mat);
  LCM_Prec ~ wishart(wishart_nu,wishart_prior_mat);

   for(t in 1:T){
    Bi_AR_innov[1:2,t] ~ multi_normal_prec(zeros,Bi_AR_Prec);
      for(n in 1:N){
    LCM_eff[1:2,t+T*(n-1)] ~ multi_normal_prec(zeros, LCM_Prec);
    Y[t+T*(n-1)] ~ poisson_log(beta0[1,n]+Bi_AR_eff[1,t]+LCM_eff[1,t+T*(n-1)]);
    Y[t+T*(n-1)+N*T] ~ poisson_log(beta0[2,n]+Bi_AR_eff[2,t]+LCM_eff[2,t+T*(n-1)]);
  }
  }
}
generated quantities{
  int Y_fitted[2*N*T];
     for(t in 1:T){
      for(n in 1:N){
    Y_fitted[t+T*(n-1)] = poisson_log_rng(beta0[1,n]+Bi_AR_eff[1,t]+LCM_eff[1,t+T*(n-1)]);
    Y_fitted[t+T*(n-1)+N*T] = poisson_log_rng(beta0[2,n]+Bi_AR_eff[2,t]+LCM_eff[2,t+T*(n-1)]);
  }
  }
}