## See Sorensen & Vasishth (under review) A tutorial on fitting Bayesian linear mixed models using Stan
# The code relies heavily on the code of Tanner Sorensen
# see: http://www.ling.uni-potsdam.de/~vasishth/statistics/BayesLMMs.html

#by saving in x, x_u, and x_w the output of model.matrix(formula), the model can be used for different lmms
functions {
  real shift_max(vector shift_u, int[] subj, vector rt) {
  
    real shift_max;
    shift_max <- positive_infinity();
    for (i in 1:num_elements(rt)){
         shift_max <- fmin(shift_max,log(rt[i]) - shift_u[subj[i]] );
        }
    return shift_max;
  }


  real shifted_lognormal_log(vector y, vector mu, real sigma, vector psi){

    if (min(psi) < 0)
      reject("Shift parameter (psi) should be bigger than 0, value found ",
        min(psi));

    if (min(y-psi) < 0)
      reject(
"Shift parameter (psi) should be smaller than y;  y-psi > 0 but here there is",
min(y-psi));
 
  return lognormal_log(y- psi,mu,sigma);
  }

}

data {
  int<lower=0> N_obs; 
  int<lower=0> N_coef;  //predictors +intercept
  int<lower=0> N_coef_u;  //predictors +intercept
  int<lower=0> N_coef_w;  //predictors +intercept
  int<lower=0> N_coef_v;  //predictors +intercept
  real effsize;
  int<lower=1> subj[N_obs];    //subject id
  int<lower=1> item[N_obs];    //item id
  int<lower=1> sentence[N_obs];    //item id

  int<lower=1> N_subj;                 //number of subjects
  int<lower=1> N_item;  
  int<lower=1> N_sentence;  
  matrix[N_obs,N_coef] x;
  matrix[N_obs,N_coef_u] x_u;
  matrix[N_obs,N_coef_w] x_w;
  matrix[N_obs,N_coef_v] x_v;
  vector[N_obs] rt;

}

transformed data {
  matrix[N_obs,N_coef-1] x_betas;
  x_betas <- block(x,1,2,N_obs,N_coef-1); # I remove the intercept here
}
parameters {
  vector[N_coef-1] delta;     //delta is size of the effect
  real<lower=0> sigma;
  real<lower=0> tau_shift;
  vector[N_subj] shift_u_raw;  // subj shift
  real<upper=shift_max(shift_u_raw * tau_shift,subj,rt)> shift;

  //subj
  vector<lower=0> [N_coef_u]  tau_u;     // subj sd
  cholesky_factor_corr[N_coef_u] L_u;      // correlation matrix for random intercepts and slopes subj
  vector[N_coef_u] z_u[N_subj];

  //items
  vector<lower=0> [N_coef_w]  tau_w;     // subj sd
  cholesky_factor_corr[N_coef_w] L_w;      // correlation matrix for random intercepts and slopes item
  vector[N_coef_w] z_w[N_item];

  //sentence
  vector<lower=0> [N_coef_v]  tau_v;     // subj sd
  cholesky_factor_corr[N_coef_v] L_v;      // correlation matrix for random intercepts and slopes item
  vector[N_coef_v] z_v[N_sentence];

  real<lower=0> alpha; 

}


transformed parameters {
     
  vector[N_coef_u]  u[N_subj];         
  vector[N_coef_w]  w[N_item];
  vector[N_coef_v]  v[N_sentence];
  matrix[N_coef_u,N_coef_u] Lambda_u;
  matrix[N_coef_w,N_coef_w] Lambda_w;
  matrix[N_coef_v,N_coef_v] Lambda_v;
  vector[N_coef-1] beta;
  vector[N_obs] psi;  //each shift
  vector[N_obs] mu; 
  vector[N_subj] shift_u;  // subj shift

  beta <- delta * sigma;  

  Lambda_u <- diag_pre_multiply(tau_u,L_u);
  for (i in 1:N_subj){
   u[i] <- Lambda_u * z_u[i];
  }

  Lambda_w <-  diag_pre_multiply(tau_w,L_w);  
  for (i in 1:N_item){
    w[i] <- Lambda_w * z_w[i]; // item random effects
  } 

  Lambda_v <-  diag_pre_multiply(tau_v,L_v);  
  for (i in 1:N_sentence){
    v[i] <- Lambda_v * z_v[i]; // item random effects
  }


  shift_u <- shift_u_raw * tau_shift; // =shift_u ~normal(0,tau_shift)
  for (i in 1:N_obs){
    mu[i] <- alpha + x_betas[i] * beta +   
          x_u[i] * u[subj[i]] + 
          x_w[i] * w[item[i]] + 
          x_v[i] * v[sentence[i]];
    psi[i] <- exp(shift + shift_u[subj[i]]); //
  }

} 

model {

  sigma ~ normal(0,1);
  tau_u ~ normal(0,1);
  tau_w ~ normal(0,1);
  tau_v ~ normal(0,1);
  tau_shift  ~ normal(0,.5);

  alpha ~ normal(0, 5);
  delta ~ normal(0, effsize);

  #for (i in 1:N_subj)
  shift ~ normal(0,1); 
  shift_u_raw ~ normal(0,1);

  L_u ~ lkj_corr_cholesky(4.0);
  L_w ~ lkj_corr_cholesky(4.0);
  L_v ~ lkj_corr_cholesky(4.0);

  for (i in 1:N_subj){
    z_u[i] ~ normal(0,1);   
  }

  for (i in 1:N_item){
    z_w[i] ~ normal(0,1);   
  }

  for (i in 1:N_sentence){
    z_v[i] ~ normal(0,1);   
  }

    rt ~ shifted_lognormal(mu, sigma,psi);
}

generated quantities {
  matrix[N_coef_u,N_coef_u] Cor_u;
  matrix[N_coef_w,N_coef_w] Cor_w;
  matrix[N_coef_v,N_coef_v] Cor_v;
  real pred_rt[N_obs];
  real log_lik[N_obs];
  real resid[N_obs];



  Cor_u <- tcrossprod(L_u);  //Correlations between random effects by subj
  Cor_w <- tcrossprod(L_w);  //Correlations between random effects by item
  Cor_v <- tcrossprod(L_v);  //Correlations between random effects by item

  for (i in 1:N_obs){
  pred_rt[i] <- lognormal_rng(mu[i],sigma) + psi[i];
  log_lik[i] <- lognormal_log(rt[i]-psi[i],mu[i], sigma);
  resid[i] <- log(rt[i]-psi[i]) - mu[i];
  }

}
