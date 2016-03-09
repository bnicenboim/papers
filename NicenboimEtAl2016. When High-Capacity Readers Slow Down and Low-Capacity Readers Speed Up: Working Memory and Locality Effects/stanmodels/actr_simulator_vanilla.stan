data {
  int<lower=0> N_obs; 
  real<lower=0> sigma;
  real expF;  
  real<lower=0> W_slope;
  real weight_slope;
  vector[N_obs] wmc;
  real<lower=0> W_int;
  real<lower=log(5)> MAS;
  real beta;
  real tau_int;  #
  vector[N_obs] decay_time;
  real<upper=1> d_int;
  real d_slope;

}

transformed data {
 real shared_cues;
 shared_cues <- 5.0;
}

model {
}

generated quantities {
  vector[N_obs]  Base_level_activation;
  vector<lower=0>[N_obs]  Spreading_activation;
  vector[N_obs]  A;
  vector<lower=0> [N_obs]  W;
  vector[N_obs]  tau;
  real s;
  vector<lower=0>[N_obs]  Pr;
  vector<lower=0,upper=1>[N_obs] weight;
  vector<upper=0>[N_obs] d;

  real pred_latency[N_obs];
  real pred_accuracy[N_obs];


  for (i in 1:N_obs){
    d[i] <- -exp(log(d_int)+d_slope* wmc[i]);
    Base_level_activation[i] <- log(decay_time[i]^d[i] ) +beta;
    W[i] <- exp(log(W_int) + W_slope * wmc[i]);
    tau[i] <- tau_int;
    weight[i] <- inv_logit(logit(1.0/3) + weight_slope * wmc[i]);
    Spreading_activation[i] <- W[i] * (weight[i] * (MAS-log(1)) + //unique cue
                        (1-weight[i]) * (MAS-log(shared_cues))) ;//non-unique cues
  }

  A <- Base_level_activation+Spreading_activation;
  s <- sqrt(3.0)*sigma/pi();
  Pr <- 1.0 ./(1.0+exp(-(A-tau)/s) );

  for (i in 1:N_obs){
    pred_accuracy[i] <- bernoulli_rng(Pr[i]);
    if(pred_accuracy[i]==1) {
      pred_latency[i] <- exp(logistic_rng( -A[i] + expF , sigma));
    } else {
      pred_latency[i] <- exp(-tau[i]+ expF)  ;
    }
  } //for
} // generated quantities
