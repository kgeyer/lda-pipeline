data {
  int <lower=2> K;                      // num topics
  int <lower=1> D;                      // num docs
  int <lower=2> J;                      // num of unique words (vocab size)
  int <lower=1> M;                      // total # of occurring words
  int <lower=1, upper=J> featw[M];      // word m ID
  int <lower=1, upper=D> docw[M];       // doc ID for each word instance
  vector <lower=0>[K] alpha;            // topic prior
  vector <lower=0>[J] mu;               // word prior
}

parameters {
  simplex[K] theta[D];   // topic dist for doc d
  simplex[J] phi[K];     // word dist for topic k
}

model {
  for (d in 1:D){
    theta[d] ~ dirichlet(alpha);    // prior
  }
  for (k in 1:K){
    phi[k] ~ dirichlet(mu);         // prior
  }
  for (m in 1:M){
    real wrd_gamma[K];
    for (k in 1:K){
      wrd_gamma[k] = log(theta[docw[m], k]) + log(phi[k, featw[m]]);
    }
    target += log_sum_exp(wrd_gamma);   // log-likelihood
  }
}
