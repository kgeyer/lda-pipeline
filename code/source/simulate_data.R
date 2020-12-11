#-------------------------------------------------------------------------------
# HEADER
#-------------------------------------------------------------------------------
# Author: Kelly Geyer, klgeyer@bu.edu
# Date: February 7, 2020
# Development: R v4.0.0 on MacOS v10.15.4
#
# Description: This script contains functions for generating test data for the 
# MM-LDA model.

#-------------------------------------------------------------------------------
# LIBRARIES
#-------------------------------------------------------------------------------
library(Matrix)
library(MCMCpack)
library(extraDistr)

#-------------------------------------------------------------------------------
# GENERATE DATA FOR MM-LDA TESTING
#-------------------------------------------------------------------------------
generate.mmlda.data <- function(K, D, J, I, lam=300, alp=NA, sig=NA, mu=NA){
  # This function generates test data for MM-LDA 
  #
  # Args:
  #   K (>= 2): Number of factors
  #   D: Number of documents
  #   J: vocabulary count
  #   I: number of distinct image regions
  #   lam: Avg. number of words/regions per document
  #   alp, sig, mu: Dirichlet hyperparameter values
  # Check parameters
  if (K < 2) {stop("'K' must be at least 2")}
  if (D < 1) {stop("'D' must be at least 1")}
  if (J < 2) {stop("'J' must be at least 2")}
  if (I < 2) {stop("'I' must be at least 2")}
  if (lam < 10) {stop("'lam' must be at least 10")}
  if (any(is.na(alp))) {alp = matrix(1, nrow = K, ncol = 1)}
  if (any(is.na(sig))) {sig = matrix(1, nrow = I, ncol = 1)}
  if (any(is.na(mu))) {mu = matrix(1, nrow = J, ncol = 1)}
  if (any(alp <= 0)) {stop("'alpha' must be greater than 0.")}
  if (any(sig <= 0)) {stop("'sigma' must be greater than 0.")}
  if (any(mu <= 0)) {stop("'mu' must be greater than 0.")}
  # Create true parameters
  theta <- drop(replicate(D, MCMCpack::rdirichlet(1, alpha=alp)))
  beta <- drop(replicate(K, MCMCpack::rdirichlet(1, alpha=sig)))
  phi <- drop(replicate(K, MCMCpack::rdirichlet(1, alpha=mu)))
  # Create document data
  R <- Matrix::Matrix(0, nrow = D, ncol = I, sparse = TRUE, byrow = TRUE)  # img
  W <- Matrix::Matrix(0, nrow = D, ncol = J, sparse = TRUE, byrow = TRUE) # word
  for (d in 1:D){
    # Number of words & image regions within document
    n.word <- rpois(1, lambda = lam)
    n.reg <- rpois(1, lambda = lam)
    # Sample labels
    vocab.topic.count <- count.vector(sample(1:K, size=n.word, replace=TRUE, prob=theta[,d]), K)
    region.topic.count <- count.vector(sample(1:K, size=n.reg, replace=TRUE, prob=theta[,d]), K)
    # Generates words/images for each topic...
    Rk <- sapply(1:K, function(k){count.vector(sample(1:I, size=region.topic.count[k], replace=TRUE, prob=beta[,k]), I)})
    Wk <- sapply(1:K, function(k){count.vector(sample(1:J, size=vocab.topic.count[k], replace=TRUE, prob=phi[,k]), J)})
    W[d,] <- rowSums(Wk)
    R[d,] <- rowSums(Rk)
  }
  return(list(theta=t(theta), beta=t(beta), phi=t(phi), W=W, R=R, lam=lam, 
              alp=alp, sig=sig, mu=mu))
}

generate.mixed.signal.data <- function(K, D, J, I, lam=300, a=0.5, b=2){
  # This function generates test data for MM-LDA where signatures have varying 
  # sparsity levels
  #
  # Args:
  #   K (>= 2): Number of factors
  #   D: Number of documents
  #   J: vocabulary count
  #   I: number of disticnt image regions
  #   lam: Avg. number of words/regions per document
  #   a: gamma shape parameter
  #   b: gamma rate parameter
  # Check parameters
  if (K < 2) {stop("'K' must be at least 2")}
  if (D < 1) {stop("'D' must be at least 1")}
  if (J < 2) {stop("'J' must be at least 2")}
  if (I < 2) {stop("'I' must be at least 2")}
  if (lam < 10) {stop("'lam' must be at least 10")}
  if (a <= 0) {stop("'a' must be greater than 0.")}
  if (b <= 0) {stop("'b' must be greater than 0.")}
  # Create true parameters
  phi = matrix(NA, nrow = K, ncol = J)
  beta = matrix(NA, nrow = K, ncol = I)
  sig = matrix(rgamma(K*I, shape = a, rate = b), nrow = K, ncol = I)
  mu = matrix(rgamma(K*J, shape = a, rate = b), nrow = K, ncol = J)
  # Generate theta
  alp = matrix(rgamma(K*D, shape = a, rate = b), nrow = K, ncol = D)
  theta = matrix(NA, nrow = D, ncol = K)
  
  # done
  alp = matrix(nrow = K, ncol = 3); colnames(alp) = c("theta", "beta", "phi")
  for (k in 1:K){
    alp[k,'theta'] = rgamma(1, shape = a, rate = b)
    theta[,k] = MCMCpack::rdirichlet(1, rep(alp[k,'theta'],D))
    while (any(is.na(theta[,k]))){
      alp[k,'theta'] = 10*alp[k,'theta']
      theta[,k] = MCMCpack::rdirichlet(1, rep(alp[k,'theta'],D))
    }
  }
  for (k in 1:K){
    alp[k,'beta'] = rgamma(1, shape = a, rate = b)
    beta[k,] = MCMCpack::rdirichlet(1, rep(alp[k,'beta'],I))
    while (any(is.na(beta[k,]))){
      alp[k,'beta'] = 10*alp[k,'beta']
      beta[k,] = MCMCpack::rdirichlet(1, rep(alp[k,'beta'],I))
    }
  }
  for (k in 1:K){
    alp[k,'phi'] = rgamma(1, shape = a, rate = b)
    phi[k,] = MCMCpack::rdirichlet(1, rep(alp[k,'phi'],J))
    while (any(is.na(phi[k,]))){
      alp[k,'phi'] = 10*alp[k,'phi']
      phi[k,] = MCMCpack::rdirichlet(1, rep(alp[k,'phi'],J))
    }
  }
  # Create document data
  R <- Matrix::Matrix(0, nrow = D, ncol = I, sparse = TRUE, byrow = TRUE)    # image
  W <- Matrix::Matrix(0, nrow = D, ncol = J, sparse = TRUE, byrow = TRUE)    # word
  for (d in 1:D){
    # Number of words & image regions within document
    n.word <- rpois(1, lambda = lam)
    n.reg <- rpois(1, lambda = lam)
    # Sample labels
    vocab.topic.count <- count.vector(sample(1:K, size=n.word, replace=TRUE, prob=theta[,d]), K)
    region.topic.count <- count.vector(sample(1:K, size=n.reg, replace=TRUE, prob=theta[,d]), K)
    # Generates words/images for each topic...
    Rk <- sapply(1:K, function(k){count.vector(sample(1:I, size=region.topic.count[k], replace=TRUE, prob=beta[,k]), I)})
    Wk <- sapply(1:K, function(k){count.vector(sample(1:J, size=vocab.topic.count[k], replace=TRUE, prob=phi[,k]), J)})
    W[d,] <- rowSums(Wk)
    R[d,] <- rowSums(Rk)
  }
  return(list(theta=t(theta), beta=t(beta), phi=t(phi), W=W, R=R, lam=lam, 
              a=a, b=b, all_alpha=alp,
              alp=mean(alp[,'theta']), mu=mean(alp[,'phi']), 
              sig=mean(alp[,'beta'])))
}

#-------------------------------------------------------------------------------
# GENERATE DATA FOR sparseMMLDA TESTING
#-------------------------------------------------------------------------------

# K = 3
# S = 50
# N1 = 30
# N2 = 30
# alp1=1
# alp2=1
# mu=1
# gam1=1
# gam2=2
# lam1=300
# lam2=300


generate.sparseMMLDA.data.v1 <- function(K, S, N1, N2, alp1=1, alp2=1, mu=1,
                                      gam1=0.5, gam2=0.5, lam1=300, lam2=300){
  # This function generates test data for testing sparse MM-LDA 
  #
  # We will assume that mu, sig, and alp all have the sample element.
  #
  # Args:
  #   K (>= 2): Number of factors
  #   S (>= 2): Number of samples
  #   D (>= 2): Number of data modalities
  #   Nd list: List of the total number of mutations for each modality
  #   lam (>= 2): Average number of mutations in data modality
  #   mu: S-dim vector of positive, non-zero elements; Hyperparameter of 
  #     sample-signature proportions
  #   alp: List of Nd-dim vectors with positive, non-zero elements; 
  #     hyperparameter of phi
  #   gam1, gam2 (>= 0): Hyperparameters of Pi
  #   sig: List of K-dim vectors with positive, non-zero elements; 
  #     hyperparameter of gamma
  # Check parameters
  if (K < 2) {stop("'K' must be at least 2.")}
  if (S < 2) {stop("'S' must be at least 2.")}
  if (N1 < 2) {stop("'N1' must be at least 2.")}
  if (N2 < 2) {stop("'N2' must be at least 2.")}
  if (alp1 <= 0) {stop("'alp1' must be greater than 0.")}
  if (alp2 <= 0) {stop("'alp2' must be greater than 0.")}
  if (mu <= 0) {stop("'mu' must be greater than 0.")}
  if (gam1 <= 0) {stop("'gam1' must be greater than 0.")}
  if (gam2 <= 0) {stop("'gam2' must be greater than 0.")}
  # Generate parameters
  theta = MCMCpack::rdirichlet(S, alpha = matrix(mu, nrow=K, ncol=1))
  phi1 = MCMCpack::rdirichlet(K, alpha = matrix(alp1, nrow=N1, ncol=1))
  phi2 = MCMCpack::rdirichlet(K, alpha = matrix(alp2, nrow=N2, ncol=1))
  b = matrix(NA, nrow=K, ncol=2)
  b[,1] = rbeta(K, shape1=gam1, shape2=gam2)
  b[,2] = 1 - b[,1]
  # Generate data
  X1 <- Matrix::Matrix(0, nrow=S, ncol=N1, sparse=TRUE, byrow=TRUE)
  X2 <- Matrix::Matrix(0, nrow=S, ncol=N2, sparse=TRUE, byrow=TRUE)
  for (s in 1:S) {
    # Generate the number of mutations within a sample
    nmut1 <- rpois(1, lambda=lam1)
    nmut2 <- rpois(1, lambda=lam2)
    # sample signature labels
    v1 = theta[s,] * b[,1]
    p1 = v1 / sqrt(sum(v1^1))
    mut.sig.count1 <- count.vector(sample(1:K, size=nmut1, replace=TRUE, prob=p1), K)
    v2 = theta[s,] * b[,2]
    p2 = v2 / sqrt(sum(v2^2))
    mut.sig.count2 <- count.vector(sample(1:K, size=nmut2, replace=TRUE, prob=p2), K)
    # Generate mutations for each signature
    Xk1 <- sapply(1:K, function(k){count.vector(sample(1:N1, size=mut.sig.count1[k], 
                                                       replace=TRUE, prob=phi1[k,]), N1)})
    X1[s,] <- rowSums(Xk1)
    Xk2 <- sapply(1:K, function(k){count.vector(sample(1:N2, size=mut.sig.count2[k], 
                                                       replace=TRUE, prob=phi2[k,]), N2)})
    X2[s,] <- rowSums(Xk2)
  }
  return(list(theta=theta, phi1=phi1, phi2=phi2, X1=X1, X2=X2, b=b))
}

generate.sparseMMLDA.data.v2 <- function(K, S, N1, N2, alp1=1, alp2=1, mu=1,
                                         gam1=3, gam2=1, lam1=300, lam2=300){
  # This function generates test data for testing sparse MM-LDA 
  #
  # We will assume that mu, sig, and alp all have the sample element.
  #
  # Args:
  #   K (>= 2): Number of factors
  #   S (>= 2): Number of samples
  #   D (>= 2): Number of data modalities
  #   Nd list: List of the total number of mutations for each modality
  #   lam (>= 2): Average number of mutations in data modality
  #   mu: S-dim vector of positive, non-zero elements; Hyperparameter of 
  #     sample-signature proportions
  #   alp: List of Nd-dim vectors with positive, non-zero elements; 
  #     hyperparameter of phi
  #   gam1, gam2 (>= 0): Hyperparameters of Pi
  #   sig: List of K-dim vectors with positive, non-zero elements; 
  #     hyperparameter of gamma
  # Check parameters
  if (K < 2) {stop("'K' must be at least 2.")}
  if (S < 2) {stop("'S' must be at least 2.")}
  if (N1 < 2) {stop("'N1' must be at least 2.")}
  if (N2 < 2) {stop("'N2' must be at least 2.")}
  if (alp1 <= 0) {stop("'alp1' must be greater than 0.")}
  if (alp2 <= 0) {stop("'alp2' must be greater than 0.")}
  if (mu <= 0) {stop("'mu' must be greater than 0.")}
  if (gam1 <= 0) {stop("'gam1' must be greater than 0.")}
  if (gam2 <= 0) {stop("'gam2' must be greater than 0.")}
  # Generate parameters
  theta = MCMCpack::rdirichlet(S, alpha = matrix(mu, nrow=K, ncol=1))
  phi1 = MCMCpack::rdirichlet(K, alpha = matrix(alp1, nrow=N1, ncol=1))
  phi2 = MCMCpack::rdirichlet(K, alpha = matrix(alp2, nrow=N2, ncol=1))
  c.sq = 1/rgamma(1, shape=gam1, rate=gam2)
  tau = extraDistr::rhcauchy(1)
  beta = extraDistr::rhcauchy(K)
  beta.tilde = calc_beta_tilde(beta, tau, c.sq)
  pp = rnorm(K, mean=0, sd = (tau*beta.tilde))   # pi
  b = matrix(NA, nrow=K, ncol=2)
  for (k in 1:K){
    b[k,1] = exp(log(1 - pnorm(pp[k], mean=0, sd=tau*beta.tilde[k])) + 
                   sum(log(pnorm(pp[-c(k)], mean=0, sd=tau*beta.tilde[-c(k)]))))
  }
  b[,1] = b[,1]/sum(b[,1])
  b[,2] = 1 - b[,1]
  # Generate data
  X1 <- Matrix::Matrix(0, nrow=S, ncol=N1, sparse=TRUE, byrow=TRUE)
  X2 <- Matrix::Matrix(0, nrow=S, ncol=N2, sparse=TRUE, byrow=TRUE)
  for (s in 1:S) {
    # Generate the number of mutations within a sample
    nmut1 <- rpois(1, lambda=lam1)
    nmut2 <- rpois(1, lambda=lam2)
    # sample signature labels
    v1 = theta[s,] * b[,1]
    p1 = v1 / sqrt(sum(v1^1))
    mut.sig.count1 <- count.vector(sample(1:K, size=nmut1, replace=TRUE, prob=p1), K)
    v2 = theta[s,] * b[,2]
    p2 = v2 / sqrt(sum(v2^2))
    mut.sig.count2 <- count.vector(sample(1:K, size=nmut2, replace=TRUE, prob=p2), K)
    # Generate mutations for each signature
    Xk1 <- sapply(1:K, function(k){count.vector(sample(1:N1, size=mut.sig.count1[k], 
                                                       replace=TRUE, prob=phi1[k,]), N1)})
    X1[s,] <- rowSums(Xk1)
    Xk2 <- sapply(1:K, function(k){count.vector(sample(1:N2, size=mut.sig.count2[k], 
                                                       replace=TRUE, prob=phi2[k,]), N2)})
    X2[s,] <- rowSums(Xk2)
  }
  return(list(theta=theta, phi1=phi1, phi2=phi2, X1=X1, X2=X2, b=b, c.sq=c.sq, 
              tau=tau, beta=beta, beta.tilde=beta.tilde, pi=pp, b=b))
}


# generate.sparseMMLDA.data <- function(K, S, N1, N2, alp1=1, alp2=1, mu=1,
#                                       gam1=0.5, gam2=0.5){
#   # This function generates test data for testing sparse MM-LDA 
#   #
#   # We will assume that mu, sig, and alp all have the sample element.
#   #
#   # Args:
#   #   K (>= 2): Number of factors
#   #   S (>= 2): Number of samples
#   #   D (>= 2): Number of data modalities
#   #   Nd list: List of the total number of mutations for each modality
#   #   lam (>= 2): Average number of mutations in data modality
#   #   mu: S-dim vector of positive, non-zero elements; Hyperparameter of 
#   #     sample-signature proportions
#   #   alp: List of Nd-dim vectors with positive, non-zero elements; 
#   #     hyperparameter of phi
#   #   gam1, gam2 (>= 0): Hyperparameters of Pi
#   #   sig: List of K-dim vectors with positive, non-zero elements; 
#   #     hyperparameter of gamma
#   # Check parameters
#   if (K < 2) {stop("'K' must be at least 2.")}
#   if (S < 2) {stop("'S' must be at least 2.")}
#   if (N1 < 2) {stop("'N1' must be at least 2.")}
#   if (N2 < 2) {stop("'N2' must be at least 2.")}
#   if (alp1 < 2) {stop("'alp1' must be greater than 0.")}
#   if (alp2 < 2) {stop("'alp2' must be greater than 0.")}
#   if (mu < 2) {stop("'mu' must be greater than 0.")}
#   if (gam1 <= 0) {stop("'gam1' must be greater than 0.")}
#   if (gam2 <= 0) {stop("'gam2' must be greater than 0.")}
#   # Generate parameters, using hyperparameters
#   theta = MCMCpack::rdirichlet(S, alpha = matrix(mu, nrow=K, ncol=1))
#   phi1 = lapply(1:K, function(k){
#     MCMCpack::rdirichlet(K, alpha=matrix(alp1, nrow=N1, ncol=1))})
#   phi2 = lapply(1:K, function(k){
#     MCMCpack::rdirichlet(K, alpha=matrix(alp2, nrow=N2, ncol=1))})
#   pi1 = matrix(rbeta(K, shape1=gam1, shape2=gam2), nrow=K, ncol=1)
#   pi2 = 1 - pi1
#   # 
#   
#   Theta = MCMCpack::rdirichlet(S, alpha=matrix(mu, nrow=K, ncol=1))
#   Phi = lapply(1:length(Nd), function(d){MCMCpack::rdirichlet(K, alpha=matrix(alp, nrow=Nd[d], ncol=1))})
#   Pi = matrix(rbeta(K*D, shape1=5, shape2=1), nrow=K, ncol=D)
#   B = matrix(0, nrow=K, ncol=D)
#   for (d in 1:D){ for (k in 1:K){ B[k,d] = rbinom(1, 1, Pi[k,d]) } }
#   if (any(rowSums(B) < 1)){print("FAILURE IN GENERATING B!!!!!!!")}
#   # Generate data
#   X = vector('list', length=length(Nd))
#   for (d in 1:length(Nd)){
#     X[[d]] <- Matrix::Matrix(0, nrow=S, ncol=Nd[d], sparse=TRUE, byrow=TRUE)
#     for (s in 1:S){
#       # generate number of mutations within data modality d
#       n.mut <- rpois(1, lambda=lam)
#       # sample signature labels
#       mut.sig.count <- count.vector(sample(1:K, size=n.mut, replace=TRUE, prob=Theta[s,]), K)
#       # generate mutations for each signature
#       Xk <- sapply(1:K, function(k){count.vector(sample(1:Nd[d], size=mut.sig.count[k], replace=TRUE, prob=Phi[[d]][k,]), Nd[d])})
#       X[[d]][s,] <- rowSums(Xk)
#     }
#   }
#   return(list(Theta=Theta, Phi=Phi, Pi=Pi, B=B, X=X))
# }

#-------------------------------------------------------------------------------
# SUPPORTING FUNCTIONS
#-------------------------------------------------------------------------------
calc_beta_tilde <- function(beta, tau, c.sq){
  return(sqrt((c.sq * (beta^2)) / (c.sq + ((tau^2)*(beta^2)))))
}

count.vector <- function(q, q.dim){
  # This function returns a categorical count vector 
  #
  # Args:
  #   q - list of factor variables
  #   q.dim - number of factor levels
  q2 = as.data.frame(table(q))
  counts <- matrix(0, nrow = q.dim)
  counts[as.numeric(as.character(q2$q))] = q2$Freq
  return(counts)
}

generate.param.majority <- function(K, val = 0.8){
  # This function generates probability parameters
  # Assume that about 1/K of docs belong to each topic
  # Generate a majority topic for the document: Assume that (val)% majority and, 
  # 100*(1-val)/(K-1) for the rest
  # Repeat this for both vocabulary & image region features
  # 
  # Input:
  # K - number of topics
  # majority.prob - topic proportion
  #
  # Returns:
  # probs - (K x DIM) vector of probabilities
  # sample dominant topic
  if ((val < 0) || (val > 1)) {stop("'majority.prob' must be in [0,1]")}
  dominant.topic <- sample(1:K, size = 1)
  probs <- rep((1 - val)/(K - 1), K)
  probs[dominant.topic] = val
  return(probs)
}



