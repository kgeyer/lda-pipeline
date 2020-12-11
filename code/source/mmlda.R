#-------------------------------------------------------------------------------
# HEADER
#-------------------------------------------------------------------------------
# Author: Kelly Geyer, klgeyer@bu.edu
# Date: February 7, 2020
# Development: R v4.0.0 on MacOS v10.15.4
#
# Description: This script contains functions for implementing the MM-LDA topic 
# model, using a collapsed Gibbs sampler.

#-------------------------------------------------------------------------------
# LIBRARIES
#-------------------------------------------------------------------------------
library(Matrix)
library(MCMCpack)

#-------------------------------------------------------------------------------
# FUNCTIONS
#-------------------------------------------------------------------------------
gibbs.sampler <- function(K, W, R, sig = 1, mu=1, alp = 1, iters = 1000, verbose=TRUE){
  # This function performs collapsed Gibbs sampling for the MM-LDA model
  #
  # Args:
  #   K (>= 2): Number of factors
  #   R (D x I): Document-region count
  #   W (D x J): Document-vocabulary count
  #   sigma (default 1's):
  #   mu (default 1's):
  #   alp (default 1's): alpha hyper parameter
  # Set parameters
  D = dim(R)[1]     # number of docs
  I = dim(R)[2]     # unique img. regions
  J = dim(W)[2]     # vocabulary
  # Check input
  if ((!is.null(alp)) & !is.numeric(alp) & (alp <= 0)){stop("'alp' must be a positive value")}
  if ((!is.null(sig)) & !is.numeric(sig) & (sig <= 0)){stop("'sig' must be a positive value")}
  if ((!is.null(mu)) & !is.numeric(mu) & (mu <= 0)){stop("'mu' must be a positive value")}
  alp = rep(alp, K)
  sig = rep(sig, I)
  mu = rep(mu, J)
  # Store each topic label of each word/region in each document (2d lists)
  Zdn = as.list(rep(NA, D*I))
  dim(Zdn) <- c(D, I)
  Vdm = as.list(rep(NA, D*J))
  dim(Vdm) = c(D, J)
  # Initialize count matrices
  doc.topic.count = Matrix::Matrix(0, nrow=D, ncol=K, sparse=TRUE)
  wrd.topic.count = Matrix::Matrix(0, nrow=J, ncol=K, sparse=TRUE)
  reg.topic.count = Matrix::Matrix(0, nrow=I, ncol=K, sparse=TRUE)
  v.topic.count = Matrix::Matrix(0, nrow=K, ncol=1, sparse=TRUE)
  z.topic.count = Matrix::Matrix(0, nrow=K, ncol=1, sparse=TRUE)
  lapply(1:D, function(d){
    j.idx = which(W[d,] > 0)
    lapply(j.idx, function(j){
      Vdm[[d, j]] <<- c(sample(1:K, size=W[d,j], replace=TRUE))
      lapply(Vdm[[d, j]], function(k){
        doc.topic.count[d,k] <<- doc.topic.count[d,k] + 1
        wrd.topic.count[j,k] <<- wrd.topic.count[j,k] + 1
        v.topic.count[k] <<- v.topic.count[k] + 1
      })
    })
    i.idx = which(R[d,] > 0)
    lapply(i.idx, function(i){
      Zdn[[d, i]] <<- c(sample(1:K, size=R[d,i], replace=TRUE))
      lapply(Zdn[[d,i]], function(k){
        doc.topic.count[d,k] <<- doc.topic.count[d,k] + 1
        reg.topic.count[i,k] <<- reg.topic.count[i,k] + 1
        z.topic.count[k] <<- z.topic.count[k] + 1
      })
    })
  })
  # Define a single Gibbs step as a function (keep the <<-!)
  gibbs.step <- function(){
    # For each document d...
    lapply(1:D, function(d){
      # For each word in doc d
      j.idx = which(W[d,] > 0)
      lapply(j.idx, function(j){
        lapply(1:W[d,j], function(w){
          # discount current topic
          cur_topic = Vdm[[d,j]][w]
          doc.topic.count[d,cur_topic] <<- doc.topic.count[d,cur_topic] - 1
          wrd.topic.count[j,cur_topic] <<- wrd.topic.count[j,cur_topic] - 1
          v.topic.count[cur_topic] <<- v.topic.count[cur_topic] - 1
          # generate new topic
          pv = as.vector((wrd.topic.count[j,] + mu[j])*(doc.topic.count[d,] + alp)/(v.topic.count + (J * mu[j])))
          new_topic = which.max(stats::rmultinom(1, size=1, prob=pv/sum(pv)))
          Vdm[[d,j]][w] <<- new_topic
          # add new topic to counts
          doc.topic.count[d,new_topic] <<- doc.topic.count[d,new_topic] + 1
          wrd.topic.count[j,new_topic] <<- wrd.topic.count[j,new_topic] + 1
          v.topic.count[new_topic] <<- v.topic.count[new_topic] + 1
        })
      })
      # For each word in doc d
      i.idx = which(R[d,] > 0)
      lapply(i.idx, function(i){
        lapply(1:R[d,i], function(r){
          # discount current topic
          cur_topic = Zdn[[d,i]][r]
          doc.topic.count[d,cur_topic] <<- doc.topic.count[d,cur_topic] - 1
          reg.topic.count[i,cur_topic] <<- reg.topic.count[i,cur_topic] - 1
          z.topic.count[cur_topic] <<- z.topic.count[cur_topic] - 1
          # generate new topic
          pz = as.vector((reg.topic.count[i,] + sig[i])*(doc.topic.count[d,] + alp)/(z.topic.count + (I * sig[i])))
          new_topic = which.max(stats::rmultinom(1, size=1, prob=pz/sum(pz)))
          Zdn[[d,i]][r] <<- new_topic
          # add new topic to counts
          doc.topic.count[d,new_topic] <<- doc.topic.count[d,new_topic] + 1
          reg.topic.count[i,new_topic] <<- reg.topic.count[i,new_topic] + 1
          z.topic.count[new_topic] <<- z.topic.count[new_topic] + 1
        })
      })
    })
  }
  # Begin Gibbs sampling
  theta = matrix(nrow = iters, ncol = D*K)
  beta = matrix(nrow = iters, ncol = I*K)
  phi = matrix(nrow = iters, ncol = J*K)
  iter.time = matrix(nrow = iters, ncol = 1)
  gibbs.step.time = matrix(nrow = iters, ncol = 1)
  lapply(1:iters, function(it){
    if (verbose) {print(paste0("Iteration ", it, " out of ", iters))}
    init.time = proc.time()[3]
    # Take a single Gibbs step
    gibbs.step()
    gibbs.step.time[it] <- proc.time()[3] - init.time
    if (verbose) {print(paste0("This Gibbs step took ", gibbs.step.time[it], " seconds."))}
    # Calcualte theta, beta, and phi
    theta[it,] <<- as.vector(calc.theta(doc.topic.count, alp))
    beta[it,] <<- as.vector(calc.beta(reg.topic.count, sig))
    phi[it,] <<- as.vector(calc.phi(wrd.topic.count, mu))
    iter.time[it] <<- proc.time()[3] - init.time
    if (verbose) {print(paste0("This iteration took ", iter.time[it], " seconds."))}
  })
  if (verbose) {print(paste0("Return results..."))}
  # Return results
  return(list(theta=theta, beta=beta, phi=phi, doc.topic.count=doc.topic.count, 
             reg.topic.count=reg.topic.count, wrd.topic.count=wrd.topic.count, 
             time=list(iter=iter.time, gibb.step=gibbs.step.time)))
}

initialize.counts <- function(K, R, W){
  # This function initializes the count matrices
  D = dim(W)[1]
  J = dim(W)[2]
  I = dim(R)[2]
  doc.topic.count = Matrix::Matrix(0, nrow=D, ncol=K, sparse=TRUE)
  wrd.topic.count = Matrix::Matrix(0, nrow=J, ncol=K, sparse=TRUE)
  reg.topic.count = Matrix::Matrix(0, nrow=I, ncol=K, sparse=TRUE)
  v.topic.count = Matrix::Matrix(0, nrow=K, ncol=1, sparse=TRUE)
  z.topic.count = Matrix::Matrix(0, nrow=K, ncol=1, sparse=TRUE)
  lapply(1:D, function(d){
    j.idx = which(W[d,] > 0)
    lapply(j.idx, function(j){
      Vdm[[d, j]] <<- c(sample(1:K, size=W[d,j], replace=TRUE))
      lapply(Vdm[[d, j]], function(k){
        doc.topic.count[d,k] <<- doc.topic.count[d,k] + 1
        wrd.topic.count[j,k] <<- wrd.topic.count[j,k] + 1
        v.topic.count[k] <<- v.topic.count[k] + 1
      })
    })
    i.idx = which(R[d,] > 0)
    lapply(i.idx, function(i){
      Zdn[[d, i]] <<- c(sample(1:K, size=R[d,i], replace=TRUE))
      lapply(Zdn[[d,i]], function(k){
        doc.topic.count[d,k] <<- doc.topic.count[d,k] + 1
        reg.topic.count[i,k] <<- reg.topic.count[i,k] + 1
        z.topic.count[k] <<- z.topic.count[k] + 1
      })
    })
  })
  out = list(doc.topic.count=doc.topic.count, )
}

calc.theta <- function(doc.topic.count, alp){
  # This function calculates theta
  # 
  # Args:
  #   doc.topic.count: (D x K) vector of document-topic counts
  #   alpha: (K x 1) hyperparameter vector
  #
  # Returns:
  #   (K x D) matrix
  D = dim(doc.topic.count)[1]
  K = dim(doc.topic.count)[2]
  s = rowSums(doc.topic.count)
  return(t(sapply(1:D, function(d){(doc.topic.count[d,] + alp)/(K*alp + s[d])})))
}

calc.beta <- function(reg.topic.count, sig){
  # This function calcualtes beta
  # 
  # Args:
  #   reg.topic.count:  (I x K) vector of image region-topic counts
  #   sig: (I x 1) hyperparamter vector
  # 
  # Returns:
  #   (I x K) matrix
  I = dim(reg.topic.count)[1]
  K = dim(reg.topic.count)[2]
  s = colSums(reg.topic.count)
  return(sapply(1:K, function(k){(reg.topic.count[,k] + sig)/(I*sig + s[k])}))
}

calc.phi <- function(wrd.topic.count, mu){
  # This function calculates phi
  # 
  # Args:
  #   wrd.topic.count : (J x K) vector of word-topic counts
  #   mu: (J x 1) hyperparamter vector
  # 
  # Returns:
  #   (J x K) matrix
  J = dim(wrd.topic.count)[1]
  K = dim(wrd.topic.count)[2]
  s = colSums(wrd.topic.count)
  return(sapply(1:K, function(k){(wrd.topic.count[,k] + mu)/(J*mu + s[k])}))
}



