# HEADER------------------------------------------------------------------------
# Author: Kelly Geyer, klgeyer@bu.edu
# Date: February 7, 2020
# Development: R v4.0.0 on MacOS v10.15.4
#
# Description: This script contains functions for implementing the MM-LDA topic 
# model, using STAN to implement a Gibbs sampler.
#
# References
# 1. LDA Model: https://mc-stan.org/docs/2_23/stan-users-guide/latent-dirichlet-allocation.html
# 2. Threading: https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html

# TODO more robust parameter checking within functions
# TODO find topic coherence
# TODO Rubin-Geman Diagnostic plots


# LIBRARIES---------------------------------------------------------------------
library(Matrix)
library(MCMCpack)
library(ggplot2)
library(ggpubr)
library(ggmcmc)
library(parallel)
library(rstan)

# FUNCTIONS TO FIT MODEL--------------------------------------------------------
mmlda.stan.gibbs.sampler <- function(K, W, R, sig=NA, mu=NA, alp=NA, iters=1000, 
                                     ncores=1, nchains=4, verbose=TRUE, 
                                     seed=123, codedir=getwd()){
  # This function performs Gibbs sampling for the MM-LDA model
  #
  # Args:
  #   K (>= 2): Number of factors
  #   W: (D x J)-sparse matrix of document-vocabulary count
  #   R: (D x I)-sparse amtrix of document-region count
  #   sigma (default 1's): concentration parameter of
  #   mu (default 1's): concentration parameters of
  #   alp (default 1's): concentration parameters of
  #   iters: Number of MCMC iterations
  #   nchains: Number of independent MCMC nchains
  #   verbose: Bool indicating wheter to print messages
  #   seed: Random seed
  #   codedir: Code directory
  #
  # Returns list containing:
  #   model: compiled stan model
  #   data: data used for sampling
  #   fit: Gibbs sampling results
  #   time: List object of time profiles
  # Set parameters
  D = dim(R)[1]     # number of docs
  I = dim(R)[2]     # unique img. regions
  J = dim(W)[2]     # vocabulary size
  # Check input
  if (K < 2){stop("'K' must be at least 2")}
  if (any(is.na(alp))) {alp = rep(1, K)}
  if (any(is.na(mu))) {mu = rep(1, J)}
  if (any(is.na(sig))) {sig = rep(1, I)}
  if (length(alp) != K){alp = rep(alp[1], K)}
  if (length(mu) != J){mu = rep(mu[1], J)}
  if (length(sig) != I){sig = rep(sig[1], I)}
  if (any(alp <= 0)) {stop("'alpha' must be greater than 0.")}
  if (any(mu <= 0)) {stop("'mu' must be greater than 0.")}
  if (any(sig <= 0)) {stop("'sig' must be greater than 0.")}
  # # Remove any complied stan files to avoid issues
  # if (file.exists(file.path(codedir, "source", "mmlda.rds"))){
  #   file.remove(file.path(codedir, "source", "mmlda.rds"))
  # }
  # Allocate cores for the stan model
  if (ncores == 0){
    options(mc.cores = parallel::detectCores())
  }else{
    options(mc.cores = ncores)
  }
  rstan::rstan_options(auto_write = TRUE)
  # Compile the stan model 
  if (verbose) {print("Compiling stan model...")}
  init.time = proc.time()[3]
  mod = rstan::stan_model(file.path(codedir, 'source','mmlda.stan'))
  compile.time = proc.time()[3] - init.time
  if (verbose) {print(paste0("That took ", compile.time, " seconds."))}
  # Store data for the stan model
  if (verbose){print("Formatting data for the stan model...")}
  init.time = proc.time()[3]
  w.vals = get.instances(W)
  r.vals = get.instances(R)
  format.time = proc.time()[3] - init.time
  if (verbose) {print(paste0("That took ", format.time, " seconds."))}
  # Sample model for posteriors
  if (verbose) {print("Estimating posterior with Gibbs sampling...")}
  init.time = proc.time()[3]
  fit = rstan::sampling(object = mod, iter = iters, chains = nchains,
                        verbose=verbose, seed=seed, 
                        data = list(K=K, D=D, J=J, I=I, M=w.vals$total.count, 
                                    N=r.vals$total.count, 
                                    featw=w.vals$feat.id, featr=r.vals$feat.id, 
                                    docw=w.vals$doc.id, docr=r.vals$doc.id,
                                    alpha=drop(alp), mu=drop(mu), 
                                    sigma=drop(sig))) 
  samp.time = proc.time()[3] - init.time
  rm(w.vals, r.vals)
  if (verbose) {print(paste0("That took ", samp.time, " seconds."))}
  # return everything
  t = list(compile=compile.time, fitting=samp.time, format=format.time)
  p = list(K=K, D=D, I=I, J=J)
  return(list(model=mod, fit=fit, time=t, params=p))
}

mmlda.stan.variational.bayes <- function(K, W, R, sig=NA, mu=NA, alp=NA, 
                                         alg="meanfield", iters=10000, 
                                         verbose=TRUE, seed=123, 
                                         codedir=getwd()){
  # This function performs variational Bayes estimation for the MM-LDA model
  #
  # Args:
  #   K (>= 2): Number of factors
  #   W: (D x J)-sparse matrix of document-vocabulary count
  #   R: (D x I)-sparse amtrix of document-region count
  #   sigma (default 1's): concentration parameter of topic-words
  #   mu (default 1's): concentration parameters of topic-images
  #   alp (default 1's): concentration parameters of topic-docs
  #   alg: VB algorithm, either 'meanfield' or 'fullrank' (see rstan doc)
  #   iters: Number of VB iterations
  #   verbose: Bool indicating wheter to print messages
  #
  # Returns list containing:
  #   model: compiled stan model
  #   data: data used for sampling
  #   fit: VB results
  #   time: List object of time profiles
  # Set parameters
  D = dim(R)[1]     # number of docs
  I = dim(R)[2]     # unique img. regions
  J = dim(W)[2]     # vocabulary size
  # Check input
  if (K < 2){stop("'K' must be at least 2")}
  if (any(is.na(alp))) {alp = rep(1, K)}
  if (any(is.na(mu))) {mu = rep(1, J)}
  if (any(is.na(sig))) {sig = rep(1, I)}
  if (length(alp) != K){alp = rep(alp[1], K)}
  if (length(mu) != J){mu = rep(mu[1], J)}
  if (length(sig) != I){sig = rep(sig[1], I)}
  if (any(alp <= 0)) {stop("'alpha' must be greater than 0.")}
  if (any(mu <= 0)) {stop("'mu' must be greater than 0.")}
  if (any(sig <= 0)) {stop("'sig' must be greater than 0.")}
  # # Remove any complied stan files to avoid issues
  # if (file.exists(file.path("source", "mmlda.rds"))){
  #   file.remove(file.path("source", "mmlda.rds"))
  # }
  # Allocate cores for the Stan model
  options(mc.cores = 1)
  rstan::rstan_options(auto_write=TRUE)
  # Compile the stan model 
  if (verbose) {print("Compiling stan model...")}
  init.time = proc.time()[3]
  mod = rstan::stan_model(file.path(codedir, 'source','mmlda.stan'))
  compile.time = proc.time()[3] - init.time
  if (verbose) {print(paste0("That took ", compile.time, " seconds."))}
  # Store data for the stan model
  init.time = proc.time()[3]
  w.vals = get.instances(W)
  r.vals = get.instances(R)
  format.time = proc.time()[3] - init.time
  # Sample model for posteriors
  if (verbose) {print(paste0("Estimating the posterior with VB + ",alg,"..."))}
  init.time = proc.time()[3]
  fit = rstan::vb(object=mod, 
                  data=list(K=K, D=D, J=J, I=I,
                            M=w.vals$total.count, N=r.vals$total.count, 
                            featw=w.vals$feat.id, featr=r.vals$feat.id, 
                            docw=w.vals$doc.id, docr=r.vals$doc.id, 
                            alpha=drop(alp), mu=drop(mu), sigma=drop(sig)),
                  seed=seed, iter=iters, algorithm=alg)  
  samp.time = proc.time()[3] - init.time
  rm(w.vals, r.vals)
  if (verbose) {print(paste0("That took ", samp.time, " seconds."))}
  # return everything
  t = list(compile=compile.time, fitting=samp.time, format=format.time)
  p = list(K=K, D=D, I=I, J=J)
  return(list(model=mod, fit=fit, time=t, params=p))
}

get.instances <- function(W){
  # This extracts coordinate and values of nonzero values in W
  #
  # Args:
  #   W - (D x J) dgCMatrix/matrix of word counts
  #
  # Return:
  #   doc.id - vector of doc IDs instances
  #   feat.id - vector of feature IDs instances
  #   total.count - total number of feature instances
  D = dim(W)[1]
  J = dim(W)[2]
  total.count = sum(W)
  doc.id = rep(0, total.count)
  feat.id = rep(0, total.count)
  extract <- Matrix::summary(W)
  nonzero.entries.count <- Matrix::nnzero(W)
  idx = 1
  for (ii in 1:nonzero.entries.count){
    doc.id[idx:(idx + extract$x[ii] - 1)] = rep(extract$i[ii], extract$x[ii])
    feat.id[idx:(idx + extract$x[ii] - 1)] = rep(extract$j[ii], extract$x[ii])
    idx = idx + extract$x[ii]
  }
  return(list(doc.id=doc.id, feat.id=feat.id, total.count=total.count))
}


# FUNCTIONS TO EXTRACT PARAMETERS-----------------------------------------------
calculate_parameter <- function(resfit, chain=1, param.name='theta', param.est='maxll'){
  # This function calculates the parameter for each MCMC chain
  #
  # Inputs
  # all.param - 3d array: with dimensions (iterations, MCMC nchains, parameters)
  # param.name - name of parameter: 'theta', 'beta' or 'phi'
  # param.est - name of method to estimate parameter: 'avg', 'maxll', or 'lastmcmc'
  all.param = rstan::extract(resfit, par = param.name, permuted = FALSE)
  niters = dim(all.param)[1]
  nchains = dim(all.param)[2]
  nparams = dim(all.param)[3]
  if (param.est == 'avg'){
    out <- apply(all.param[,chain,], 2, mean)
  } else if (param.est == 'maxll') {
    ll = extract_loglikelihood(resfit)
    max.idx = apply(ll, 2, which.max)
    out <-all.param[max.idx[chain],chain,]
  } else if (param.est == 'lastmcmc') {
    out <- all.param[niters,chain,]
  } else {
    stop(paste0("The value param.est= '", param.est, "' is not recognized."))
  }
  return(out)
}

mmlda_greedy_topic_permute <- function(true.theta, true.beta, true.phi, 
                                       est.theta, est.beta, est.phi, 
                                       codedir=getwd()){
  # This function returns indexes to align the estimated topics with the truth.
  # For each parameter, cosine similarity is computed by topic
  # The most commonly selected topic becomes the new assignment
  # In the case of a tie, only the ambiguous topics are randomly assigned 
  #  
  # Inputs:
  #   true.theta, true.beta, true.phi - the true parameters
  #   est.theta, est.beta, est.phi - the estimated parameters
  #   resfit - the result$fit (rstan::stanfit object) object from the function
  #            stan.gibbs.sampler
  #   param.est - method of calculating the posterior
  #
  # Returns
  source(file.path(code_dir, "source", "cosine_similarity.R"))
  #   (nchains x K) matrix of topic reassigments, to plug into estimated params
  D = dim(true.theta)[1]
  K = dim(true.theta)[2]
  J = dim(true.phi)[2]
  I = dim(true.beta)[2]
  # Calculate cosine similarity scores
  theta.CS = sapply(1:K, function(x){apply(true.theta, 2, calculate_cosine_similarity, vec2=est.theta[,x])})
  beta.CS = sapply(1:K, function(x){apply(true.beta, 1, calculate_cosine_similarity, vec2=est.beta[x,])})
  phi.CS = sapply(1:K, function(x){apply(true.phi, 1, calculate_cosine_similarity, vec2=est.phi[x,])})
  mean.CS = apply(array(unlist(c(theta.CS, beta.CS, phi.CS)), c(K, K, 3)), 1:2, mean)
  newlabs = apply(mean.CS, 1, which.max)
  # Deal with duplicate values, set bad values to NA
  dup_topics = newlabs[duplicated(newlabs)]
  for (ii in dup_topics){
    jj = which(newlabs == ii)
    jj.val = mean.CS[newlabs==ii, ii]
    win.idx = jj[which.max(jj.val)]
    newlabs[jj[which.max(jj.val)]] = ii
    lose.idx = jj[jj != win.idx]
    newlabs[lose.idx] = NA
  }
  # Deal with unassigned labels
  unassined.vals = setdiff(1:K, newlabs)
  while(length(unassined.vals) > 0){
    undef.newlabs = which(is.na(newlabs))
    vv = max(as.vector(mean.CS[undef.newlabs, unassined.vals]))
    vv.idx = which(mean.CS == vv, arr.ind=TRUE)
    newlabs[vv.idx[1]] = vv.idx[2]
    unassined.vals = setdiff(1:K, newlabs)
  }
  return(newlabs)
}

# FUNCTIONS TO EVALUATE MODEL---------------------------------------------------
calc_convergence_measure <- function(resfit, fun = rstan::Rhat){
  # This function calcualtes either Rhat, Bulk Effective Sample Size (ESS), or 
  # Tail ESS, for each parameter estimated by the result of stan Gibbs sampler 
  # of MM-LDA.
  #
  # Input:
  #   resfit - the result$fit (rstan::stanfit object) object from the function 
  #            stan.gibbs.sampler
  #   fun - Measure function. May be 'rstan::Rhat' (default), 
  #           'rstan::ess_buld', or 'rstan::ess_tail'
  all.theta = rstan::extract(resfit, par = 'theta', permuted = FALSE)
  theta.meas = matrix(nrow=dim(all.theta)[3], ncol=1)
  theta.meas = sapply(1:dim(all.theta)[3], 
                          function(i){return(fun(all.theta[,,i]))})
  names(theta.meas) = dimnames(all.theta)$parameters
  all.beta = rstan::extract(resfit, par = 'beta', permuted = FALSE)
  beta.meas = matrix(nrow=dim(all.beta)[3], ncol=1)
  beta.meas = sapply(1:dim(all.beta)[3], 
                         function(i){return(fun(all.beta[,,i]))})
  names(beta.meas) = dimnames(all.beta)$parameters
  all.phi = rstan::extract(resfit, par = 'phi', permuted = FALSE)
  phi.meas = matrix(nrow=dim(all.phi)[3], ncol=1)
  phi.meas = sapply(1:dim(all.phi)[3], 
                        function(i){return(fun(all.phi[,,i]))})
  names(phi.meas) = dimnames(all.phi)$parameters
  return(list(theta=theta.meas, phi=phi.meas, beta=beta.meas))
}

make_conv_stat_plot <- function(resfit, plt_fn, plt_title){
  # This function creates and saves a figure summarizing the convergence 
  # statistics Rhat,Bulk Effective Sample Size (ESS), and Tail ESS. These 
  # statistics are calcualted for each parameter, and each subfigure displays 
  # one statistic.
  #
  # Inputs:
  #   resfit - the result$fit (rstan::stanfit object) object from the function 
  #            stan.gibbs.sampler
  #   plt_fn - filename of plot
  #   plt_title - main title of plot
  # Calculate Rhat stats
  Rhat.vals = calc_convergence_measure(resfit, fun=rstan::Rhat)
  # Calcualte bulk ESS values
  bulk.ess.vals = calc_convergence_measure(resfit, fun=rstan::ess_bulk)
  # Calculate tail ESS values
  tail.ess.vals = calc_convergence_measure(resfit, fun=rstan::ess_tail)
  # Create plots
  X = data.frame(Rhat=c(Rhat.vals$theta, Rhat.vals$beta, Rhat.vals$phi), 
                 bulkESS=c(bulk.ess.vals$theta, bulk.ess.vals$beta, bulk.ess.vals$phi),
                 tailESS=c(tail.ess.vals$theta, tail.ess.vals$beta, tail.ess.vals$phi),
                 lab=c(rep("Theta", length(Rhat.vals$theta)), 
                       rep("Beta", length(Rhat.vals$beta)), 
                       rep("Phi", length(Rhat.vals$phi))))
  p1 <- ggplot2::ggplot(X, aes(x=Rhat, col=lab, fill=lab)) + 
    ggplot2::geom_density(alpha = 0.2) + 
    ggplot2::ggtitle("R-hat Statistics") + 
    ggplot2::xlab("R-hat") + ggplot2::ylab("Density")
  p2 <- ggplot2::ggplot(X, aes(x=bulkESS, col=lab, fill=lab)) + 
    ggplot2::geom_density(alpha = 0.2) + 
    ggplot2::ggtitle("Bulk Effective Sample Size") + 
    ggplot2::xlab("Bulk Effective Sample Size") + ggplot2::ylab("Density")
  p3 <- ggplot2::ggplot(X, aes(x=tailESS, col=lab, fill=lab)) + 
    ggplot2::geom_density(alpha = 0.2) + 
    ggplot2::ggtitle("Tail Effective Sample Size") + 
    ggplot2::xlab("Tail Effective Sample Size") + ggplot2::ylab("Density")
  plt <- ggpubr::ggarrange(p1, p2, p3, nrow=1, ncol=3, common.legend = TRUE, legend = "bottom")
  plt <- ggpubr::annotate_figure(plt, top=plt_title)
  ggplot2::ggsave(plt_fn, plt, width=11, height=11, units="in")
}

make_geweke_plot <- function(resfit, plt_fn, plt_title, width=11, height=30){
  # This function creates and saves a figure summarizing the Geweke scores for 
  # each parameter
  #
  # Inputs:
  #   resfit - the result$fit (rstan::stanfit object) object from the function 
  #            stan.gibbs.sampler
  #   plt_fn - filename of plot
  #   plt_title - main title of plot
  # Convert 'stanfit' object to be compatible with the package ggmcmc
  plt <- ggmcmc::ggs_geweke(ggmcmc::ggs(resfit))
  plt <- ggpubr::annotate_figure(plt, top=plt_title)
  # Save Geweke diganostic plot
  ggplot2::ggsave(plt_fn, plt, width=width, height=height, units="in", 
                  limitsize=FALSE)
}

extract_loglikelihood <- function(resfit){
  # This function extracts the log-likelihood values from model output
  #
  # Inputs:
  #   resfit - the result$fit (rstan::stanfit object) object from the function 
  #            stan.gibbs.sampler
  # Returns:
  #   (niter x nchains) matrix of log-likelihood probabilities
  return(drop(rstan::extract(resfit, par = "lp__", permuted = FALSE)))
}

make_vb_conv_stat_plot <- function(resfit, plt_fn, plt_title){
  # This function creates and saves a figure summarizing the convergence 
  # statistics Pareto-k (good values are < 0.7). These statistics are calcualted 
  # for each parameter.
  #
  # Inputs:
  #   resfit - the result$fit (rstan::stanfit object) object from the function 
  #            stan.variational.bayes()
  #   plt_fn - filename of plot
  #   plt_title - main title of plot
  theta.all = rstan::summary(resfit, pars='theta')
  beta.all = rstan::summary(resfit, pars='beta')
  phi.all = rstan::summary(resfit, pars='phi')
  X = data.frame(khat=c(theta.all$summary[,'khat'], beta.all$summary[,'khat'], 
                        phi.all$summary[,'khat']),
                 parameter=c(rep("Theta", dim(theta.all$summary)[1]), 
                             rep("Beta", dim(beta.all$summary)[1]), 
                             rep("Phi", dim(phi.all$summary)[1])))
  plt <- ggplot2::ggplot(X, aes(x=khat, col=parameter, fill=parameter)) + 
    ggplot2::geom_density(alpha = 0.2) + 
    ggplot2::ggtitle("Pareto-k Diagnostic") + 
    ggplot2::xlab("khat") + ggplot2::ylab("Density")
  plt <- ggpubr::annotate_figure(plt, top=plt_title)
  ggplot2::ggsave(plt_fn, plt, width=7, height=9, units="in")
}

make_mmlda_cosine_plots <- function(true.theta, true.beta, true.phi, est.theta, 
                                    est.beta, est.phi, plt_fn, plt_title, 
                                    codedir=getwd()){
  # 6 plots in one
  # Assumes that the parameters have been matched and sorted by topic
  # load source files
  source(file.path(code_dir, "source", "cosine_similarity.R"))
  # Make topic-CS plots
  theta_title = "Theta: Doc-Topic Probabilities"
  p1 <- make_topic_cs_barplot(true.param=true.theta, est.param=est.theta, plt_title=theta_title, param.name='theta')
  beta_title = "Beta: Img-Topic Probabilities"
  p2 <- make_topic_cs_barplot(true.param=true.beta, est.param=est.beta, plt_title=beta_title, param.name='beta')
  phi_title = "Phi: Word-Topic Probabilities"
  p3 <- make_topic_cs_barplot(true.param=true.phi, est.param=est.phi, plt_title=phi_title, param.name='phi')
  topic_plt <- ggpubr::ggarrange(p1, p2, p3, nrow=1, ncol=3)
  topic_plt <- ggpubr::annotate_figure(topic_plt, left="Topic Cosine Similarity", fig.lab.face="bold")
  # Make feat-CS plots
  p4 <- make_feat_cs_barplot(true.param=true.theta, est.param=est.theta, plt_title=theta_title, param.name='theta')
  p5 <- make_feat_cs_barplot(true.param=true.beta, est.param=est.beta, plt_title=beta_title, param.name='beta')
  p6 <- make_feat_cs_barplot(true.param=true.phi, est.param=est.phi, plt_title=phi_title, param.name='phi')
  feat_plt <- ggpubr::ggarrange(p4, p5, p6, nrow=1, ncol=3)
  feat_plt <- ggpubr::annotate_figure(feat_plt, left="Feature Cosine Similarity", fig.lab.face="bold")
  # Combine 
  plt <- ggpubr::ggarrange(topic_plt, feat_plt, nrow=2, ncol=1)
  plt <- ggpubr::annotate_figure(plt, top=plt_title, fig.lab.face = "bold")
  ggplot2::ggsave(plt_fn, plt, width=11, height=11, units="in")
}

# OLD FUNCTIONS-----------------------------------------------------------------
mmlda.stan.optimizer <- function(K, W, R, seed, alg, sig = 1, mu=1, alp = 1, iters = 2000, verbose=TRUE){
  # This function gets a point estimates of the posteriors by maximizing the 
  # joint posterior from the MM-LDA model.
  # It uses the algorithm Limited-memory BFGS.
  #
  # Args:
  #   K (>= 2): Number of factors
  #   W: (D x J)-sparse matrix of document-vocabulary count
  #   R: (D x I)-sparse amtrix of document-region count
  #   seed: Random seed
  #   alg: optimization algorithm ("LBFGS", "BFGS", or "Newton")
  #   sigma (default 1's): concentration parameter of
  #   mu (default 1's): concentration parameters of
  #   alp (default 1's): concentration parameters of
  #   iters: Number of MCMC iterations
  #   verbose: Bool indicating wheter to print messages
  #
  # Returns list containing:
  #   model: compiled stan model
  #   data: data used for sampling
  #   fit: Optimizer results
  #   time: List object of time profiles
  # Set parameters
  D = dim(R)[1]     # number of docs
  I = dim(R)[2]     # unique img. regions
  J = dim(W)[2]     # vocabulary size
  # Check input
  if (K < 2){stop("'K' must be at least 2")}
  if ((!is.null(alp)) & !is.numeric(alp) & (alp <= 0)){stop("'alp' must be a positive value")}
  if ((!is.null(sig)) & !is.numeric(sig) & (sig <= 0)){stop("'sig' must be a positive value")}
  if ((!is.null(mu)) & !is.numeric(mu) & (mu <= 0)){stop("'mu' must be a positive value")}
  alp = rep(alp, K)
  sig = rep(sig, I)
  mu = rep(mu, J)
  # Remove any complied stan files to avoid issues
  if (file.exists(file.path("source", "mmlda.rds"))){
    file.remove(file.path("source", "mmlda.rds"))
  }
  # Compile the stan model 
  if (verbose) {print("Compiling stan model...")}
  init.time = proc.time()[3]
  mod = rstan::stan_model(file.path('source','mmlda.stan'))
  compile.time = proc.time()[3] - init.time
  if (verbose) {print(paste0("That took ", compile.time, " seconds."))}
  # Store data for the stan model
  w.vals = get.instances(W)
  r.vals = get.instances(R)
  dat = list(K=K, D=D, J=J, I=I, M=w.vals$M, N=r.vals$M, featw=w.vals$feat.id, 
             featr=r.vals$feat.id, docw=w.vals$doc.id, docr=r.vals$doc.id, 
             countw=w.vals$count, countr=r.vals$count, 
             alpha=alp, mu=mu, sigma=sig)
  # Sample model for posteriors
  if (verbose) {print(paste0("Estimating the posterior with ",alg,"..."))}
  init.time = proc.time()[3]
  fit = rstan::optimizing(object=mod, data=dat, iter=iters, verbose=verbose, seed=seed, algorithm=alg)
  samp.time = proc.time()[3] - init.time
  if (verbose) {print(paste0("That took ", samp.time, " seconds."))}
  # return everything
  t = list(compile=compile.time, fitting=samp.time, total=(compile.time+samp.time))
  return(list(model=mod, data=dat, fit=fit, time=t))
}

Mode <- function(x){
  # Find the mode of a set of data
  #
  # Input:
  #   x - vector of numeric values
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

greedy_topic_permute.allchains <- function(true.theta, true.beta, true.phi, resfit, param.est="maxll"){
  # This function returns indexes to align the estimated topics with the truth.
  # For each parameter, cosine similarity is computed by topic
  # The most seleted topic becomes the new assignment
  # In the case of a tie, only the ambigous topics are randomly assigned 
  #  
  # Inputs:
  #   true.theta, true.beta, true.phi - the true parameters
  #   resfit - the result$fit (rstan::stanfit object) object from the function
  #            stan.gibbs.sampler
  #   param.est - method of calculating the posterior
  #
  # Returns
  #   (nchains x K) matrix of topic reassigments, to plug into estimated params
  D = dim(true.theta)[1]
  K = dim(true.theta)[2]
  J = dim(true.phi)[2]
  I = dim(true.beta)[2]
  # Estimate the parameters
  est.theta.allnchains = calculate_parameter(resfit, 'theta', param.est)
  est.beta.allnchains = calculate_parameter(resfit, 'beta', param.est)
  est.phi.allnchains = calculate_parameter(resfit, 'phi', param.est)
  nchains = dim(est.theta.allnchains)[2]
  # Find new topic assignments
  new.topic.idx = matrix(nrow=nchains, ncol=K)
  for (cc in 1:nchains){  # for each chain...
    theta.CS = matrix(0, nrow = K, ncol = K)
    beta.CS = matrix(0, nrow = K, ncol = K)
    phi.CS = matrix(0, nrow = K, ncol = K)
    est.theta <- matrix(est.theta.allnchains[,cc], nrow = D, ncol = K)
    est.beta <- matrix(est.beta.allnchains[,cc], nrow = K, ncol = I)
    est.phi <- matrix(est.phi.allnchains[,cc], nrow=K, ncol = J)
    for (ii in 1:K){    # for each topic in the true param...
      for (jj in 1:K){    # for each topic in the est. param...
        # Calc CS
        theta.CS[ii,jj] = calculate_cosine_similarity(true.theta[,ii], est.theta[,jj])
        beta.CS[ii,jj] = calculate_cosine_similarity(true.beta[ii,], est.beta[jj,])
        phi.CS[ii,jj] = calculate_cosine_similarity(true.phi[ii,], est.phi[jj,])
      }
    }
    # Pick new topics
    cs.topics = matrix(nrow=3, ncol=K)
    cs.topics[1,] = apply(theta.CS, 1, which.max)
    cs.topics[2,] = apply(beta.CS, 1, which.max)
    cs.topics[3,] = apply(phi.CS, 1, which.max)
    newlabs = apply(cs.topics, 2, Mode)
    # if a topic is assigned more than once...random assignment
    if (length(unique(newlabs)) < K){
      a = table(newlabs)
      given.labs = strtoi(rownames(a))
      lab.to.reassign = setdiff(1:K, given.labs)              # unassigned topics
      lab.to.reassign = c(lab.to.reassign, given.labs[a > 1]) # topics assigned more than once
      bad.topic.idx = which(newlabs %in% lab.to.reassign)
      newlabs[bad.topic.idx] = sample(lab.to.reassign, size=length(lab.to.reassign), replace=FALSE)
    }
    new.topic.idx[cc,] = newlabs
    if (length(unique(newlabs)) < K){stop("BAD LABELS!!!!")}
  }
  return(new.topic.idx)
}

