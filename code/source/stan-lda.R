# HEADER------------------------------------------------------------------------
# Author: Kelly Geyer, klgeyer@bu.edu
# Date: September 30, 2020
# Development: R v4.0.2 on MacOS v10.15.6
#
# Description: This script contains functions implementing and evaluating LDA

# LIBRARIES---------------------------------------------------------------------
library(Matrix)
library(MCMCpack)
library(rstan)
library(ggplot2)
library(ggpubr)
library(ggmcmc)
library(parallel)

# FUNCTIONS TO FIT MODEL--------------------------------------------------------
lda.stan.gibbs.sampler <- function(K, W, mu=NA, alp=NA, iters=1000, ncores=1, 
                                   nchains=4, verbose=TRUE, seed=123, 
                                   codedir=getwd()){
  # This function performs Gibbs sampling for the LDA model
  #
  # Args:
  #   K (>= 2): Number of factors
  #   W: (D x J)-sparse matrix of document-vocabulary count
  #   mu (default 1's): concentration parameters of
  #   alp (default 1's): concentration parameters of
  #   iters: Number of MCMC iterations
  #   nchains: Number of independent MCMC nchains
  #   verbose: Bool indicating wheter to print messages
  #   seed: Random seed
  #
  # Returns list containing:
  #   model: compiled stan model
  #   data: data used for sampling
  #   fit: Gibbs sampling results
  #   time: List object of time profiles
  # Set parameters
  D = dim(W)[1]     # number of docs
  J = dim(W)[2]     # vocabulary size
  # Check input
  if (K < 2){stop("'K' must be at least 2")}
  if (any(is.na(alp))) {alp = rep(1, K)}
  if (any(is.na(mu))) {mu = rep(1, J)}
  if (length(alp) != K){alp = rep(alp[1], K)}
  if (length(mu) != J){mu = rep(mu[1], J)}
  if (any(alp <= 0)) {stop("'alpha' must be greater than 0.")}
  if (any(mu <= 0)) {stop("'mu' must be greater than 0.")}
  # Allocate cores for the Stan model
  if (ncores == 0){
    options(mc.cores = parallel::detectCores())
  }else{
    options(mc.cores = ncores)
  }
  rstan::rstan_options(auto_write = TRUE)
  # Compile the stan model 
  if (verbose) {print("Compiling stan model...")}
  init.time = proc.time()[3]
  mod = rstan::stan_model(file.path(codedir, 'source', 'lda.stan'))
  compile.time = proc.time()[3] - init.time
  if (verbose) {print(paste0("That took ", compile.time, " seconds."))}
  # Store data for the stan model
  if (verbose){print("Formatting data for the stan model...")}
  init.time = proc.time()[3]
  w.vals = lda.get.instances(W)
  format.time = proc.time()[3] - init.time
  if (verbose) {print(paste0("That took ", format.time, " seconds."))}
  # Sample model for posteriors
  if (verbose) {print("Estimating the posterior with Gibbs sampling...")}
  init.time = proc.time()[3]
  fit = rstan::sampling(object = mod, iter = iters, chains = nchains,
                        verbose=verbose, seed=seed, 
                        data=list(K=K, D=D, J=J, M=w.vals$total.count, 
                                  featw=w.vals$feat.id, docw=w.vals$doc.id, 
                                  alpha=drop(alp), mu=drop(mu))) 
  rm(w.vals)
  samp.time = proc.time()[3] - init.time
  if (verbose) {print(paste0("That took ", samp.time, " seconds."))}
  # return everything
  t = list(compile=compile.time, fitting=samp.time, format=format.time)
  return(list(model=mod, fit=fit, time=t))
}

lda.stan.variational.bayes <- function(K, W, mu=NA, alp=NA, alg="meanfield", 
                                       iters=10000, verbose=TRUE, seed=123, 
                                       codedir=getwd()){
  # This function performs variation Bayes estimation for the MM-LDA model
  #
  # Args:
  #   K (>= 2): Number of factors
  #   W: (D x J)-sparse matrix of document-vocabulary count
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
  D = dim(W)[1]     # number of docs
  J = dim(W)[2]     # vocabulary size
  # Check input
  if (K < 2){stop("'K' must be at least 2")}
  if (any(is.na(alp))) {alp = rep(1, K)}
  if (any(is.na(mu))) {mu = rep(1, J)}
  if (length(alp) != K){alp = rep(alp[1], K)}
  if (length(mu) != J){mu = rep(mu[1], J)}
  if (any(alp <= 0)) {stop("'alpha' must be greater than 0.")}
  if (any(mu <= 0)) {stop("'mu' must be greater than 0.")}
  # # Remove any complied stan files to avoid issues
  # if (file.exists(file.path("source", "lda.rds"))){
  #  file.remove(file.path("source", "lda.rds"))
  # }
  # Allocate cores for the Stan model
  options(mc.cores = 1)
  rstan::rstan_options(auto_write=TRUE)
  # Compile the stan model 
  if (verbose) {print("Compiling stan model...")}
  init.time = proc.time()[3]
  mod = rstan::stan_model(file.path(codedir,'source','lda.stan'))
  compile.time = proc.time()[3] - init.time
  if (verbose) {print(paste0("That took ",compile.time," seconds."))}
  # Store data for the stan model
  if (verbose){print("Formatting data for the stan model...")}
  init.time = proc.time()[3]
  w.vals = lda.get.instances(W)
  format.time = proc.time()[3] - init.time
  if (verbose) {print(paste0("That took ",format.time," seconds."))}
  # Sample model for posteriors
  if (verbose) {print(paste0("Estimating the posterior with VB + ",alg,"..."))}
  init.time = proc.time()[3]
  fit = rstan::vb(object = mod, 
                  data = list(K=K, D=D, J=J, M=w.vals$total.count, 
                              featw=w.vals$feat.id, docw=w.vals$doc.id, 
                              alpha=drop(alp), mu=drop(mu)), 
                  seed=seed, iter=iters, algorithm=alg)  
  samp.time = proc.time()[3] - init.time
  rm(w.vals)
  if (verbose) {print(paste0("That took ",samp.time," seconds."))}
  # return everything
  t = list(compile=compile.time, fitting=samp.time, format=format.time)
  return(list(model=mod, fit=fit, time=t))
}

lda.get.instances <- function(W){
  # This extracts coordinate and values of nonzero values in W
  #
  # Args:
  #   W - (D x J) dgCMatrix/matrix of word counts
  #
  # Return:
  #   doc.id - vector of doc IDs instances
  #   feat.id - vector of feature IDs instances
  #   total.count - total number of feature instances
  #print("ENTERS GET.INSTANCES :)")
  #print(paste("The class of W is ",class(W)))
  #print(W)
  D = dim(W)[1]
  #print(paste("D =",D))
  J = dim(W)[2]
  #print(paste("J =",J))
  total.count = sum(W)
  #print(paste("total count is", total.count))
  doc.id = rep(0, total.count)
  feat.id = rep(0, total.count)
  extract <- Matrix::summary(W)
  #print(extract)
  #nonzero.entries.count <- nrow(extract)
  #nonzero.entries.count1 <- dim(extract)[1]
  #nonzero.entries.count2 <- length(extract)
  #nonzero.entries.count3 <- nrow(extract[,1])
  nonzero.entries.count <- Matrix::nnzero(W)
  #print(paste("num. of nonzero entries:", nonzero.entries.count))
  #print(paste("num. of nonzero entries1:", nonzero.entries.count1))
  #print(paste("num. of nonzero entries2:", nonzero.entries.count2))
  #print(paste("num. of nonzero entries3:", nonzero.entries.count3))
  #print(paste("num. of nonzero entries4:", nonzero.entries.count4))
  idx = 1
  for (ii in 1:nonzero.entries.count){
    doc.id[idx:(idx + extract$x[ii] - 1)] = rep(extract$i[ii], extract$x[ii])
    feat.id[idx:(idx + extract$x[ii] - 1)] = rep(extract$j[ii], extract$x[ii])
    idx = idx + extract$x[ii]
  }
  #print(doc.id)
  #print(feat.id)
  return(list(doc.id=doc.id, feat.id=feat.id, total.count=total.count))
}

# FUNCTIONS TO EVALUATE MODEL---------------------------------------------------
lda_greedy_topic_permute <- function(true.theta, true.beta, est.theta, est.beta, codedir=getwd()){
  # This function returns indexes to align the estimated topics with the truth.
  # For each parameter, cosine similarity is computed by topic
  # The most commonly selected topic becomes the new assignment
  #
  # theta is DxK
  # beta is KxI
  #  
  # Inputs:
  #   true.theta (DxK), true.beta - the true parameters
  #   est.theta, est.beta - the estimated parameters
  source(file.path(code_dir, "source", "cosine_similarity.R"))
  K = dim(true.theta)[2]
  # Calculate cosine similarity scores
  theta.CS = sapply(1:K, function(x){apply(true.theta, 2, calculate_cosine_similarity, vec2=est.theta[,x])})
  beta.CS = sapply(1:K, function(x){apply(true.beta, 1, calculate_cosine_similarity, vec2=est.beta[x,])})
  mean.CS = apply(array(unlist(c(theta.CS, beta.CS)), c(K, K, 2)), 1:2, mean)
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

lda_calc_convergence_measure <- function(resfit, fun = rstan::Rhat){
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
  all.phi = rstan::extract(resfit, par = 'phi', permuted = FALSE)
  phi.meas = matrix(nrow=dim(all.phi)[3], ncol=1)
  phi.meas = sapply(1:dim(all.phi)[3], 
                    function(i){return(fun(all.phi[,,i]))})
  names(phi.meas) = dimnames(all.phi)$parameters
  return(list(theta=theta.meas, phi=phi.meas))
}

lda_make_conv_stat_plot <- function(resfit, plt_fn, plt_title){
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
  Rhat.vals = lda_calc_convergence_measure(resfit, fun=rstan::Rhat)
  # Calcualte bulk ESS values
  bulk.ess.vals = lda_calc_convergence_measure(resfit, fun=rstan::ess_bulk)
  # Calculate tail ESS values
  tail.ess.vals = lda_calc_convergence_measure(resfit, fun=rstan::ess_tail)
  # Create plots
  X = data.frame(Rhat=c(Rhat.vals$theta, Rhat.vals$phi), 
                 bulkESS=c(bulk.ess.vals$theta, bulk.ess.vals$phi),
                 tailESS=c(tail.ess.vals$theta, tail.ess.vals$phi),
                 lab=c(rep("Theta", length(Rhat.vals$theta)), 
                       rep("Beta", length(Rhat.vals$phi))))
  p1 <- ggplot2::ggplot(X, aes(x=Rhat, col=lab, fill=lab)) + 
    ggplot2::geom_density(alpha = 0.2) + 
    ggplot2::ggtitle("R-hat Statistics") + 
    xlab("R-hat") + ylab("Density")
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

lda_make_vb_conv_plot <- function(resfit, plt_fn, plt_title){
  # This function creates and saves a figure summarizing the convergence 
  # statistics Pareto-k (good values are < 0.7). These statistics are calcualted 
  # for each parameter.
  #
  # Inputs:
  #   resfit - the result$fit (rstan::stanfit object) object from the function 
  #            stan.variational.bayes()
  #   plt_fn - filename of plot
  #   plt_title - main title of plot
  theta.all = summary(resfit, pars='theta')
  phi.all = summary(resfit, pars='phi')
  X = data.frame(khat=c(theta.all$summary[,'khat'], phi.all$summary[,'khat']),
                 parameter=c(rep("Theta", dim(theta.all$summary)[1]), 
                             rep("Beta", dim(phi.all$summary)[1])))
  plt <- ggplot2::ggplot(X, aes(x=khat, col=parameter, fill=parameter)) + 
    ggplot2::geom_density(alpha = 0.2) + 
    ggplot2::ggtitle("Pareto-k Diagnostic") + 
    ggplot2::xlab("khat") + ggplot2::ylab("Density")
  plt <- ggpubr::annotate_figure(plt, top=plt_title)
  ggplt2::ggsave(plt_fn, plt, width=7, height=9, units="in")
}

make_lda_cosine_plots <- function(true.theta, true.beta, est.theta, est.beta,
                                  plt_fn, plt_title, codedir=getwd()){
  # 4 plots in one
  # Assumes that the parameters have been matched and sorted by topic
  # load source files
  source(file.path(code_dir, "source", "cosine_similarity.R"))
  # Make topic-CS plots
  theta_title = "Theta: Doc-Topic Probabilities"
  p1 <- make_topic_cs_barplot(true.param=true.theta, est.param=est.theta, plt_title=theta_title, param.name='theta')
  beta_title = "Beta: Img-Topic Probabilities"
  p2 <- make_topic_cs_barplot(true.param=true.beta, est.param=est.beta, plt_title=beta_title, param.name='beta')
  topic_plt <- ggpubr::ggarrange(p1, p2, nrow=1, ncol=2)
  topic_plt <- ggpubr::annotate_figure(topic_plt, left="Topic Cosine Similarity", fig.lab.face="bold")
  # Make feat-CS plots
  p4 <- make_feat_cs_barplot(true.param=true.theta, est.param=est.theta, plt_title=theta_title, param.name='theta')
  p5 <- make_feat_cs_barplot(true.param=true.beta, est.param=est.beta, plt_title=beta_title, param.name='beta')
  feat_plt <- ggpubr::ggarrange(p4, p5, nrow=1, ncol=2)
  feat_plt <- ggpubr::annotate_figure(feat_plt, left="Feature Cosine Similarity", fig.lab.face="bold")
  # Combine 
  plt <- ggpubr::ggarrange(topic_plt, feat_plt, nrow=2, ncol=1)
  plt <- ggpubr::annotate_figure(plt, top=plt_title, fig.lab.face = "bold")
  ggplot2::ggsave(plt_fn, plt, width=11, height=11, units="in")
}
