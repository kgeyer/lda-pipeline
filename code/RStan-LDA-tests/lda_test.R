#!/usr/bin/env Rscript

# HEADER------------------------------------------------------------------------
# Author: Kelly Geyer, klgeyer@bu.edu
# Date: September 26, 2020
# Development: R v4.0.2 on MacOS v10.15.6
# Description: This script runs topicmodel's implementation of LDA.

# TODO: Rstan nested parallelization

# SET UP SCRIPT-----------------------------------------------------------------
print("FIT TOPICMODELS'S LDA")
args <- commandArgs()
code_dir = args[6]
data_dir = args[7]
model_dir = args[8]
nReps = strtoi(args[9])
ncores = strtoi(args[10])
nchains = strtoi(args[11])
vb_iters = strtoi(args[12])
mcmc_iters = strtoi(args[13])
out_fn = args[14]
print(paste("The code directory is", code_dir))
print(paste("The data directory is", data_dir))
print(paste("The model directory is", model_dir))
print(paste("We will have", nReps, "replications."))
print(paste("We use", ncores, "cores."))
print(paste("We will run", nchains, "MCMC chains."))
print(paste("We will run", mcmc_iters, "for each MCMC chain."))
print(paste("We will use", vb_iters, "VI iterations."))
# Check input
if (nReps <= 0){stop(paste("'nReps' must be at least 1"))}
if (ncores <= 0){stop(paste("'ncores' must be at least 1"))}
if (!dir.exists(code_dir)){stop(paste("The directory", code_dir, "doesn't exist."))}
if (!dir.exists(data_dir)){stop(paste("The directory", data_dir, "doesn't exist."))}
if (!dir.exists(model_dir)){dir.create(model_dir, showWarnings = TRUE)}
# Load packages & set seed
seed.val = 123
set.seed(seed.val)
library(Matrix)
library(rhdf5)
library(topicmodels)
library(foreach)
library(doParallel)
suppressMessages(suppressWarnings(source(file.path(code_dir, "source", "stan-lda.R"))))
suppressMessages(suppressWarnings(source(file.path(code_dir, "source", "stan-mmlda.R"))))
# Generate all data file names
data_types = c('symmdata', 'spardata')
data_fns = list()
for (dd in 1:length(data_types)){
  data_fns = append(data_fns, lapply(1:nReps, function(x){return(file.path(data_dir, paste0(data_types[dd], x, ".h5")))}))
}

# FUNCTION-----------------------------------------------------------------------
save_params <- function(model_h5, theta, beta, rt, data_fn, trace_fn, K, D, I){
  # We want theta to be DxK
  # We want the rows of theta to each sum to 1
  # We want beta to be KxI
  # We want the rows of beta to each sum to 1
  # Check parameters
  #print(paste("K =", K))
  #print(paste("D =", D))
  #print(paste("I =", I))
  #print(dim(theta))
  #print(dim(beta))
  if (any(dim(theta) != c(D, K))){
    stop(paste0("The dim of theta is bad: ", dim(theta)[1], "x", dim(theta)[2]))
  }
  if (any(dim(beta) != c(K, I))){
    stop(paste0("The dim of beta is bad: ", dim(beta)[1], "x", dim(beta)[2]))
  }
  if (!all.equal(rowSums(theta), rep(1,D))){
    stop("The rows of theta should sum to 1.")
  }
  if (!all.equal(rowSums(beta), rep(1,K))){
    stop("The rows of phi should sum to 1.")
  }
  # Save parameters
  rhdf5::h5createFile(model_h5)
  rhdf5::h5write(theta, model_h5, "theta")
  rhdf5::h5write(beta, model_h5, "beta")
  rhdf5::h5write(rt, model_h5, "rt")
  rhdf5::h5write(data_fn, model_h5, "data_fn")
  rhdf5::h5write(trace_fn, model_h5, "trace_fn")
}

fit_lda <- function(data_fn, nchains, mcmc_iters, seed.val, code_dir, 
                    model_dir){
  print(paste0("PROCESS DATA FILE ", data_fn))
  data_set = tools::file_path_sans_ext(basename(data_fn))
  data_type = gsub('+[[:digit:]]+', '', data_set)
  # Read data file
  data.h5f = rhdf5::H5Fopen(data_fn)
  K = as.integer(data.h5f$params$K)
  D = as.integer(data.h5f$params$D)
  I = as.integer(data.h5f$params$I)
  # 3. Run RStan
  model_h5 = file.path(model_dir, paste0('lda_', data_set, '_rstan_gibbs.h5'))
  model_rds = file.path(model_dir, paste0('lda_', data_set, '_rstan_gibbs.rds'))
  if (!file.exists(model_h5)){
    if (!file.exists(model_rds)){
      # Run LDA
      mod = lda.stan.gibbs.sampler(K=K, iters=mcmc_iters,
                                     W=Matrix::Matrix(data.h5f$data$R, sparse=TRUE),
                                     alp=data.h5f$params$alp,
                                     mu=data.h5f$params$sig,
                                     ncores=nchains, nchains=nchains, 
                                     seed=seed.val, codedir=code_dir)
      saveRDS(mod, file=model_rds)
    } else {
      # Read RDS file
      mod = readRDS(model_rds)
    }
    # Here, run time is the average time (seconds) for each chain
    rt = sum(colMeans(rstan::get_elapsed_time(mod$fit)))
    # Extract parameters
    est.theta = calculate_parameter(resfit=mod$fit, param.name='theta')
    est.theta = matrix(est.theta, nrow = D, ncol = K)
    est.beta = calculate_parameter(resfit=mod$fit, param.name='phi')
    est.beta = matrix(est.beta, nrow = K, ncol = I)
    # Save parameters
    save_params(model_h5=model_h5, theta=est.theta, beta=est.beta, rt=rt, 
                data_fn=data_fn, trace_fn=model_rds, K=K, D=D, I=I)
  }
}


# RUN LDA ----------------------------------------------------------------------
# Run all of the implementations of LDA

# B. Nested parallel loop
cl <- parallel::makeCluster(ncores, outfile=out_fn)
parallel::clusterExport(cl, varlist=ls(envir=.GlobalEnv), envir = .GlobalEnv)
out <- parallel::parLapply(cl, 1:length(data_fns),
                           function(x){fit_lda(data_fn=data_fns[[x]],
                                               nchains=nchains,
                                               mcmc_iters=mcmc_iters,
                                               seed.val=seed.val,
                                               code_dir=code_dir,
                                               model_dir=model_dir)})
# Close all HDF5 files
rhdf5::h5closeAll()
# Stop the cluster
parallel::stopCluster(cl)


# A. Simple parallel loop
# Set up cluster
cl <- parallel::makeCluster(ncores, outfile="")
mcoptions <- list(set.seed=seed.val)
doParallel::registerDoParallel(cl)
out = foreach (dd = 1:length(data_fns), .combine='c', .options.multicore=mcoptions,
               .verbose=TRUE, .export=ls(envir=.GlobalEnv)) %dopar% {
                 #for (dd in 1:length(data_fns)){
                 print(paste0("PROCESS DATA FILE #", dd, " OUT OF ", length(data_fns)))
                 # Configure file names
                 data_fn = data_fns[[dd]]
                 data_set = tools::file_path_sans_ext(basename(data_fn))
                 data_type = gsub('+[[:digit:]]+', '', data_set)
                 # Read data file
                 data.h5f = rhdf5::H5Fopen(data_fn)
                 K = as.integer(data.h5f$params$K)
                 D = as.integer(data.h5f$params$D)
                 I = as.integer(data.h5f$params$J)
                 # 1. Run topicmodel's implementation of LDA
                 model_h5 = file.path(model_dir, paste0('lda_', data_set, '_topicmodels_VEM.h5'))
                 if (!file.exists(model_h5)){
                   print(paste("Estimate the posterior: ", model_h5))
                   # Run LDA with VEM
                   start_time = proc.time()
                   Xmod <- topicmodels::LDA(data.h5f$data$R, k=K, method="VEM",
                                            alpha=as.integer(data.h5f$params$alp[1]),
                                            estimate.beta=TRUE)
                   rt = proc.time() - start_time
                   rt = rt[3]
                   # Transform parameters
                   est.theta = as.matrix(Xmod@gamma)
                   est.beta = exp(as.matrix(Xmod@beta))
                   # Save parameters
                   save_params(model_h5=model_h5, theta=est.theta, beta=est.beta, rt=rt,
                               data_fn=data_fn, trace_fn="", K=K, D=D, I=I)
                 }
                 # 2. Run RStan's implementation of LDA with mean-field VI
                 model_h5 = file.path(model_dir, paste0('lda_', data_set, '_rstan_mfVI.h5'))
                 model_rds = file.path(model_dir, paste0('lda_', data_set, '_rstan_mfVI.rds'))
                 if (!file.exists(model_h5)){
                   if (!file.exists(model_rds)){
                     good.mod = FALSE
                     print(paste("Run RStan mean-field VI for", model_rds))
                     # Try statement in case model can't be fitted, save outcome as NA
                     try({
                       mod <- lda.stan.variational.bayes(K=K,
                                                         W=Matrix::Matrix(data.h5f$data$R, sparse=TRUE),
                                                         alg="meanfield", iters=vb_iters,
                                                         mu=data.h5f$params$mu,
                                                         alp=data.h5f$params$alp,
                                                         codedir=code_dir, seed=seed.val)
                       saveRDS(mod, model_rds)
                       good.mod = TRUE
                     })
                     if (!good.mod){
                       mod = NULL
                       saveRDS(mod, model_rds)
                     }
                   } else {
                     # Load the model
                     mod = readRDS(model_rds)
                   }
                   if (!is.null(mod)){
                     # Extract run time
                     rt = mod$time$format + mod$time$compile
                     # Extract parameters
                     est.theta = calculate_parameter(resfit=mod$fit, param.name='theta',
                                                     param.est="lastmcmc")
                     est.theta = matrix(est.theta, nrow = D, ncol = K)
                     est.beta = calculate_parameter(resfit=mod$fit, param.name='phi',
                                                    param.est="lastmcmc")
                     est.beta = matrix(est.beta, nrow = K, ncol = I)
                     # Normalize row sums
                     for (dd in 1:D){
                       est.theta[dd,] = est.theta[dd,]/sum(est.theta[dd,])
                     }
                     for (kk in 1:K){
                       est.beta[kk,] = est.beta[kk,]/sum(est.beta[kk,])
                     }
                     # Save parameters
                     save_params(model_h5=model_h5, theta=est.theta, beta=est.beta, rt=rt,
                                 data_fn=data_fn, trace_fn=model_rds, K=K, D=D, I=I)
                   }
                 }
                 # 3. Run RStan's implementation of LDA with full-rank VI
                 model_h5 = file.path(model_dir, paste0('lda_', data_set, '_rstan_frVI.h5'))
                 model_rds = file.path(model_dir, paste0('lda_', data_set, '_rstan_frVI.rds'))
                 if (!file.exists(model_h5)){
                   if (!file.exists(model_rds)){
                     good.mod = FALSE
                     print(paste("Run RStan full-rank VI for", model_rds))
                     # Try statement in case model can't be fitted, save outcome as NA
                     try({
                       mod <- lda.stan.variational.bayes(K=K, iters=vb_iters,
                                                         W=Matrix::Matrix(data.h5f$data$R, sparse=TRUE),
                                                         alg="fullrank", mu=data.h5f$params$mu,
                                                         alp=data.h5f$params$alp,
                                                         codedir=code_dir, seed=seed.val)
                       # Save results
                       saveRDS(mod, model_rds)
                       good.mod = TRUE
                     })
                     if (!good.mod){
                       mod = NULL
                       saveRDS(mod, model_rds)
                     }
                   } else {
                     # Load the model
                     mod = readRDS(model_rds)
                   }
                   if (!is.null(mod)){
                     # Extract run time
                     rt = mod$time$format + mod$time$compile
                     # Extract parameters
                     est.theta = calculate_parameter(resfit=mod$fit, param.name='theta',
                                                     param.est="lastmcmc")
                     est.theta = matrix(est.theta, nrow = D, ncol = K)
                     est.beta = calculate_parameter(resfit=mod$fit, param.name='phi',
                                                    param.est="lastmcmc")
                     est.beta = matrix(est.beta, nrow = K, ncol = I)
                     # Normalize row sums
                     for (dd in 1:D){
                       est.theta[dd,] = est.theta[dd,]/sum(est.theta[dd,])
                     }
                     for (kk in 1:K){
                       est.beta[kk,] = est.beta[kk,]/sum(est.beta[kk,])
                     }
                     # Save parameters
                     save_params(model_h5=model_h5, theta=est.theta, beta=est.beta, rt=rt,
                                 data_fn=data_fn, trace_fn=model_rds, K=K, D=D, I=I)
                   }
                 }
                 # Close the data file
                 rhdf5::H5Fclose(data.h5f)
               }
# Stop the cluster
parallel::stopCluster(cl)
# Close all HDF5 files
rhdf5::h5closeAll()

print("COMPLETED LDA TESTING SUCCESSFULLY.")
