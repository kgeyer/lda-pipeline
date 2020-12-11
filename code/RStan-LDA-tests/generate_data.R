#!/usr/bin/env Rscript


# HEADER------------------------------------------------------------------------
# Author: Kelly Geyer, klgeyer@bu.edu
# Date: September 24, 2020
# Development: R v4.0.2 on MacOS v10.15.4
# Description: Generate data sets for the LDA and MM-LDA models. In, all there 
# are three methods of simulating data:
# (i) Assume symmetric hyperparameters. That is, all alpha_k's, mu_j's,and 
#     sigma_i's are equal to 1. This implies that the topics occur equally among 
#     documents and features.
# (ii) Assume sparse topics. That is, all alpha_k's, mu_j's,and 
#     sigma_i's are equal to 0.01. This implies that a single topic will have 
#     very a high probability among documents and features while the others are 
#     low.
# (iii) We can simulate signatures with varying sparsity levels. We sample each 
#     topic i ~ Dir(alpha_i), where alpha_i ~ Gamma(0.5, 2).


# SET UP SCRIPT-----------------------------------------------------------------
# Read parameters
print("SIMULATE TEST DATA")
args <- commandArgs()
code_dir = args[6]
data_dir = args[7]
nReps = strtoi(args[8])
ncores = strtoi(args[9])
K = strtoi(args[10])
D = strtoi(args[11])
I = strtoi(args[12])
J = strtoi(args[13])
lam = strtoi(args[14])
print(paste("The code directory is", code_dir))
print(paste("The data directory is", data_dir))
print(paste("We will generate", nReps, "dataset replications."))
print(paste("We will use", ncores, "cores."))
print(paste("We will set K =", K, "."))
print(paste("We will set D =", D, "."))
print(paste("We will set I =", I, "."))
print(paste("We will set J =", J, "."))
print(paste("We will set lam =", lam, "."))
# Check input
if (nReps <= 0){stop(paste("'nReps' must be at least 1"))}
if (ncores <= 0){stop(paste("'ncores' must be at least 1"))}
if (!dir.exists(code_dir)){stop(paste("The directory", code_dir, "doesn't exist."))}
if (!dir.exists(data_dir)){dir.create(data_dir, showWarnings = TRUE)}
# Load dependencies
set.seed(123)
source(file.path(code_dir, "source", "simulate_data.R"))
library(rhdf5)
library(foreach)
library(doParallel)
# Set up cluster
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

# DEFINE DATA PARAMETERS--------------------------------------------------------
# Recommended parameters (from JC on 5.14.20 and 5.22.20)
# K = 5
# 200-1000 documents
# 100 vocab size
# 100 distinct image regions
# 500 average number of words/regions per document
# K = 10
# D = 500
# J = 100
# I = 100?
# lam = 500

# mini-test
K = 3
D = 50
J = 20
I = 20
lam = 100

# GENERATE SYMMETRIC DATA-------------------------------------------------------
# Generate several instances of symmetric data sets. That is, each topic is
# equally likely among documents and features.
out = foreach (rr = 1:nReps) %dopar% {
  data_fn = file.path(data_dir, paste0("symmdata", toString(rr), ".h5"))
  if (!file.exists(data_fn)){
    print(paste("Creating the dataset:", data_fn))
    alp = matrix(1, nrow = K, ncol = 1)
    sig = matrix(1, nrow = I, ncol = 1)
    mu = matrix(1, nrow = J, ncol = 1)
    data = generate.mmlda.data(K=K, D=D, J=J, I=I, lam=lam, alp=alp, sig=sig, mu=mu)
    rhdf5::h5createFile(data_fn)
    rhdf5::h5createGroup(data_fn,"data")
    rhdf5::h5createGroup(data_fn, "params")
    rhdf5::h5write(as.matrix(data$W), data_fn, "data/W")
    rhdf5::h5write(as.matrix(data$R), data_fn, "data/R")
    rhdf5::h5write(data$theta, data_fn, "params/theta")
    rhdf5::h5write(data$beta, data_fn, "params/beta")
    rhdf5::h5write(data$phi, data_fn, "params/phi")
    rhdf5::h5write(K, data_fn, "params/K")
    rhdf5::h5write(D, data_fn, "params/D")
    rhdf5::h5write(J, data_fn, "params/J")
    rhdf5::h5write(I, data_fn, "params/I")
    rhdf5::h5write(lam, data_fn, "params/lam")
    rhdf5::h5write(alp, data_fn, "params/alp")
    rhdf5::h5write(sig, data_fn, "params/sig")
    rhdf5::h5write(mu, data_fn, "params/mu")
    rhdf5::h5closeAll()
  }
}

# GENERATE SPARSE DATA----------------------------------------------------------
# Generate several instances of sparse data sets by setting low Dirichlet
# parameters. That is, a single topic will  have very high probability among
# documents and features while the others are low.
out = foreach (rr = 1:nReps) %dopar% {
  data_fn = file.path(data_dir, paste0("spardata", toString(rr), ".h5"))
  if (!file.exists(data_fn)){
    print(paste("Creating the dataset:", data_fn))
    alp = matrix(0.01, nrow = K, ncol = 1)
    sig = matrix(0.01, nrow = I, ncol = 1)
    mu = matrix(0.01, nrow = J, ncol = 1)
    data = generate.mmlda.data(K=K, D=D, J=J, I=I, lam=lam, alp=alp, sig=sig, mu=mu)
    rhdf5::h5createFile(data_fn)
    rhdf5::h5createGroup(data_fn,"data")
    rhdf5::h5createGroup(data_fn, "params")
    rhdf5::h5write(as.matrix(data$W), data_fn, "data/W")
    rhdf5::h5write(as.matrix(data$R), data_fn, "data/R")
    rhdf5::h5write(data$theta, data_fn, "params/theta")
    rhdf5::h5write(data$beta, data_fn, "params/beta")
    rhdf5::h5write(data$phi, data_fn, "params/phi")
    rhdf5::h5write(K, data_fn, "params/K")
    rhdf5::h5write(D, data_fn, "params/D")
    rhdf5::h5write(J, data_fn, "params/J")
    rhdf5::h5write(I, data_fn, "params/I")
    rhdf5::h5write(lam, data_fn, "params/lam")
    rhdf5::h5write(alp, data_fn, "params/alp")
    rhdf5::h5write(sig, data_fn, "params/sig")
    rhdf5::h5write(mu, data_fn, "params/mu")
    rhdf5::h5closeAll()
  }
}

# Stop the cluster
parallel::stopCluster(cl)
# Close all HDF5 files
rhdf5::h5closeAll()

print("COMPLETED GENERATING DATA SUCCESSFULLY.")


