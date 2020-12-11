#!/bin/bash

# Set directories
datadir=/Users/kelly/GoogleDriveBU/MM-LDA/lda-testing/simulation-tests/sim_12112020/data
codedir=/Users/kelly/Documents/github/lda-pipeline/code
modeldir=/Users/kelly/GoogleDriveBU/MM-LDA/lda-testing/simulation-tests/sim_12112020/models
resultdir=/Users/kelly/GoogleDriveBU/MM-LDA/lda-testing/simulation-tests/sim_12112020/results

# Define the number of relications or simulations for test
nreps=1
# Define number of cores
ncores=2
# Data generation parameters
K=3
D=20
I=10
J=10
lam=50
# R test parameters
mcmciters=500
viiters=1000
nchain=4

# Set environment variables
export HDF5_USE_FILE_LOCKING=FALSE

# Generate data
Rscript generate_data.R $codedir $datadir $nreps $ncores $K $D $I $J $lam

# Run models in R
outfn=$codedir/RStan-LDA-tests/lda_test_output.out
Rscript lda_test.R $codedir $datadir $modeldir $nreps $ncores $nchain $viiters $mcmciters $outfn

# Evaluate the tests
Rscript evaluate_lda_results.R $codedir $datadir $modeldir $resultdir $nreps
