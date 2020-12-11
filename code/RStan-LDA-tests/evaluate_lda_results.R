#!/usr/bin/env Rscript

# HEADER------------------------------------------------------------------------
# Author: Kelly Geyer, klgeyer@bu.edu
# Date: September 24, 2020
# Development: R v4.0.2 on MacOS v10.15.6
# Description: This script evaluates the results of the tests for the pyMC3 
# models.

# SET UP SCRIPT-----------------------------------------------------------------
print("EVALUATE RESULTS")
args <- commandArgs()
code_dir = args[6]
data_dir = args[7]
model_dir = args[8]
results_dir = args[9]
nReps = strtoi(args[10])

print(paste("The code directory is", code_dir))
print(paste("The data directory is", data_dir))
print(paste("The model directory is", model_dir))
print(paste("The results direcotry is", results_dir))
print(paste("We will have", nReps, "replications."))
# Check input
if (nReps <= 0){stop(paste("'nReps' must be at least 1"))}
if (!dir.exists(code_dir)){stop(paste("The directory", code_dir, "doesn't exist."))}
if (!dir.exists(data_dir)){stop(paste("The directory", data_dir, "doesn't exist."))}
if (!dir.exists(model_dir)){stop(paste("The directory", model_dir, "doesn't exist."))}
if (!dir.exists(results_dir)){dir.create(results_dir, showWarnings = TRUE)}
#
set.seed(123)
source(file.path(code_dir, "source", "bagel_plts.R"))
source(file.path(code_dir, "source", "stan-mmlda.R"))
source(file.path(code_dir, "source", "stan-lda.R"))
library(rhdf5)
# Generate all model file names
nReps = c(1,2,3)
lda_est_types = c('pyMC3_mfADVI', 'pyMC3_frADVI', 'pyMC3_SVGD', 'pyMC3_NUTS', 'topicmodels_VEM', 'rstan_gibbs', 'rstan_mfVI', 'rstan_frVI')
model_fns = list()
for (rr in 1:length(nReps)){
  model_fns = append(model_fns, lapply(lda_est_types, function(x){return(file.path(model_dir, paste0("lda_symmdata",nReps[rr],"_",x,'.h5')))}))
  model_fns = append(model_fns, lapply(lda_est_types, function(x){return(file.path(model_dir, paste0("lda_spardata",nReps[rr],"_",x,'.h5')))}))
}
# Close all HDF5 files
rhdf5::h5closeAll()


# FUNCTION-----------------------------------------------------------------------
check_params <- function(theta, beta, K, D, I){
  # We want theta to be DxK
  # We want the rows of theta to each sum to 1
  # We want beta to be KxI
  # We want the rows of beta to each sum to 1
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
}

# LDA RESULTS-------------------------------------------------------------------
# Data frame of statistics
lda.summary <- data.frame(ds=character(0), dt=character(0), mod=character(0), 
                          feat=character(0),rt=numeric(0), theta.cs=numeric(0), 
                          beta.cs=numeric(0))
colnames(lda.summary) <- c("Data Set", "Data Type", "Model Type", "Feature", "Run Time", 
                           "Theta CS", "Beta CS")
# Extract results for each model file
for (fn in 1:length(model_fns)){
  model_fn = model_fns[[fn]]
  x = strsplit(tools::file_path_sans_ext(basename(model_fn)), '_')
  data_set = x[[1]][2]
  data_type = gsub('+[[:digit:]]+', '', data_set)
  model_type = paste(x[[1]][3], x[[1]][4], sep="_")
  package_name = x[[1]][3]
  data_fn = file.path(data_dir, paste0(data_set, '.h5'))
  if (file.exists(model_fn)){
    print(paste0("Generating results for ", model_fn, "..."))
    # Read data file
    data.h5f = rhdf5::H5Fopen(data_fn)
    true.theta = data.h5f$params$theta
    true.beta = data.h5f$params$beta
    K = as.integer(data.h5f$params$K)
    D = as.integer(data.h5f$params$D)
    I = as.integer(data.h5f$params$I)
    # Read model file
    model.h5f = rhdf5::H5Fopen(model_fn)
    if (x[[1]][4] == "NUTS"){
      run.time = abs(model.h5f$rt)/4
    } else {
      run.time = abs(model.h5f$rt)
    }
    est.theta = model.h5f$theta
    est.beta = model.h5f$beta
    h5closeAll()
    check_params(est.theta, est.beta, K, D, I)
    # Match topics
    new.idx = lda_greedy_topic_permute(true.theta, true.beta, est.theta, 
                                       est.beta, code_dir)
    est.theta = est.theta[,new.idx]
    est.beta = est.beta[new.idx,]
    # Store summary statistics
    # # Find generalized CS
    # all.cs <- data.frame(ds=data_set, dt=data_type, mod=model_type, feat="all",
    #                      rt=run.time,
    #                      theta.cs=calculate_cosine_similarity(est.theta, true.theta),
    #                      beta.cs=calculate_cosine_similarity(est.beta, true.beta))
    # # Find CS by topic
    # topic.cs <- data.frame(ds=rep(data_set,K), dt=rep(data_type,K),
    #                        mod=rep(model_type,K), feat=rep("topic",K),
    #                        rt=rep(NA,K),
    #                        theta.cs=sapply(1:K, function(ii){calculate_cosine_similarity(est.theta[,ii], true.theta[,ii])}),
    #                        beta.cs=sapply(1:K, function(ii){calculate_cosine_similarity(est.beta[ii,], true.beta[ii,])}))
    # # Find CS by feature
    # theta.feat.cs <- data.frame(ds=rep(data_set,D), dt=rep(data_type,D),
    #                             mod=rep(model_type,D), feat=rep("feat",D),
    #                             rt=rep(NA,D), beta.cs=rep(NA,D),
    #                             theta.cs=sapply(1:D, function(ii){calculate_cosine_similarity(est.theta[ii,], true.theta[ii,])}))
    # beta.feat.cs <- data.frame(ds=rep(data_set,I), dt=rep(data_type,I),
    #                            mod=rep(model_type,I), feat=rep("feat",I),
    #                            rt=rep(NA,I), theta.cs=rep(NA,I),
    #                            beta.cs=sapply(1:I, function(ii){calculate_cosine_similarity(est.beta[,ii], true.beta[,ii])}))
    # lda.summary <- rbind(lda.summary, all.cs, topic.cs, theta.feat.cs, beta.feat.cs)
    # Create BAGEL plot
    true_bagel_fn = file.path(results_dir, paste0('lda_', data_set, '_true_bagel.pdf'))
    est_bagel_fn = file.path(results_dir, paste0('lda_', data_set, '_', model_type, '_bagel.pdf'))
    data_lab = "Sparse"
    if (data_type == "symmdata"){data_lab = "Dense"}
    plt_title = paste0("LDA Estimation: ", model_type, "\n",data_lab," Topic Data: K = ", K, ", D = ", D, ", I = ", I)
    create_test_lda_bagel_plt(true.theta=true.theta, true.beta=true.beta, 
                              est.theta=est.theta, est.beta=est.beta, 
                              est_bagel_fn=est_bagel_fn, 
                              true_bagel_fn=true_bagel_fn, plt_title=plt_title)  # write this!
    # Create CS plot
    cs_plt_fn = file.path(results_dir, paste0('lda_', data_set, '_', model_type, '_cs.pdf'))
    make_lda_cosine_plots(true.theta=true.theta, true.beta=true.beta, 
                          est.theta=est.theta, est.beta=true.beta, 
                          plt_fn=cs_plt_fn, plt_title=plt_title, 
                          codedir=code_dir)   
    base_fn = tools::file_path_sans_ext(basename(model_fn))
    fn_rds = file.path(model_dir, paste0(base_fn, '.rds'))
    if (file.exists(fn_rds)){
      mod = readRDS(fn_rds)
      # Create RStan Geweke plot
      geweke_plt_fn = file.path(results_dir, paste0(base_fn, '_geweke.pdf'))
      make_geweke_plot(resfit=mod$fit, plt_fn=geweke_plt_fn, plt_title=plt_title, height=80)
      # Create convergence plots
      conv_plt_fn = file.path(results_dir, paste0(base_fn, '_convstat.pdf'))
      lda_make_conv_stat_plot(resfit=mod$fit, plt_fn=conv_plt_fn, plt_title=plt_title)
    }
  }
}
# Close all HDF5 files
rhdf5::h5closeAll()

# SUMMARY RESULTS---------------------------------------------------------------
print("CREATE LDA TABLES")
summary_rds = file.path(results_dir, "lda-summary.rds")
lda.summary[which(!is.finite(lda.summary))] <- NA
saveRDS(lda.summary, summary_rds)
# TABLE 1: means
summary_fn = file.path(results_dir, "lda-summary-mean.csv")
X.all = lda.summary[lda.summary$feat=="all", c(2,3,5,6,7)]
agg.all.mean = aggregate(X.all, by=list(X.all$dt, X.all$mod), FUN=function(x){round(mean(x, na.rm=TRUE),2)})
agg.all.mean = agg.all.mean[-c(3,4)]
X.topic = lda.summary[lda.summary$feat=="topic", c(2,3,6,7)]
agg.topic.mean = aggregate(X.topic, by=list(X.topic$dt, X.topic$mod), FUN=function(x){round(mean(x, na.rm=TRUE),2)})
agg.topic.mean = agg.topic.mean[-c(3,4)]
X.feat = lda.summary[lda.summary$feat=="feat", c(2,3,6,7)]
agg.feat.mean = aggregate(X.feat, by=list(X.feat$dt, X.feat$mod), FUN=function(x){round(mean(x, na.rm=TRUE),2)})
agg.feat.mean  = agg.feat.mean[-c(3,4)]
# merge data
X = data.frame(dt=agg.all.mean$Group.1, mod=agg.all.mean$Group.2, 
               rt=agg.all.mean$rt, 
               theta.cs=agg.all.mean$theta.cs, 
               beta.cs=agg.all.mean$beta.cs, 
               theta.topic.cs=agg.topic.mean$theta.cs, 
               beta.topic.cs=agg.topic.mean$beta.cs,
               theta.feat.cs=agg.feat.mean$theta.cs,
               beta.feat.cs=agg.feat.mean$beta.cs)
colnames(X) <- c("Data Type", "Estimation", "Run Time (seconds)", "Theta CS", 
                 "Beta CS", paste0("Topic Theta CS: K=",K), 
                 paste0("Topic Beta CS: K=",K), 
                 paste0("Doc Theta CS: D=",D),
                 paste0("Feat Beta CS: I=", I))
write.csv(X, summary_fn, row.names=FALSE)

# TABLE 2: sd
summary_fn = file.path(results_dir, "lda-summary-sd.csv")
X.all = lda.summary[lda.summary$feat=="all", c(2,3,5,6,7)]
agg.all.sd = aggregate(X.all, by=list(X.all$dt, X.all$mod), FUN=function(x){round(sd(x, na.rm=TRUE),2)})
agg.all.sd = agg.all.sd[-c(3,4)]
X.topic = lda.summary[lda.summary$feat=="topic", c(2,3,6,7)]
agg.topic.sd = aggregate(X.topic, by=list(X.topic$dt, X.topic$mod), FUN=function(x){round(sd(x, na.rm=TRUE),2)})
agg.topic.sd = agg.topic.sd[-c(3,4)]
X.feat = lda.summary[lda.summary$feat=="feat", c(2,3,6,7)]
agg.feat.sd = aggregate(X.feat, by=list(X.feat$dt, X.feat$mod), FUN=function(x){round(sd(x, na.rm=TRUE),2)})
agg.feat.sd  = agg.feat.sd[-c(3,4)]
# merge data
X = data.frame(dt=agg.all.sd$Group.1, mod=agg.all.sd$Group.2, 
               rt=agg.all.sd$rt, 
               theta.cs=agg.all.sd$theta.cs, 
               beta.cs=agg.all.sd$beta.cs, 
               theta.topic.cs=agg.topic.sd$theta.cs, 
               beta.topic.cs=agg.topic.sd$beta.cs,
               theta.feat.cs=agg.feat.sd$theta.cs,
               beta.feat.cs=agg.feat.sd$beta.cs)
colnames(X) <- c("Data Type", "Estimation", "Run Time (seconds)", "Theta CS", 
                 "Beta CS", paste0("Topic Theta CS: K=",K), 
                 paste0("Topic Beta CS: K=",K), 
                 paste0("Doc Theta CS: D=",D),
                 paste0("Feat Beta CS: I=", I))
write.csv(X, summary_fn, row.names=FALSE)

