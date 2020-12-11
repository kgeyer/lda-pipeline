# HEADER------------------------------------------------------------------------
# Author: Kelly Geyer, klgeyer@bu.edu
# Date: October 10, 2020
# Development: R v4.0.2 on MacOS v10.15.6
#
# Description: This script contains functions for evaluating cosine similarity 
# between true and estimated parameters.

# LIBRARIES---------------------------------------------------------------------
library(ggplot2)
library(ggpubr)

# FUNCTIONS---------------------------------------------------------------------
make_topic_cs_barplot <- function(true.param, est.param, plt_title, param.name='theta'){
  # Make a single topic CS plot for a parameter
  trans.param = TRUE
  if (param.name == 'theta'){trans.param = FALSE}
  if (trans.param){
    true.param = t(true.param)
    est.param = t(est.param)
  }
  # Check params
  if (any(dim(true.param) != dim(est.param))){
    stop(paste0("Dimension missmatch! true is ", dim(true.param)[1], "x", dim(true.param)[2], " and ", dim(est.param)[1], "x", dim(est.param)[2]))
  }
  # get dims
  NROW = dim(true.param)[1]  # feats
  K = dim(true.param)[2]  # topics
  # Calc the cosine similarity
  cs = matrix(0, nrow=K, ncol=1)
  for (kk in 1:K){
    cs[kk] = calculate_cosine_similarity(est.param[,kk], true.param[,kk])
  }
  # Make plot
  cs.df = data.frame(cs=cs, topic=factor(1:K))
  plt <- ggplot2::ggplot(cs.df, aes(x=topic, y=cs, fill=topic)) + 
    ggplot2::geom_bar(stat="identity", color="black", position = position_dodge()) +
    ggplot2::xlab("Topic") +  ggplot2::coord_cartesian(ylim=c(0,1)) +
    ggplot2::ggtitle(plt_title) + ggplot2::theme(legend.position = "none") 
  return(plt)
}

make_feat_cs_barplot <- function(true.param, est.param, plt_title, param.name='theta'){
  # Make a single topic CS plot for a parameter
  trans.param = TRUE
  if (param.name == 'theta'){trans.param = FALSE}
  if (trans.param){
    true.param <- t(true.param)
    est.param <- t(est.param)
  }
  # Check params
  if (any(dim(true.param) != dim(est.param))){
    stop(paste0("Dimension missmatch! true is ", dim(true.param)[1], "x", dim(true.param)[2], " and ", dim(est.param)[1], "x", dim(est.param)[2]))
  }
  # get dims
  N = dim(true.param)[1]  # feats
  K = dim(true.param)[2]  # topics
  # Calc the cosine similarity
  cs = matrix(0, nrow=N, ncol=1)
  for (kk in 1:N){
    cs[kk] = calculate_cosine_similarity(est.param[kk,], true.param[kk,])
  }
  # Make plot
  cs.df = data.frame(cs=cs, feat=factor(1:N))
  plt <- ggplot2::ggplot(cs.df, aes(x=feat, y=cs, fill=feat)) + 
    ggplot2::geom_bar(stat="identity", color="black", position = position_dodge()) +
    ggplot2::xlab("Feature") + ggplot2::coord_cartesian(ylim=c(0,1)) +
    ggplot2::ggtitle(plt_title) + ggplot2::theme(legend.position = "none") 
  return(plt)
}

calculate_cosine_similarity <- function(vec1, vec2){
  # This function calculates the cosine similarity of two vectors
  # 
  # Inputs:
  #   vec1, vec2: Vectors/lists/matrices of equal size
  vec1 = as.vector(vec1)
  vec2 = as.vector(vec2)
  if (length(vec1) != length(vec2)){
    stop("'vec1' and 'vec2' must have the same size.")
  }
  sim.score = 0
  norm1 = norm(vec1, type = "2")
  norm2 = norm(vec2, type = "2")
  denom = norm1*norm2
  if (denom > 0){sim.score = sum(vec1*vec2)/denom}
  return(sim.score)
}

