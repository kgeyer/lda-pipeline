#-------------------------------------------------------------------------------
# HEADER
#-------------------------------------------------------------------------------
# Author: Kelly Geyer, klgeyer@bu.edu
# Date: May 21, 2020
# Development: R v4.0.0 on MacOS v10.15.4
#
# Description: Create BAGEL plots

#-------------------------------------------------------------------------------
# LIBRARIES
library(BAGEL)
library(ggplot2)
library(ggpubr)

#-------------------------------------------------------------------------------
# FUNCTIONS
create_test_smmlda_bagel_plt_fixb <- function(true.theta, true.beta, true.phi, 
                                         est.theta, est.beta, est.phi, 
                                         true_bagel_fn, est_bagel_fn, bagel_fn,
                                         plt_title){
  # This function creates BAGEL plots for MM-LDA parameters
  #
  # Args:
  #   theta: (D x K)
  #   beta: (K x I)
  #   phi: (K x J)
  # Create true parameter plots
  D = dim(true.theta)[1]
  K = dim(true.theta)[2]
  I = dim(true.beta)[2]
  J = dim(true.phi)[2]
  true.beta <- t(true.beta)
  colnames(true.beta) <- paste("Signature", 1:K, sep = "")
  rownames(true.beta) <- lapply(1:I, toString)
  true.theta <- t(true.theta)
  colnames(true.theta) <- lapply(1:D, toString)
  rownames(true.theta) <- paste("Signature", 1:K, sep = "")
  true.phi <- t(true.phi)
  colnames(true.phi) <- paste("Signature", 1:K, sep = "")
  rownames(true.phi) <- lapply(1:J, toString)
  # Convert matrices to BAGEL objects
  setClass("Result", representation(exposures="matrix", signatures="matrix"))
  true.bag <- methods::new("Result", exposures=true.theta, signatures=true.beta)
  p1 <- BAGEL::plot_exposures(true.bag, proportional = TRUE, 
                              label_samples = FALSE, no_legend = TRUE)
  p1 <- p1 + ggplot2::ggtitle('Theta: Signature-Sample Proportions')
  p2 <- BAGEL::plot_signatures(true.bag, no_legend = TRUE, text_size = 10, 
                               facet_size = 10)
  p2 <- p2 + ggtitle("Phi1: Modality #1 Signature-Mutation Proportion")  
  true.bag <- methods::new("Result", exposures=true.theta, signatures=true.phi)
  p3 <- BAGEL::plot_signatures(true.bag, no_legend = TRUE, text_size = 10, 
                               facet_size = 10)
  p3 <- p3 + ggplot2::ggtitle("Phi2: Modality #2 Signature-Mutation Proportions")
  true.plt <- ggpubr::ggarrange(p1, p2, p3, nrow=3, ncol=1, common.legend = FALSE)
  true.plt <- ggpubr::annotate_figure(true.plt, top="True Parameters", fig.lab.face="bold")
  ggplot2::ggsave(true_bagel_fn, true.plt, width=10, height=15, units="in")
  # Create estimated parameter plots
  est.beta <- t(est.beta)
  colnames(est.beta) <- paste("Signature", 1:K, sep = "")
  rownames(est.beta) <- lapply(1:I, toString)
  est.theta <- t(est.theta)
  colnames(est.theta) <- lapply(1:D, toString)
  rownames(est.theta) <- paste("Signature", 1:K, sep = "")
  est.phi <- t(est.phi)
  colnames(est.phi) <- paste("Signature", 1:K, sep = "")
  rownames(est.phi) <- lapply(1:J, toString)
  # Convert matrices to BAGEL objects
  est.bag <- methods::new("Result", exposures=est.theta, signatures=est.beta)
  p4 <- BAGEL::plot_exposures(est.bag, proportional = TRUE, 
                              label_samples = FALSE, no_legend = TRUE)
  p4 <- p4 + ggplot2::ggtitle('Theta: Signature-Sample Proportions')
  p5 <- BAGEL::plot_signatures(est.bag, no_legend = TRUE, text_size = 10, 
                               facet_size = 10)
  p5 <- p5 + ggplot2::ggtitle("Phi1: Modaltiy #1 Signature-Mutation Proportion")  
  est.bag <- methods::new("Result", exposures=est.theta, signatures=est.phi)
  p6 <- BAGEL::plot_signatures(est.bag, no_legend = TRUE, text_size = 10, 
                               facet_size = 10)
  p6 <- p6 + ggplot2::ggtitle("Phi2: Modality #2 Signature-Mutation Proportions") 
  est.plt <- ggpubr::ggarrange(p4, p5, p6, nrow=3, ncol=1, common.legend = FALSE)
  est.plt <- ggpubr::annotate_figure(est.plt, top="Estimated Parameters", fig.lab.face="bold")
  ggplot2::ggsave(est_bagel_fn, est.plt, width=10, height=15, units="in")
  # Combine plots
  plt <- ggpubr::ggarrange(true.plt, est.plt, nrow=1, ncol=2, common.legend = TRUE)
  plt <- ggpubr::annotate_figure(plt, top=plt_title, fig.lab.face = "bold")
  ggplot2::ggsave(bagel_fn, plt, width=15, height=35, units="in")
}


create_test_smmlda_bagel_plt <- function(true.theta, true.beta, true.phi, 
                                         true.b1,
                                         est.theta, est.beta, est.phi, est.b1,
                                         true_bagel_fn, est_bagel_fn, bagel_fn,
                                         plt_title){
  # This function creates BAGEL plots for MM-LDA parameters
  #
  # Args:
  #   theta: (D x K)
  #   beta: (K x I)
  #   phi: (K x J)
  # Create true parameter plots
  D = dim(true.theta)[1]
  K = dim(true.theta)[2]
  I = dim(true.beta)[2]
  J = dim(true.phi)[2]
  true.beta <- t(true.beta)
  colnames(true.beta) <- paste("Signature", 1:K, sep = "")
  rownames(true.beta) <- lapply(1:I, toString)
  true.theta <- t(true.theta)
  colnames(true.theta) <- lapply(1:D, toString)
  rownames(true.theta) <- paste("Signature", 1:K, sep = "")
  true.phi <- t(true.phi)
  colnames(true.phi) <- paste("Signature", 1:K, sep = "")
  rownames(true.phi) <- lapply(1:J, toString)
  # Convert matrices to BAGEL objects
  setClass("Result", representation(exposures="matrix", signatures="matrix"))
  true.bag <- methods::new("Result", exposures=true.theta, signatures=true.beta)
  p1 <- BAGEL::plot_exposures(true.bag, proportional = TRUE, 
                              label_samples = FALSE, no_legend = TRUE)
  p1 <- p1 + ggplot2::ggtitle('Theta: Signature-Sample Proportions')
  p2 <- BAGEL::plot_signatures(true.bag, no_legend = TRUE, text_size = 10, 
                               facet_size = 10)
  p2 <- p2 + ggtitle("Phi1: Modality #1 Signature-Mutation Proportion")  
  true.bag <- methods::new("Result", exposures=true.theta, signatures=true.phi)
  p3 <- BAGEL::plot_signatures(true.bag, no_legend = TRUE, text_size = 10, 
                               facet_size = 10)
  p3 <- p3 + ggplot2::ggtitle("Phi2: Modality #2 Signature-Mutation Proportions")
  true.df <- data.frame(sig=as.factor(1:K), prob=true.b1)
  p7 <- ggplot2::ggplot(data=true.df, aes(x=sig, y=prob)) + 
    geom_bar(stat="identity") + ylim(0,1) + xlab("Signature") + ylab("Proportion")
  p7 <- p7 + ggplot2::ggtitle("b1: Signature Proportions for Modality #1")
  true.plt <- ggpubr::ggarrange(p1, p2, p3, p7, nrow=4, ncol=1, common.legend = FALSE)
  true.plt <- ggpubr::annotate_figure(true.plt, top="True Parameters", fig.lab.face="bold")
  ggplot2::ggsave(true_bagel_fn, true.plt, width=10, height=15, units="in")
  # Create estimated parameter plots
  est.beta <- t(est.beta)
  colnames(est.beta) <- paste("Signature", 1:K, sep = "")
  rownames(est.beta) <- lapply(1:I, toString)
  est.theta <- t(est.theta)
  colnames(est.theta) <- lapply(1:D, toString)
  rownames(est.theta) <- paste("Signature", 1:K, sep = "")
  est.phi <- t(est.phi)
  colnames(est.phi) <- paste("Signature", 1:K, sep = "")
  rownames(est.phi) <- lapply(1:J, toString)
  # Convert matrices to BAGEL objects
  est.bag <- methods::new("Result", exposures=est.theta, signatures=est.beta)
  p4 <- BAGEL::plot_exposures(est.bag, proportional = TRUE, 
                              label_samples = FALSE, no_legend = TRUE)
  p4 <- p4 + ggplot2::ggtitle('Theta: Signature-Sample Proportions')
  p5 <- BAGEL::plot_signatures(est.bag, no_legend = TRUE, text_size = 10, 
                               facet_size = 10)
  p5 <- p5 + ggplot2::ggtitle("Phi1: Modaltiy #1 Signature-Mutation Proportion")  
  est.bag <- methods::new("Result", exposures=est.theta, signatures=est.phi)
  p6 <- BAGEL::plot_signatures(est.bag, no_legend = TRUE, text_size = 10, 
                               facet_size = 10)
  p6 <- p6 + ggplot2::ggtitle("Phi2: Modality #2 Signature-Mutation Proportions") 
  est.df <- data.frame(sig=as.factor(1:K), prob=est.b1)
  p8 <- ggplot2::ggplot(data=est.df, aes(x=sig, y=prob)) + 
    geom_bar(stat="identity") + ylim(0,1) + xlab("Signature") + ylab("Proportion")
  p8 <- p8 + ggplot2::ggtitle("b1: Signature Proportions for Modality #1")
  est.plt <- ggpubr::ggarrange(p4, p5, p6, p8, nrow=4, ncol=1, common.legend = FALSE)
  est.plt <- ggpubr::annotate_figure(est.plt, top="Estimated Parameters", fig.lab.face="bold")
  ggplot2::ggsave(est_bagel_fn, est.plt, width=10, height=15, units="in")
  # Combine plots
  plt <- ggpubr::ggarrange(true.plt, est.plt, nrow=1, ncol=2, common.legend = TRUE)
  plt <- ggpubr::annotate_figure(plt, top=plt_title, fig.lab.face = "bold")
  ggplot2::ggsave(bagel_fn, plt, width=15, height=35, units="in")
}

create_test_mmlda_bagel_plt <- function(true.theta, true.beta, true.phi, 
                                        est.theta, est.beta, est.phi, 
                                        true_bagel_fn, est_bagel_fn, 
                                        plt_title){
  # This function creates BAGEL plots for MM-LDA parameters
  #
  # Args:
  #   theta: (D x K)
  #   beta: (K x I)
  #   phi: (K x J)
  # Create true parameter plots
  D = dim(true.theta)[1]
  K = dim(true.theta)[2]
  I = dim(true.beta)[2]
  J = dim(true.phi)[2]
  true.beta <- t(true.beta)
  colnames(true.beta) <- paste("Topic", 1:K, sep = "")
  rownames(true.beta) <- lapply(1:I, toString)
  true.theta <- t(true.theta)
  colnames(true.theta) <- lapply(1:D, toString)
  rownames(true.theta) <- paste("Topic", 1:K, sep = "")
  true.phi <- t(true.phi)
  colnames(true.phi) <- paste("Topic", 1:K, sep = "")
  rownames(true.phi) <- lapply(1:J, toString)
  # Convert matrices to BAGEL objects
  setClass("Result", representation(exposures="matrix", signatures="matrix"))
  true.bag <- methods::new("Result", exposures=true.theta, signatures=true.beta)
  p1 <- BAGEL::plot_exposures(true.bag, proportional = TRUE, 
                              label_samples = FALSE, no_legend = TRUE)
  p1 <- p1 + ggplot2::ggtitle('Theta: Topic-Document Proportions')
  p2 <- BAGEL::plot_signatures(true.bag, no_legend = TRUE, text_size = 15, 
                               facet_size = 10)
  p2 <- p2 + ggtitle("Beta: Topic-Image Proportion")  
  true.bag <- methods::new("Result", exposures=true.theta, signatures=true.phi)
  p3 <- BAGEL::plot_signatures(true.bag, no_legend = TRUE, text_size = 10, 
                               facet_size = 10)
  p3 <- p3 + ggplot2::ggtitle("Phi: Topic-Word Proportions") 
  true.plt <- ggpubr::ggarrange(p1, p2, p3, nrow=3, ncol=1, common.legend = FALSE)
  true.plt <- ggpubr::annotate_figure(true.plt, top="True Parameters", fig.lab.face="bold")
  ggplot2::ggsave(true_bagel_fn, true.plt, width=10, height=15, units="in")
  # Create estimated parameter plots
  est.beta <- t(est.beta)
  colnames(est.beta) <- paste("Topic", 1:K, sep = "")
  rownames(est.beta) <- lapply(1:I, toString)
  est.theta <- t(est.theta)
  colnames(est.theta) <- lapply(1:D, toString)
  rownames(est.theta) <- paste("Topic", 1:K, sep = "")
  est.phi <- t(est.phi)
  colnames(est.phi) <- paste("Topic", 1:K, sep = "")
  rownames(est.phi) <- lapply(1:J, toString)
  # Convert matrices to BAGEL objects
  est.bag <- methods::new("Result", exposures=est.theta, signatures=est.beta)
  p4 <- BAGEL::plot_exposures(est.bag, proportional = TRUE, 
                              label_samples = FALSE, no_legend = TRUE)
  p4 <- p4 + ggplot2::ggtitle('Theta: Topic-Document Proportions')
  p5 <- BAGEL::plot_signatures(est.bag, no_legend = TRUE, text_size = 15, 
                               facet_size = 10)
  p5 <- p5 + ggplot2::ggtitle("Beta: Topic-Image Proportion")  
  est.bag <- methods::new("Result", exposures=est.theta, signatures=est.phi)
  p6 <- BAGEL::plot_signatures(est.bag, no_legend = TRUE, text_size = 10, 
                               facet_size = 10)
  p6 <- p6 + ggplot2::ggtitle("Phi: Topic-Word Proportions") 
  est.plt <- ggpubr::ggarrange(p4, p5, p6, nrow=3, ncol=1, common.legend = FALSE)
  est.plt <- ggpubr::annotate_figure(est.plt, top="Estimated Parametrs", fig.lab.face="bold")
  ggplot2::ggsave(est_bagel_fn, est.plt, width=10, height=15, units="in")
  # # Combine plots
  # plt <- ggpubr::ggarrange(true.plt, est.plt, nrow=1, ncol=2, common.legend = TRUE)
  # plt <- ggpubr::annotate_figure(plt, top=plt_title, fig.lab.face = "bold")
  # ggplot2::ggsave(bagel_fn, plt, width=15, height=35, units="in")
}

create_test_lda_bagel_plt <- function(true.theta, true.beta, est.theta, 
                                      est.beta, true_bagel_fn, est_bagel_fn, 
                                      plt_title){
  # This function creates BAGEL plots for LDA parameters. It creates two plots, 
  # one containing the true parameters and ther other with estimated parameters.
  #
  # Args:
  #   theta: (D x K) doc-topic proportion matrix
  #   beta:  (K x I) topic-word proportion matrix
  # Create true parameter plots
  D = dim(true.theta)[1]
  K = dim(true.theta)[2]
  I = dim(true.beta)[2]
  true.beta <- t(true.beta)
  colnames(true.beta) <- paste("Topic", 1:K, sep = "")
  rownames(true.beta) <- lapply(1:I, toString)
  true.theta <- t(true.theta)
  colnames(true.theta) <- lapply(1:D, toString)
  rownames(true.theta) <- paste("Topic", 1:K, sep = "")
  # Convert matrices to BAGEL objects
  setClass("Result", representation(exposures="matrix", signatures="matrix"))
  true.bag <- methods::new("Result", exposures=true.theta, signatures=true.beta)
  p1 <- BAGEL::plot_exposures(true.bag, proportional=TRUE, label_samples=FALSE, 
                              no_legend=TRUE)
  p1 <- p1 + ggplot2::ggtitle('Theta: Topic-Document Proportions')
  p2 <- BAGEL::plot_signatures(true.bag, no_legend = TRUE, text_size = 15, 
                               facet_size = 10)
  p2 <- p2 + ggplot2::ggtitle("Beta: Topic-Image Proportion")  
  true.plt <- ggpubr::ggarrange(p1, p2, nrow=2, ncol=1, common.legend = FALSE)
  true.plt <- ggpubr::annotate_figure(true.plt, top="True Parameters", fig.lab.face="bold")
  ggplot2::ggsave(true_bagel_fn, true.plt, width=10, height=15, units="in")
  # Create estimated parameter plots
  est.beta <- t(est.beta)
  colnames(est.beta) <- paste("Topic", 1:K, sep = "")
  rownames(est.beta) <- lapply(1:I, toString)
  est.theta <- t(est.theta)
  colnames(est.theta) <- lapply(1:D, toString)
  rownames(est.theta) <- paste("Topic", 1:K, sep = "")
  # Convert matrices to BAGEL objects
  est.bag <- methods::new("Result", exposures=est.theta, signatures=est.beta)
  p4 <- BAGEL::plot_exposures(est.bag, proportional = TRUE, 
                              label_samples = FALSE, no_legend = TRUE)
  p4 <- p4 + ggplot2::ggtitle('Theta: Topic-Document Proportions')
  p5 <- BAGEL::plot_signatures(est.bag, no_legend = TRUE, text_size = 15, 
                               facet_size = 10)
  p5 <- p5 + ggplot2::ggtitle("Beta: Topic-Image Proportion")  
  est.plt <- ggpubr::ggarrange(p4, p5, nrow=2, ncol=1, common.legend = FALSE)
  est.plt <- ggpubr::annotate_figure(est.plt, top="Estimated Parameters", fig.lab.face="bold")
  ggplot2::ggsave(est_bagel_fn, est.plt, width=10, height=15, units="in")
  # Combine plots
  #plt <- ggpubr::ggarrange(true.plt, est.plt, nrow=1, ncol=2, common.legend = TRUE)
  #plt <- ggpubr::annotate_figure(plt, top=plt_title, fig.lab.face = "bold")
  #ggplot2::ggsave(bagel_fn, plt, width=15, height=35, units="in")
}

create_panlung_mmlda_bagel_plt <- function(theta, beta, phi, Rcolnames, 
                                           Rrownames, Wcolnames, K, bagel_fn, 
                                           plt_title){
  # This function creates the bagel plots for results from the panlung data
  # Assumes that theta is ?x?
  # Assumes that beta is ?x?
  # Assumes that phi is ?x?
  # Transform parameters
  beta <- t(beta)
  rownames(beta) <- Rcolnames
  colnames(beta) <- paste("Signature", 1:K, sep = "")
  theta <- t(theta)
  rownames(theta) <- paste("Signature", 1:K, sep = "")
  colnames(theta) <- Rrownames
  phi <- t(phi)
  rownames(phi) <- Wcolnames
  colnames(phi) <- paste("Signature", 1:K, sep = "")
  # Convert matrices to BAGEL objects
  setClass("Result", representation(exposures="matrix", signatures="matrix"))
  bag <- methods::new("Result", exposures=theta, signatures=beta)
  # Create BAGEL plots
  p1 <- BAGEL::plot_exposures(bag, proportional = TRUE, 
                              label_samples = FALSE, 
                              sort_samples = paste("Signature", K:1, sep=""))
  p1 <- p1 + ggtitle('Theta - Exposures')
  p2 <- BAGEL::plot_signatures(bag, no_legend = FALSE, text_size = 15, 
                               facet_size = 10)
  p2 <- p2 + ggtitle("Beta - DBS Signatures")  
  bag <- methods::new("Result", exposures=theta, signatures=phi)
  p3 <- BAGEL::plot_signatures(bag, no_legend = FALSE, text_size = 10, 
                               facet_size = 10)
  p3 <- p3 + ggtitle("Phi - SBS Signatures") 
  plt <- ggpubr::ggarrange(p1, p2, p3, nrow=3, ncol=1, common.legend = FALSE)
  plt <- ggpubr::annotate_figure(plt, top=plt_title, fig.lab.face="bold")
  ggsave(bagel_fn, plt, width=15, height=35, units="in")
}


# # ------
# make_final_bagel_plots <- function(data, result, data_type_long='Sparse', param.est='maxll'){
#   # Make BAGEL plots
#   # Check the dimensionality of the true parameters
#   if (!all(dim(data$theta) == c(result$data$D,result$data$K))){stop("theta dim incorrect")}
#   if (!all(dim(data$phi) == c(result$data$K,result$data$J))){stop("phi dim incorrect")}
#   if (!all(dim(data$beta) == c(result$data$K,result$data$I))){stop("beta dim incorrect")}
#   # Create bagel plots for true parameters
#   plt_left = "True Parameter"
#   truth_plt = make_bagel_plots(data$theta, data$beta, data$phi, plt_left)
#   # Find new topic assignmnets
#   new.topic = greedy_topic_permute(data$theta, data$beta, data$phi, resfit=result$fit, param.est=param.est)
#   # Create bagel plots for estimated parameters
#   theta.allchains = calculate_parameter(result$fit, param.name='theta', param.est=param.est)
#   beta.allchains = calculate_parameter(result$fit, param.name='beta', param.est=param.est)
#   phi.allchains = calculate_parameter(result$fit, param.name='phi', param.est=param.est)
#   nchains = dim(theta.allchains)[2]
#   plt_left = "Estimated Parameter"
#   chain_plts = vector("list", length = nchains)
#   for (ii in 1:nchains){
#     # Make plt for estimated parameter
#     theta = matrix(theta.allchains[,ii], nrow = data$D, ncol = data$K)
#     theta <- theta[,new.topic[ii,]]
#     beta = matrix(beta.allchains[,ii], nrow = data$K, ncol = data$I)
#     beta <- beta[new.topic[ii,],]
#     phi = matrix(phi.allchains[,ii], nrow = data$K, ncol = data$J)
#     phi <- phi[new.topic[ii,],]
#     est_plt = make_bagel_plots(theta, beta, phi, plt_left)
#     # Make combined plt
#     #plt <- ggpubr::ggarrange(truth_plt, est_plt, nrow = 1, ncol = 2)
#     #plt_title = paste0(data_type_long," Data: K=", data$K, ", D=", data$D, ", J=", data$J, ", I=", data$I, "\nParameter Estimation Method: ", param.est,"\nchain ", ii)
#     #plt <- ggpubr::annotate_figure(plt, top=plt_title, fig.lab.face = "bold")
#     #chain_plts[[ii]] = plt
#     chain_plts[[ii]] <- ggpubr::ggarrange(truth_plt, est_plt, nrow = 1, ncol = 2)
#   }
#   return(chain_plts)
# }
# 
# make_bagel_plots <- function(theta, beta, phi, plt_left){
#   # This function creates and saves 'exposures' and 'signatures' plots from BAGEL
#   #
#   # Inputs:
#   #   theta - (K x D) matrix of topic-document probabilities
#   #   beta - (I x K) matrix of img.region-topic probabilities
#   #   phi - (J x K) matrix of word-topic probabilities
#   #   plt_fn - filename of plot
#   #   plt_title - main title of plot
#   # Extract dimensions
#   K = dim(theta)[2]
#   D = dim(theta)[1]
#   J = dim(phi)[2]
#   plt_wrd_labs = J > 20
#   I = dim(beta)[2]
#   plt_img_labs = I > 20
#   # Make doc-topic plt
#   X = methods::new("Result", samples=t(theta), signatures=t(phi))
#   rownames(X@samples) = paste("Topic", 1:K, sep = "")
#   colnames(X@samples) = paste("Doc", 1:D, sep="")
#   rownames(X@signatures) = paste("Word ", 1:J, sep="")
#   colnames(X@signatures) = paste("Topic", 1:K, sep="")
#   # testing
#   #print(dim(X@samples))
#   #print(rownames(X@samples))
#   #print(colnames(X@samples))
#   #print(X)
#   #saveRDS(X, "/Users/kelly/Desktop/example.rds")
#   # end of testing
#   p1 <- BAGEL::plot_exposures(X, sort_samples = paste("Topic", 1:K, sep=""))
#   #p1 <- ggplot() + theme_void()
#   p2 <- BAGEL::plot_signatures(X, no_legend=plt_wrd_labs)
#   p2 <- ggpubr::annotate_figure(p2, top="Word-Topic Proportion")
#   # Make img plts
#   X = methods::new("Result", samples=t(theta), signatures=t(beta))
#   colnames(X@samples) = paste("Doc", 1:D, sep = "")
#   rownames(X@samples) = paste("Topic", 1:K, sep="")
#   colnames(X@signatures) = paste("Topic ", 1:K, sep="")
#   rownames(X@signatures) = paste("Img", 1:I, sep="")
#   p4 <- BAGEL::plot_signatures(X, no_legend=plt_img_labs)
#   p4 <- ggpubr::annotate_figure(p4, top = "Img-Topic Proportion")
#   # Combine and save
#   plt <- ggpubr::ggarrange(p1, p2, p4, nrow=3, ncol=1, common.legend = FALSE)
#   plt <- ggpubr::annotate_figure(plt, top=plt_left, fig.lab.face = "bold")
#   #plt <- ggpubr::annotate_figure(plt, top=plt_title, fig.lab.face="bold")
#   #ggsave(plt_fn, plt, width=11, height=11, units="in")
#   return(plt)
# }
# 
# make_panlung_bagel <- function(result, phi.names, beta.names, param.est='avg'){
#   # Make BAGEL plot using only estimated data
#   # estimate parameters
#   theta.maxll.allchains = calculate_parameter(result$fit, param.name='theta', param.est=param.est)
#   beta.maxll.allchains = calculate_parameter(result$fit, param.name='beta', param.est=param.est)
#   phi.maxll.allchains = calculate_parameter(result$fit, param.name='phi', param.est=param.est)
#   nchains = dim(theta.maxll.allchains)[2]
#   plt_left = ""
#   chain_plts = vector("list", length = nchains)
#   print(paste0("nchains = ", nchains))
#   K = result$data$K
#   D = result$data$D
#   J = result$data$J
#   I = result$data$I
#   for (ii in 1:nchains){
#     # Make plt for estimated parameter
#     theta = matrix(theta.maxll.allchains[,ii], nrow = D, ncol = K)
#     beta = matrix(beta.maxll.allchains[,ii], nrow = K, ncol = I)
#     phi = matrix(phi.maxll.allchains[,ii], nrow = K, ncol = J)
#     # Make doc-topic plt
#     X = methods::new("Result", samples=t(theta), signatures=t(phi))
#     rownames(X@samples) = paste("Topic", 1:K, sep = "")
#     colnames(X@samples) = paste("Doc", 1:D, sep="")
#     rownames(X@signatures) = phi.names
#     colnames(X@signatures) = paste("Topic", 1:K, sep="")
#     p1 <- BAGEL::plot_exposures(X, label_samples = FALSE, 
#                                 sort_samples = paste("Topic", 1:K, sep=""))
#     p2 <- BAGEL::plot_signatures(X) + ggtitle("Double-base Variants")
#     # Make img plts
#     X = methods::new("Result", samples=t(theta), signatures=t(beta))
#     colnames(X@samples) = paste("Gene ", 1:D, sep = "")
#     rownames(X@samples) = paste("Topic ", 1:K, sep="")
#     colnames(X@signatures) = paste("Topic ", 1:K, sep="")
#     rownames(X@signatures) = beta.names
#     p4 <- BAGEL::plot_signatures(X) + ggtitle("Single-base Variants")
#     p5 <- ggpubr::ggarrange(p2, p4, nrow = 1, ncol = 2)
#     # Combine and save
#     plt <- ggpubr::ggarrange(p1, p5, nrow=2, ncol=1, common.legend = FALSE)
#     #plt_fn = file.path(panlung_dir, paste0("panlung_bagel_chain",toString(ii),"_",param.est,".pdf"))
#     #plt_title = paste0("PanLung Data: ", K, " Topics\nPosterior Estimation Method: ", param.est,"\nchain ", ii)
#     plt_title = paste0("PanLung Data: K=",result$data$K,", D=",result$data$D,", J=",result$data$J,", I=",result$data$I,"\nParameter Estimation Method: ",param.est,"\nchain ",ii)
#     plt <- ggpubr::annotate_figure(plt, top=plt_title, fig.lab.face="bold")
#     plt_fn = file.path(result_dir, paste0("panlung_K",toString(K),"_gibbs",toString(niters),"_chain",toString(ii),"_bagel.pdf"))
#     ggsave(plt_fn, plt, width=15, height=17, units="in")
#   }
# }
