########################################################
####           List of Bootstrap Resampling         ####
########################################################

#' Generate Bootstrapped Datasets
#' Constructs a list of resampled datasets
#'
#' @param data A numeric matrix of real-world dataset
#' @param num  The length of the list
#' @param num_sample The number of samples to be selected
#' @param Replace Sample with or without replacement?
#' 
#' @return List of resampled datasets

List_Bootstrap_Sample <- function(data, num, num_sample, Replace){
  lapply(1:num, function(i) {as.matrix(dplyr::sample_n(as.data.frame(data), size = num_sample, replace = Replace))})
}

#' @example
amgut1_bootstrap <- List_Bootstrap_Sample(amgut1.filt, nrow(amgut1.filt), nrow(amgut1.filt), TRUE)