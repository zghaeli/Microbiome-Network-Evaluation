#######################################################################
####         Noisy data from truncated normal distribution         ####
#######################################################################

#` Constructs a list of noisy datasets from truncated normal distribution
#
#' @param Data A numeric matrix
#' @param percent Percentage of elements receiving noise from each column
#' @param num The length of the list
#' @return The length of the list

List_TRUNC_NoisyData <- function(Data, percent, num) {
  INC <- Index__NN__Cols(Data)
  replicate(num, TRUNC_NoisyData(Data, PERCENT(INC, percent)), simplify = FALSE)
}

#' Add truncated normal noise to selected elements
#' @param Data A numeric matrix
#' @param SubsetList List of row indices per column
#' @return Data with noise added to specified elements
#' @noRd

TRUNC_NoisyData <- function(Data, SubsetList) {
  noisydata <- Data
  for (i in seq_len(ncol(Data))) {
    idx <- SubsetList[[i]]
    if (length(idx) > 0) {
      noise <- truncnorm::rtruncnorm(length(idx), a = min(Data[, i]), b = max(Data[, i]), 
                                     mean = mean(Data[, i]), sd = sd(Data[, i]))
      noisydata[idx, i] <- noisydata[idx, i] + noise
    }
  }
  round_noisydata <- round(noisydata)
  return(round_noisydata)
}

#' Sample a percentage of indices from a list
#' @param List List of index vectors
#' @param perc Numeric between 0 and 1
#' @return List of randomly sampled indices
#' @noRd
PERCENT <- function(List, perc) { lapply(List, function(x) sample(x, size = perc * length(x))) }

#' Get row indices of non-zero elements in each column
#' @param Data A numeric matrix
#' @return List of row indices for non-zero values per column
#' @noRd
Index__NN__Cols <- function(Data){ lapply(1:ncol(Data), function(i) {which(Data[, i] != 0)}) } 


#' @example
amgut1_noise05_100 <- List_TRUNC_NoisyData(amgut1.filt, 0.05, 100)
amgut1_noise20_100 <- List_TRUNC_NoisyData(amgut1.filt, 0.2, 100)


