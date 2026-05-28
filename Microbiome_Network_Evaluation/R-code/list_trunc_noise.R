#######################################################################
####         Noisy data from truncated normal distribution         ####
#######################################################################

#' Constructs a list of noisy datasets from truncated normal distribution
#'
#' @param Data A data frame or numeric matrix
#' @param percent Percentage of elements receiving noise from each column (between 0 and 1)
#' @param num The length of the list
#' @return A list of length \code{num}, each element being a noisy version of \code{Data}
#' @export
#'
List_TRUNC_NoisyData <- function(Data, percent, num) {
  replicate(num, TRUNC_NoisyData(Data, percent), simplify = FALSE)
}

#' Add truncated normal noise to selected elements
#'
#' For each column with positive standard deviation, a random subset of rows
#' (size = ceiling(percent * nrow)) is selected. Noise is drawn from a truncated
#' normal distribution with mean = 0 and sd = sd(column)/2, truncated between
#' -original_value and the column's maximum value. The result is rounded to integers.
#'
#' @param Data A data frame or numeric matrix
#' @param percent Numeric between 0 and 1. Proportion of rows to perturb per column.
#' @return Data frame with noise added to specified elements (rounded integers)
#' @export
#'
TRUNC_NoisyData <- function(Data, percent) {
  n <- nrow(Data)
  col_sds <- apply(Data, 2, sd, na.rm = TRUE)
  noisy <- Data 
  
  for (j in seq_len(ncol(Data))) {
    if (col_sds[j] > 0) {
      idx <- sample.int(n, size = ceiling(percent * n))
      a_vec <- -Data[idx, j]
      b_vec <- max(Data[, j], na.rm = TRUE)
      noise <- truncnorm::rtruncnorm(length(idx), a = a_vec, b = b_vec, mean = 0, sd = col_sds[j] / 2)
      noisy[idx, j] <- noisy[idx, j] + noise
    }
  }
  return(round(noisy))
}

#' @example
amgut1_truncnoise05_100 <- List_TRUNC_NoisyData(amgut1.filt, 0.05, 100)
amgut1_truncnoise20_100 <- List_TRUNC_NoisyData(amgut1.filt, 0.2, 100)
