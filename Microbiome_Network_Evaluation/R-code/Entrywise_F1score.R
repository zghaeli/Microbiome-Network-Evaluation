#====================================
#         calculate_fscore          #
#====================================

#' Calculate F1-Score with Automatic Binary Conversion
#'
#' Computes the F1-score for a predicted matrix compared to a true matrix, after
#' converting both to binary (non-zero to 1, zero to 0). The F1-score is the harmonic
#' mean of precision and recall. Returns 0 if both precision and recall are zero.
#'
#' @param predicted_matrix A numeric matrix of predicted values (e.g., from a model).
#' @param true_matrix A numeric matrix of true values (ground truth).
#'
#' @return A numeric value representing the F1-score (between 0 and 1).
#' @export
#'
#'
calculate_fscore <- function(predicted_matrix, true_matrix) {
  # Check dimensions
  if (!identical(dim(true_matrix), dim(predicted_matrix))) {
    stop("Predicted and true matrices must have identical dimensions")
  }
  
  # Convert to binary (TRUE for non-zero, FALSE for zero)
  true_binary <- true_matrix != 0
  pred_binary <- predicted_matrix != 0
  
  # Calculate confusion matrix components
  true_positives  <- sum(true_binary & pred_binary, na.rm = TRUE)
  false_positives <- sum(!true_binary & pred_binary, na.rm = TRUE)
  false_negatives <- sum(true_binary & !pred_binary, na.rm = TRUE)
  
  # Precision and recall
  precision <- ifelse(true_positives + false_positives == 0, 0,
                      true_positives / (true_positives + false_positives))
  recall <- ifelse(true_positives + false_negatives == 0, 0,
                   true_positives / (true_positives + false_negatives))
  
  # F1-score
  fscore <- ifelse(precision + recall == 0, 0,
                   2 * (precision * recall) / (precision + recall))
  
  return(fscore)
}



#' @example
library(vegan)
source(Entry_Fscore.R)

# Bootstrap
amgut1_bootstrap_fscore <- lapply(amgut1_bootstrap, calculate_fscore, amgut1.filt)

# Noisy
amgut1_truncnoise05_fscore <- lapply(amgut1_truncnoise05_100, calculate_fscore, amgut1.filt)
amgut1_truncnoise20_fscore <- lapply(amgut1_truncnoise20_100, calculate_fscore, amgut1.filt)

# Synthetic
amgut1_gcoda_SE_fscore <- lapply(amgut1_gcoda_SE_100, calculate_fscore, amgut1.filt)
amgut1_gcoda_SP_fscore <- lapply(amgut1_gcoda_SP_100, calculate_fscore, amgut1.filt)
amgut1_cluster_SE_fscore <- lapply(amgut1_cluster_SE_100, calculate_fscore, amgut1.filt)
amgut1_cluster_SP_fscore <- lapply(amgut1_cluster_SP_100, calculate_fscore, amgut1.filt)


# par(mfrow = c(1,1), mar = c(2.5, 2, 2.5, 2), oma = c(10, 23, 10, 22)) #
boxplot(as.numeric(amgut1_cluster_SE_fscore),
        as.numeric(amgut1_gcoda_SE_fscore),
        as.numeric(amgut1_cluster_SP_fscore),
        as.numeric(amgut1_gcoda_SP_fscore),
        as.numeric(amgut1_truncnoise05_fscore),
        as.numeric(amgut1_truncnoise20_fscore),
        as.numeric(amgut1_bootstrap_fscore),
        names = c("cluster", "gCoda", "cluster", "gCoda", "5%", "20%", "Bootstrap"),
        col = c("red3", "red3", "blue3", "blue3", "green3", "green3", "darkgreen"),
        border = c("red3", "red3", "blue3", "blue3", "green3", "green3", "darkgreen"),
        ylim = c(0,1.2),
        cex = .7)
abline(v = 2.5, col = "gray", lty = 1, lwd = 1.5)
abline(v = 4.5, col = "gray", lty = 1, lwd = 1.5)
abline(v = 6.5, col = "gray", lty = 1, lwd = 1.5)
mtext(expression(bold("SE")), side = 3, at = .4, line = .1, cex = .7)
mtext(expression(bold("SP")), side = 3, at = 2.7, line = .1, cex = .7)
mtext(expression(bold("Noise")), side = 3, at = 4.8, line = .1, cex = .7)
mtext(expression(bold("Bootstrap")), side = 3, at = 7, line = .1, cex = .7)
mtext(expression(bold("F-score")), side = 2, line = 2.5, cex = 1, las = 0)
mtext(expression(bold("Matrix entry similarity")), side = 3, line = 1.7, cex = 1.1, las = 0)

# dev.off()

