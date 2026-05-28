#####################################################
#            KS test comparison for lists           #
#####################################################

#' Compare OTU distributions for a list of matrices against real data
#'
#' For each matrix in a list (e.g., noisy replicates), this function performs
#' Kolmogorov-Smirnov tests between each column of the simulated matrix and the
#' corresponding column of the real data. P-values are adjusted using the
#' Benjamini-Hochberg (BH) method. The function returns the proportion of OTUs
#' (columns) with adjusted p-value below a significance threshold for each matrix.
#'
#' @param list_of_matrices A list of numeric matrices (or data frames) with the same
#'        dimensions as \code{real_data}. Each matrix represents a simulated dataset.
#' @param real_data A numeric matrix or data frame of real data (rows = samples,
#'        columns = OTUs).
#' @param alpha Significance level (default: 0.05). Adjusted p-values below this
#'        threshold are considered significantly different.
#'
#' @return A list with one element: \code{percentages}, a numeric vector of the same
#'         length as \code{list_of_matrices}, where each entry is the proportion of
#'         OTUs with significant KS test for that simulated matrix.
#' @export
#'
#'
compare_OTU_distributions_to_real <- function(list_of_matrices, real_data, alpha = 0.05) {
  n_otus <- ncol(real_data)
  percentages <- numeric(length(list_of_matrices))
  
  for (i in seq_along(list_of_matrices)) {
    sim_data <- as.matrix(list_of_matrices[[i]])
    p_values <- sapply(1:n_otus, function(j) {
      ks_test <- ks.test(real_data[, j], sim_data[, j])
      return(ks_test$p.value)
    })
    p_adj <- p.adjust(p_values, method = "BH")
    n_significant <- sum(p_adj < alpha)
    percentages[i] <- n_significant / n_otus
  }
  
  result <- list(percentages = percentages)
  return(result)
}


#' @examples
amgut1_bootstrap_ks <- lapply(amgut1_bootstrap,  compare_OTU_distributions_to_real, amgut1.filt, 0.05)

amgut1_truncnoise05_ks <- lapply(amgut1_truncnoise05_100, compare_OTU_distributions_to_real, amgut1.filt, 0.05)
amgut1_truncnoise20_ks <- lapply(amgut1_truncnoise20_100, compare_OTU_distributions_to_real, amgut1.filt, 0.05)

amgut1_gcoda_SE_ks <- lapply(amgut1_gcoda_SE_100, compare_OTU_distributions_to_real, amgut1.filt, 0.05)
amgut1_gcoda_SP_ks <- lapply(amgut1_gcoda_SP_100, compare_OTU_distributions_to_real, amgut1.filt, 0.05)

amgut1_cluster_SE_ks <- lapply(amgut1_cluster_SE_100, compare_OTU_distributions_to_real, amgut1.filt, 0.05)
amgut1_cluster_SP_ks <- lapply(amgut1_cluster_SP_100, compare_OTU_distributions_to_real, amgut1.filt, 0.05)

# Plot
# ##*************       KS.TEST        ****************
# par(mar = c(2.5, 2, 2.5, 2), mfrow = c(2,2), oma = c(6, 19, 5, 28))

vioplot::vioplot(as.numeric(amgut1_cluster_SE_ks),
                 as.numeric(amgut1_gcoda_SE_ks),
                 col = c(rep("deepskyblue1", 2)),
                 cex.axis = .8,
                 ylim = c(.8,1),
                 names = c("cluster", "gCoda"))
mtext(expression(bold("% OTUs")),  side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Topology")), side = 1, line = 2.5, cex = .9)
mtext(expression(bold("SE")), side = 3, line = .1, cex = .7, at = .5)
#mtext(expression(bold("A")),        side = 3, line = .6,  cex = 2, at = -.5)

vioplot::vioplot(as.numeric(amgut1_truncnoise05_ks), 
                 as.numeric(amgut1_truncnoise20_ks),
                 col = c("tan", "salmon1"),
                 border = c("tan", "salmon1"),
                 cex.axis = .8,
                 ylim = c(0,1),
                 names = c("5%", "20%"))
mtext(expression(bold("% OTUs")), side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Noise")), side = 1, line = 2.2, cex = .9)
mtext(expression(bold("Noise")), side = 3, line = .1, cex = .7, at = .57)
#mtext(expression(bold("C")),     side = 3,    line = .6, cex = 2, at = -.5)

vioplot::vioplot(as.numeric(amgut1_cluster_SP_ks),
                 as.numeric(amgut1_gcoda_SP_ks),
                 col = c(rep("mediumvioletred", 2)),
                 cex.axis = .8,
                 names = c("cluster", "gCoda"))
mtext(expression(bold("% OTUs")), side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Topology")), side = 1, line = 2.5, cex = .9)
mtext(expression(bold("SP")), side = 3, line = .1, cex = .7, at = .5)
#mtext(expression(bold("A")),        side = 3, line = .6,  cex = 2, at = -.5)

vioplot::vioplot(as.numeric(amgut1_bootstrap_ks), cex.axis = .8, col = "mediumseagreen")
mtext(expression(bold("% OTUs")), side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Bootstrap")), side = 1, line = 2.5, cex = .9)
mtext(expression(bold("Bootstrap")), side = 3, line = .1, cex = .7, at = .57)
#mtext(expression(bold("D")),         side = 3, line = .6, cex = 2, at = .32)

mtext(expression(bold("KS__amgut1")), side = 3, line = -.7, outer = TRUE, cex = 1.2)

# dev.off()