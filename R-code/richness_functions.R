# richness_functions.R
#
# Function to compare species richness between real data and a list of datasets using Wilcoxon tests.
# Requires: vegan package for specnumber, significance_functions.R for wilcoxon_pvalue.

#' Compare Species Richness
#'
#' Computes species richness for a real dataset and a list of datasets using
#' \code{vegan::specnumber}, then compares each with the real dataset using
#' a Wilcoxon rank-sum test via \code{wilcoxon_pvalue}.
#'
#' @param real_data A numeric matrix or data frame with non-negative values (rows as samples, columns as species).
#' @param list_data A list of numeric matrices or data frames with the same structure as \code{real_data}.
#' @return A named list of p-values from Wilcoxon tests comparing the richness of each dataset in \code{list_data} to \code{real_data}.
#' 
compare_richness <- function(real_data, list_data) {
  # Calculate richness for real_data
  real_richness <- vegan::specnumber(real_data)
  # Calculate richness for each dataset in list_data
  list_richness <- lapply(list_data, vegan::specnumber)
  # Perform Wilcoxon tests and name the results
  p_values <- lapply(seq_along(list_richness), function(i) {
    wilcoxon_pvalue(list_richness[[i]], real_richness)
  })
  return(p_values)
}

#' Perform Wilcoxon Test for Two Datasets
#'
#' Performs a Wilcoxon rank-sum test to compare two numeric datasets and returns
#' the p-value to assess statistical significance.
#'
#' @param synthetic_data A numeric vector of synthetic data.
#' @param real_data A numeric vector of real data.
#' @return A numeric value representing the p-value from the Wilcoxon test.
#' 
#' @noRd
wilcoxon_pvalue <- function(synthetic_data, real_data) {
  wilcox_result <- stats::wilcox.test(synthetic_data, real_data)
  return(wilcox_result$p.value)
}

#' @example
#' library(vegan)
#' source(richness_functions.R)
amgut1_bootstrap_richness_pvalues <- compare_richness(amgut1.filt, amgut1_bootstrap)

amgut1_noise05_richness_pvalues <- compare_richness(amgut1.filt, amgut1_noise05_100)
amgut1_noise20_richness_pvalues <- compare_richness(amgut1.filt, amgut1_noise20_100)

amgut1_gcoda_SE_richness_pvalues <- compare_richness(amgut1.filt, amgut1_gcoda_SE_100)
amgut1_gcoda_SP_richness_pvalues <- compare_richness(amgut1.filt, amgut1_gcoda_SP_100)

amgut1_cluster_SE_richness_pvalues <- compare_richness(amgut1.filt, amgut1_cluster_SE_100)
amgut1_cluster_SP_richness_pvalues <- compare_richness(amgut1.filt, amgut1_cluster_SP_100)

 
# ##*****************            Richness             ****************
par(mar = c(2.5, 2, 2.5, 2), mfrow = c(2,2), oma = c(6, 19, 5, 28))

vioplot::vioplot(as.numeric(amgut1_cluster_SE_richness_pvalues),
                 as.numeric(amgut1_gcoda_SE_richness_pvalues),
                 col = c(rep("deepskyblue1", 2)),
                 cex.axis = .8,
                 names = c("cluster", "gCoda"))
mtext(expression(bold("P-value")),  side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Topology")), side = 1, line = 2.5, cex = .9)
mtext(expression(bold("SE")), side = 3, line = .1, cex = .7, at = .5)
#mtext(expression(bold("A")),        side = 3, line = .6,  cex = 2, at = -.5)

boxplot(as.numeric(amgut1_noise05_richness_pvalues), 
        as.numeric(amgut1_noise20_richness_pvalues),
        col = c("tan", "salmon1"),
        border = c("tan", "salmon1"),
        cex.axis = .8,
        names = c("5%", "20%"))
mtext(expression(bold("P-value")), side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Noise")), side = 1, line = 2.2, cex = .9)
mtext(expression(bold("Noise")), side = 3, line = .1, cex = .7, at = .57)
#mtext(expression(bold("C")),     side = 3,    line = .6, cex = 2, at = -.5)


vioplot::vioplot(as.numeric(amgut1_cluster_SP_richness_pvalues),
                 as.numeric(amgut1_gcoda_SP_richness_pvalues),
                 col = c(rep("mediumvioletred", 2)),
                 cex.axis = .8,
                 names = c("cluster", "gCoda"))
mtext(expression(bold("P-value")), side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Topology")), side = 1, line = 2.5, cex = .9)
mtext(expression(bold("SP")), side = 3, line = .1, cex = .7, at = .5)
#mtext(expression(bold("A")),        side = 3, line = .6,  cex = 2, at = -.5)

vioplot::vioplot(as.numeric(amgut1_bootstrap_richness_pvalues), cex.axis = .8, col = "mediumseagreen")
mtext(expression(bold("P-value")), side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Bootstrap")), side = 1, line = 2.5, cex = .9)
mtext(expression(bold("Bootstrap")), side = 3, line = .1, cex = .7, at = .57)
#mtext(expression(bold("D")),         side = 3, line = .6, cex = 2, at = .32)

mtext(expression(bold("Richness__amgut1")), side = 3, line = -.7, outer = TRUE, cex = 1.2)

# dev.off()
