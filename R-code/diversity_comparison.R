######################################################################################
####          Compare Diversity Indices (Shannon, InvSimpson, Richness)           ####
######################################################################################

#' Compare Multiple Diversity Indices Between Real Data and a List of Simulated Datasets
#'
#' Computes Shannon, Inverse Simpson, and Richness for a real matrix and for each
#' simulated matrix in a list. Then performs paired or unpaired Wilcoxon tests
#' between each simulated replicate and the real data, adjusts p-values for multiple
#' comparisons (FDR), and returns detailed results.
#'
#' @param list_data A list of numeric matrices or data frames (rows = samples, columns = species).
#' @param real_matrix A numeric matrix or data frame with the same structure as each element of \code{list_data}.
#' @param paired Logical. If \code{TRUE} (default), uses paired Wilcoxon test (requires same number of rows). 
#'        If \code{FALSE}, uses unpaired test.
#' @param adjust_method Method for p-value adjustment (see \code{p.adjust}). Default: "BH" (Benjamini-Hochberg).
#' @param seed Integer seed for reproducibility (default: 123). Set to \code{NULL} to avoid setting seed.
#'
#' @return A list with five components:
#' \item{real_values}{List of three vectors: Shannon, InvSimpson, Richness for the real data.}
#' \item{bootstrap_mats}{List of three matrices (each row = replicate, column = sample) for the three indices.}
#' \item{pvalues}{List of three data frames (Shannon, InvSimpson, Richness) with columns: replicate, p_raw, p_adj, significant.}
#' \item{summary}{Data frame with columns: Index (diversity measure), N_replicates, Significant_after_FDR, Percent_significant.}
#' \item{settings}{List of input parameters.}
#' @export
#'
#'
compare_diversity_list <- function(list_data, real_matrix, paired = TRUE, 
                                   adjust_method = "BH", seed = 123) {
  if (!is.null(seed)) set.seed(seed)
  
  # Real diversity indices
  real_div <- list(
    Shannon    = vegan::diversity(real_matrix, index = "shannon"),
    InvSimpson = vegan::diversity(real_matrix, index = "invsimpson"),
    Richness   = vegan::specnumber(real_matrix)
  )
  
  # Diversity matrices for each simulated dataset
  div_mats <- list_data %>%
    purrr::map(~ list(
      shannon    = vegan::diversity(.x, index = "shannon"),
      invsimpson = vegan::diversity(.x, index = "invsimpson"),
      specnumber = vegan::specnumber(.x)
    )) %>%
    purrr::transpose() %>%
    purrr::map(~ do.call(cbind, .x)) %>%
    rlang::set_names(c("shannon_mat", "invsimpson_mat", "specnumber_mat"))
  
  # Wilcoxon tests and p-value adjustment
  pvalues_list <- purrr::map2(real_div, div_mats, ~ {
    p_raw <- apply(.y, 2, function(col) {
      stats::wilcox.test(.x, col, paired = paired)$p.value
    })
    tibble::tibble(
      replicate = seq_along(p_raw),
      p_raw = p_raw,
      p_adj = stats::p.adjust(p_raw, method = adjust_method),
      significant = p_adj < 0.05
    )
  }) %>%
    rlang::set_names(c("Shannon", "InvSimpson", "Richness"))
  
  # Summary table
  summary <- purrr::map2_dfr(pvalues_list, names(pvalues_list), ~ {
    tibble::tibble(
      Index = .y,
      N_replicates = nrow(.x),
      Significant_after_FDR = sum(.x$significant),
      Percent_significant = mean(.x$significant) * 100
    )
  })
  
  list(
    real_values    = real_div,
    bootstrap_mats = div_mats,
    pvalues        = pvalues_list,
    summary        = summary,
    settings       = list(paired = paired, adjust_method = adjust_method, 
                          n_replicates = length(list_data), seed = seed)
  )
}

#' @examples
amgut1_bootstrap_diversity <- compare_diversity_list(amgut1_bootstrap, amgut1.filt, paired = TRUE, adjust_method = "BH")

amgut1_truncnoise05_diversity <- compare_diversity_list(amgut1_truncnoise05_100, amgut1.filt)
amgut1_truncnoise20_diversity <- compare_diversity_list(amgut1_truncnoise20_100, amgut1.filt)

amgut1_gcoda_SE_diversity <- compare_diversity_list(amgut1_gcoda_SE_100, amgut1.filt)
amgut1_gcoda_SP_diversity <- compare_diversity_list(amgut1_gcoda_SP_100, amgut1.filt)

amgut1_cluster_SE_diversity<- compare_diversity_list(amgut1_cluster_SE_100, amgut1.filt)
amgut1_cluster_SP_diversity<- compare_diversity_list(amgut1_cluster_SP_100, amgut1.filt)


# ============================================================
# Plots for Richness p-values using new compare_diversity_list
# ============================================================

par(mar = c(2.5, 2, 2.5, 2), mfrow = c(2,2), oma = c(6, 19, 5, 28))

# Panel 1: cluster vs gCoda (SE)
vioplot::vioplot(amgut1_cluster_SE_diversity$pvalues$Richness$p_raw,
                 amgut1_gcoda_SE_diversity$pvalues$Richness$p_raw,
                 col = c(rep("deepskyblue1", 2)),
                 cex.axis = 0.8,
                 names = c("cluster", "gCoda"))
mtext(expression(bold("P-value")), side = 2, line = 2, cex = 0.9, las = 0)
mtext(expression(bold("Topology")), side = 1, line = 2.5, cex = 0.9)
mtext(expression(bold("SE")), side = 3, line = 0.1, cex = 0.7, at = 0.5)

# Panel 2: Noise levels 5% vs 20%
boxplot(amgut1_truncnoise05_diversity$pvalues$Richness$p_raw,
        amgut1_truncnoise20_diversity$pvalues$Richness$p_raw,
        col = c("tan", "salmon1"),
        border = c("tan", "salmon1"),
        cex.axis = 0.8,
        names = c("5%", "20%"))
mtext(expression(bold("P-value")), side = 2, line = 2, cex = 0.9, las = 0)
mtext(expression(bold("Noise")), side = 1, line = 2.2, cex = 0.9)
mtext(expression(bold("Noise")), side = 3, line = 0.1, cex = 0.7, at = 0.57)

# Panel 3: cluster vs gCoda (SP)
vioplot::vioplot(amgut1_cluster_SP_diversity$pvalues$Richness$p_raw,
                 amgut1_gcoda_SP_diversity$pvalues$Richness$p_raw,
                 col = c(rep("mediumvioletred", 2)),
                 cex.axis = 0.8,
                 names = c("cluster", "gCoda"))
mtext(expression(bold("P-value")), side = 2, line = 2, cex = 0.9, las = 0)
mtext(expression(bold("Topology")), side = 1, line = 2.5, cex = 0.9)
mtext(expression(bold("SP")), side = 3, line = 0.1, cex = 0.7, at = 0.5)

# Panel 4: Bootstrap
vioplot::vioplot(amgut1_bootstrap_diversity$pvalues$Richness$p_raw, cex.axis = 0.8, col = "mediumseagreen")
mtext(expression(bold("P-value")), side = 2, line = 2, cex = 0.9, las = 0)
mtext(expression(bold("Bootstrap")), side = 1, line = 2.5, cex = 0.9)
mtext(expression(bold("Bootstrap")), side = 3, line = 0.1, cex = 0.7, at = 0.57)

# Outer title
mtext(expression(bold("Richness__amgut1")), side = 3, line = -0.7, outer = TRUE, cex = 1.2)


