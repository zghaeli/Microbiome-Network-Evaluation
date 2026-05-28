####################################################################
####           gCoda Network Inference and Evaluation           ####
####################################################################

# Dependencies:
# - The script "gcoda.R" must be sourced before using these functions.
#   It should contain the function gcoda() (e.g., from https://github.com/...)
# - The function MY_Diagnostic() must be defined elsewhere (user‑supplied).

#' Run gCoda Analysis on a Dataset (Fixed Parameters)
#'
#' Performs gCoda analysis on input count data using fixed parameters:
#' \code{pseudo = 0.5}, \code{lambda.min.ratio = 1e-4}, \code{nlambda = 15},
#' and \code{ebic.gamma = 0.5}. Returns the raw gcoda output object.
#'
#' @param data A matrix or data frame of count data (rows = samples, columns = OTUs).
#'
#' @return The full output from \code{gcoda()} (a list containing at least \code{path}).
#' @export
#'
gC_NEW <- function(data) { gcoda(data, counts = TRUE, pseudo = 0.5, lambda.min.ratio = 1e-4, nlambda = 15, ebic.gamma = 0.5) }

#' Process gCoda on a List of Sample Datasets and Compute Diagnostics
#'
#' Applies \code{gC_NEW} to each dataset in a list (e.g., bootstrap replicates),
#' extracts a specific adjacency matrix from the regularization path (using a
#' pre‑determined index), computes diagnostic metrics by comparing each extracted
#' network to a gold standard network, and saves all results as RDS files.
#'
#' The function uses \code{MY_Diagnostic(gold, adj, ncol(data))} to obtain
#' diagnostics (e.g., F1, precision, recall). You must define \code{MY_Diagnostic}
#' yourself; it should return a list whose 10th element (e.g., \code{d[[1]][10]})
#' is taken as the F‑score.
#'
#' @param data Original data matrix (used for dimension info, e.g., number of OTUs).
#' @param sample_list A list of matrices or data frames (each with same dimensions as \code{data}).
#' @param gold A gold standard adjacency matrix (e.g., inferred from the original data).
#' @param adj_index Integer; the index of the adjacency matrix to extract from the
#'        \code{path} component of each gcoda output (e.g., obtained from
#'        \code{select_adj_by_edge_fraction} on the original data).
#' @param save_dir Directory to save RDS files (default: current directory ".").
#' @param timed Logical; if \code{TRUE}, prints total processing time.
#'
#' @return A list with four components:
#' \item{gc_data}{List of full gcoda outputs for each sample.}
#' \item{adj_matrices}{List of extracted adjacency matrices (at index \code{adj_index}).}
#' \item{diagnostics}{List of diagnostic results from \code{MY_Diagnostic()}.}
#' \item{f_scores}{Numeric vector of F‑scores (the 10th element of each diagnostic).}
#' @export
#'
#'
process_gcoda <- function(data, sample_list, gold, adj_index, save_dir = ".", timed = FALSE) {
  start <- if (timed) Sys.time() else NULL
  sample_name <- deparse(substitute(sample_list))
  
  # Run gC_NEW on each sample in the list
  gc_data <- lapply(sample_list, gC_NEW)
  if (timed) print(Sys.time() - start)
  saveRDS(gc_data, file.path(save_dir, paste0("gcoda_", sample_name, ".rds")))
  
  # Extract adjacency matrices from the chosen index
  adj_matrices <- lapply(gc_data, function(x) x$path[[adj_index]])
  saveRDS(adj_matrices, file.path(save_dir, paste0("adj_gcoda_", sample_name, ".rds")))
  
  # Compute diagnostics and F-scores (requires MY_Diagnostic)
  diagnostics <- lapply(adj_matrices, function(adj) MY_Diagnostic(gold, adj, ncol(data)))
  f_scores <- vapply(diagnostics, function(d) d[[1]][10], numeric(1))
  saveRDS(diagnostics, file.path(save_dir, paste0("diag_gcoda_", sample_name, ".rds")))
  saveRDS(f_scores, file.path(save_dir, paste0("fscore_gcoda_", sample_name, ".rds")))
  
  list(gc_data = gc_data, adj_matrices = adj_matrices, diagnostics = diagnostics, f_scores = f_scores)
}

#' @example
amgut1_bootstrap_gcoda <- process_gcoda(amgut1.filt, amgut1_bootstrap, adj_gcoda_amgut1, index_adj_gcoda_amgut1, save_dir = ".", timed = TRUE)

amgut1_truncnoise05_gcoda <- process_gcoda(amgut1.filt, amgut1_truncnoise05_100, adj_gcoda_amgut1, index_adj_gcoda_amgut1, save_dir = ".", timed = TRUE)
amgut1_truncnoise20_gcoda <- process_gcoda(amgut1.filt, amgut1_truncnoise20_100, adj_gcoda_amgut1, index_adj_gcoda_amgut1, save_dir = ".", timed = TRUE)

amgut1_gcoda_SE_gcoda <- process_gcoda(amgut1.filt, amgut1_gcoda_SE_100, adj_gcoda_amgut1, index_adj_gcoda_amgut1, save_dir = ".", timed = TRUE)
amgut1_gcoda_SP_gcoda <- process_gcoda(amgut1.filt, amgut1_gcoda_SP_100, adj_gcoda_amgut1, index_adj_gcoda_amgut1, save_dir = ".", timed = TRUE)

amgut1_cluster_SE_gcoda <- process_gcoda(amgut1.filt, amgut1_cluster_SE_100, adj_cluster_amgut1, index_adj_gcoda_amgut1, save_dir = ".", timed = TRUE)
amgut1_cluster_SP_gcoda <- process_gcoda(amgut1.filt, amgut1_cluster_SP_100, adj_cluster_amgut1, index_adj_gcoda_amgut1, save_dir = ".", timed = TRUE)


# Plot
vioplot::vioplot(amgut1_cluster_SE_gcoda$f_scores,
                 amgut1_gcoda_SE_gcoda$f_scores,
                 amgut1_cluster_SP_gcoda$f_scores,
                 amgut1_gcoda_SP_gcoda$f_scores,
                 amgut1_truncnoise05_gcoda$f_scores,
                 amgut1_truncnoise20_gcoda$f_scores,
                 amgut1_bootstrap_gcoda$f_scores,
                 names = c("cluster", "gCoda", "cluster", "gCoda", "5%", "20%", "Bootstrap"),
                 col = c("gold", "#C2DF23FF", "gold", "#C2DF23FF", "#AEFFAE", "#63C163", "slateblue2"),
                 border = c("gray50", "gray50", "gray50", "gray50", "gray50", "gray50", "gray50"),
                 cex.axis = 1,
                 ylim = c(0,1))
abline(v = 2.5, col = "gray", lty = 1, lwd = 1.5)
abline(v = 4.5, col = "gray", lty = 1, lwd = 1.5)
abline(v = 6.5, col = "gray", lty = 1, lwd = 1.5)
mtext(expression(bold("SE")),        side = 3, at = .5,  line = .1,  cex = .7)
mtext(expression(bold("SP")),        side = 3, at = 2.7, line = .1,  cex = .7)
mtext(expression(bold("Noise")),     side = 3, at = 4.9, line = .1,  cex = .7)
mtext(expression(bold("Bootstrap")), side = 3, at = 6.9, line = .1,  cex = .7)
mtext(expression(bold("F-score")),   side = 2, las = 0,  line = 2.3, cex = 1.1)
title(paste0("gCoda"), cex.main = 1)
mtext(expression(bold("amgut1")), side = 3, at = .1, line = -.8, outer = TRUE, cex = 1.3) 

