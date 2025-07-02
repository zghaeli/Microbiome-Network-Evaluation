#' Run gcoda Analysis on a Dataset
#'
#' Performs gcoda analysis on input data and returns results.
#'
#' @param data Matrix or data frame for gcoda analysis.
#' @param Thre Numeric value for ebic.gamma parameter.
#' @return A list with `gC` (gcoda output) and `gC_adj` (refit adjacency matrix).

source("/Microbiome-Network-Evaluation/R-code/gcoda.R")

gC_NEW <- function(data, Thre) {
  # if (!is.matrix(data) && !is.data.frame(data)) stop("data must be a matrix or data frame")
  # if (!is.numeric(Thre) || Thre < 0) stop("Thre must be a non-negative numeric value")
  start <- Sys.time()
  gC <- gcoda(data, counts = TRUE, pseudo = 1, lambda.min.ratio = 1e-4, nlambda = 50, ebic.gamma = Thre)
  cat("gC_NEW execution time:", Sys.time() - start, "\n")
  list(gC = gC, gC_adj = gC$refit)
}

#' @example
gcoda_amgut1 <- gC_NEW(amgut1.filt, .85)
adj_gcoda_amgut1 <- gcoda_amgut1$gC_adj





#' Process gCoda Data and Compute Diagnostics
#'
#' Applies gCoda processing to sample lists, extracts adjusted matrices, computes
#' diagnostics, and saves results. Returns a list of processed data and metrics.
#'
#' @param data Numeric matrix of input data.
#' @param sample_list List of sample datasets.
#' @param gold Numeric matrix for diagnostic comparison.
#' @param threshold Numeric threshold for gCoda processing (default: 0.85).
#' @param save_dir Directory to save RDS files (default: current directory).
#' @param timed Logical to print processing time (default: FALSE).
#' @return List containing gCoda data, adjusted matrices, diagnostics, and F-scores (as a numeric vector).
#' @examples
#' \dontrun{
#'   data <- matrix(rnorm(100), nrow=10)
#'   samples <- replicate(5, data, simplify=FALSE)
#'   gold <- matrix(rnorm(100), 10, 10)
#'   process_gcoda(data, samples, gold)
#' }
process_gcoda <- function(data, sample_list, gold, threshold = 0.85, save_dir = ".", timed = FALSE) {
  start <- if (timed) Sys.time() else NULL
  sample_name <- deparse(substitute(sample_list))
  
  # Process gCoda data
  gc_data <- lapply(sample_list, gC_NEW, Thre = threshold)
  if (timed) print(Sys.time() - start)
  saveRDS(gc_data, file.path(save_dir, paste0("gcoda_", sample_name, ".rds")))
  
  # Extract adjusted matrices
  adj_matrices <- lapply(gc_data, `[[`, 2)
  saveRDS(adj_matrices, file.path(save_dir, paste0("adj_gcoda_", sample_name, ".rds")))
  
  # Compute diagnostics and F-scores
  diagnostics <- lapply(adj_matrices, function(adj) MY_Diagnostic(gold, adj, ncol(data)))
  f_scores <- vapply(diagnostics, function(d) d[[1]][10], numeric(1))
  saveRDS(diagnostics, file.path(save_dir, paste0("diag_gcoda_", sample_name, ".rds")))
  saveRDS(f_scores, file.path(save_dir, paste0("fscore_gcoda_", sample_name, ".rds")))
  
  list(gc_data = gc_data, adj_matrices = adj_matrices, diagnostics = diagnostics, f_scores = f_scores)
}

#' @example
amgut1_bootstrap_gcoda <- process_gcoda(amgut1.filt, amgut1_bootstrap, adj_gcoda_amgut1, threshold = 0.85, save_dir = ".", timed = TRUE)

amgut1_noise05_gcoda <- process_gcoda(amgut1.filt, amgut1_noise05_100, adj_gcoda_amgut1, threshold = 0.85, save_dir = ".", timed = TRUE)
amgut1_noise20_gcoda <- process_gcoda(amgut1.filt, amgut1_noise20_100, adj_gcoda_amgut1, threshold = 0.85, save_dir = ".", timed = TRUE)

amgut1_gcoda_SE_gcoda <- process_gcoda(amgut1.filt, amgut1_gcoda_SE_100, adj_gcoda_amgut1, threshold = 0.85, save_dir = ".", timed = TRUE)
amgut1_gcoda_SP_gcoda <- process_gcoda(amgut1.filt, amgut1_gcoda_SP_100, adj_gcoda_amgut1, threshold = 0.85, save_dir = ".", timed = TRUE)

amgut1_cluster_SE_gcoda <- process_gcoda(amgut1.filt, amgut1_cluster_SE_100[1:5], adj_cluster_amgut1, threshold = 0.85, save_dir = ".", timed = TRUE)
amgut1_cluster_SP_gcoda <- process_gcoda(amgut1.filt, amgut1_cluster_SP_100, adj_cluster_amgut1, threshold = 0.85, save_dir = ".", timed = TRUE)


# Plot
vioplot::vioplot(amgut1_cluster_SE_gcoda$f_scores,
                 amgut1_gcoda_SE_gcoda$f_scores,
                 amgut1_cluster_SP_gcoda$f_scores,
                 amgut1_gcoda_SP_gcoda$f_scores,
                 amgut1_noise05_gcoda$f_scores,
                 amgut1_noise20_gcoda$f_scores,
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





