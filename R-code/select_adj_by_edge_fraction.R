#' Select adjacency matrix from a list based on target edge fraction
#'
#' Given a list of square adjacency matrices (e.g., from a regularization path),
#' this function selects the matrix whose number of off‑diagonal non‑zero entries
#' (i.e., edges) is closest to a target fraction of the maximum possible edges.
#' Returns both the selected matrix and its index in the input list.
#'
#' @param adj_list A list of square matrices (numeric or logical) of the same dimension.
#' @param target_frac Numeric between 0 and 1; desired proportion of edges relative
#'        to the total possible edges \eqn{n(n-1)/2}.
#'
#' @return A list with two components:
#' \item{adj}{The selected adjacency matrix.}
#' \item{index}{The index (position) of that matrix in `adj_list`.}
#' @export
#'
#' @examples
#' \dontrun{
#' gC_result <- gC_NEW(amgut1.filt)
#' selected <- select_adj_by_edge_fraction(gC_result$path, target_frac = 0.02)
#' gold <- selected$adj
#' idx <- selected$index
#' }
select_adj_by_edge_fraction <- function(adj_list, target_frac) {
  if (!is.list(adj_list) || length(adj_list) == 0) {
    stop("adj_list must be a non-empty list")
  }
  n <- ncol(adj_list[[1]])
  max_edges <- n * (n - 1) / 2
  target <- target_frac * max_edges
  
  edge_counts <- vapply(adj_list, function(adj) {
    (sum(adj != 0, na.rm = TRUE) - sum(diag(adj) != 0, na.rm = TRUE)) / 2
  }, numeric(1))
  
  idx <- which.min(abs(edge_counts - target))
  list(adj = adj_list[[idx]], index = idx)
}

#' @examples
gcoda_amgut1 <- gC_NEW(amgut1.filt)
gold_gcoda_amgut1 <- select_adj_by_edge_fraction(gcoda_amgut1$path, target_frac = 0.02)
adj_gcoda_amgut1 <- gold_gcoda_amgut1$adj
index_adj_gcoda_amgut1 <- gold_gcoda_amgut1$index

