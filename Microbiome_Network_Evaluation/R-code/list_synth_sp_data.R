####################################################################
####        SPRING-Based Synthetic Data Generation             ####
####################################################################

#' Generate a List of Synthetic Datasets Using SPRING (SP)
#'
#' Creates multiple synthetic datasets that preserve the topological structure
#' of a given adjacency matrix. The function uses SPRING's synthetic data generator
#' based on the empirical cumulative distribution function (ECDF) of the original data.
#'
#' @param data A numeric matrix or data frame of real count data (rows = samples, columns = OTUs).
#' @param ADJ An adjacency matrix (binary or weighted) representing the underlying network topology.
#' @param num Integer. Number of synthetic datasets to generate.
#' @param seed Integer seed for reproducibility (optional). If provided, each replicate
#'        uses a different seed (`seed + i`) to ensure independent datasets.
#'
#' @return A list of length \code{num}, each element being a synthetic count matrix
#'         with the same dimensions as \code{data}.
#' @export
#'
#' @importFrom SPRING synthData_from_ecdf
#' @importFrom SpiecEasi graph2prec prec2cov
#'
#' @examples
#' \dontrun{
#' # Assuming 'amgut1.filt' is a count matrix and 'adj_gcoda_amgut1' is an adjacency matrix
#' sp_synth_list <- List_SP_Data(amgut1.filt, adj_gcoda_amgut1, num = 10, seed = 123)
#' }
#'
List_SP_Data <- function(data, ADJ, num, seed = NULL) {
  if (!is.matrix(data) && !is.data.frame(data)) stop("data must be a matrix or data frame")
  if (!is.matrix(ADJ)) stop("ADJ must be a matrix (adjacency)")
  if (num <= 0) stop("num must be a positive integer")
  
  # Precompute correlation matrix from adjacency (once)
  Cor <- Cor_AdjGraph(ADJ)
  n <- nrow(data)
  
  # Generate replicates
  result <- lapply(seq_len(num), function(i) {
    if (!is.null(seed)) set.seed(seed + i)
    X <- SPRING::synthData_from_ecdf(data, Sigma = Cor, n = n, seed = NULL)  # seed handled by set.seed
    return(X)
  })
  
  return(result)
}


#' Convert Adjacency Matrix to Correlation Matrix via Precision Matrix
#'
#' Internal function that transforms an adjacency matrix (graph) into a
#' correlation matrix suitable for SPRING's synthetic data generator.
#' Uses SpiecEasi's graph2prec and prec2cov utilities.
#'
#' @param ADJ Adjacency matrix (binary or weighted).
#' @return A correlation matrix (positive definite).
#' @noRd
Cor_AdjGraph <- function(ADJ) {
  class(ADJ) <- "graph"
  Prec_ADJ <- SpiecEasi::graph2prec(ADJ)
  Cor_ADJ <- cov2cor(SpiecEasi::prec2cov(Prec_ADJ))
  return(Cor_ADJ)
}


# For backward compatibility (keep old name, but mark as deprecated)
#' @rdname List_SP_Data
#' @export
ListData_SP <- function(data, ADJ, num) {
  .Deprecated("List_SP_Data")
  List_SP_Data(data, ADJ, num)
}


