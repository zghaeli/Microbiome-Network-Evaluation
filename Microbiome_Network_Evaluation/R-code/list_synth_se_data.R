####################################################################
####         SpiecEasi-Based Synthetic Data Generation          ####
####################################################################

#' Generate a List of Synthetic Datasets Using SpiecEasi (SE)
#'
#' Creates multiple synthetic datasets that preserve the topological structure
#' of a given adjacency matrix. The function uses the SpiecEasi framework
#' (synthetic community from counts) to generate new count data.
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
#' @importFrom SpiecEasi norm_to_total synth_comm_from_counts
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom huge huge.generator
#'
#'
List_SE_Data <- function(data, ADJ, num, seed = NULL) {
  if (!is.matrix(data) && !is.data.frame(data)) stop("data must be a matrix or data frame")
  if (!is.matrix(ADJ)) stop("ADJ must be a matrix (adjacency)")
  if (num <= 0) stop("num must be a positive integer")
  
  # Precompute once: depths, normalized data, correlation matrix from ADJ
  depths <- rowSums(data)
  data_norm <- t(apply(data, 1, SpiecEasi::norm_to_total))
  data_cs <- round(data_norm * min(depths))
  Cor <- Cor_AdjGraph(ADJ)
  
  # Generate replicates
  result <- lapply(seq_len(num), function(i) {
    if (!is.null(seed)) set.seed(seed + i)
    X <- SpiecEasi::synth_comm_from_counts(data_cs, mar = 2, distr = 'zinegbin',
                                           Sigma = Cor, n = nrow(data_cs))
    return(X)
  })
  
  return(result)
}


#' Convert Adjacency Matrix to Correlation Matrix via Precision Matrix
#'
#' Internal function that transforms an adjacency matrix (graph) into a
#' correlation matrix suitable for SpiecEasi's synthetic community generator.
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


# For backward compatibility (optional, keep old name if needed)
#' @rdname List_SE_Data
#' @export
ListData_SE <- function(data, ADJ, num) {
  .Deprecated("List_SE_Data")
  List_SE_Data(data, ADJ, num)
}


#####################################################################
####    Example: Generating Synthetic Datasets from Topologies    ###
#####################################################################

# Load necessary functions
source("gcoda_functions.R") 
source("list_synth_se_data.R")
source("list_synth_sp_data.R")
source("select_adj_by_edge_fraction.R")
library(SpiecEasi)

# 1. Create topology using the gCoda algorithm
gcoda_amgut1 <- gC_NEW(amgut1.filt)
gold_gcoda_amgut1 <- select_adj_by_edge_fraction(gcoda_amgut1$path, target_frac = 0.02)
adj_gcoda_amgut1 <- gold_gcoda_amgut1$adj
index_adj_gcoda_amgut1 <- gold_gcoda_amgut1$index


# 2. Generate synthetic datasets from gCoda topology
#    SpiecEasi (SE) method
amgut1_gcoda_SE_100 <- List_SE_Data(amgut1.filt, adj_gcoda_amgut1, num = 100, seed = 123)
#    SPRING (SP) method
amgut1_gcoda_SP_100 <- List_SP_Data(amgut1.filt, adj_gcoda_amgut1, num = 100, seed = 456)


# 3. Create a cluster topology using SpiecEasi
graph_cluster_amgut1 <- SpiecEasi::make_graph('cluster', ncol(amgut1.filt), 0.02 * (ncol(amgut1.filt)(ncol(amgut1.filt) -1)/2))
adj_cluster_amgut1 <- as.matrix(graph_cluster_amgut1)

# 4. Generate synthetic datasets from cluster topology
#    SpiecEasi (SE) method
amgut1_cluster_SE_100 <- List_SE_Data(amgut1.filt, adj_cluster_amgut1, num = 100, seed = 789)
#    SPRING (SP) method
amgut1_cluster_SP_100 <- List_SP_Data(amgut1.filt, adj_cluster_amgut1, num = 100, seed = 101)



