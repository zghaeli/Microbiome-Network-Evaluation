#########################################################
####    SpiecEasi-Based Synthetic Data Generation    ####
#########################################################

#' Generate Synthetic Data via SpiecEasi (SE) Method
#' Creates a list of synthetic datasets preserving topological structure
#'
#' @param data A numeric matrix of real-world dataset
#' @param ADJ Adjacency matrix of a topology
#'`@param num The length of the list
#' 
#' @return List of synthetic SE-based datasets

ListData_SE <- function(data, ADJ, num){
  lapply(1:num, function(i) { Synth_SE(data, ADJ, i) } )
}

#' Internal: Generate Single Synthetic Dataset
#' @noRd
Synth_SE <- function(data, ADJ, numberseed){
  depths <- rowSums(data)
  data.n  <- t(apply(data, 1, SpiecEasi::norm_to_total))
  data.cs <- round(data.n * min(depths))
  n <- nrow(data.cs)
  Cor <- Cor_AdjGraph(ADJ)
  set.seed(numberseed)
  X <- SpiecEasi::synth_comm_from_counts(data.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)
  return(X)
}

#' Convert Adjacency Matrix to Correlation Matrix
#' @noRd
Cor_AdjGraph <- function(ADJ){
  class(ADJ) <- "graph"
  Prec_ADJ  <- graph2prec(ADJ)
  Cor_ADJ <- cov2cor(prec2cov(Prec_ADJ))
  return(Cor_ADJ)
}




#' @example
# Create topology using the gCoda algorithm.
source(gcoda_functions.R)
gcoda_amgut1 <- gC_NEW(amgut1.filt, .85)
adj_gcoda_amgut1 <- gcoda_amgut1$gC_adj

# Create 100 synthetic datasets using the gCoda topology and SpiecEasi (SE) algorithm
amgut1_gcoda_SE_100 <- ListData_SE(amgut1.filt, adj_gcoda_amgut1, 100)
# Create 100 synthetic datasets using the gCoda topology and SPRING (SP) algorithm
amgut1_gcoda_SP_100 <- ListData_SP(amgut1.filt, adj_gcoda_amgut1, 100)

# Create a cluster topology using the SpiecEasi (SE) algorithm
graph_cluster_amgut1 <- SpiecEasi::make_graph('cluster', ncol(amgut1.filt), 175)
adj_cluster_amgut1 <- as.matrix(graph_cluster_amgut1)

# Create 100 synthetic datasets using the cluster topology and SpiecEasi (SE) algorithm
amgut1_cluster_SE_100 <- ListData_SE(amgut1.filt, adj_cluster_amgut1, 100)
amgut1_cluster_SP_100 <- ListData_SP(amgut1.filt, adj_cluster_amgut1, 100)
