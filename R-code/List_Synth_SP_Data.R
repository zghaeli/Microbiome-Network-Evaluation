#########################################################
####     SPRING-Based Synthetic Data Generation      ####
#########################################################

#' Generate Synthetic Data via SPRING (SP) Method
#' Creates a list of synthetic datasets preserving topological structure
#
#' @param data A numeric matrix of real-world dataset
#' @param ADJ Adjacency matrix of a topology
#'`@param num The length of the list
#' 
#' @return List of synthetic SP-based datasets

ListData_SP <- function(data, ADJ, num){
  lapply(1:num, function(i) {Synth_SP(data, ADJ, i)})
}

#' Internal: Generate Single Synthetic Dataset
#' @noRd
Synth_SP <- function(data, ADJ, numberseed){
  n <- nrow(data)
  Cor1 <- Cor_AdjGraph(ADJ)
  X1_count <- SPRING::synthData_from_ecdf(data, Sigma = Cor1, n = n, seed = numberseed)
  return(X1_count)
}

#' Convert Adjacency Matrix to Correlation Matrix
#' @noRd
Cor_AdjGraph <- function(ADJ){
  class(ADJ) <- "graph"
  Prec_ADJ  <- SpiecEasi::graph2prec(ADJ)
  Cor_ADJ <- cov2cor(SpiecEasi::prec2cov(Prec_ADJ))
  return(Cor_ADJ)
}


