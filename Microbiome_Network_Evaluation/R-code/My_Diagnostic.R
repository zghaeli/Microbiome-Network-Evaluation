############################################################
####         Metrics for compairing two networks        ####
############################################################

MY_Diagnostic <- function(true_ADJ, predicted_ADJ, number_of_genes){
  
  TrueNet1 <- true_ADJ
  
  for ( i in 1:number_of_genes) {
    for ( j in 1:number_of_genes) {
      if (TrueNet1[i,j] > 1) { TrueNet1[i,j] = 1} else {TrueNet1[i,j]=TrueNet1[i,j]}
    }}
  
  #Adjacent matrix for undirected graph
  #the value in each cell : the number of edges between two nodes.
  
  AD <- predicted_ADJ
  
  for ( i in 1:number_of_genes) {
    for ( j in 1:number_of_genes) {
      if (AD[i,j] > 1) { AD[i,j] = 1} else {AD[i,j]=AD[i,j]}
    }}
  
  #####Diagnostic Measures #####
  
  matches <- matrix(0,number_of_genes,number_of_genes)
  for ( i in 1:number_of_genes) {
    for ( j in i:number_of_genes) {
      if (AD[i,j] == 1 && TrueNet1[i,j]==1) {matches [i,j]=1} else {
        if (AD[i,j] == 0 && TrueNet1[i,j]==0) {matches [i,j]=0} else {
          if (AD[i,j] == 1 && TrueNet1[i,j]==0) {matches [i,j]=2} else {
            matches [i,j]=3}}}
    }}
  
  matches_UP <-  matches[upper.tri(matches,diag=FALSE)]
  
  TP <- length (which(matches_UP == 1))
  TN <- length (which(matches_UP == 0))
  FP <- length (which(matches_UP == 2))
  FN <- length (which(matches_UP == 3))
  
  #Diagnostic_measures_undir <- cbind(TP,TN,FP,FN)
  
  ###################
  
  Recall      <- TP / (TP+FN)
  FPR         <- FP / (FP+TN)
  Precision   <- TP / (TP+FP)
  Specificity <- TN / (TN+FP)
  Accuracy    <- (TP+TN) / (TP+FP+TN+FN)
  F_score     <- 2*(Recall*Precision)/(Recall+Precision)
  
  Accuracy_measures_undir <- cbind(TP,TN,FP,FN,Recall,FPR,Specificity,Precision,Accuracy,F_score)
  
  ################################################################################
  ################################################################################
  
  TrueNet <- true_ADJ
  
  #Adjacent Matrix for directed Network
  
  AD.dir <- predicted_ADJ
  AD.dir[number_of_genes,number_of_genes] <- 0
  
  #####Diagnostic Measures #####
  
  matches.dir <- matrix(0,number_of_genes,number_of_genes)
  for ( i in 1:number_of_genes) {
    for ( j in 1:number_of_genes) {
      if (AD.dir[i,j] == 1 && TrueNet[i,j]==1) {matches.dir [i,j]=1} else {
        if (AD.dir[i,j] == 0 && TrueNet[i,j]==0) {matches.dir [i,j]=0} else {
          if (AD.dir[i,j] == 1 && TrueNet[i,j]==0) {matches.dir [i,j]=2} else {
            matches.dir [i,j]=3}}}
    }}
  
  TP.dir <- length (which(matches.dir == 1))
  TN.dir <- length (which(matches.dir == 0))-(number_of_genes)
  FP.dir <- length (which(matches.dir == 2))
  FN.dir <- length (which(matches.dir == 3))
  
  #Diagnostic_measures_dir <- cbind(TP.dir,TN.dir,FP.dir,FN.dir)
  
  ###################
  Recall.dir      <- TP.dir / (TP.dir+FN.dir)
  FPR.dir<- FP.dir / (FP.dir+TN.dir)
  Precision.dir   <- TP.dir / (TP.dir+FP.dir)
  Specificity.dir <- TN.dir / (TN.dir+FP.dir)
  Accuracy.dir    <- TP.dir / (TP.dir+FP.dir+TN.dir+FN.dir)
  F_score.dir     <- 2*(Recall.dir*Precision.dir)/(Recall.dir+Precision.dir)
  
  Accuracy_measures_dir <- cbind(TP.dir,TN.dir,FP.dir,FN.dir,Recall.dir,FPR.dir,Specificity.dir,Precision.dir,Accuracy.dir,F_score.dir)
  
  Diag_results <- list(Accuracy_measures_undir,Accuracy_measures_dir)
  
  return(Diag_results)
}
