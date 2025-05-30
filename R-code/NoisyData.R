#############################################################################
####         Noisy data from truncated normal distribution               ####
#############################################################################
library(truncnorm)

TOTAL_TRUNC_NoisyData <- function(Data, percent, num){
  INC <- Index__NN__Cols(Data)
  L_NoisyData <- list(1:num)
  for (s in 1:num) {
    L_per <- PERCENT(INC, percent)
    L_NoisyData[[s]] <- TRUNC_NoisyData(Data, L_per)
  }
  return(L_NoisyData)
}

TRUNC_NoisyData <- function(Data, SubsetList){
  noisydata <- Data
  for (i in 1:dim(Data)[2]) {
    rtruncnorm_noise <- rtruncnorm(length(SubsetList[[i]]), a = min(Data[,i]), b = max(Data[,i]), mean = mean(Data[,i]), sd = sd(Data[,i]))
    ll <- noisydata[,i][SubsetList[[i]]] + rtruncnorm_noise
    t <- 1
    for (j in SubsetList[[i]]) {
      noisydata[,i][j] <- ll[t]
      t <- t + 1
    }
  }
  return(noisydata)
}

PERCENT <- function(List, perc){
  L_percent <- list(1:length(List))
  for (i in 1:length(List)) {
    L_percent[[i]] <- sample(List[[i]], perc*length(List[[i]]))
  }
  return(L_percent)
}

Index__NN__Cols <- function(Data){
  index__nn__zero <- list(1:ncol(Data))
  for (i in 1:ncol(Data)) {
    index__nn__zero[[i]] <- which(Data[,i] != 0)
  }
  return(index__nn__zero)
}

