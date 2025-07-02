# Microbiome-Network-Evaluation
A Pipeline for Evaluating New Algorithms for Constructing Microbial Co-occurrence Networks




## Description

In this work evaluating microbial network inference algorithms—gCoda, OIPCQ, S-E(glasso), S-E(mb), SPRING, and SparCC—under varied conditions, we generated three types of datasets derived from real data: synthetic, noisy, and bootstrap datasets. To assess how closely these datasets resemble real microbiome data, we conducted three types of comparisons: diversity analysis, matrix entry similarity, and distributional alignment. To evaluate their performance, we compared the inferred networks to reference networks using the F-score metric. We adjusted the parameters of each method so that the number of edges in the generated networks fell within a comparable range for a fair comparison.

## Table of Contents
- [Algorithms used for evaluation comparison](#algorithms-used-for-comparing)
- [Data generation methods](#data-generation-methods)
    - [Bootstrap](#bootstrap)
    - [Noisy](#noisy)
    - [Synthetic](#synthetic)
- [Comparing generated data with real data](#comparing-generated-data-with-real-data)
  - [Diversity indices](#diversity-indices)
  - [Matrix entry similarity](#matrix-entry-similarity)
  - [Distributional similarity (KS test)](#distributional-similarity-ks-test)
- [performance on constructing microbiome networks based on generated data](#performance-on-constructing-microbiome-networks-based-on-generated-data)
## Algorithms used for evaluation comparison

+ gCoda ([R code on GitHub](https://github.com/huayingfang/gCoda))
+ OIPCQ ([CMIMN package](https://github.com/rosaaghdam/CMiNet/tree/main))
+ S-E(glasso) ([SpiecEasi package](https://github.com/zdk123/SpiecEasi))
+ S-E(mb) ([SpiecEasi package](https://github.com/zdk123/SpiecEasi))
+ SPRING ([SPRING package](https://github.com/GraceYoon/SPRING))
+ SparCC ([SpiecEasi package](https://github.com/zdk123/SpiecEasi))


## Data generation methods

#### loading the Data
We use the American Gut data from [SpiecEasi package](https://github.com/zdk123/SpiecEasi).
```R
library(SpiecEasi)
data(amgut1.filt)
```

##### Bootstrap
```R
library(dplyr)
library(SpiecEasi)
source(List_Bootstrap_Sample_Data.R)
data(amgut1.filt)

# Create bootstrap datasets for amgut1.filt
amgut1_bootstrap <- List_Bootstrap_Sample(amgut1.filt, nrow(amgut1.filt), nrow(amgut1.filt), TRUE)
```

#### Noisy
```R
library(rtruncnorm)
library(SpiecEasi)
source(List_NoisyData.R)
data(amgut1.filt)

# Create 100 noisy datasets with 5% and 20% noise levels
amgut1_noise05_100 <- List_TRUNC_NoisyData(amgut1.filt, 0.05, 100)
amgut1_noise20_100 <- List_TRUNC_NoisyData(amgut1.filt, 0.2, 100)
```

#### Synthetic
```R
library(SpiecEasi)
library(SPRING)
source(gcoda.R)
source(gcoda_functions.R)
source(List_Synth_SE_Data.R)
source(List_Synth_SP_Data.R)

data(amgut1.filt)

# Create topology using the gCoda algorithm.
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

```

## Comparing generated data with real data

### Diversity indices
```R
library(vegan)
source(richness_functions.R)

amgut1_bootstrap_richness_pvalues <- compare_richness(amgut1.filt, amgut1_bootstrap)

amgut1_noise05_richness_pvalues <- compare_richness(amgut1.filt, amgut1_noise05_100)
amgut1_noise20_richness_pvalues <- compare_richness(amgut1.filt, amgut1_noise20_100)

amgut1_gcoda_SE_richness_pvalues <- compare_richness(amgut1.filt, amgut1_gcoda_SE_100)
amgut1_gcoda_SP_richness_pvalues <- compare_richness(amgut1.filt, amgut1_gcoda_SP_100)

amgut1_cluster_SE_richness_pvalues <- compare_richness(amgut1.filt, amgut1_cluster_SE_100)
amgut1_cluster_SP_richness_pvalues <- compare_richness(amgut1.filt, amgut1_cluster_SP_100)
```
```R
# Plots
par(mar = c(2.5, 2, 2.5, 2), mfrow = c(2,2), oma = c(6, 19, 5, 28))

vioplot::vioplot(as.numeric(amgut1_cluster_SE_richness_pvalues),
                 as.numeric(amgut1_gcoda_SE_richness_pvalues),
                 col = c(rep("deepskyblue1", 2)),
                 cex.axis = .8,
                 names = c("cluster", "gCoda"))
mtext(expression(bold("P-value")),  side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Topology")), side = 1, line = 2.5, cex = .9)
mtext(expression(bold("SE")), side = 3, line = .1, cex = .7, at = .5)
#mtext(expression(bold("A")),        side = 3, line = .6,  cex = 2, at = -.5)

boxplot(as.numeric(amgut1_noise05_richness_pvalues), 
        as.numeric(amgut1_noise20_richness_pvalues),
        col = c("tan", "salmon1"),
        border = c("tan", "salmon1"),
        cex.axis = .8,
        names = c("5%", "20%"))
mtext(expression(bold("P-value")), side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Noise")), side = 1, line = 2.2, cex = .9)
mtext(expression(bold("Noise")), side = 3, line = .1, cex = .7, at = .57)
#mtext(expression(bold("C")),     side = 3,    line = .6, cex = 2, at = -.5)

vioplot::vioplot(as.numeric(amgut1_cluster_SP_richness_pvalues),
                 as.numeric(amgut1_gcoda_SP_richness_pvalues),
                 col = c(rep("mediumvioletred", 2)),
                 cex.axis = .8,
                 names = c("cluster", "gCoda"))
mtext(expression(bold("P-value")), side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Topology")), side = 1, line = 2.5, cex = .9)
mtext(expression(bold("SP")), side = 3, line = .1, cex = .7, at = .5)
#mtext(expression(bold("A")),        side = 3, line = .6,  cex = 2, at = -.5)

vioplot::vioplot(as.numeric(amgut1_bootstrap_richness_pvalues), cex.axis = .8, col = "mediumseagreen")
mtext(expression(bold("P-value")), side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Bootstrap")), side = 1, line = 2.5, cex = .9)
mtext(expression(bold("Bootstrap")), side = 3, line = .1, cex = .7, at = .57)
#mtext(expression(bold("D")),         side = 3, line = .6, cex = 2, at = .32)

mtext(expression(bold("Richness__amgut1")), side = 3, line = -.7, outer = TRUE, cex = 1.2)
```

![richness](/image/richness.png)

### Matrix entry similarity

```R
source(Entry_Fscore.R)

# Bootstrap
amgut1_bootstrap_fscore <- lapply(amgut1_bootstrap, calculate_fscore, amgut1.filt)

# Noisy
amgut1_noise05_fscore <- lapply(amgut1_noise05_100, calculate_fscore, amgut1.filt)
amgut1_noise20_fscore <- lapply(amgut1_noise20_100, calculate_fscore, amgut1.filt)

# Synthetic
amgut1_gcoda_SE_fscore <- lapply(amgut1_gcoda_SE_100, calculate_fscore, amgut1.filt)
amgut1_gcoda_SP_fscore <- lapply(amgut1_gcoda_SP_100, calculate_fscore, amgut1.filt)
amgut1_cluster_SE_fscore <- lapply(amgut1_cluster_SE_100, calculate_fscore, amgut1.filt)
amgut1_cluster_SP_fscore <- lapply(amgut1_cluster_SP_100, calculate_fscore, amgut1.filt)
```
```R
# Plots
boxplot(as.numeric(amgut1_cluster_SE_fscore),
        as.numeric(amgut1_gcoda_SE_fscore),
        as.numeric(amgut1_cluster_SP_fscore),
        as.numeric(amgut1_gcoda_SP_fscore),
        as.numeric(amgut1_noise05_fscore),
        as.numeric(amgut1_noise20_fscore),
        as.numeric(amgut1_bootstrap_fscore),
        names = c("cluster", "gCoda", "cluster", "gCoda", "5%", "20%", "Bootstrap"),
        col = c("red3", "red3", "blue3", "blue3", "green3", "green3", "darkgreen"),
        border = c("red3", "red3", "blue3", "blue3", "green3", "green3", "darkgreen"),
        ylim = c(0,1),
        cex = 1.2)
abline(v = 2.5, col = "gray", lty = 1, lwd = 1.5)
abline(v = 4.5, col = "gray", lty = 1, lwd = 1.5)
abline(v = 6.5, col = "gray", lty = 1, lwd = 1.5)
mtext(expression(bold("SE")), side = 3, at = .4, line = .1, cex = .8)
mtext(expression(bold("SP")), side = 3, at = 2.7, line = .1, cex = .8)
mtext(expression(bold("Noise")), side = 3, at = 4.8, line = .1, cex = .8)
mtext(expression(bold("Bootstrap")), side = 3, at = 7, line = .1, cex = .8)
mtext(expression(bold("F-score")), side = 2, line = 2.5, cex = 1.3, las = 0)
mtext(expression(bold("Matrix entry similarity")), side = 3, line = 1.7, cex = 1.4, las = 0)

```
![Matrix Entry Similarity](/image/entry_fscore.png)

### Distributional similarity (KS test)

```R
source(ks_signif_ratio.R)

amgut1_bootstrap_ks <- lapply(amgut1_bootstrap,  ks_signif_ratio, amgut1.filt, 0.05)

amgut1_noise05_ks <- lapply(amgut1_noise05_100, ks_signif_ratio, amgut1.filt, 0.05)
amgut1_noise20_ks <- lapply(amgut1_noise20_100, ks_signif_ratio, amgut1.filt, 0.05)

amgut1_gcoda_SE_ks <- lapply(amgut1_gcoda_SE_100, ks_signif_ratio, amgut1.filt, 0.05)
amgut1_gcoda_SP_ks <- lapply(amgut1_gcoda_SP_100, ks_signif_ratio, amgut1.filt, 0.05)

amgut1_cluster_SE_ks <- lapply(amgut1_cluster_SE_100, ks_signif_ratio, amgut1.filt, 0.05)
amgut1_cluster_SP_ks <- lapply(amgut1_cluster_SP_100, ks_signif_ratio, amgut1.filt, 0.05)
```
```R
# Plots
par(mar = c(2.5, 2, 2.5, 2), mfrow = c(2,2), oma = c(6, 19, 5, 28))

vioplot::vioplot(as.numeric(amgut1_cluster_SE_ks),
                 as.numeric(amgut1_gcoda_SE_ks),
                 col = c(rep("deepskyblue1", 2)),
                 cex.axis = .8,
                 ylim = c(.8,1),
                 names = c("cluster", "gCoda"))
mtext(expression(bold("% OTUs")),  side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Topology")), side = 1, line = 2.5, cex = .9)
mtext(expression(bold("SE")), side = 3, line = .1, cex = .7, at = .5)
#mtext(expression(bold("A")),        side = 3, line = .6,  cex = 2, at = -.5)

vioplot::vioplot(as.numeric(amgut1_noise05_ks), 
                 as.numeric(amgut1_noise20_ks),
                 col = c("tan", "salmon1"),
                 border = c("tan", "salmon1"),
                 cex.axis = .8,
                 ylim = c(0,1),
                 names = c("5%", "20%"))
mtext(expression(bold("% OTUs")), side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Noise")), side = 1, line = 2.2, cex = .9)
mtext(expression(bold("Noise")), side = 3, line = .1, cex = .7, at = .57)
#mtext(expression(bold("C")),     side = 3,    line = .6, cex = 2, at = -.5)


vioplot::vioplot(as.numeric(amgut1_cluster_SP_ks),
                 as.numeric(amgut1_gcoda_SP_ks),
                 col = c(rep("mediumvioletred", 2)),
                 cex.axis = .8,
                 names = c("cluster", "gCoda"))
mtext(expression(bold("% OTUs")), side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Topology")), side = 1, line = 2.5, cex = .9)
mtext(expression(bold("SP")), side = 3, line = .1, cex = .7, at = .5)
#mtext(expression(bold("A")),        side = 3, line = .6,  cex = 2, at = -.5)

vioplot::vioplot(as.numeric(amgut1_bootstrap_ks), cex.axis = .8, col = "mediumseagreen")
mtext(expression(bold("% OTUs")), side = 2, line = 2, cex = .9, las = 0)
mtext(expression(bold("Bootstrap")), side = 1, line = 2.5, cex = .9)
mtext(expression(bold("Bootstrap")), side = 3, line = .1, cex = .7, at = .57)
#mtext(expression(bold("D")),         side = 3, line = .6, cex = 2, at = .32)

mtext(expression(bold("KS__amgut1")), side = 3, line = -.7, outer = TRUE, cex = 1.2)
```
![ks](/image/ks.png)


## performance on constructing microbiome networks based on generated data

```R
source("gcoda.R")
source("gcoda_functions.R")
source("My_Diagnostic.R")

amgut1_bootstrap_gcoda <- process_gcoda(amgut1.filt, amgut1_bootstrap, adj_gcoda_amgut1, threshold = 0.85, save_dir = ".", timed = TRUE)

amgut1_noise05_gcoda <- process_gcoda(amgut1.filt, amgut1_noise05_100, adj_gcoda_amgut1, threshold = 0.85, save_dir = ".", timed = TRUE)
amgut1_noise20_gcoda <- process_gcoda(amgut1.filt, amgut1_noise20_100, adj_gcoda_amgut1, threshold = 0.85, save_dir = ".", timed = TRUE)

amgut1_gcoda_SE_gcoda <- process_gcoda(amgut1.filt, amgut1_gcoda_SE_100, adj_gcoda_amgut1, threshold = 0.85, save_dir = ".", timed = TRUE)
amgut1_gcoda_SP_gcoda <- process_gcoda(amgut1.filt, amgut1_gcoda_SP_100, adj_gcoda_amgut1, threshold = 0.85, save_dir = ".", timed = TRUE)

amgut1_cluster_SE_gcoda <- process_gcoda(amgut1.filt, amgut1_cluster_SE_100[1:5], adj_cluster_amgut1, threshold = 0.85, save_dir = ".", timed = TRUE)
amgut1_cluster_SP_gcoda <- process_gcoda(amgut1.filt, amgut1_cluster_SP_100, adj_cluster_amgut1, threshold = 0.85, save_dir = ".", timed = TRUE)
```
```R
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
```
![gcoda fscores](/image/gcoda_fscores.png)