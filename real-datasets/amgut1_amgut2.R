##########################################################
####                load amgut1 & amgut2              ####
##########################################################
library(SpiecEasi)

data("amgut1.filt")

library(phyloseq)
data('amgut2.filt.phy')
amgut2.filt <- t(amgut2.filt.phy@otu_table@.Data)
saveRDS(amgut2.filt, "amgut2.filt.RData")

# or
amgut2.filt <- readRDS("amgut2.filt.RData")
