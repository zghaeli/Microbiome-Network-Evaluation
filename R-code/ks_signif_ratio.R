#' Calculate KS Test Significance Ratio
#'
#' Computes the ratio of columns in synthetic data where the Kolmogorov-Smirnov
#' test p-value is below a specified threshold compared to real data.
#'
#' @param synth_data Numeric matrix of synthetic data.
#' @param real_data Numeric matrix of real data with same dimensions as synth_data.
#' @param alpha Significance level (default: 0.05).
#' @return Numeric ratio of significant p-values, rounded to 4 decimal places.
#' 
ks_signif_ratio <- function(synth_data, real_data, alpha = 0.05) {
  p_vals <- vapply(seq_len(ncol(synth_data)), 
                   function(i) { ks.test(synth_data[, i], real_data[, i])$p.value }, 
                   numeric(1))
  round(mean(p_vals < alpha, na.rm = TRUE), digits = 4)
}

#' @examples
amgut1_bootstrap_ks <- lapply(amgut1_bootstrap,  ks_signif_ratio, amgut1.filt, 0.05)

amgut1_noise05_ks <- lapply(amgut1_noise05_100, ks_signif_ratio, amgut1.filt, 0.05)
amgut1_noise20_ks <- lapply(amgut1_noise20_100, ks_signif_ratio, amgut1.filt, 0.05)

amgut1_gcoda_SE_ks <- lapply(amgut1_gcoda_SE_100, ks_signif_ratio, amgut1.filt, 0.05)
amgut1_gcoda_SP_ks <- lapply(amgut1_gcoda_SP_100, ks_signif_ratio, amgut1.filt, 0.05)

amgut1_cluster_SE_ks <- lapply(amgut1_cluster_SE_100, ks_signif_ratio, amgut1.filt, 0.05)
amgut1_cluster_SP_ks <- lapply(amgut1_cluster_SP_100, ks_signif_ratio, amgut1.filt, 0.05)

# Plot
# ##*************       KS.TEST        ****************
# par(mar = c(2.5, 2, 2.5, 2), mfrow = c(2,2), oma = c(6, 19, 5, 28))

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

# dev.off()
