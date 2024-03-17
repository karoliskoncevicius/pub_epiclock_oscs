#--- init ----------------------------------------------------------------------

# read epigenetic age predictions from the supplement
dat <- read.csv("input/PBMCtimeofday_clocks_IEAA.csv", skip = 1)


#--- obtain differences between 12:45 and 16:15 --------------------------------

# extract time and donor information from ids
dat$time  <- sub("_.*", "", dat$id)
dat$donor <- sub(".*_", "", dat$id)

# make sure pairs of donors are matched
stopifnot(all.equal(dat$donor[dat$time == "12:45"], dat$donor[dat$time == "16:15"]))


#--- plot ----------------------------------------------------------------------

plotclock <- function(x, time, title) {
  pval <- t.test(x ~ time, paired = TRUE)$p.value

  boxplot(x ~ time, axes = FALSE, ylim = extendrange(pretty(x), f = c(0,0.05)), medlty = ifelse(pval <= 0.05, 1, 2), medlwd = 2, col = NA, xlab = "", ylab = "")

  legend("topleft", legend = title, bty = "n", inset = c(-0.03,0))
  legend("topright", legend = paste0("p = ", signif(pval, 2)), bty="n")

  axis(1, at = 1:length(unique(time)), lwd = 0, labels = unique(time))
  axis(2, las = 2)
  box()

  title(xlab = "time of day")
  title(ylab = "prediction", line = 3)
}


pdf("fig_s9_abc.pdf", width = 18/2.54, height = 5/2.54, pointsize = 8)

par(mar = c(4,4,4,4), oma = c(0,2,0,0), mgp = c(2.5,0.5,0), cex.main = 1)
par(mfrow = c(1,3))

plotclock(dat$Horvath.pan.tissue.2013, dat$time, "Horvath Pan-tissue 2013")
title("e", adj = 0, cex.main = 1.5)

plotclock(dat$Teschendorff.epiTOC2.2020, dat$time, "Teschendorff epiTOC2 2020")
title("f", adj = 0, cex.main = 1.5)

plotclock(dat$Lu.GrimAge2.2022, dat$time, "Lu GrimAge2 2022")
title("g", adj = 0, cex.main = 1.5)

invisible(dev.off())
