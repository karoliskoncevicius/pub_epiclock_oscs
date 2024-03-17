#--- init ----------------------------------------------------------------------

# read epigenetic age predictions from the supplement
dat <- read.csv("input/Neu30yo_clocks.csv", skip = 1)


#--- obtain oscillation results ------------------------------------------------

# replace id with time
dat$id <- sub(".*?_", "", dat$id)
dat$id <- as.numeric(gsub("_R1|_R2", "", dat$id))
colnames(dat)[1] <- "time"

# combine technical replicates
dat <- apply(dat, 2, tapply, dat$time, mean)


#--- plot ----------------------------------------------------------------------

# function to calculate oscillation and plot the results
plotclock <- function(x, zt, title) {
  zt <- zt %% 24

  osc  <- lm(x ~ sin(2*pi/24*zt) + cos(2*pi/24*zt))
  null <- lm(x ~ 1)
  pval <- anova(osc, null)[2,6]

  plot.new()
  plot.window(xlim = c(0,24), ylim = extendrange(x, f = c(0, 0.175)))

  points(zt, x)

  xs <- seq(-1, 25, 0.1)
  lines(xs, coef(osc)[1] + coef(osc)[2] * sin(2*pi/24*xs) + coef(osc)[3] * cos(2*pi/24*xs), lty = ifelse(pval <= 0.05, 1, 2))
  legend("topleft", legend = title, bty = "n", inset = c(-0.03,0))
  legend("topright", legend = paste0("p = ", signif(pval, 2)), bty="n")

  axis(1, at = seq(0,24,4), lwd = 0)
  axis(2, lwd = 0, las = 2)
  box()

  title(xlab = "time of day")
  title(ylab = "prediction", line = 3)
}


pdf("fig_s12_abc.pdf", width = 18/2.54, height = 5/2.54, pointsize = 8)

par(mar = c(4,4,4,4), oma = c(0,2,0,0), mgp = c(2.5,0.5,0), cex.main = 1)
par(mfrow = c(1,3))

plotclock(dat[,"Horvath.pan.tissue.2013"], dat[,1], "Horvath Pan-tissue 2013")
title("a", adj = 0, cex.main = 1.5)

plotclock(dat[,"Teschendorff.epiTOC2.2020"], dat[,1], "Teschendorff epiTOC2 2020")
title("b", adj = 0, cex.main = 1.5)

plotclock(dat[,"Lu.GrimAge2.2022"], dat[,1], "Lu GrimAge2 2022")
title("c", adj = 0, cex.main = 1.5)

invisible(dev.off())
