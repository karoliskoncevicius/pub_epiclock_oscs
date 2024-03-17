#--- init ----------------------------------------------------------------------

# read epigenetic age predictions from the supplement
dat1 <- read.csv("input/Neu52yo_clocks.csv", skip = 1)
dat2 <- read.csv("input/Neu54yo_clocks.csv", skip = 1)
dat3 <- read.csv("input/Neu30yo_clocks.csv", skip = 1)


#--- obtain oscillation results ------------------------------------------------

# replace id with time
dat1$id <- as.numeric(gsub(".*_CT|_R1", "", dat1$id))
colnames(dat1)[1] <- "time"

dat2$id <- sub(".*?_", "", dat2$id)
dat2$id <- as.numeric(gsub("_R1|_R2", "", dat2$id))
colnames(dat2)[1] <- "time"

dat3$id <- sub(".*?_", "", dat3$id)
dat3$id <- as.numeric(gsub("_R1|_R2", "", dat3$id))
colnames(dat3)[1] <- "time"

# combine technical replicates
dat1 <- as.data.frame(apply(dat1, 2, tapply, dat1$time, mean))
dat2 <- as.data.frame(apply(dat2, 2, tapply, dat2$time, mean))
dat3 <- as.data.frame(apply(dat3, 2, tapply, dat3$time, mean))

# combine datasets
dat1 <- data.frame(donor = "A", time = dat1$time, dat1[,-1])
dat2 <- data.frame(donor = "B", time = dat2$time, dat2[,-1])
dat3 <- data.frame(donor = "C", time = dat3$time, dat3[,-1])

dat <- rbind(dat1, dat2, dat3)


#--- plot ----------------------------------------------------------------------

# function to calculate oscillation and plot the results
plotclock <- function(x, donor, zt, title) {
  # NOTE: here we allow each donor to have separate oscillation curve
  #       this is achieved by adding an interaction term between sin/cos and donor
  osc  <- lm(x ~ sin(2*pi/24*zt) * donor + cos(2*pi/24*zt) * donor)
  null <- lm(x ~ donor)
  pval <- anova(osc, null)[2,6]

  zt <- zt %% 24

  plot.new()
  plot.window(xlim = c(0,24), ylim = extendrange(x, f = c(0, 0.175)))

  points(zt, x, col=basetheme::lab2col(donor, pal=c("cornflowerblue", "blueviolet", "magenta2"), ref=c("A", "B", "C")))

  xs <- seq(-1, 25, 0.1)
  fitA <- lm(x[donor == "A"] ~ sin(2*pi/24*zt)[donor == "A"] + cos(2*pi/24*zt)[donor == "A"])
  lines(xs, coef(fitA)[1] + coef(fitA)[2] * sin(2*pi/24*xs) + coef(fitA)[3] * cos(2*pi/24*xs), col="cornflowerblue", lty=ifelse(pval <= 0.05, 1, 2))
  fitB <- lm(x[donor == "B"] ~ sin(2*pi/24*zt)[donor == "B"] + cos(2*pi/24*zt)[donor == "B"])
  lines(xs, coef(fitB)[1] + coef(fitB)[2] * sin(2*pi/24*xs) + coef(fitB)[3] * cos(2*pi/24*xs), col="blueviolet", lty=ifelse(pval <= 0.05, 1, 2))
  fitC <- lm(x[donor == "C"] ~ sin(2*pi/24*zt)[donor == "C"] + cos(2*pi/24*zt)[donor == "C"])
  lines(xs, coef(fitC)[1] + coef(fitC)[2] * sin(2*pi/24*xs) + coef(fitC)[3] * cos(2*pi/24*xs), col="magenta2", lty=ifelse(pval <= 0.05, 1, 2))

  legend("topleft", legend = title, bty = "n", inset = c(-0.03,0))
  legend("topright", legend = paste0("p = ", signif(pval, 2)), bty="n")

  axis(1, at = seq(0,24,4), lwd = 0)
  axis(2, lwd = 0, las = 2)
  box()

  title(xlab = "time of day")
  title(ylab = "prediction", line = 3)
}


pdf("fig_2_klm.pdf", width = 18/2.54, height = 5/2.54, pointsize = 8)

par(mar = c(4,4,4,4), oma = c(0,2,0,0), mgp = c(2.5,0.5,0), cex.main = 1)
par(mfrow = c(1,3))

plotclock(dat[,"Horvath.pan.tissue.2013"], dat[,1], dat[,2], "Horvath Pan-tissue 2013")
title("k", adj = 0, cex.main = 1.5)
legend("bottomright", c("30 yr", "52 yr", "54 yr"), col = c("magenta2", "cornflowerblue", "blueviolet"), lty = 1, bty = "n", cex = 0.9, inset = c(0, 1), horiz = TRUE, xpd = TRUE)

plotclock(dat[,"Teschendorff.epiTOC2.2020"], dat[,1], dat[,2], "Teschendorff epiTOC2 2020")
title("l", adj = 0, cex.main = 1.5)

plotclock(dat[,"Lu.GrimAge2.2022"], dat[,1], dat[,2], "Lu GrimAge2 2022")
title("m", adj = 0, cex.main = 1.5)

invisible(dev.off())
