#--- init ----------------------------------------------------------------------

# read epigenetic age predictions from the supplement
dat <- read.csv("input/WBCjohansson_clocks.csv", skip = 1)

# read cell count data from the supplement
cells <- read.csv("input/WBCjohansson_cellcounts_Reinius.csv", skip = 1)

# check if samples are matched
stopifnot(all.equal(dat$id, cells$id))


#--- obtain ages ---------------------------------------------------------------

# replace id age
dat$id <- as.numeric(sub(".*_", "", dat$id))
colnames(dat)[1] <- "age"

cells$id <- as.numeric(sub(".*_", "", cells$id))
colnames(cells)[1] <- "age"


#--- plot ----------------------------------------------------------------------

plotcor <- function(x, y, col) {
  cor   <- cor.test(x,y)
  issig <- cor$p.value <= 0.05

  plot.new()
  plot.window(xlim = range(x), ylim = extendrange(y, f = c(0,0.25)))
  points(x, y, col = col, pch = 19, cex = 0.25)
  abline(lm(y ~ x), col = col, lty = ifelse(issig, 1, 2), lwd = 1.5)
  axis(2, lwd = 0)
  axis(1, lwd = 0)
  box()
  legend("topright", legend = c(paste0("cor = ", signif(cor$estimate, 2)), paste0("p    = ", signif(cor$p.value, 2))), bty = 'n', inset = c(0, -0.025))
}

pdf("fig_2_efghij.pdf", width = 18/2.54, height = 4/2.54, pointsize = 8)

par(mfrow=c(1,6), mar=c(4,0,2,2), oma=c(0,6,2,1), mgp=c(2,0.25,0), las = 1)

plotcor(cells[,"nk"] * 100,    dat$Horvath.pan.tissue.2013 - dat$age, "#F15854")
title(ylab = expression(paste(Delta, " (predicted age - real age)")), line = 2.5, xpd = NA)
title(xlab = "NK-cell %")
title("e", adj = 0, cex.main = 1.5, line = 1.5, xpd = NA)

plotcor(cells[,"bcell"] * 100, dat$Horvath.pan.tissue.2013 - dat$age, "limegreen")
title(xlab = "B-cell %")
title("f", adj = 0, cex.main = 1.5, line = 1.5, xpd = NA)

plotcor(cells[,"tcd4"] * 100,  dat$Horvath.pan.tissue.2013 - dat$age, "deepskyblue")
title(xlab = "CD4+ T-cell %")
title("g", adj = 0, cex.main = 1.5, line = 1.5, xpd = NA)

plotcor(cells[,"tcd8"] * 100,  dat$Horvath.pan.tissue.2013 - dat$age, "dodgerblue")
title(xlab = "CD8+ T-cell %")
title("h", adj = 0, cex.main = 1.5, line = 1.5, xpd = NA)

plotcor(cells[,"mono"] * 100,  dat$Horvath.pan.tissue.2013 - dat$age, "#DECF3F")
title(xlab = "Monocyte %")
title("i", adj = 0, cex.main = 1.5, line = 1.5, xpd = NA)

plotcor(cells[,"gran"] * 100,  dat$Horvath.pan.tissue.2013 - dat$age, "#FAA43A")
title(xlab = "Granulocyte %")
title("j", adj = 0, cex.main = 1.5, line = 1.5, xpd = NA)

invisible(dev.off())
