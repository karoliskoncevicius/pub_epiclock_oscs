#--- init ----------------------------------------------------------------------

# read epigenetic age predictions from the supplement
dat1 <- read.csv("input/CELLSwang_clocks.csv", skip = 1)
dat2 <- read.csv("input/CELLSreinius_clocks.csv", skip = 1)


#--- subtract whole blood ------------------------------------------------------

# extract donor and celltype from ids
dat1$donor    <- sub(".*_", "", dat1$id)
dat1$celltype <- sub("_.*", "", dat1$id)

dat2$donor    <- sub(".*_", "", dat2$id)
dat2$celltype <- sub("_.*", "", dat2$id)

# reorder cell types
dat1$celltype <- factor(dat1$celltype, levels = c("NK-cells", "B-cells", "CD4+ T-cells", "CD8+ T-cells", "Monocytes", "Neutrophils", "Whole blood"))
dat1 <- dat1[order(dat1$celltype, dat1$donor),]

dat2$celltype <- factor(dat2$celltype, levels = c("NK-cells", "B-cells", "CD4+ T-cells", "CD8+ T-cells", "Monocytes", "Neutrophils", "Whole blood"))
dat2 <- dat2[order(dat2$celltype, dat2$donor),]

# make sure donors are matched
stopifnot(length(unique(split(dat1$donor, dat1$celltype))) == 1)
stopifnot(length(unique(split(dat2$donor, dat2$celltype))) == 1)

# select only chronological clocks
# NOTE: skip Horvath pan-tissue since it is in Figure 2.
dat1 <- dat1[,c(1,19,20,2,4:8)]
dat2 <- dat2[,c(1,19,20,2,4:8)]

# subtract whole blood from all other celltypes
tmp1 <- split(dat1[,-(1:3)], dat1$celltype)
tmp1 <- Map(`-`, tmp1, tmp1["Whole blood"])
dat1[,-(1:3)] <- do.call(rbind, tmp1)

tmp2 <- split(dat2[,-(1:3)], dat2$celltype)
tmp2 <- Map(`-`, tmp2, tmp2["Whole blood"])
dat2[,-(1:3)] <- do.call(rbind, tmp2)

# remove whole blood
dat1 <- dat1[dat1$celltype != "Whole blood",]
dat1$celltype <- droplevels(dat1$celltype)

dat2 <- dat2[dat2$celltype != "Whole blood",]
dat2$celltype <- droplevels(dat2$celltype)


#--- plot ----------------------------------------------------------------------

violin <- function(x, g, cols, ylim) {
  dens <- tapply(x, g, density, na.rm=TRUE)
  meds <- tapply(x, g, median, na.rm=TRUE)
  ps   <- tapply(x, g, function(x) t.test(x)$p.value)

  xs <- Map(getElement, dens, "y")
  ys <- Map(getElement, dens, "x")

  xs <- Map(function(x) (x-min(x)) / max(x-min(x))/ 2.5, xs)

  xs <- Map(c, xs, Map(rev, Map(`*`, -1, xs)))
  ys <- Map(c, ys, Map(rev, ys))

  xs <- Map(`+`, xs, 1:length(xs))

  plot.new()
  plot.window(xlim = range(xs)+c(-0.15,0.15), ylim = ylim, xaxs = "i")
  abline(h = 0, lwd = 1, lty = 2, col = "grey")

  Map(polygon, xs, ys, col = adjustcolor(cols, 0.75))
  points(1:length(unique(g)), meds, cex = 3, pch = 18)
  points(as.factor(g), x, col = "white")
  ptxt <- paste0("p = ", ifelse(ps < 0.001, formatC(ps, format="e", digits=1), signif(ps, digits=2)))
  text(1:length(unique(g)), sapply(ys, max), ptxt, col = ifelse(ps < 0.05, "black", "gray50"), font = ifelse(ps < 0.05, 2, 1), pos = 3)

  axis(2, las = 1, lwd = 0)
  box()
}


pdf("fig_s5.pdf", width = 18/2.54, height = 21/2.54, pointsize = 4)

par(mfrow=c(6,2), mar=c(1,1,2,3), oma=c(2,5,0,0), tck=-0.02, mgp=c(2.5,0.5,0), cex = 0.9)

cols <- c("#F15854", "limegreen", "deepskyblue", "dodgerblue", "#DECF3F", "#FAA43A")

for (i in 4:ncol(dat1)) {
  ylim <- extendrange(pretty(c(dat1[,i], dat2[,i], -dat1[,i], -dat2[,i])), f=c(0.25,0.25))

  violin(dat1[,i], dat1$celltype, cols, ylim)
  if (i == ncol(dat1)) axis(1, at = 1:length(unique(dat1$celltype)), labels = unique(dat1$celltype), lwd = 0, font = 1)

  legend("bottomleft", legend = colnames(dat1)[i], bty = "n", inset = c(-0.03,1), xpd = NA)

  violin(dat2[,i], dat2$celltype, cols, ylim)
  if (i == ncol(dat2)) axis(1, at = 1:length(unique(dat2$celltype)), labels = unique(dat2$celltype), lwd = 0, font = 1)
}

title(ylab = "predicted age", line = 3, outer = TRUE)
title(ylab = expression(paste(Delta, " (celltype - whole blood)")), line = 1.5, outer = TRUE)

invisible(dev.off())

