#--- init ----------------------------------------------------------------------

# read cell count data from the supplement
cells <- read.csv("input/WBCminusNeu_cellcounts_Reinius.csv", skip = 1)


#--- obtain cellcount oscillations ---------------------------------------------

# replace id with time
cells$id <- as.numeric(gsub(".*_CT|_R1", "", cells$id))
colnames(cells)[1] <- "time"

cells <- apply(cells, 2, tapply, cells$time, mean)
cells <- as.data.frame(cells)

colnames(cells) <- c("time", "CD4+ T-cells", "CD8+ T-cells", "B-cells", "NK-cells", "Monocytes", "Granulocytes")


#--- plot ----------------------------------------------------------------------

times  <- cells[,1]
counts <- data.matrix(cells[,-1]) * 100
colors <- c("deepskyblue", "dodgerblue", "limegreen", "#F15854", "#DECF3F", "#FAA43A")


pdf("fig_2_a.pdf", width = 9/2.54, height = 5/2.54, pointsize = 4)

par(mar = c(1,1,3,2), oma = c(3,4,0,0), mgp = c(2.5,0.5,0), cex.main = 1)

plot.new()
plot.window(xlim = c(-6,98), ylim = c(0,35), xaxs = "i", yaxs = "i")
grid(nx=NA, ny=NULL)

points(rep(times, ncol(counts)), counts, col = rep(colors, each = nrow(counts)), cex = 0.8)

for(i in 1:ncol(counts)) {
  fit0 <- lm(counts[,i] ~ 1)
  fit1 <- lm(counts[,i] ~ sin(2*pi/24*times) + cos(2*pi/24*times))
  pval <- anova(fit1, fit0)[2,6]

  xs <- seq(12, 96, 0.1)
  lines(xs, coef(fit1)[1] + coef(fit1)[2] * sin(2*pi/24*xs) + coef(fit1)[3] * cos(2*pi/24*xs), col = colors[i], lty = ifelse(pval <= 0.05, 1, 2))
  offset <- ifelse(colnames(counts)[i] == "NK-cells", 0.9, 0)            # NOTE: to better position celltype labels
  offset <- ifelse(colnames(counts)[i] == "CD4+ T-cells", -0.75, offset) # NOTE: to better position celltype labels
  text(12, coef(fit1)[1] + coef(fit1)[2] * sin(2*pi/24*12) + coef(fit1)[3] * cos(2*pi/24*12) + offset, colnames(counts)[i], col = colors[i], pos = 2, cex = 0.8)
}

axis(1, at=seq(0,96,12), lwd=0)
axis(2, las=1, lwd=0)
box()

title("a", adj = 0, cex.main = 1.5)
title(xlab = "CT", xpd = NA)
title(ylab = "predicted cell type %", line = 3, xpd = NA)

invisible(dev.off())

