#--- init ----------------------------------------------------------------------

# for cosinor test
library(matrixTests)

# read cell count data from the supplement
cells <- read.csv("input/WBCminusNeu_cellcounts_Reinius.csv", skip = 1)


#--- obtain cellcount oscillations ---------------------------------------------

# replace id with time
cells$id <- as.numeric(gsub(".*_CT|_R1", "", cells$id))
colnames(cells)[1] <- "time"

# combine technical replicates
cells <- apply(cells, 2, tapply, cells$time, mean)
cells <- as.data.frame(cells)

# add full celltype names
colnames(cells) <- c("time", "CD4+ T-cells", "CD8+ T-cells", "B-cells", "NK-cells", "Monocytes", "Granulocytes")

# get oscillation results
res <- col_cosinor(cells[,-1], cells[,"time"])

# reorder cells and select result columns
res <- res[c(4,3,1,2,5,6), c(9,5,2,3,4)]


#--- plot ----------------------------------------------------------------------

# function to format p-value to text
fpval <- function(x) {
  res <- ifelse(x < 0.001, formatC(x, format = "e", digits = 1), signif(x, 2))
  ifelse(p.adjust(x) <= 0.05, paste0(res, "*"), res)
}

# function to format digits to text
fdigit <- function(x) {
  x <- round(x, 2)
  x[abs(x) > 10]  <- round(x[abs(x) > 10], 1)
  x[abs(x) > 100] <- round(x[abs(x) > 100])
  formatC(x)
}

# function to format hours to text
fhour <- function(x) {
  res <- round(60 * x %% 1)
  res <- ifelse(nchar(res) == 1, paste0("0", res), res)
  paste(floor(x), res, sep=":")
}


pdf("fig_2_b.pdf", width = 9/2.54, height = 5/2.54, pointsize = 4)

par(mar = c(1,3,3,3), oma = c(3,5,0,0), mgp = c(2.5,0.5,0), cex.main = 1)

ys <- col(res)
xs <- row(res)

plot.new()
plot.window(xlim = c(0, nrow(res)), ylim = c(0, 7), xaxs = "i", yaxs = "i")

rect(xs-1, ys-1, xs, ys)

text(xs-0.5, 4.5, fpval(res$pvalue),     col = ifelse(res$pvalue <= 0.05, "black", "gray50"), font = ifelse(res$pvalue <= 0.05, 2, 1))
text(xs-0.5, 3.5, fdigit(res$rsquared),  col = ifelse(res$pvalue <= 0.05, "black", "gray50"), font = ifelse(res$pvalue <= 0.05, 2, 1))
text(xs-0.5, 2.5, fdigit(res$mesor),     col = ifelse(res$pvalue <= 0.05, "black", "gray50"), font = ifelse(res$pvalue <= 0.05, 2, 1))
text(xs-0.5, 1.5, fdigit(res$amplitude), col = ifelse(res$pvalue <= 0.05, "black", "gray50"), font = ifelse(res$pvalue <= 0.05, 2, 1))
text(xs-0.5, 0.5, fhour(res$acrophase),  col = ifelse(res$pvalue <= 0.05, "black", "gray50"), font = ifelse(res$pvalue <= 0.05, 2, 1))

text(0, ys[1,-4]-0.5, c("acrophase", "amplitude", "mesor", "p-value"), pos = 2, xpd = NA)
text(0, ys[1,4]-0.5, expression(paste("R"^"2")), pos = 2, xpd = NA)
text(xs-0.75, 5.35, rownames(res), pos = 4, srt = 45, xpd = NA)

title("b", adj = 0, cex.main = 1.5)

invisible(dev.off())

