#--- init ----------------------------------------------------------------------

# for t-tests
library(matrixTests)

# read epigenetic age predictions from the supplement
dat <- read.csv("input/PBMCtimeofday_clocks_IEAA.csv", skip = 1)


#--- obtain differences between 12:45 and 16:15 --------------------------------

# extract time and donor information from ids
dat$time  <- sub("_.*", "", dat$id)
dat$donor <- sub(".*_", "", dat$id)

# make sure pairs of donors are matched
stopifnot(all.equal(dat$donor[dat$time == "12:45"], dat$donor[dat$time == "16:15"]))

res <- col_t_paired(dat[dat$time == "16:15", 2:18], dat[dat$time == "12:45", 2:18])


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


pdf("fig_s9_d.pdf", width = 18/2.54, height = 5/2.54, pointsize = 8)

par(mar = c(0,4,8,4), oma = c(0,2,0,0), mgp = c(2.5,0.5,0), cex.main = 1)
layout(matrix(c(1,1,1), nrow = 1, byrow = TRUE))

plot.new()
plot.window(xlim = c(0, nrow(res)+1), ylim = c(-1, 4), xaxs = "i")

res <- res[,c("mean.diff", "pvalue")]

# get table x and y placements and add offsets
ys <- col(res) + 2
xs <- row(res)
xs[xs > 7]  <- xs[xs > 7]  + 0.5
xs[xs > 12] <- xs[xs > 12] + 0.5

rect(xs-1, ys-1, xs, ys)

text(xs-0.5, 3.5, fpval(res$pvalue),     col = ifelse(res$pvalue <= 0.05, "black", "gray50"), font = ifelse(res$pvalue <= 0.05, 2, 1))
text(xs-0.5, 2.5, fdigit(res$mean.diff),  col = ifelse(res$pvalue <= 0.05, "black", "gray50"), font = ifelse(res$pvalue <= 0.05, 2, 1))

text(0, ys[1,2]-0.5, "p-value", pos = 2, xpd = NA)
text(0, ys[1,1]-0.5, expression(paste(Delta, " 16:15-12:45")), pos = 2, xpd=NA)
text(xs-0.75, 4.35, rownames(res), pos = 4, srt = 45, xpd = NA)

text(c(3.5, 9.5, 15), 1.25, c("Chronological clocks", "Mitotic clocks", "Biological clocks"))

title("h", adj = 0, cex.main = 1.5, line = 7)

invisible(dev.off())
