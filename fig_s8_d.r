#--- init ----------------------------------------------------------------------

# for cosinor test
library(matrixTests)

# read epigenetic age predictions from the supplement
dat <- read.csv("input/WBCminusNeu_clocks_IEAA.csv", skip = 1)


#--- obtain oscillation results ------------------------------------------------

# replace id with time
dat$id <- as.numeric(gsub(".*_CT|_R1", "", dat$id))
colnames(dat)[1] <- "time"

# combine technical replicates
dat <- apply(dat, 2, tapply, dat$time, mean)

# get oscillation results
res <- col_cosinor(dat[,-1], dat[,"time"])


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


pdf("fig_s8_d.pdf", width = 18/2.54, height = 5/2.54, pointsize = 8)

par(mar = c(0,4,8,4), oma = c(0,2,0,0), mgp = c(2.5,0.5,0), cex.main = 1)
layout(matrix(c(1,1,1), nrow = 1))

plot.new()
plot.window(xlim = c(0, nrow(res)+1), ylim = c(-1, 5), xaxs = "i")

res  <- res[,c("pvalue", "rsquared", "mesor", "amplitude", "acrophase")]

# get table x and y placements and add offsets
ys <- col(res)
xs <- row(res)
xs[xs > 7]  <- xs[xs > 7]  + 0.5
xs[xs > 12] <- xs[xs > 12] + 0.5

rect(xs-1, ys-1, xs, ys)

text(xs-0.5, 4.5, fpval(res$pvalue),     col = ifelse(res$pvalue <= 0.05, "black", "gray50"), font = ifelse(res$pvalue <= 0.05, 2, 1))
text(xs-0.5, 3.5, fdigit(res$rsquared),  col = ifelse(res$pvalue <= 0.05, "black", "gray50"), font = ifelse(res$pvalue <= 0.05, 2, 1))
text(xs-0.5, 2.5, fdigit(res$mesor),     col = ifelse(res$pvalue <= 0.05, "black", "gray50"), font = ifelse(res$pvalue <= 0.05, 2, 1))
text(xs-0.5, 1.5, fdigit(res$amplitude), col = ifelse(res$pvalue <= 0.05, "black", "gray50"), font = ifelse(res$pvalue <= 0.05, 2, 1))
text(xs-0.5, 0.5, fhour(res$acrophase),  col = ifelse(res$pvalue <= 0.05, "black", "gray50"), font = ifelse(res$pvalue <= 0.05, 2, 1))

text(0, ys[1,-4]-0.5, c("acrophase", "amplitude", "mesor", "p-value"), pos = 2, xpd = NA)
text(0, ys[1,4]-0.5, expression(paste("R"^"2")), pos = 2, xpd = NA)
text(xs-0.75, 5.35, rownames(res), pos = 4, srt = 45, xpd = NA)

text(c(3.5, 9.5, 15), -0.75, c("Chronological clocks", "Mitotic clocks", "Biological clocks"))

title("d", adj = 0, cex.main = 1.5, line = 7)

invisible(dev.off())
