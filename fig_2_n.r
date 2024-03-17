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

# get oscillations
# NOTE: here we allow each donor to have separate oscillation curve
#       this is achieved by adding an interaction term between sin/cos and donor
res <- data.frame(row.names = colnames(dat)[-(1:2)])
res$pvalue <- NA
for(cl in rownames(res)) {
  osc  <- lm(dat[,cl] ~ sin(2*pi/24*dat$time) * dat$donor + cos(2*pi/24*dat$time) * dat$donor)
  null <- lm(dat[,cl] ~ dat$donor)
  res[cl, "pvalue"] <- anova(osc, null)[2,6]
}


#--- plot ----------------------------------------------------------------------

# function to format p-value to text
fpval <- function(x) {
  res <- ifelse(x < 0.001, formatC(x, format = "e", digits = 1), signif(x, 2))
  ifelse(p.adjust(x) <= 0.05, paste0(res, "*"), res)
}

pdf("fig_2_n.pdf", width = 18/2.54, height = 5/2.54, pointsize = 8)

par(mar = c(0,4,8,4), oma = c(0,2,0,0), mgp = c(2.5,0.5,0), cex.main = 1)
layout(matrix(c(1,1,1), nrow = 1))

plot.new()
plot.window(xlim = c(0, nrow(res)+1), ylim = c(-1, 5), xaxs = "i")

# get table x and y placements and add offsets
ys <- 5
xs <- row(res)
xs[xs > 7]  <- xs[xs > 7]  + 0.5
xs[xs > 12] <- xs[xs > 12] + 0.5

rect(xs-1, ys-1, xs, ys)

text(xs-0.5, 4.5, fpval(res$pvalue), col = ifelse(res$pvalue <= 0.05, "black", "gray50"), font = ifelse(res$pvalue <= 0.05, 2, 1))

text(0, 4.5, "p-value", pos = 2, xpd=NA)
text(xs-0.75, 5.35, rownames(res), pos = 4, srt = 45, xpd = NA)

text(c(3.5, 9.5, 15), 3.25, c("Chronological clocks", "Mitotic clocks", "Biological clocks"))

title("n", adj = 0, cex.main = 1.5, line = 7)

invisible(dev.off())
