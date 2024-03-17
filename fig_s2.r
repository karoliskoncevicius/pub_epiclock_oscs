#--- init ----------------------------------------------------------------------

library(basetheme)

basetheme("clean")


#--- obtain cellcount results --------------------------------------------------

# read cell count data from the supplement
dat <- read.csv("input/WBCminusNeu_cellcounts_Salas.csv", row.names = 1, skip = 1)

# shorten row names
rownames(dat) <- sub(".*?_", "", rownames(dat))
rownames(dat) <- sub("_R1", " replicate", rownames(dat))


#--- plot ----------------------------------------------------------------------

# coloring celltypes
cols <- c(tcd4n = "deepskyblue", tcd4m = "dodgerblue", tcd8n = "lightblue2", tcd8m = "skyblue3",
          treg  = "#B276B2",        bn = "lightgreen",    bm = "#60BD68",       nk = "#F15854",
          bas   = "#C4A484",       eos = "#FFF804",     mono = "#DECF3F",      neu = "#FAA43A")

cell2col <- lab2col(pal = cols, ref = names(cols))


pdf("fig_s2.pdf", width = 18/2.54, height = 20/2.54, pointsize = 8)

par(mar = c(4,6,5,3))

barplot(t(dat[nrow(dat):1,]) * 100, horiz = TRUE, xlim = c(0,100), col = cell2col(colnames(dat)), xlab = "predicted cell type %", yaxs = "i", xaxs = "i")

legend("bottom",
       c("CD4+ naive T-cells", "CD4+ memory T-cells", "CD8+ naive T-cells", "CD8+ memory T-cells", "Regulatory T-cells", "Naive B-cells", "Memory B-cells",  "NK-cells", "Basofils", "Eosonofils", "Neutrophils", "Monocytes"),
       fill = cell2col(c("tcd4n", "tcd4m", "tcd8n", "tcd8m", "treg", "bn", "bm", "nk", "bas", "eos", "neu", "mono")),
       inset = c(0, 1.01), xpd = TRUE, ncol = 3, bty = "n", cex = 0.8)

invisible(dev.off())
