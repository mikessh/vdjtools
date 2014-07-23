
filein = "intersection2_aa.txt"
id_col1 = "X.1_sample_id"
id_col2 = "X2_sample_id"
measure_col = "freq12"
factor_col1 = "X1_age"
factor_col2 = "X2_age"

require(ape); require(reshape)

tbl <- read.delim(file_in)
df <- data.frame(tbl)
df.2 <- df
df.2[c(id_col1, id_col2)] <- df.2[c(id_col2, id_col1)]
df.sym <- rbind(df, df.2)
df.m <- cast(df.sym, id_col1 ~ id_col2, mean, value = measure_col)
df.d <- dist(df.m)
hcl <- hclust(df.d)

# trick as we need both columns to get whole data on a factor
# due to the fact that the output is upper triangle of intersection matrix
factor <- unique(df[c(id_col2, factor_col2)])
names(factor) <- names(df[c(id_col1, factor_col1)])
factor <- unique(rbind(factor, df[c(id_col1, factor_col1)]))

# switch to numeric values and sort by factor
factor[, factor_col1] <- sapply(factor[, factor_col1], as.character)
factor[, factor_col1] <- sapply(factor[, factor_col1], as.numeric)
factor <- factor[order(factor[, factor_col1]), ]

# design a palette to color by unique factor levels

pal <- colorRampPalette(c("red", "black", "green"))
fu <- unique(factor[, factor_col1])
col <- pal(length(fu))

# for nan
if(any(is.na(fu))) {
    col[length(fu)]="grey"
}

#
ind1 <- match(factor[, factor_col1], fu)
ind2 <- match(hcl$labels, factor[, id_col1])


plot(as.phylo(hcl), type = "fan", tip.color = col[ind1[ind2]])

dev.off()

mds <- isoMDS(df.d, k = 2)

apply.bounds <- function(x, y, lp, hp) {
   qx <- quantile(x, c(lp, hp))
   qy <- quantile(y, c(lp, hp))

   oxl <- x <= qx[1]; oxh <- x >= qx[2]
   oyl <- y <= qy[1]; oyh <- y >= qy[2]

   xx <- x
   xx[oxl] = qx[1]; xx[oxh] = qx[2]
   yy <- y
   yy[oyl] = qy[1]; yy[oyh] = qy[2]

   o <- (oxl | oxh | oyl | oyh)

   list(x = xx, y = yy, o = o)
}

x <- mds$points[,1]
y <- mds$points[,2]
xy <- apply.bounds(x, y, 0.1, 0.9)

ind3 <- match(row.names(mds$points), factor$X.1_sample_id)

plot(xy.x, xy.y, xlab="mds1", ylab="mds2", main="MDS", type = "n")#, col=ifelse(xy.o, "red", "black"), pch="+")
text(xy.x, xy.y, labels = row.names(as.matrix(df.d)), col = col[ind1[ind3]], cex=.5)

