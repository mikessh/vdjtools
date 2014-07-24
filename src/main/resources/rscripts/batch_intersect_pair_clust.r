args<-commandArgs(TRUE)

# read in

file_in           = args[1]              #"intersection2_aa.txt"
id_col1_index     = as.integer(args[2])  #1#"X.1_sample_id"
id_col2_index     = as.integer(args[3])  #9#"X2_sample_id"
measure_col_index = as.integer(args[4])  #22#"freq12"
factor_col1_index = as.integer(args[5])  #4#"X1_age"
factor_col2_index = as.integer(args[6])  #12#"X2_age"
lbl_col1_index    = as.integer(args[7])  #4  #
lbl_col2_index    = as.integer(args[8])  #12 #
file_out_hc       = args[9]
file_out_mds      = args[10]

color_by_factor <- TRUE

if(factor_col1_index < 1) {
   factor_col1_index = id_col1_index
   factor_col2_index = id_col2_index
   color_by_factor   = FALSE
}

if(lbl_col1_index < 1) {
   lbl_col1_index = id_col1_index
   lbl_col2_index = id_col2_index
}

require(ape); require(reshape); require(MASS)

# read data

tbl  <- read.delim(file_in)
df   <- data.frame(tbl)

# Split tables to protect from column overlap
# Also rename them
# all this stuff is done as R re-formats column names
# e.g. 1_sample_id to X1_sample_id, ...

df.aux <- data.frame(
    id_col1 = df[, id_col1_index], id_col2 = df[, id_col2_index],
    factor_col1 = df[, factor_col1_index], factor_col2 = df[, factor_col2_index],
    lbl_col1 = df[, lbl_col1_index], lbl_col2 = df[, lbl_col2_index]
    )

df <- data.frame(
    id_col1 = df[, id_col1_index], id_col2 = df[, id_col2_index],
    measure_col = df[, measure_col_index]
    )

## Auxillary table

# trick as we need both columns to get all the data on a factor
# due to the fact that the output is upper triangle of intersection matrix
aux        <- unique(df.aux[c("id_col2", "factor_col2", "lbl_col2")])
names(aux) <- names(df.aux[c("id_col1", "factor_col1", "lbl_col1")])
aux        <- unique(rbind(aux, df.aux[c("id_col1", "factor_col1", "lbl_col1")]))

## Factor & coloring

if (color_by_factor) {
    # switch to numeric values and sort by factor
    aux[, "factor_col1"] <- sapply(aux[, "factor_col1"], as.character)
    aux[, "factor_col1"] <- sapply(aux[, "factor_col1"], as.numeric)
    aux <- aux[order(aux[, "factor_col1"]), ]

    # design a palette to color by unique factor levels

    pal <- colorRampPalette(c("red", "black", "green"))
    fu <- unique(aux[, "factor_col1"])
    cc <- pal(length(fu))

    # for nan
    if(any(is.na(fu))) {
        cc[length(fu)] <- "grey"
    }
}

## Distance

# Symmetrize matrix & compute distance measure

df.2 <- df
df.2[c("id_col1", "id_col2")] <- df.2[c("id_col2", "id_col1")]
df.sym <- rbind(df, df.2)
df.m <- cast(df.sym, id_col1 ~ id_col2, mean, value = "measure_col")

df.d <- dist(df.m)

## HCL

hcl <- hclust(df.d)

## plotting

# for matching colors and labels
cc_final <- "black"

if (color_by_factor) {
   ind1 <- match(aux[, "factor_col1"], fu)
   ind2 <- match(hcl$labels, aux[, "id_col1"])
   cc_final <- cc[ind1[ind2]]
}

lbl  <- sapply(aux[match(hcl$labels, aux[, "id_col1"]), "lbl_col1"], as.character)

# set lables
phylo <- as.phylo(hcl)
phylo$tip.label <- lbl

# plot
pdf(file_out_hc)

plot(phylo, type = "fan", tip.color = cc_final)

dev.off()

## MDS

mds <- isoMDS(df.d, k = 2)

# move outliers to plot boundaries

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

# re-match label and color

if (color_by_factor) {
   ind3 <- match(row.names(mds$points), aux[, "id_col1"])
   cc_final <-  cc[ind1[ind3]]
}

lbl  <- sapply(aux[match(row.names(as.matrix(df.d)), aux[, "id_col1"]), "lbl_col1"], as.character)

pdf(file_out_mds)

plot(xy$x, xy$y, xlab="mds1", ylab="mds2", main="MDS", type = "n")
text(xy$x, xy$y, labels = lbl, col = cc_final, cex=.5)

dev.off()
