args<-commandArgs(TRUE)

# read in

file_in           = args[1] #"mmu_sep13_aa.txt"
id_col1_index     = as.integer(args[2])  #1
id_col2_index     = as.integer(args[3])  #2
measure_col_index = as.integer(args[4])  #18
factor_col1_index = as.integer(args[5])  #21
factor_col2_index = as.integer(args[6])  #22
lbl_col1_index    = as.integer(args[7])  #1#
lbl_col2_index    = as.integer(args[8])  #2
factor_name       = args[9]  #"Group"
cont_factor       = args[10] #TRUE
file_out_hc       = args[11]
file_out_mds      = args[12]
k_clust           = args[13]
file_out_clust    = args[14]
file_out_coord    = args[15]

color_by_factor <- TRUE

if(factor_col1_index < 1) {
   factor_col1_index = id_col1_index
   factor_col2_index = id_col2_index
   color_by_factor   = FALSE
   cont_factor   = FALSE
}

if(lbl_col1_index < 1) {
   lbl_col1_index = id_col1_index
   lbl_col2_index = id_col2_index
}

require(ape); require(reshape2); require(ggplot2); require(MASS); require(plotrix) #require(SDMTools)

# read data
tbl  <- read.delim(file_in)
df   <- data.frame(tbl)

# convert factor columns depending on if continuous coloring is desired or not

if (cont_factor) {
    df[, factor_col1_index] <- as.numeric(as.character(df[, factor_col1_index]))
    df[, factor_col2_index] <- as.numeric(as.character(df[, factor_col2_index]))
} else {
    df[, factor_col1_index] <- as.factor(df[, factor_col1_index])
    df[, factor_col2_index] <- as.factor(df[, factor_col2_index])
}

# Split tables to protect from column overlap, also rename them
# all this stuff is done as R re-formats column names
# e.g. 1_sample_id to X1_sample_id, ...

df.aux <- data.frame(
    id_col1 = df[, id_col1_index], id_col2 = df[, id_col2_index],
    factor_col1 = df[, factor_col1_index], factor_col2 = df[, factor_col2_index],
    lbl_col1 = df[, lbl_col1_index], lbl_col2 = df[, lbl_col2_index]
    )

df <- data.frame(
    id_col1 = df[, id_col1_index], id_col2 = df[, id_col2_index],
    measure_col = as.numeric(as.character(df[, measure_col_index]))
    )

## Auxillary table

# trick as we need both columns to get all the data on a factor
# due to the fact that the output is upper triangle of intersection matrix

aux        <- unique(df.aux[c("id_col2", "factor_col2", "lbl_col2")])
names(aux) <- names(df.aux[c("id_col1", "factor_col1", "lbl_col1")])
aux        <- unique(rbind(aux, df.aux[c("id_col1", "factor_col1", "lbl_col1")]))

## Factor & coloring

if (color_by_factor) {
    if (cont_factor) {
       # switch to numeric values and sort by factor
       aux[, "factor_col1"] <- as.numeric(as.character(aux[, "factor_col1"]))
       aux <- aux[order(aux[, "factor_col1"]), ]
    }

    # design a palette to color by unique factor levels
    pal <- colorRampPalette(c("#feb24c", "#31a354", "#2b8cbe"))
    fu <- unique(aux[, "factor_col1"])
    cc <- pal(length(fu))

    # for nan
    cc[is.na(fu)] <- "grey"
}

## Distance

# Symmetrize matrix & compute distance measure

df.1 <- data.frame(id_col1 = df$id_col1, id_col2 = df$id_col2, measure_col = df$measure_col)
df.2 <- df.1
df.2[c("id_col1", "id_col2")] <- df.2[c("id_col2", "id_col1")]
df.sym <- rbind(df.1, df.2)

# cast back to square matrix
df.m <- as.matrix(dcast(df.sym, id_col1 ~ id_col2, value.var = "measure_col"))
rownames(df.m) <- df.m[,1]
df.m <- df.m[,2:ncol(df.m)]
df.m <- df.m[order(rownames(df.m)), order(colnames(df.m))]

# compute distance
df.m <- apply(df.m, 2, as.numeric) # Don't ask me why
diag(df.m) <- 0 # replace diag NAs with 0
df.d <- as.dist(df.m)

## HCL

hcl <- hclust(df.d)

k_clust <- as.integer(k_clust)
if (k_clust > 0) {
   clusters <- cutree(hcl, k = k_clust)
   clusters <- data.frame(id_col1 = names(clusters), cluster = clusters)
   ctbl <- merge(aux, clusters, by = "id_col1")
   write.table(ctbl, file_out_clust, sep = "\t", quote = FALSE, row.names = FALSE)
}

## Plotting

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

# plotting functions

my.plot <- function(hcl, ...) {
   if (color_by_factor) {
      if (hcl) {
         fig <- c(0.05, 0.75, 0, 1.0)
         mar <- c(0, 0, 0, 0)
      } else {
         fig <- c(0.05, 0.9, 0.1, 0.9)
         mar <- c(2, 4, 2, 4)
      }
   } else {
      if (hcl) {
         fig <- c(0.05, 0.95, 0, 1.0)
         mar <- c(0, 0, 0, 0)
      } else {
         fig <- c(0, 1.0, 0, 1.0)
         mar <- c(4, 4, 4, 4)
      }
   }

   layout(matrix(1:2, ncol=2), width = c(2, 0.1), height = c(1, 1))

   par(fig = fig, mar = mar, xpd = NA) # this ensures labels are not cut

   plot(...)
}

my.legend <- function(hcl) {
   if (color_by_factor) {
      if (hcl) {
         f <- 2
         fig <- c(0.85, 0.95, 0, 1.0)
      } else {
         f <- 1.2
         fig <- c(0.8, 0.95, 0, 1.0)
      }
      par(fig = fig, mar = c(0, 0, 0, 0), xpd = NA, new=TRUE)
      if (cont_factor) {
         # custom legend.gradient
         #px = c(-0.075, 0.075, 0.075, -0.075)
         #py = c(-0.1,   -0.1,  0.1,    0.1)

         # order & get rid of NAs
         fu1 <- fu[ind1]
         cc1 <- cc[ind1]
         cc1 <- cc1[!is.na(fu1)]
         fu1 <- fu1[!is.na(fu1)]

         ux <- grconvertX(c(0.4 - 0.05 * f, 0.4 + 0.05 * f), from = "npc", to = "user")
         uy <- grconvertY(c(0.5 - 0.15, 0.5 + 0.15), from = "npc", to = "user")

         color.legend( ux[1], uy[1], ux[2], uy[2],#-0.39, -0.4,  + 0.39, 0.4,#0.07 + f, 0.1,
                      legend = fu1[c(1, length(fu1)/2 + 1, length(fu1))],
                      rect.col = cc1, gradient = "y", align="rb")

         uy <- grconvertY(0.5 + 0.15 + 0.02, from = "npc", to = "user")
         text(0.5, uy, factor_name, adj = c(0.5, 0.0))

         # old impl
         # rect(-0.075, -0.1, 0.075, 0.1)
         #legend.gradient(cbind(x = px, y = py), cols = cc[ind1], title = factor_name, limits = c(fu[1], fu[length(fu)]))
      } else {
         # vanilla legend for discrete factor
         legend("center", "(x,y)", fill = cc, pch = '', box.lwd = 0, title = factor_name, legend = fu)
      }
   }
}

# draw dendrogram
pdf(file_out_hc)

my.plot(TRUE, phylo, type = "fan", tip.color = cc_final)
my.legend(TRUE)

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
fac  <- sapply(aux[match(row.names(as.matrix(df.d)), aux[, "id_col1"]), "factor_col1"], as.character)

pdf(file_out_mds)

my.plot(FALSE, xy$x, xy$y, xlab="mds1", ylab="mds2", type = "n")
text(xy$x, xy$y, labels = lbl, col = cc_final, cex=.5)
my.legend(FALSE)

dev.off()

write.table(data.frame(id = aux[, "id_col1"], lbl = lbl, factor = fac, x=xy$x, y=xy$y),
   file_out_coord, sep = "\t", quote = FALSE, row.names = FALSE)

#df.c <- data.frame(x = xy$x, y = xy$y, f = aux[match(row.names(as.matrix(df.d)), aux[, "id_col1"]), "factor_col1"])

#write.table(df.c, "tmp.txt")

#pdf("tmp.pdf")

#ggplot() + stat_density2d(data = df.c, aes(x, y, weight = f), geom="tile", contour = FALSE)

#dev.off()
