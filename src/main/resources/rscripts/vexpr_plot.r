require(gplots); require(RColorBrewer); require(ggplot2)

args<-commandArgs(TRUE)

input_file  <- args[1] #"ms.segments.V2.txt"
v_col_count <- as.numeric(args[2]) #48
lbl_col     <- as.numeric(args[3]) #1
fac_col     <- as.numeric(args[4]) #2
cont_factor <- as.logical(args[5])
output_file <- args[6]

if (lbl_col < 1) {
   lbl_col = 1 # use sample id if not specified
}

color_by_factor <- fac_col > 0

# read-in data and pre-process
df<-read.table(input_file, header=T, sep="\t",comment="")
vcols<-(ncol(df)-v_col_count + 1):ncol(df)
df[, vcols] <- apply(df[, vcols], 2, as.numeric)
x <- as.matrix(df[, vcols])

# prepare palletes and factor levels
if (color_by_factor) {
    factor_name <- colnames(df)[fac_col]
    if (cont_factor) {
       df[,fac_col] <- as.numeric(df[,fac_col])
       fu   <- as.factor(as.numeric(cut(df[,fac_col],10)))
       scol <- colorRampPalette(c("#feb24c", "#31a354", "#2b8cbe"))(10)[fu]
    } else {
       df[, fac_col] <- as.factor(df[, fac_col])
       fu   <- levels(df[, fac_col])
       nlevels <- length(fu)
       scol <- colorRampPalette(brewer.pal(nlevels, "Set2"))(nlevels)[df[,fac_col]]
    }
}

# plotting function for simplicity
my.plot <- function(...) {
   heatmap.2(t(x),  na.rm=F, labCol = df[,lbl_col],
      na.col="grey50",
      col=colorRampPalette(c("#2f98ce", "#e0f3db", "#f47104")),
      density.info="none", trace="none",scale="column",
      ...)
}

if (color_by_factor) {
   pdf(output_file)

   # layout plot
   fig <- c(0.05, 0.75, 0, 1.0)
   mar <- c(0, 0, 0, 0)

   layout(matrix(1:2, ncol=2), width = c(2, 0.1), height = c(1, 1))

   par(fig = fig, mar = mar, xpd = NA) # this ensures labels are not cut

   # plot with coloring columns
   my.plot(ColSideColors=scol)

   # draw a separate legend
   fig <- c(0.85, 0.95, 0, 1.0)

   par(fig = fig, mar = c(0, 0, 0, 0), xpd = NA, new=TRUE)

   if (cont_factor) {
      # set colorbar pos
      f <- 2
      ux <- grconvertX(c(0.4 - 0.05 * f, 0.4 + 0.05 * f), from = "npc", to = "user")
      uy <- grconvertY(c(0.8 - 0.15, 0.8 + 0.15), from = "npc", to = "user")

      color.legend(ux[1], uy[1], ux[2], uy[2],
                   legend = fu[c(1, length(fu)/2 + 1, length(fu))],
                   rect.col = scol, gradient = "y", align="rb")

      # position text
      uy <- grconvertY(0.5 + 0.15 + 0.02, from = "npc", to = "user")
      text(0.5, uy, factor_name, adj = c(0.5, 0.0))
   } else {
      # vanilla legend for discrete factor
      legend("top", "(x,y)", fill = scol, pch = '', box.lwd = 0, title = factor_name, legend = fu)
   }

   dev.off()
} else {
   pdf(output_file)

   # layout plot
   fig <- c(0.05, 0.95, 0, 1.0)
   mar <- c(0, 0, 0, 0)

   layout(matrix(1:2, ncol=2), width = c(2, 0.1), height = c(1, 1))

   par(fig = fig, mar = mar, xpd = NA) # this ensures labels are not cut

   my.plot()

   dev.off()
}