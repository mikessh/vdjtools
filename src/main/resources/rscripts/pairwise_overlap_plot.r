require(circlize); require(reshape2); require(RColorBrewer)

args<-commandArgs(TRUE)

file_in  <- args[1]
file_out <- args[2]

df <- read.table(file_in, header=T, comment="", quote="")
df[,] <- apply(df[,], 2, as.character)
df <- data.frame(sample1 = c(df$X.1_sample_id, df$X2_sample_id),
                 sample2 = c(df$X2_sample_id, df$X.1_sample_id),
                 count = as.numeric(c(df$count12, df$count21)),
                 frequency = as.numeric(c(df$freq12, df$freq21)),
                 diversity = as.numeric(c(df$div12, df$div21)))

# equate the number of characters for visualizaiton
adjust_labels <- function(mat) {
  rn <- rownames(mat)
  cn <- colnames(mat)

  maxrn <- max(nchar(rn))
  maxcn <- max(nchar(cn))

  for(i in seq_len(length(rn))) {
    rn[i] <- paste(rn[i], paste(rep(" ", maxrn - nchar(rn[i])), collapse = ''))
  }

  for(i in seq_len(length(cn))) {
    cn[i] <- paste(cn[i], paste(rep(" ", maxcn - nchar(cn[i])), collapse = ''))
  }

  rownames(mat) <- rn
  colnames(mat) <- cn
}

# plotting function
render_plot <- function(var) {
  mat <- acast(df, sample1 ~ sample2, value.var=var)
  
  diag(mat) <- 0
  #diag(mat) <- df[df$sample1==rownames(mat),"f1"] * (nrow(df) / 2 - 1) - # total possible overlap freq
  #  rowSums(mat, na.rm=T) # overlapped freq, note that row -> col is default direction for circlize
  adjust_labels(mat)

  # viz using circlize
  #circos.par(gap.degree = c(rep(3, nrow(mat)-1), 10, rep(3, ncol(mat)-1), 10), start.degree = 5)

  cols <- colorRampPalette(brewer.pal(12, "Paired"))(nrow(mat))

  chordDiagram(mat, annotationTrack = "grid", grid.col = cols,
               preAllocateTracks = list(track.height = 0.2), transparency = 0.5)

  circos.trackPlotRegion(track.index = 1, bg.border = NA,
                         panel.fun = function(x, y) {
                           sector.name = get.cell.meta.data("sector.index")
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           circos.text(mean(xlim), ylim[1], cex = 0.5, sector.name, facing = "clockwise", adj = c(0, 0.5))
                         }
  )
}

pdf(file_out)
vars <- c("count", "frequency", "diversity")
layout(t(matrix(c(0,0,0,1,2,3,0,0,0),3,3)))
for(var in vars) {
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  render_plot(var)
  title(var)
  circos.clear()
}
dev.off()