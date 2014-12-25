args<-commandArgs(TRUE)

#file_in  <- "4.fancyvj.txt"
#file_out <- "4.fancyvj.pdf"
file_in  <- args[1]
file_out <- args[2]

require(circlize)


# load data and preproc to fit formats

temp <- read.table(file_in, sep="\t", comment="")
n <- nrow(temp)
m <- ncol(temp)
rn = as.character(temp[2:n,1])
cn = apply(temp[1,2:m], 2 , as.character) #as.character(colnames(temp)[2:m])
mat <- matrix(apply(temp[2:n, 2:m], 1:2, as.numeric), n - 1, m-1) * 100

n <- nrow(temp)
m <- ncol(temp)

rownames(mat) <- rn
colnames(mat) <- cn

# sort

col_sum = apply(mat, 2, sum)
row_sum = apply(mat, 1, sum)

mat <- mat[order(row_sum), order(col_sum)]


# equal number of characters for visualizaiton

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


# viz using circlize

pdf(file_out)

circos.par(gap.degree = c(rep(3, nrow(mat)-1), 10, rep(3, ncol(mat)-1), 15), start.degree = 5)

chordDiagram(mat, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.2), transparency = 0.5)

circos.trackPlotRegion(track.index = 1, bg.border = NA,
       panel.fun = function(x, y) {
                   sector.name = get.cell.meta.data("sector.index")
                   xlim = get.cell.meta.data("xlim")
                   ylim = get.cell.meta.data("ylim")
                   circos.text(mean(xlim), ylim[1], cex = 0.5, sector.name, facing = "clockwise", adj = c(0, 0.5))
                   }
       )

circos.clear()

dev.off()