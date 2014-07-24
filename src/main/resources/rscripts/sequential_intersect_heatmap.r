args <- commandArgs(TRUE)

require(ggplot2); require(reshape)

inputF  <- args[1]
inputD  <- args[2]
inputR  <- args[3]
points  <- args[4]
file_out  <- args[5]

#fileConn<-file("debugF.txt")
#writeLines(inputF, fileConn)
#close(fileConn)

#fileConn<-file("debugD.txt")
#writeLines(inputD, fileConn)
#close(fileConn)

#fileConn<-file("debugR.txt")
#writeLines(inputR, fileConn)
#close(fileConn)

# preprocess data

proc_input <- function(input, type) {
   tbl <- read.table(text = input, sep="\t")#text = input)
   tbl <- t(apply(tbl, 1:2, as.numeric))

   # append type for facet
   tbl.m <- melt(tbl)
   tbl.m$T <- rep(type, nrow(tbl))

   tbl.m
}

dF.m <- proc_input(inputF, "Frequency")
dD.m <- proc_input(inputD, "Diversity")
dR.m <- proc_input(inputR, "Correlation")

p <- apply(as.vector(read.table(text = points, sep = ";")), 1, as.character)
n <- length(p)

pdf(file_out)

ggplot(rbind(dF.m, dD.m, dR.m), aes(X1, X2)) +
    geom_tile(aes(fill = value), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    coord_fixed(ratio = 1) +
    facet_wrap( ~ T) +
    scale_x_discrete(expand = c(0, 0), labels = p) +
    scale_y_discrete(limits = 1:n, breaks = 1:n, expand = c(0, 0), labels = p) +
    xlab("") + ylab("") +
    guides(fill = guide_legend(title = "value")) +
    theme(
             legend.position = "bottom",
             axis.ticks = element_blank(),
             axis.text.x = element_text(angle = 90, hjust = 0)
    )

dev.off()