args<-commandArgs(TRUE)

label    <- args[1]
input    <- args[2]
file_out <- args[3]

require(ggplot2); require(reshape)

# preprocess data

tbl <- read.delim(input, sep = "\t")
tbl <- apply(tbl, 1:2, as.numeric)

df <- data.frame(tbl)
df <- melt(df, id = 'time', variable_name = 'series')

# plot series as facets

if (grepl("\\.pdf$",file_out)){
   pdf(file_out)
} else if (grepl("\\.png$",file_out)) {
   png(file_out, width     = 4.25,
                 height    = 2.25,
                 units     = "in",
                 res       = 1200,
                 pointsize = 4)
} else {
   stop('Unknown plotting format')
}

ggplot(df, aes(time, value)) +
    geom_point(aes(colour = series)) +
    geom_step(aes(colour = series)) +
    scale_x_continuous(limit = c(min(df$time), max(df$time)),
                        breaks = df$time) +
    scale_y_continuous(limit = c(0, 1)) +
    facet_grid(series ~ .) +
    xlab(label) +
    ylab("") +
    guides(colour=FALSE)

dev.off()