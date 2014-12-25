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

pdf(file_out)

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