#
# Prints a stacked histogram
#

args<-commandArgs(TRUE)

require(ggplot2); require(reshape); require(RColorBrewer)

# data input

file_in    <- args[1]
file_out   <- args[2]
label      <- args[3]
gradient   <- args[4]

# transform data

df <- read.table(file_in, sep ="\t", header = TRUE, comment="", quote="")

# Don't ask me why
df[, 1:ncol(df)] <- apply(df[, 1:ncol(df)], 2, as.numeric)

df.m <- melt(df, id = "Len")

# trick to add clonotype labels with geom_text, as suggested at StackExchange

#dd <- subset(df, sample == "s2")

#dd <- dd[with(dd, order(cdr3nt, levels(dd$cdr3aa))), ]
#dd$cum <- cumsum(dd$expr)

# custom palette (color blind)

if (gradient) {
   pal <- colorRampPalette(c("#2b8cbe", "#e0f3db", "#fdbb84"))(ncol(df) - 2)
} else {
   pal <- brewer.pal(ncol(df) - 2, "Paired")
}

# plotting

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

ggplot(df.m, aes(x = Len, y = value, fill = variable)) +
  geom_bar(width = 1, stat = "identity") +
  xlab("CDR3 length, bp") +
  labs(fill=label) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=c("grey75", pal)) +
  theme_bw() +
  theme(legend.text=element_text(size=8), axis.title.y=element_blank()) +
  guides(fill = guide_legend(reverse = TRUE))

dev.off()