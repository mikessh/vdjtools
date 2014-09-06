#
# Prints a stacked histogram
#

args<-commandArgs(TRUE)

require(ggplot2); require(reshape)

# data input

table      <- args[1]
file_out   <- args[2]

# transform data

df <- read.table(text = table, sep ="\t", header = TRUE)

# Oh really?
# We have to use it otherwise the table will never be treated as completely numeric
# Moreover, df <- apply(df, as.numeric) generates lots of fun
# Stop using R guys, this is pathetic
df[, 1:ncol(df)] <- apply(df[, 1:ncol(df)], 2, as.numeric)

df.m <- melt(df, id = "Len")

# trick to add clonotype labels with geom_text, as suggested at StackExchange

#dd <- subset(df, sample == "s2")

#dd <- dd[with(dd, order(cdr3nt, levels(dd$cdr3aa))), ]
#dd$cum <- cumsum(dd$expr)

# custom palette (color blind)

gs.pal <- colorRampPalette(c("#2b8cbe", "#e0f3db", "#fdbb84"))

# plotting

pdf(file_out)

ggplot(df.m, aes(x = Len, y = value, fill = variable)) +
  geom_bar(width = 1, stat = "identity") +
  xlab("CDR3 length, bp") +
  labs(fill="Clonotype") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=c("grey75", gs.pal(ncol(df) - 2))) +
  theme(legend.text=element_text(size=8), axis.title.y=element_blank())

dev.off()