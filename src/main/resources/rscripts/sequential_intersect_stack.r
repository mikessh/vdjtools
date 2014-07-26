# data input

args <- commandArgs(TRUE)

require(ggplot2); require(reshape)

label    <- args[1] #"time since HSCT, months"
points   <- args[2] #"-48;-0.5;4;10;25;37;72"
file_in  <- args[3] #"luc_table_collapsed.txt"
file_out <- args[4] #"out.pdf"

# load time points, create some auxillary variables
x <- apply(as.vector(read.table(text=points, sep=";")), 1, as.numeric)
xcols <- (ncol(df)-n+1):ncol(df)
fcols <- 1:(ncol(df)-n)
xlbls <- colnames(df)[xcols]
n <- length(x)

# load data
df <- data.frame(read.delim(file_in))

# convert abundance columns to numeric
df[, xcols] <- apply(df[, xcols], 2, as.numeric)

# reshape data
df.m <- melt(df, id = fcols)

# replace sample ids (factor) by time (numeric)
ind <- match(df.m$variable, xlbls)
df.m$variable <- x[ind]

pal <- colorRampPalette(c("blue", "white", "red"))

# plot

pdf(file_out)

ggplot(df.m, aes(variable, value, group = cdr3nt, fill = factor(peak))) +
   geom_area(colour = "black", position = 'stack', size = 0.1) +
   scale_fill_manual(
      name   = "Peak position",
      breaks = c(-2, -1, 0:(n-1)),
      labels = c("Non-overlapping", "Not-shown", x),
      values = c("grey20", "grey50", pal(n))
      ) +
   scale_x_continuous(expand = c(0,0), limit = c(min(x), max(x)), breaks = x) +
   scale_y_continuous(expand = c(0,0)) +
   xlab(label) +
   ylab("abundance") +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dev.off()