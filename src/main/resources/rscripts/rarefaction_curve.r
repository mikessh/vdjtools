args<-commandArgs(TRUE)

datasets <- args[1] # "A4-i127\tA4-i127\tA4-i127\tA4-i118\tA4-i118\tA4-i118\tA4-i122\tA4-i122\tA4-i122\tA3-i145\tA3-i145\tA3-i145\tA3-i150\tA3-i150\tA3-i150"
points   <- args[2] # "0\t100000\t200000\t300000\t400000\t500000\t600000\t700000\t800000\t900000\t1000000"
input    <- args[3] # "4.rarefaction_aa.txt"
file_out <- args[4]

require(ggplot2); require(reshape)

# read data & do some pre-processing
df <- read.table(input, sep ="\t")
labels <- sapply(read.table(text=datasets, sep=","), as.character)
x <- apply(as.vector(read.table(text=points, sep=",")), 1, as.numeric)

# strip unnescessary columns and conver to numeric
df.s <- df[,(ncol(df) - length(x) + 1):ncol(df)]
df.s[, 1:ncol(df.s)] <- apply(df.s[, 1:ncol(df.s)], 2, as.numeric)

# format to a suitable data frame
df.f <- data.frame(x=x, series=t(df.s))
names(df.f) <- c("ss", labels)
df.m <- melt(df.f, id="ss")

# plot

pdf(file_out)

ggplot(df.m,aes(x=ss, y=value, colour=variable, group=variable)) +
	xlab("Sample size") +
	ylab("Diveristy") +
	labs(colour="Dataset") +
	scale_x_continuous(expand = c(0, 0)) +
	scale_y_continuous(expand = c(0, 0)) +
	geom_smooth()

dev.off()