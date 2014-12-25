args<-commandArgs(TRUE)

require(reshape); require(ggplot2); require(gridExtra); require(FField)

file_in  <- args[1] #"rarefaction.strict.txt"
lbl_col  <- as.integer(args[2]) #5
fac_col  <- as.integer(args[3]) #6
num_fac  <- as.logical(args[4]) #"F"
add_lbl  <- as.logical(args[5]) #"T"
wide     <- as.logical(args[6]) #"F"
file_out <- args[7] # "rarefaction.strict.pdf"

df <- read.table(file_in,header=T,comment="",sep="\t")

# get name of label and factor
fac_name<-names(df)[fac_col]

# make sure that we parse numeric factor as needed
if (num_fac) {
   df[,fac_col] <- as.numeric(as.character(df[,fac_col]))
}

# rename columns so we can access them
# this is done as some column names could be quite complex to pass by command line,
# e.g. #1_label will be converted to X1.label, etc
# we also protect from column name collisions
df <- data.frame(dataset = df[,1],
                 fac  =  df[,fac_col],
                 lbl  =  df[,lbl_col],
                 type =  df[,ncol(df)],
                 ciU  =  df[,ncol(df)-1],
                 ciL  =  df[,ncol(df)-2],
                 mean =  df[,ncol(df)-3],
                 x    =  df[,ncol(df)-4])

# those are the last points with observed diversity
# we'll highlight them with point and (if required) a label
df.p <- subset(df, type==2 & x==max(df$x))

if (add_lbl) {
   df.l <- df.p
   # We absolutely sure need this trick as labels are going to overlap leading to huge mess..
   x.fact <- 100 / max(df.l$x)
   y.fact <- 100 / max(df.l$mean)
   coords <- FFieldPtRep(coords = cbind(df.l$x * x.fact, df.l$mean * y.fact), iter.max = 1000)
   df.l$x <- coords$x / x.fact
   df.l$mean <- coords$y / y.fact
} else {
   df.l <- df.l[0, ]
}

# interpolated & extrapolated
df.i <- subset(df, type==0)
df.o <- subset(df, type==1)
df.e <- subset(df, type==2)

g <- ggplot(df, aes(x=x, y=mean, group=dataset)) +
     geom_point(data=df.o, aes(colour=fac)) +
     geom_line(data=df.i, aes(colour=fac), linetype="solid") +
     geom_line(data=df.e, aes(colour=fac), linetype="dashed") +
     geom_ribbon(aes(ymin=ciL, ymax=ciU), alpha=0.2) +
	 xlab("Sample size") + ylab("Diveristy") +
	 labs(colour=fac_name) +
	 scale_x_continuous(oob=scales::rescale_none) +
	 scale_y_continuous(oob=scales::rescale_none) +
	 theme_bw()

# add corresponding fancy axis
if (num_fac) {
   g <- g + scale_colour_gradient2(low="#feb24c", mid="#31a354", high="#2b8cbe", midpoint=(max(df$fac) + min(df$fac))/2)
} else {
   g <- g + scale_colour_brewer(palette="Set2")
}

g <- g + geom_text(data=df.l, aes(label=lbl), color="black", fontface = "plain", hjust=1.0, vjust=0.5, cex=2.5)

pdf(file_out)

if (wide) {
   # wide plot layout
   grid.arrange(g, ncol=1, nrow=2)
} else {
   g
}

dev.off()