args<-commandArgs(TRUE)

require(reshape); require(ggplot2)

file_in  <- args[1]
lbl_col  <- as.integer(args[2])#"1"
fac_col  <- as.integer(args[3])#"3"
num_fac  <- as.logical(args[4])#"T"
add_lbl  <- as.logical(args[5])#"T"
file_out <- args[6]

df<-read.table(file_in, comment="",header=F,sep="\t",stringsAsFactor=F)

# get name of label and factor
lbl_name<-df[1,lbl_col]
fac_name<-df[1,fac_col]

nums <- sapply(df, is.numeric) # sampling points
x<-as.integer(df[1,nums])
df<-df[2:nrow(df),]
m<-nrow(df)

df.m <- melt(df, id.vars=paste("V",(1:ncol(df))[!nums],sep=""))
df.m$variable <- rep(x,1,each=m)

df.m <- data.frame(value = df.m$value, variable = df.m$variable,
                   dataset = df.m[,1], fac = df.m[,fac_col], lbl = df.m[,lbl_col])

if (num_fac) {
   df.m$fac <- as.numeric(as.character(df.m$fac))
}

# selecting last points only using a skipping variable
df.m[ ,"sk"] <- NA
for (d in unique(df.m$dataset)) {
    ind <- df.m[,"dataset"]==d
    df.m[ind,"sk"][which.max(df.m[ind,"variable"]+df.m[ind,"value"])] <- 1
}

g<-ggplot(df.m,aes(x=variable, y=value, colour=fac, group=dataset)) +
     geom_point(aes(x = sk * variable)) +
     geom_smooth(size = 1) +
	 xlab("Sample size") + ylab("Diveristy") +
	 labs(colour=fac_name) +
	 scale_x_continuous(expand = c(0, 0)) +
	 scale_y_continuous(expand = c(0, 0)) +
	 theme_bw()

if (add_lbl) {
   g <- g + geom_text(aes(x = sk * variable, label=lbl), color="black", vjust=1.0, hjust=1.0)
}

pdf(file_out)

if (num_fac) {
   g + scale_colour_gradient2(low="#feb24c", mid="#31a354", high="#2b8cbe", midpoint=(max(df.m$fac) + min(df.m$fac))/2)
} else {
   g
}

dev.off()