args<-commandArgs(TRUE)

require(reshape); require(ggplot2); require(gridExtra); require(FField)

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

# last point
count_col<-which(df[1,]=="count")
div_col<-which(df[1,]=="diversity")

# sampling points
nums <- sapply(df, is.numeric) # sampling points
x<-as.integer(df[1,nums])
df<-df[2:nrow(df),]
m<-nrow(df)

# melt data

## firsty select last points
df.L <- data.frame(value = as.integer(df[, div_col]), variable = as.integer(df[, count_col]),
                   dataset = df[,1], fac = df[,fac_col], lbl = df[,lbl_col])
df.L[ ,"sk"] <- 1
df.L[ ,"sk1"] <- 1e100

## select sampling points
df.m <- melt(df, id.vars=paste("V",(1:ncol(df))[!nums],sep=""))
df.m$variable <- rep(x,1,each=m)

df.m <- data.frame(value = df.m$value, variable = df.m$variable,
                   dataset = df.m[,1], fac = df.m[,fac_col], lbl = df.m[,lbl_col])
df.m[ ,"sk"] <- NA
df.m[ ,"sk1"] <- -1e100

df.m <- rbind(df.m, df.L)

if (num_fac) {
   df.m$fac <- as.numeric(as.character(df.m$fac))
}

if (add_lbl) {
   x.fact <- 100 / max(df.m$variable)
   yy <- df.m$value
   yy[which(is.na(yy))] <- -1e100
   y.fact <- 100 / max(yy)
   coords <- FFieldPtRep(coords = cbind(pmin(df.m$variable, df.m$sk1) * x.fact, pmin(yy, df.m$sk1) * y.fact))
   df.m <- cbind(df.m, data.frame(lbl.x = df.m$sk * coords$x / x.fact, lbl.y = df.m$sk * coords$y / y.fact))
}

# selecting last points only using a skipping variable
#df.m[ ,"sk"] <- NA
#for (d in unique(df.m$dataset)) {
#    ind <- df.m[,"dataset"]==d
#    df.m[ind,"sk"][which.max(df.m[ind,"variable"]+df.m[ind,"value"])] <- 1
#}

g<-ggplot(df.m,aes(x=variable, y=value, colour=fac, group=dataset)) +
     geom_point(aes(x = sk * variable)) +
     stat_smooth(level = 0.99, fullrange = TRUE) +
	 xlab("Sample size") + ylab("Diveristy") +
	 labs(colour=fac_name) +
	 scale_x_continuous(expand = c(0, 0)) +
	 scale_y_continuous(expand = c(0, 0)) +
	 theme_bw()

if (add_lbl) {
   g <- g + geom_text(aes(x = lbl.x, y = lbl.y, label=lbl), color="black", fontface = "plain", vjust=1.0, hjust=1.0, cex = 5)
}

pdf(file_out)

if (num_fac) {
   g <- g + scale_colour_gradient2(low="#feb24c", mid="#31a354", high="#2b8cbe", midpoint=(max(df.m$fac) + min(df.m$fac))/2)
   grid.arrange(g, ncol=1, nrow=2)
} else {
   g <- g + scale_colour_brewer(palette="Set1")
   grid.arrange(g, ncol=1, nrow=2)
}

dev.off()