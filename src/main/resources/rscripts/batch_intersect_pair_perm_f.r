require(ggplot2)

args<-commandArgs(TRUE)

file_in  <- args[1]
file_out <- args[2]

df <- read.table(file_in, sep="\t", comment="", header=T, stringsAsFactors=F)
df[,3:5]<-apply(df[,3:5],2,as.numeric)

df.l<-data.frame(factor=df$factor,type=df$type,obs=df$obs,p=df$p)
df.l<-df.l[!duplicated(df.l), ]

g<-ggplot(df,aes(x=perm))+geom_histogram(colour = "darkgrey", fill = "grey") +
geom_vline(aes(xintercept = obs), colour="red", size=1, linetype = "longdash") +
theme_bw() + xlab("") + scale_y_continuous(name="permutations") + facet_grid(type~factor)

df.l$perm<-rep(max(ggplot_build(g)$data[[1]]$count),nrow(df.l))

pdf(file_out)

g+geom_text(data = df.l, aes(x=obs, label=as.character(p), y=perm), hjust=-0.1)

dev.off()