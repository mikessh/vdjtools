require(ggplot2);require(reshape)

df<-read.table("~/Documents/healthy70.txt", comment="",header=T)

df<-read.table("pwmgrid.MS.txt", comment="",header=T)
df.m<-melt(df, id.vars=c("X.v","len","pos","div","freq"))
df.m$value<-as.numeric(as.character(df.m$value))/log2(20)
df.s<-subset(df.m, X.v == "TRBV5-6")
ggplot(df.s, aes(x=pos, y=value, fill=variable))+geom_bar(stat="identity")+
scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0), limits=c(0,1))+
facet_grid(X.v ~ len, space="free_x", scales="free_x")+theme_bw()