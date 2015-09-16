require(ggplot2);require(reshape)

df<-read.table("~/Documents/healthy70.txt", comment="",header=T)

df<-read.table("pwmgrid.MS.txt", comment="",header=T, sep="\t")
df.m<-melt(df, id.vars=c("X.v","len","pos","div","freq"))
df.m$value<-as.numeric(as.character(df.m$value))/log2(20)
df.s<-subset(df.m, X.v == "TRBV5-6")
ggplot(df.s, aes(x=pos, y=value, fill=variable))+geom_bar(stat="identity")+
scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0), limits=c(0,1))+
facet_grid(X.v ~ len, space="free_x", scales="free_x")+theme_bw()
require(ggplot2);require(reshape);

colors.aa <- c(R="#E60606",
K="#C64200",
Q="#FF6600",
N="#FF9900",
E="#FFCC00",
D="#FFCC99",
H="#FFFF99",
P="#FFFF00",
Y="#CCFFCC",
W="#CC99FF",
S="#CCFF99",
T="#00FF99",
G="#00FF00",
A="#CCFFFF",
M="#99CCFF",
C="#00FFFF",
F="#00CCFF",
L="#3366FF",
V="#0000FF",
I="#000080")


pdf("ms+hsct.pdf")
df<-read.table("pwmgrid.MS+HSCT.txt", comment="",header=T)
df.m<-melt(df, id.vars=c("X.v","len","pos","div","freq"))
df.m$value<-as.numeric(as.character(df.m$value))/log2(20)
df.m$variable<-factor(df.m$variable, levels=c("R","K","Q","N","E","D","H","P","Y","W","S","T","G","A","M","C","F","L","V","I"))

df.s<-subset(df.m, X.v %in% c("TRBV5-1","TRBV5-6","TRBV20-1"))
ggplot(df.s, aes(x=pos, y=value, group=pos, fill=factor(variable)))+
geom_bar(stat="identity")+
scale_y_continuous(expand=c(0,0), limits=c(0,1))+
scale_fill_manual(values=colors.aa)+
facet_grid(X.v ~ len, space="free_x", scales="free_x")#+theme_bw()
dev.off()

pdf("ms.pdf")
df<-read.table("pwmgrid.MS.txt", comment="",header=T)
df.m<-melt(df, id.vars=c("X.v","len","pos","div","freq"))
df.m$value<-as.numeric(as.character(df.m$value))/log2(20)
df.m$variable<-factor(df.m$variable, levels=c("R","K","Q","N","E","D","H","P","Y","W","S","T","G","A","M","C","F","L","V","I"))

df.s<-subset(df.m, X.v %in% c("TRBV5-1","TRBV5-6","TRBV20-1"))
ggplot(df.s, aes(x=pos, y=value, group=pos, fill=factor(variable)))+
geom_bar(stat="identity")+
scale_y_continuous(expand=c(0,0), limits=c(0,1))+
scale_fill_manual(values=colors.aa)+
facet_grid(X.v ~ len, space="free_x", scales="free_x")#+theme_bw()
dev.off()

pdf("c.pdf")
df<-read.table("pwmgrid.C.txt", comment="",header=T)
df.m<-melt(df, id.vars=c("X.v","len","pos","div","freq"))
df.m$value<-as.numeric(as.character(df.m$value))/log2(20)
df.m$variable<-factor(df.m$variable, levels=c("R","K","Q","N","E","D","H","P","Y","W","S","T","G","A","M","C","F","L","V","I"))

df.s<-subset(df.m, X.v %in% c("TRBV5-1","TRBV5-6","TRBV20-1"))
ggplot(df.s, aes(x=pos, y=value, group=pos, fill=factor(variable)))+
geom_bar(stat="identity")+
scale_y_continuous(expand=c(0,0), limits=c(0,1))+
scale_fill_manual(values=colors.aa)+
facet_grid(X.v ~ len, space="free_x", scales="free_x")#+theme_bw()
dev.off()

