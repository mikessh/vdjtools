require(gplots); require(RColorBrewer); require(ggplot2)

args<-commandArgs(TRUE)

input_file  <- args[1] #"ms.segments.V2.txt"
v_col_count <- args[2] #48
lbl_col     <- args[3] #1
fac_col     <- args[4] #2
num_fac     <- args[5]
output_file <- args[6]

v_col_count <- as.numeric(v_col_count)
lbl_col <- as.numeric(lbl_col)
fac_col <- as.numeric(fac_col)
num_fac <- num_fac == "TRUE"

if (lbl_col < 1) {
   lbl_col = 1
}

df<-read.table(input_file, header=T, sep="\t",comment="")
#df$age<-NULL
vcols<-(ncol(df)-v_col_count + 1):ncol(df)
df[, vcols] <- apply(df[, vcols], 2, as.numeric)
x <- as.matrix(df[, vcols])

if (fac_col > 0) {
    if (num_fac) {
       df[,fac_col] <- as.numeric(df[,fac_col])
       scol <- colorRampPalette(c("#feb24c", "#31a354", "#2b8cbe"))(10)[as.factor(as.numeric(cut(df[,fac_col],10)))]
    } else {
       df[, fac_col] <- as.factor(df[, fac_col])
       nlevels <- length(levels(df[, fac_col]))
       scol <- colorRampPalette(brewer.pal(nlevels, "Set1"))(nlevels)[df[,fac_col]]
    }
}

if (fac_col > 0) {
pdf(output_file)
heatmap.2(t(x),  na.rm=F, labCol = df[,lbl_col],
 na.col="grey50",
 col=colorRampPalette(c("#2f98ce", "#e0f3db", "#f47104")),
 density.info="none", trace="none",scale="column",
 ColSideColors=scol)
dev.off()
} else {
pdf(output_file)
heatmap.2(t(x),  na.rm=F, labCol = df[,lbl_col],
 na.col="grey50",
 col=colorRampPalette(c("#2f98ce", "#e0f3db", "#f47104")),
 density.info="none", trace="none",scale="column")
dev.off()
}