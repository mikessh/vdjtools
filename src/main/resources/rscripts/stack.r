# data input

args <- commandArgs(TRUE)

label    <- args[1] #"time since HSCT, months"
points   <- args[2] #"-48;-0.5;4;10;25;37"
file_in  <- args[3] #"luc_table_collapsed.txt"
file_out <- args[4] #"out.pdf"

# data processing

x <- apply(as.vector(read.table(text=points, sep=";")),1,as.numeric)
table <- read.delim(file_in)
n <- length(x)
nn <- ncol(table)
mm <- nrow(table)
cols <- (nn-n+1):nn

y <- t(apply(as.matrix(table[-1,][1:(mm-2),cols]),1:2,as.numeric))

#xx <- x
#x <- unique(c(x, seq(from = min(x), to = max(x), by = (max(x)-min(x)) / 30)))
#y <- apply(y, 2, function(z) {
#       ss <- smooth.spline(xx, z, spar=0.05)
#       p <- predict(ss, x)$y       
#    }) 
#y <- apply(y, 1:2, function(z) { 
#    max(1e-7, z) + runif(1,0,1e-6) 
#    }) 


#Implementation from RBloggers
#
#plot.stacked makes a stacked plot where each y series is plotted on top
#of the each other using filled polygons
#
#Arguments include:
#'x' - a vector of values
#'y' - a matrix of data series (columns) corresponding to x
#'order.method' = c("as.is", "max", "first") 
#  "as.is" - plot in order of y column
#  "max" - plot in order of when each y series reaches maximum value
#  "first" - plot in order of when each y series first value > 0
#'...' - other plot arguments
 
plot.stacked <- function(
	x, y, 
	order.method = "as.is",
	ylab="", xlab="",
	ylim=NULL,
	...
){
    ord <- order(apply(y, 2, function(r) min(which(r>0))))
    y2 <- y[, ord]
    pal <- colorRampPalette(c("blue", "cyan", "yellow", "red"))
    col <- pal(ncol(y2))
 
	col <- as.vector(matrix(col, nrow=ncol(y), ncol=1))
	lwd <- as.vector(matrix(1, nrow=ncol(y), ncol=1))
 
	if(order.method == "max") {
		ord <- order(apply(y, 2, which.max))
		y <- y[, ord]
		col <- col[ord]
	}
 
	if(order.method == "first") {
		ord <- order(apply(y, 2, function(r) min(which(r>0))))
		y <- y[, ord]
		col <- col[ord]
	}
 
	top.old <- x*0
	polys <- vector(mode="list", ncol(y))
	for(i in seq(polys)){
		top.new <- top.old + y[,i]
		polys[[i]] <- list(x=c(x, rev(x)), y=c(top.old, rev(top.new)))
		top.old <- top.new
	}
 
	if(is.null(ylim)) ylim <- range(sapply(polys, function(x) range(x$y, na.rm=TRUE)), na.rm=TRUE)
	plot(x,y[,1], ylab=ylab, xlab=xlab, ylim=ylim, t="n", ...)
	for(i in seq(polys)){
		polygon(polys[[i]],  col=col[i], lwd=lwd[i], border=NA)
	}
 
}

# draw

pdf(file_out)

plot.stacked(x, y, xlab=label, ylab="abundance")

dev.off()