#
# Prints a scatterplot with marginal distributions
# Note that while scatterplot is based only on vector of species present in sample intersection,
# marginal plots are built based on complete datasets
#
# Adapted from R-bloggers (http://www.r-bloggers.com/ggplot2-cheatsheet-for-visualizing-distributions/)
# With minor additions based on various StackOverflow replies
#

require(ggplot2); require(grid); require(gridExtra)

# data input

args<-commandArgs(TRUE)

sample1_id <- args[1]
sample2_id <- args[2]
file_xy    <- args[3]
file_xx    <- args[4]
file_yy    <- args[5]
file_out   <- args[6]

# transform data

getcol = function(lst, col){
	apply(lst[col], 2, as.numeric)
}

xy <- read.delim(file_xy)
xy <- data.frame(getcol(xy,"x"),getcol(xy,"y"))

xx <- getcol(read.delim(file_xx), "xx")
xx <- data.frame(xx)

yy <- getcol(read.delim(file_yy), "yy")
yy <- data.frame(yy)

# plotting
 
pdf(file_out)

# placeholder (top right)

empty <- ggplot() +
     geom_point(aes(1,1), colour="white") +
     theme(                              
       plot.background  = element_blank(), 
       panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank(), 
       panel.border     = element_blank(), 
       panel.background = element_blank(),
       axis.title.x     = element_blank(),
       axis.title.y     = element_blank(),
       axis.text.x      = element_blank(),
       axis.text.y      = element_blank(),
       axis.ticks       = element_blank()
       #plot.margin      = unit(c(3, -5.5, 4, 3), "mm")
     )

# scatterplot

scatter <- ggplot(xy, aes(x, y, size = (x + y) / 2)) + 
  theme_bw() +
  geom_point(
     fill   = "grey50", 
     colour = "black",
     alpha  = 0.3, 
     pch    = 21
     ) +     
  scale_x_continuous(limit = c(-7, 0), breaks = -7:0) +     
  scale_y_continuous(limit = c(-7, 0), breaks = -7:0) + 
  scale_size_continuous(guide = "none", range = c(0.1, 5)) +
  xlab(sample1_id) +
  ylab(sample2_id) +  
  theme(legend.position = c(1, 1), legend.justification = c(1, 1))

# marginal density of x

plot_top <- ggplot(xx, aes(xx, weight=10^xx)) +  
  theme_bw() +
  stat_density(
       aes(ymax = ..density..,  ymin = -..density..),
       fill     = "grey50", 
       colour   = "grey25",      
       alpha    = 0.3, 
       position = "identity",
       adjust   = 5
     ) + 
  ylab("") +
  scale_x_continuous(limit = c(-7, 0), breaks = -7:0) +
  scale_y_continuous(limit = c(0, 4.5), breaks = c(0, 1.5, 3, 4.5)) +
  theme(legend.position = "none", axis.title.x = element_blank())


# marginal density of y

plot_right <- ggplot(yy, aes(yy, weight=10^yy)) + 
  theme_bw() +
  stat_density(
       aes(ymax = ..density..,  ymin = -..density..),
       fill     = "grey50",
       colour   = "grey25",   
       alpha    = 0.3, 
       position = "identity",
       adjust   = 5
     ) + 
  scale_x_continuous(limit = c(-7, 0), breaks = -7:0) +
  scale_y_continuous(limit = c(0, 4.5), breaks = c(0, 1.5, 3, 4.5)) +
  coord_flip() + 
  ylab("") +
  theme(legend.position = "none", axis.title.y = element_blank()) 

# arrange the plots together, with appropriate height and width for each row and column

grid.arrange(plot_top, empty, scatter, plot_right, ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

dev.off()