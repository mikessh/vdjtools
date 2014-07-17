args<-commandArgs(TRUE)

sample1_id<-args[1]
sample2_id<-args[2]
file_xy<-args[3]
file_xx<-args[4]
file_yy<-args[5]
file_out<-args[6]

getcol=function(lst, col){
	apply(lst[col], 2, as.numeric)
}

xy<-read.delim(file_xy)
xy<-data.frame(getcol(xy,"x"),getcol(xy,"y"))

n<-dim(xy)[1]

xx<-getcol(read.delim(file_xx), "xx")
xx<-sample(xx, n)
xx<-data.frame(xx)

yy<-getcol(read.delim(file_yy), "yy")
yy<-sample(yy, n)
yy<-data.frame(yy)

require(ggplot2); require(grid); require(gridExtra)

lm_eqn = function(df){
    m = lm(y ~ x, df)
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(
              #a = format(coef(m)[1], digits = 2), 
              #b = format(coef(m)[2], digits = 2), 
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq))                 
}

pdf(file_out)

#placeholder plot - prints nothing at all
empty <- ggplot() +
     geom_point(aes(1,1), colour="white") +
     theme(                              
       plot.background = element_blank(), 
       panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank(), 
       panel.border = element_blank(), 
       panel.background = element_blank(),
       axis.title.x = element_blank(),
       axis.title.y = element_blank(),
       axis.text.x = element_blank(),
       axis.text.y = element_blank(),
       axis.ticks = element_blank()
     )

#scatterplot of x and y variables
scatter <- ggplot(xy, aes(x, y, size=(x+y)/2, weight=10^(x+y))) + 
  geom_point(alpha = 0.5) +
  #geom_smooth(method = "lm", formula = y ~ x) +
  scale_x_continuous(limit=c(-7, 0), breaks=-7:0) +     
  scale_y_continuous(limit=c(-7, 0), breaks=-7:0) + 
  scale_size_continuous(guide="none", range = c(0.1, 5)) +
  xlab(sample1_id) +
  ylab(sample2_id) +  
  theme(legend.position=c(1,1),legend.justification=c(1,1)) 

#marginal density of x - plot on top
plot_top <- ggplot(xx, aes(xx, weight=10^xx)) +  
stat_density(aes(ymax = ..density..,  ymin = -..density..),
    fill = "grey50", colour = "grey25",
    geom = "ribbon", position = "identity",
    adjust=5) + 
  scale_x_continuous(limit=c(-7, 0), breaks=-7:0) +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())

#marginal density of y - plot on the right
plot_right <- ggplot(yy, aes(yy, weight=10^yy)) + 
  stat_density(aes(ymax = ..density..,  ymin = -..density..),
    fill = "grey50", colour = "grey25",
    geom = "ribbon", position = "identity",
    adjust=5) + 
  scale_x_continuous(limit=c(-7, 0), breaks=-7:0) +
 coord_flip() + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) 

#arrange the plots together, with appropriate height and width for each row and column
grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

dev.off()