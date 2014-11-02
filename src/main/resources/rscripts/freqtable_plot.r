args<-commandArgs(TRUE)

require(ggplot2);

# data input

file_in    <- args[1]
file_out   <- args[2]
require(ggplot2)

df <- read.table(file_in, sep ="\t", header = TRUE, comment = "")

# mandatory transform data
df[, ] <- apply(df[, ], 2, as.numeric)
df$clonotype_size_l[df$clonotype_size_l<1] <- 1
df <- df[1:(nrow(df)-1), ]

# For label plotting, source: http://goo.gl/K4yh

x <- log10(df$X.clonotype_size)
y <- log10(df$compl_cdf)
w <- df$number_of_clonotypes

dd <- data.frame(x = x, y = y, w = w)

m <- lm(y ~ x, dd, weights = w);
eq <- substitute(italic(y) == a %.% italic(x)^b*","~~italic(r)^2~"="~r2,
     list(a = format(10 ^ coef(m)[1], digits = 2),
          b = format(coef(m)[2], digits = 2),
         r2 = format(summary(m)$r.squared, digits = 3)))
lbl<-as.character(as.expression(eq))

#plot

pdf(file_out)

ggplot() +
    geom_ribbon(data = df, aes(x = compl_cdf, y = X.clonotype_size, ymin = clonotype_size_l, ymax = clonotype_size_u), alpha = 0.3, fill="blue") +
    geom_line(data = df, aes(x = compl_cdf, y = X.clonotype_size), color="blue") +
    scale_y_log10(expand = c(0,0), limits=c(1, max(df$X.clonotype_size))) + scale_x_log10(expand = c(0,0)) +
    theme_bw() + coord_flip() +
    xlab("1-CDF") + ylab("clonotype size") +
    geom_line(stat="smooth",data = df, aes(x = compl_cdf, y = X.clonotype_size, weight = number_of_clonotypes), method = 'lm', formula = y~x, se=FALSE, linetype ="dashed", size = 1.0, color="black") +
    geom_text(aes(x = max(df$compl_cdf), y = max(df$X.clonotype_size), label = lbl), hjust=1.1, vjust=1.2, parse = TRUE)

dev.off()