#
# Prints a stacked area plot for two samples
#

args<-commandArgs(TRUE)

require(ggplot2)

# data input

sample1_id <- args[1]
sample2_id <- args[2]
file_xy    <- args[3] # should contain cdr3nt, cdr3aa cols; two last columns should correspond to clonotype freq in sample pair
file_out   <- args[4]

# transform data

getcol_n = function(lst, col){
	apply(lst[col], 2, as.numeric)
}

getcol_c = function(lst, col){
	apply(lst[col], 2, as.character)
}

table <- read.delim(file_xy)
n <- nrow(table)
m <- ncol(table)

x1 <- getcol_n(table, m - 1)
x2 <- getcol_n(table, m)

# here we do change sample labels and further refer to them as s1 and s2
# otherwise R will re-order samples alphabetically by name ...

df <- data.frame(cdr3nt = rep(getcol_c(table, "cdr3nt"), 2),
                 cdr3aa = rep(getcol_c(table, "cdr3aa"), 2),
                 sample = rep(c("s1", "s2"), each = n),
                 expr   = c(x1, x2))

# trick to add clonotype labels with geom_text, as suggested at StackExchange

dd <- subset(df, sample == "s2")

#dd <- dd[with(dd, order(cdr3nt, levels(dd$cdr3aa))), ]
dd$cum <- cumsum(dd$expr)

# custom palette (color blind)

gs.pal <- colorRampPalette(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

# plotting

pdf(file_out)

ggplot(df, aes(x = sample, y = expr)) +
    geom_area(aes(fill = cdr3nt, group = cdr3nt), guide = "none", position="stack") +
    ylab("cumulative abundance") +
    xlab("") +
    scale_x_discrete(expand = c(0,0), labels=c(sample1_id, sample2_id)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = gs.pal(n)) +
    theme(legend.position = "none") +
    geom_text(data = dd, aes(x = "s2", y = cum, label = cdr3aa, size = expr), hjust = 1, vjust = 1)

dev.off()