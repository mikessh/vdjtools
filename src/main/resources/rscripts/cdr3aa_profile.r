require(ggplot2)

args <- commandArgs(T)

file_in  <- args[1] 
file_out <- args[2] 
group_id <- as.integer(args[3])

if (group_id < 1) {
   group_id = 1 # use sample id if not specified
}

df <- read.table(file_in, header=T, comment="", quote="", sep="\t", stringsAsFactors = F)

df$value <- as.numeric(df$value)
df$total <- as.numeric(df$total)

# collect required columns and select grouping column
groupName <- colnames(df)[group_id]

df <- data.frame(bin = factor(df$bin + 1), 
                 freq = ifelse(df$total == 0, 0, df$value / df$total),
                 property = df$property,
                 cdr3.segment = df$cdr3.segment,
                 group = factor(df[,group_id]))

# set order of segments
df$cdr3.segment <- factor(df$cdr3.segment, levels = c("CDR3-full",
  "V-germ",
  "VD-junc", "D-germ", "DJ-junc", "VJ-junc",
  "J-germ"))


if (grepl("\\.pdf$",file_out)){
   pdf(file_out)
} else if (grepl("\\.png$",file_out)) {
   png(file_out, width     = 3.25,
                 height    = 3.25,
                 units     = "in",
                 res       = 1200,
                 pointsize = 4)
} else {
   stop('Unknown plotting format')
}

ggplot(df, aes(x=bin, y=freq, color=group)) +
  geom_boxplot() + facet_grid(property ~ cdr3.segment, scales="free") +
  xlab("") + ylab("") +
  theme_bw() + scale_color_brewer(groupName, palette="Set1")

dev.off()


# you can further proceed with T-tests if you want..
#t.test(freq ~ genotype, subset(df, cdr3.segment == "VJ-junc" & property == "disorder"))