require(ggplot2)

args <- commandArgs(T)

file_in  <- args[1] #"cdr3aa.profile.wt.txt"
file_out <- args[2] #"cdr3aa.profile.wt.pdf"
group_id <- as.integer(args[3]) #4
normalized <- as.logical(args[4]) #F

if (group_id < 1) {
   group_id = 1 # use sample id if not specified
}

df <- read.table(file_in, header=T, comment="", quote="", sep="\t", stringsAsFactors = F)

df$value <- as.numeric(df$value)
df$total <- as.numeric(df$total)

if (normalized) {
    df$value <- ifelse(df$total == 0, 0, df$value / df$total)
}

# collect required columns and select grouping column
groupName <- colnames(df)[group_id]

df <- data.frame(bin = factor(df$bin + 1), 
                 value = df$value,
                 property = df$property,
                 cdr3.segment = df$cdr3_segment,
                 group = factor(df[,group_id]))


# set order of segments
df$cdr3.segment <- factor(df$cdr3.segment, levels = c("CDR3-full",
  "V-germ",
  "VD-junc", "D-germ", "DJ-junc", "VJ-junc",
  "J-germ"))

get_facet_formula = function() {
  if (length(unique(df$property)) > 1) { 
    property ~ cdr3.segment
  } else {
    cdr3.segment ~ property
  }
}

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

ggplot(df, aes(x=bin, y=value, color=group)) +
  geom_boxplot() + 
  facet_grid(get_facet_formula(), scales="free") +
  xlab("") + ylab("") +
  theme_bw() + scale_color_brewer(groupName, palette="Set2")

dev.off()


# you can further proceed with T-tests if you want..
#t.test(freq ~ genotype, subset(df, cdr3.segment == "VJ-junc" & property == "disorder"))