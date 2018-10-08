# Load required packages -----------------------------------------------------
library(ggplot2)
library(reshape2)
library(plyr)
library(multcompView)
library(gridExtra)
library(grid)

# Load data table ------------------------------------------------------------
pas_data <- read.csv("pas-zotus-merged-with-meta.csv",
                     header =T,
                     na.strings = c("", " ", "NA"),
                     colClasses = c(rep("factor",11),
                                    rep("numeric", 154)))
colnames(pas_data)

# Data subgrouping -----------------------------------------------------------
# Split controls into test conditions.
har_ratio <- pas_data[grep("hartismeri ratios", pas_data$subgroup),]
ple_ratio <- pas_data[grep("plectid ratios", pas_data$subgroup),]
pen_ratio <- pas_data[grep("penetrans ratios", pas_data$subgroup),]
tsb_ratio <- pas_data[grep("tsb ratios", pas_data$subgroup),]

# Add total assembled read pair column.
pas_data$total <- rowSums(pas_data[13:165])

# Data re-formatting ---------------------------------------------------------
# For each group data must be converted to represent the total assembled 
# read pairs as a proportion of the total.

# Penetrans ratio samples.
# Create a reduced data frame with the sample and ZOTU columns only.
penetrans_ratio <- pen_ratio[, c(8,13:165)]
# Reduce the data frame further to only the plasmid input sequences.
penetrans_ratio <- penetrans_ratio[,c(1:4)]
# Transpose the ZOTU dataframe.
penetrans_ratio <- t(penetrans_ratio[2:4])
# Convert assembled read pair numbers to proportions of
# the total in each sample
penetrans_ratio <- scale(penetrans_ratio,
                         center = F,
                         scale = colSums(penetrans_ratio))
# Convert proportional table to long format. 
pen_melt <- melt(penetrans_ratio, id.vars = "sample")
# Give the variables sensible names. 
pen_melt <- rename(pen_melt, c("Var2" = "colid",
                         "Var1" = "zotu",
                         "value" = "proportion_of_sample_reads"))
# Re-merge proportional table with PCR metadata.
pen_ratio$colid <- rownames(pen_ratio)
pen_remerge <- merge(pen_melt, pen_ratio, by = "colid")
# Reduce to only the needed columns.
pen_melt <- pen_remerge[,1:15]
pen_melt$pcrid <- as.numeric(pen_melt$pcrid)

# Luffness ratio samples.
luffness_ratio <- tsb_ratio[, c(8,13:93)]
luffness_ratio <- luffness_ratio[,c(1:4)]
luffness_ratio <- t(luffness_ratio[2:4])
luffness_ratio <- scale(luffness_ratio,
                        center = F,
                        scale= colSums(luffness_ratio))
luff_melt <- melt(luffness_ratio, id.vars = "sample")
luff_melt <- rename(luff_melt, c("Var2" = "colid",
                                 "Var1" = "zotu",
                                 "value" = "proportion_of_sample_reads"))
tsb_ratio$colid <- rownames(tsb_ratio)
luff_remerge <- merge(luff_melt, tsb_ratio, by = "colid")
luff_melt <- luff_remerge[,1:15]
luff_melt$pcrid <- as.numeric(luff_melt$pcrid)

# Hartismeri ratio samples.
hartismeri_ratio <- hart_ratio[, c(8,13:93)]
hartismeri_ratio <- hartismeri_ratio[,c(1:4)]
hartismeri_ratio <- t(hartismeri_ratio[2:4])
hartismeri_ratio <- scale(hartismeri_ratio,
                          center = F,
                          scale= colSums(hartismeri_ratio))
hart_melt <- melt(hartismeri_ratio, id.vars = "sample")
hart_melt <- rename(hart_melt, c("Var2" = "colid",
                                 "Var1" = "zotu",
                                 "value" = "proportion_of_sample_reads"))
hart_ratio$colid <- rownames(hart_ratio)
hart_remerge <- merge(hart_melt, hart_ratio, by = "colid")
hart_melt <- hart_remerge[,1:15]
hart_melt$pcrid <- as.numeric(hart_melt$pcrid)

# Save the converted tables.
write.csv(pen_melt, "penetrans_plasmid_ratios.csv")
write.csv(luff_melt, "luffness_plasmid_ratios.csv")
write.csv(hart_melt, "hartismeri_plasmid_ratios.csv")

# Manual input of plot.no varibale which is the logical order of samples
# with a gap between replicate sets.
# TO DO: replace manual plot number input with automatic
# numbering from sample ID.

# Re-load ratio tables -------------------------------------------------------
pen_dat <- read.csv("penetrans_plasmid_ratios.csv")
print(levels(pen_dat$zotu))
tsb_dat <- read.csv("luffness_plasmid_ratios.csv")
print(levels(tsb_dat$zotu))
hart_dat <- read.csv("hartismeri_plasmid_ratios.csv")
print(levels(hart_dat$zotu))

# Plot ratio consistency -----------------------------------------------------
# Penetrans.
penratio <- ggplot(pen_dat, aes(x = plot.no,
                                y = proportion_of_sample_reads,
                                fill=zotu, order=zotu)) +
                   geom_bar(stat='identity')
penratio

ggsave("penetrans_ratio_plot.svg",
        dpi=600,
        width = 12,
        height = 9)

# Luffness.
tsbratio <- ggplot(tsb_dat, aes(x = plot.no,
                                y = proportion_of_sample_reads,
                                fill=zotu,
                                order=zotu)) + 
                   geom_bar(stat='identity')
tsbratio

ggsave("luffness_ratio_plot.svg",
        dpi=600,
        width = 12,
        height = 9)

# Hartismeri.
# Bar chart.
hartratio <- ggplot(hart_dat, aes(x = plot.no,
                                  y = proportion_of_sample_reads,
                                  fill=zotu, order=zotu)) + 
                    geom_bar(stat='identity')
hartratio

# Box plot.
hartratio <- ggplot(hart_dat, aes(x = sample,
                                  y = proportion_of_sample_reads,
                                  fill = ZOTU,
                                  order = zotu)) + 
                    geom_boxplot(aes(fill = zotu)) +
                    xlab("Sample") +
                    ylab("Proportion of assembled reads")
hartratio

ggsave("hartismeri_ratio_plot.svg",
       dpi=600,
       width = 8,
       height = 5)

# Seperate legend plot -------------------------------------------------------
# As each plot has the same three ZOTUs the legend is identical in each, 
# so the legend can be extracted from any one of the above plots. 

# Define a function to separate the legend from the main plot. 
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

# Extract the legend from the plot
legend <- g_legend(hartratio)
legend
g <- grid.arrange(arrangeGrob(legend))

# Save as scalable vector graphic.
ggsave("pas_ratio_legend.svg",
        g,
        dpi = 600,
        width = 4,
        height = 5)

