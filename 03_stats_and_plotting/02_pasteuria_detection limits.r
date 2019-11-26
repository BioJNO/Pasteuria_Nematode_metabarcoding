# Load required packages -----------------------------------------------------
library(ggplot2)
library(reshape2)
library(plyr)
library(multcompView)
library(RColorBrewer)
library(ggsci)

# Load data table ------------------------------------------------------------
pas_data <- read.csv("pas-zotus-merged-with-meta-minabd-10.csv",
                     header =T,
                     na.strings = c("", " ", "NA"),
                     colClasses = c(rep("factor",12),
                                    rep("numeric", 77)))

colnames(pas_data)

# Data subgrouping -----------------------------------------------------------
# Split controls into test conditions.
pen_det <- pas_data[grep("penetrans detection limits", pas_data$subgroup),]
har_det <- pas_data[grep("hartismeri detection limits", pas_data$subgroup),]

# Statistical analyses -------------------------------------------------------
# Look at the average number of reads per ul for each PCR band score.
aggregate(pas_data[, 89], list(pas_data$band), mean)

# Is the number of reads per microlitre statistically
# significant between band scores?
pasband_anova <- aov(perul ~ band, data = pas_data)
summary(pasband_anova)
# Yes

# Is the number of reads per microlitre statistically
# significant between ALL band scores?
pasband_tukey <- TukeyHSD(pasband_anova)

# Define a function to extract labels from TukeyHSD test.
generate_label_df <- function(tukey, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- tukey[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  # Need to put the labels in the same order as in the boxplot:
  Tukey.labels$perul=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$perul) , ]
  return(Tukey.labels)
}

generate_label_df(pasband_tukey, "band")
# Yes

# Get the Spearmans rank correlation of band score
# to assembled read pairs per ul.
pas_data$band <- as.numeric(pas_data$band)
cor.test(pas_data$band,
         pas_data$perul,
         method = "spearman",
         conf.level = 0.95)

# Remove string (text) from sample ID to give numeric copy number value.
har_det$template <- gsub("hartismeri ", "", har_det$sample)
pen_det$template <- gsub("penetrans ", "", pen_det$sample)
# Set the column class as numeric.
har_det$template <- as.numeric(har_det$template)
pen_det$template <- as.numeric(pen_det$template)

# Get the Spearmans rank correlation of hartismeri template copy number and
# total assembled read pairs/assembled pairs per ul.
cor.test(har_det$abundance,
         har_det$template,
         method = "spearman",
         conf.level = 0.95)

cor.test(har_det$perul,
          har_det$template,
          method = "spearman",
          conf.level = 0.95)

# Repeat for penetrans.          
cor.test(pen_det$abundance,
         pen_det$template,
         method = "spearman",
         conf.level = 0.95)
         
cor.test(pen_det$perul,
         pen_det$template,
         method = "spearman",
         conf.level = 0.95)          

# Plot the results -----------------------------------------------------------
# Band score vs assembled read pairs
pas_data$band <- as.numeric(pas_data$band)
bandscore_box <- ggplot(data = pas_data, 
                        aes(band, perul, group=band)) +
                        geom_boxplot(aes(fill = band)) + 
                        xlab("PCR band score") + 
                        ylab("Assembled reads per ul") + 
                        geom_point(aes(colour = band)) + 
                        geom_smooth(method = "auto")
bandscore_box

# Save the plot as a scalable vector graphic.
#ggsave("pas_bandscore_scatter.svg",
#       dpi=600,
#       width = 6,
#       height = 4)

# Hartismeri copy number vs assembled reads per ul.
hartdet_perulplot <- ggplot(data = har_det, aes(log(template), perul)) + 
                        geom_boxplot(aes(fill = sample)) + 
                        xlab("Log Template Copy Number") + 
                        ylab("Amplicons per ul") + 
                        geom_smooth(method="auto") + 
                        geom_point(aes(colour = sample)) +
                        scale_fill_brewer(palette = "PuOr") +
                        scale_color_brewer(palette = "PuOr")
hartdet_perulplot

ggsave("Figures/Fig2B_hartismeri_detection_scatter_curve_and_box_perul.svg",
       dpi=600,
       width = 6,
       height = 4)

# Penetrans copy number vs assembled read pairs perul.
pen_det_box	<- ggplot(data = pen_det, aes(log(template), perul)) +
                      geom_boxplot(aes(fill = sample)) +
                      xlab("Log Template Copy Number") +
                      ylab("Amplicons per ul") +
                      geom_smooth(method="auto") +
                      geom_point(aes(colour = sample)) +
                      scale_fill_brewer(palette = "BrBG") +
                      scale_color_brewer(palette = "BrBG")
pen_det_box

ggsave("Figures/Fig2A_penetrans_detection_scatter_curve_and_box_perul.svg",
        dpi= 600,
        width = 6,
        height = 4)
