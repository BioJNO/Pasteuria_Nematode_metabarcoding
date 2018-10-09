# Load required packages -----------------------------------------------------
library(ggplot2)
library(reshape2)
library(plyr)
library(multcompView)

# Load data table ------------------------------------------------------------
pas_data <- read.csv("pas-zotus-merged-with-meta.csv",
                     header =T,
                     na.strings = c("", " ", "NA"),
                     colClasses = c(rep("factor",11),
                                    rep("numeric", 154)))

colnames(pas_data)

# Add total assembled read pair and perul columns.
pas_data$total <- rowSums(pas_data[13:165])
pas_data$perul <- pas_data$total/pas_data$vol

# Data subgrouping -----------------------------------------------------------
# Split controls into test conditions.
pen_det <- pas_data[grep("penetrans detection limits", pas_data$subgroup),]
har_det <- pas_data[grep("hartismeri detection limits", pas_data$subgroup),]

# Statistical analyses -------------------------------------------------------
# Look at the average number of reads per PCR band score.
aggregate(pas_data[, 167], list(pas_data$band), mean)

# Is the number of reads per microlitre statistically
# significant between band scores?
pasband_anova <- aov(perul ~ band, data = pas_data)
summary(pasband_anova)

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
cor.test((har_det$total),
         har_det$template,
         method = "spearman",
         conf.level = 0.95)

cor.test((har_det$total/har_det$vol),
          har_det$template,
          method = "spearman",
          conf.level = 0.95)

# Repeat for penetrans.          
cor.test((pen_det$total),
         pen_det$template,
         method = "spearman",
         conf.level = 0.95)
         
cor.test((pen_det$total/pen_det$vol),
         pen_det$template,
         method = "spearman",
         conf.level = 0.95)          

# Plot the results -----------------------------------------------------------
# Band score vs assembled read pairs
bandscore_box <- ggplot(data = pas_data, 
                        aes(band, perul)) +
                        geom_boxplot(aes(fill = band)) + 
                        xlab("PCR band score") + 
                        ylab("Assembled reads per ul") + 
                        geom_point(aes(colour = band)) + 
                        geom_smooth(method = "auto")
bandscore_box

# Save the plot as a scalable vector graphic.
ggsave("pas_bandscore_scatter.svg",
       dpi=600,
       width = 6,
       height = 4)

# Hartismeri copy number vs total assembled reads.
hartdet_totalplot <- ggplot(data = hart_det,
                        aes(log(template), (total))) + 
                        geom_boxplot(aes(fill = sample)) + 
                        xlab("Log template copy number") + 
                        ylab("Total assembled reads") + 
                        geom_smooth(method="auto") + 
                        geom_point(aes(colour = sample))        
hartdet_totalplot

ggsave("hartismeri_detection_scatter_curve_and_box.svg",
       dpi=600,
       width = 6,
       height = 4)

# Hartismeri copy number vs assembled reads per ul.
hartdet_perulplot <- ggplot(data = hart_det, aes(log(template), per.ul)) + 
                        geom_boxplot(aes(fill = sample)) + 
                        xlab("Log Template Copy Number") + 
                        ylab("Amplicons per ul") + 
                        geom_smooth(method="auto") + 
                        geom_point(aes(colour = sample))
hartdet_perulplot

ggsave("hartismeri_detection_scatter_curve_and_box_perul.svg",
       dpi=600,
       width = 6,
       height = 4)

# Penetrans copy number vs total assembled read pairs.
pen_det_box <- ggplot(data = pen_det, aes(log(template), total)) +
                      geom_boxplot(aes(fill = sample)) +
                      xlab("Log Template Copy Number") +
                      ylab("Total Amplicons") +
                      geom_smooth(method="auto") +
                      geom_point(aes(colour = sample))
pen_det_box

ggsave("penetrans_detection__scatter_curve_and_box_total.svg",
       dpi= 600,
       width = 6,
       height = 4)

# Penetrans copy number vs assembled read pairs perul.
pen_det_box	<- ggplot(data = pen_det, aes(log(template), (total/vol))) +
                      geom_boxplot(aes(fill = sample)) +
                      xlab("Log Template Copy Number") +
                      ylab("Amplicons per ul") +
                      geom_smooth(method="auto") +
                      geom_point(aes(colour = sample))

pen_det_box

ggsave("penetrans_detection__scatter_curve_and_box_perul.svg",
        dpi= 600,
        width = 6,
        height = 4)
