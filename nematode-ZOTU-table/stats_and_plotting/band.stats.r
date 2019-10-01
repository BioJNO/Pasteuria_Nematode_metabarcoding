library(ggplot2)
library(reshape2)
library(plyr)
library(multcompView)

## load data ##
nem_data 						   <- read.csv("nem-zotus-merged-with-meta-minabd-10.csv", header =T, na.strings = c("", " ", "NA"), colClasses = c(rep("factor",10), rep("numeric", 680)))

colnames(nem_data)

# boxplot bands vs reads all data 
nem_data$total <- rowSums(nem_data[12:690])
nem_data$perul <- nem_data$total/nem_data$volume

# look at the average number of reads per band group
aggregate(nem_data[, 692], list(nem_data$band), mean)

# plot 
bandscore_box					<- ggplot(data = nem_data, aes(band.strength, perul)) +geom_boxplot(aes(fill = band.strength)) + xlab("PCR band score") + ylab("Assembled reads per ul") + geom_point(aes(colour = band.strength)) 
bandscore_box

# is the number of reads per microlitre statistically significant between band scores? #
ANOVA <- aov(perul ~ band.strength, data = nem_data)
summary(ANOVA)
# yes #

# is the number of reads per microlitre statistically significant between ALL band scores? #
TUKEY <- TukeyHSD(ANOVA)

generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$perul=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$perul) , ]
  return(Tukey.labels)
}

generate_label_df(TUKEY, "band.strength")
# yes #

# Define the correlation 
nem_data$band <- as.numeric(nem_data$band)
cor.test(as.numeric(nem_data$band.strength), nem_data$perul, method = "spearman", conf.level = 0.95)

# save the plot and write the above on it #
ggsave("nem_bandscore_boxplot.svg", dpi=600, width = 6, height = 4)


# Split into study groups ---------------------------------------------------
pas_control 					<- pas_data[grep("control", pas_data$group),]
pas_negatives 				<- pas_data[grep("neg", pas_data$sample),]
pas_esfn 						  <- pas_data[grep("ESFN", pas_data$group),]
pas_nsis 						  <- pas_data[grep("NSIS", pas_data$group),]

# Write out control and nsis groups -----------------------------------------
#write.csv(pas_control, "pas_control_zotu.csv")
#write.csv(pas_negatives, "pas_negative_zotus.csv")
#write.csv(pas_nsis, "pas_nsis_zotu.csv")

# Split into test conditions ------------------------------------------------
pen_det <- pas_data[grep("penetrans detection limits", pas_data$subgroup),]
hart_det <- pas_data[grep("hartismeri detection limits", pas_data$subgroup),]
hart_ratio <- pas_data[grep("hartismeri ratios", pas_data$subgroup),]
plec_ratio <- pas_data[grep("plectid ratios", pas_data$subgroup),]
pen_ratio <- pas_data[grep("penetrans ratios", pas_data$subgroup),]
tsb_ratio <- pas_data[grep("tsb ratios", pas_data$subgroup),]
control_group_negatives <- pas_data[grep("negatives", pas_data$subgroup),]


# Data filtering -------------------------------------------------------------
# Convert volume to numeric for per ul conversions
pen_det$vol <- as.numeric(pen_det$vol)
hart_det$vol <- as.numeric(hart_det$vol)

# remove string from sample ID to give numeric copy number value #
hart_det$template <- gsub("hartismeri ", "", hart_det$sample)
hart_det$template <- as.numeric(hart_det$template)
pen_det$template <- gsub("penetrans ", "", pen_det$sample)
pen_det$template <- as.numeric(pen_det$template)

#write.csv(pen_det, "pen_det_total_amps.csv")

# hartismeri
# Total reads
hart_det_plot <- ggplot(data = hart_det, aes(log(template), (total))) + geom_boxplot(aes(fill = sample)) + xlab("Log template copy number") + ylab("Total assembled reads") + geom_smooth(method="auto")+ geom_point(aes(colour = sample))

hart_det_plot

ggsave("hartismeri_detection_scatter_curve_and_box.svg", dpi=600, width = 6, height = 4)

hart_det$per.ul <- hart_det$total/hart_det$vol

# Reads per ul
hart_det_plot <- ggplot(data = hart_det, aes(log(template), per.ul)) + geom_boxplot(aes(fill = sample)) + xlab("Log Template Copy Number") + ylab("Amplicons per ul") + geom_smooth(method="auto")+ geom_point(aes(colour = sample))

hart_det_plot

ggsave("hartismeri_detection_scatter_curve_and_box_perul.svg", dpi=600, width = 6, height = 4)

# test Spearmans rank correlation #
cor.test((hart_det$total), hart_det$template, method = "spearman", conf.level = 0.95)
cor.test((hart_det$total/hart_det$vol), hart_det$template, method = "spearman", conf.level = 0.95)

# penetrans #
# total # 
pen_det_box	<- ggplot(data = pen_det, aes(log(template), total)) +geom_boxplot(aes(fill = sample)) + xlab("Log Template Copy Number") + ylab("Total Amplicons") + geom_smooth(method="auto")+ geom_point(aes(colour = sample))
pen_det_box

ggsave("penetrans_detection__scatter_curve_and_box_total.svg", dpi= 600, width = 6, height = 4)


# perul #
pen_det_box	<- ggplot(data = pen_det, aes(log(template), (total/vol))) + geom_boxplot(aes(fill = sample)) + xlab("Log Template Copy Number") + ylab("Amplicons per ul") + geom_smooth(method="auto")+ geom_point(aes(colour = sample))

pen_det_box

ggsave("penetrans_detection__scatter_curve_and_box_perul.svg", dpi= 600, width = 6, height = 4)

cor.test((pen_det$total), pen_det$template, method = "spearman", conf.level = 0.95)
cor.test((pen_det$total/pen_det$vol), pen_det$template, method = "spearman", conf.level = 0.95)

