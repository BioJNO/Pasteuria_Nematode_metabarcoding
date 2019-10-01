# barcode stats nem

library(ggplot2)
library(reshape2)
library(plyr)
library(multcompView)
library(multcomp)

nem_data <- read.csv("nem-zotus-merged-with-meta-minabd-10.csv", header =T, na.strings = c("", " ", "NA"), colClasses = c(rep("factor",10), rep("numeric", 680)))

# boxplot bands vs reads all data
nem_data$total <- rowSums(nem_data[12:690])
nem_data$perul <- nem_data$total/nem_data$volume

# look at the average number of reads per band group
aggregate(nem_data[, 692], list(nem_data$band), mean)
# plot 
ftag_box <- ggplot(data = nem_data, aes(ftag, perul)) +geom_boxplot(aes(fill = ftag)) + xlab("Forward PCR barcode") + ylab("Assembled reads per ul") + geom_point(aes(colour = ftag)) + theme(axis.text.x = element_text(angle=90,hjust=1), legend.position = "none")
ftag_box

# Test for normal distribution
shapiro.test(nem_data$perul)

# is the number of reads per microlitre statistically significant between band scores?
ANOVA <- aov(perul ~ ftag, data = nem_data)
summary(ANOVA)

# is the number of reads per microlitre statistically significant between ALL band scores?
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

generate_label_df(TUKEY, "ftag")

ggsave("ftags_nem.svg", dpi= 600, width = 6, height = 4)

# rtags
# plot 
rtag_box <- ggplot(data = nem_data, aes(rtag, perul)) +geom_boxplot(aes(fill = rtag)) + xlab("Reverse PCR barcode") + ylab("Assembled reads per ul") + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") 
rtag_box

# is the number of reads per microlitre statistically significant between band scores?
ANOVA <- aov(perul ~ rtag, data = nem_data)
summary(ANOVA)

# is the number of reads per microlitre statistically significant between ALL band scores?
TUKEY <- TukeyHSD(ANOVA)

generate_label_df(TUKEY, "rtag")

ggsave("rtags_nem.svg", dpi= 600, width = 6, height = 4)

# groups 
# plot 
group_box <- ggplot(data = nem_data, aes(plate, perul)) +geom_boxplot(aes(fill = plate)) + xlab("Template group") + ylab("Assembled reads per ul") + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") 
group_box

# is the number of reads per microlitre statistically significant between band scores?
ANOVA <- aov(perul ~ group, data = nem_data)
summary(ANOVA)

fit <- aov(perul ~ ftag + group, data=nem_data)
summary(fit)

TUKEY <- TukeyHSD(fit)

plot(TUKEY)

# is the number of reads per microlitre statistically significant between ALL band scores?
TUKEY <- TukeyHSD(ANOVA)

control <- nem_data[grep("control", nem_data$group), ]
ESFN <- nem_data[grep("ESFN", nem_data$group), ]

fit <- aov(perul ~ ftag, data=ESFN)
summary(fit)
  
TUKEY <- TukeyHSD(fit)
 
generate_label_df(TUKEY, "ftag")

rtag_box <- ggplot(data = ESFN, aes(ftag, perul)) +geom_boxplot(aes(fill = ftag)) + xlab("Forward PCR barcode") + ylab("Assembled reads per ul") + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") 
rtag_box
