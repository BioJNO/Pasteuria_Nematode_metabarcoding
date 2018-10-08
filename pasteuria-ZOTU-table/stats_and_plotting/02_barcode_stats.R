# Load required packages -----------------------------------------------------
library(ggplot2)
library(reshape2)
library(plyr)
library(multcompView)
library(multcomp)

# Load in Pasteuria ZOTU table merged with PCR metadata ----------------------
pas_data <- read.csv("pas-zotus-merged-with-meta.csv",
                     header =T,
                     na.strings = c("", " ", "NA"),
                     colClasses = c(rep("factor",11),
                     rep("numeric", 154)))
colnames(pas_data)

# Create total assembled read pair and per ul columns.
pas_data$total <- rowSums(pas_data[13:165])
pas_data$perul <- pas_data$total/pas_data$vol

# Statistical analyses -------------------------------------------------------

# Average number of reads per barcode seq.
aggregate(pas_data[, 167], list(pas_data$ftag), mean)
aggregate(pas_data[, 167], list(pas_data$rtag), mean)

# Forward barcodes -----------------------------------------------------------
# Is the number of reads per microlitre statistically significant
# between band scores?
fbar_anova <- aov(perul ~ ftag, data = pas_data)
summary(fbar_anova)

# Is the number of reads per microlitre statistically significant
# between ALL band scores?
fbar_tukey <- TukeyHSD(fbar_anova)


# Define a function which generates labels for TukeyHSD.
generate_label_df <- function(tukey, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- tukey[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  # Labels must be in the same order as in the boxplot:
  Tukey.labels$variable=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$variable) , ]
  return(Tukey.labels)
}

generate_label_df(fbar_tukey, "ftag")

# Reverse Barcodes -----------------------------------------------------------

rbar_anova <- aov(perul ~ rtag, data = pas_data)
summary(rbar_anova)
rbar_tukey <- TukeyHSD(ANOVA)
generate_label_df(rbar_tukey, "rtag")

# Barcodes vs read numbers in each sample set --------------------------------
# Split into sample groups (ESFN, NSIS, Control).
cont_samples <- pas_data[grep("control", pas_data$group), ]
esfn_samples <- pas_data[grep("ESFN", pas_data$group), ]
nsis_samples <- pas_data[grep("NSIS", pas_data$group), ]

# ANOVA rtag vs assembled read pairs per ul each set.
cont_fit <- aov(perul ~ rtag, data = cont_samples)
summary(cont_fit)
esfn_fit <- aov(perul ~ rtag, data = esfn_samples)
summary(esfn_fit)
nsis_fit <- aov(perul ~ rtag, data = nsis_samples)
summary(nsis_fit)

# Tukey HSD of ANOVA fit for each set.
cont_hsd <- TukeyHSD(cont_fit)
esfn_hsd <- TukeyHSD(esfn_fit)
nsis_hsd <- TukeyHSD(nsis_fit)

# Get HSD labels for each barcode in each set.
generate_label_df(cont_hsd, "rtag")
generate_label_df(esfn_hsd, "rtag")
generate_label_df(nsis_hsd, "rtag")


# Plot results ---------------------------------------------------------------

# Forward barcode boxplot.
fbar_box <- ggplot(data = pas_data,
                   aes(ftag, perul)) +
                   geom_boxplot(aes(fill = ftag)) + 
                   xlab("PCR band score") + 
                   ylab("Assembled reads per ul") + 
                   geom_point(aes(colour = ftag)) + 
                   theme(axis.text.x = element_text(angle=90,hjust=1),
                         legend.position = "none")
fbar_box

# Save as scalable vector graphic.
ggsave("ftags_pas.svg",
       dpi= 600,
       width = 6,
       height = 4)

# Reverse barcode boxplot.
rtag_box <- ggplot(data = pas_data,
                   aes(rtag, perul)) +
                   geom_boxplot(aes(fill = rtag)) + 
                   xlab("Reverse PCR barcode") + 
                   ylab("Assembled reads per ul") + 
                   theme(axis.text.x = element_text(angle = 90, hjust = 1),
                         legend.position = "none") 
rtag_box

ggsave("rtags_pas.svg",
       dpi= 600,
       width = 6,
       height = 4)

# Reverse barcodes by group.
# Controls.
cont_box <- ggplot(data = cont_samples,
                    aes(copy.number, perul, group=rtag)) +
                    geom_boxplot(aes(fill = rtag)) + 
                    xlab("Template DNA group") + 
                    ylab("Assembled reads per ul") + 
                    theme(axis.text.x = element_text(angle = 90, hjust = 1))
cont_box

ggsave("rtags_pas_control.svg",
       dpi= 600,
       width = 6,
       height = 4)

# ESFN.
esfn_box <- ggplot(data = esfn_samples,
                    aes(copy.number, perul, group=rtag)) +
                    geom_boxplot(aes(fill = rtag)) + 
                    xlab("Template DNA group") + 
                    ylab("Assembled reads per ul") + 
                    theme(axis.text.x = element_text(angle = 90, hjust = 1))
esfn_box

ggsave("rtags_pas_esfn.svg",
       dpi= 600,
       width = 6,
       height = 4)

# NSIS.
nsis_box <- ggplot(data = nsis_samples,
                    aes(copy.number, perul, group=rtag)) +
                    geom_boxplot(aes(fill = rtag)) + 
                    xlab("Template DNA group") + 
                    ylab("Assembled reads per ul") + 
                    theme(axis.text.x = element_text(angle = 90, hjust = 1))
nsis_box

ggsave("rtags_pas_nsis.svg",
       dpi= 600,
       width = 6,
       height = 4)
