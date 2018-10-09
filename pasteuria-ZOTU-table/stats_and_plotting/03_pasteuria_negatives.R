# Load required packages -----------------------------------------------------
library(multcompView)

# Load filtered pas-zotu-table -----------------------------------------------

pas_data <- read.csv("pas-zotus-merged-with-meta.csv",
                     header =T,
                     na.strings = c("", " ", "NA"),
                     colClasses = c(rep("factor",11),
                                    rep("numeric", 154)))
colnames(pas_data)

# Analyses -------------------------------------------------------------------

# Add total assembled read pair and perul columns.
pas_data$total <- rowSums(pas_data[13:165])
pas_data$perul <- pas_data$total/pas_data$vol

# Sub-group negatives from all samples.
negatives <- pas_data[grep("neg", pas_data$sample), ]

# Discard low abundance ZOTUs.
drop <- c(rep(TRUE,12), colSums(negatives[13:165]) >=10, rep(TRUE, 2))
negatives <- negatives[,drop]
# See which ZOTUs turn up in negative controls.
colnames(negatives)

# Get the average number of reads in negatives in each sample group.
aggregate(negatives[, 20], list(negatives$group), mean)
