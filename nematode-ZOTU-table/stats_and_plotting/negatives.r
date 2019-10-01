library(ggplot2)
library(multcompView)


nemdat <- read.csv("nem-zotus-merged-with-meta-minabd-10.csv")

negatives <- nemdat[grep("neg", nemdat$sample), ]

# convert to reads per microlitre #
colnames(negatives)
negatives[,12:690] <- negatives[,12:690]/negatives$volume

# Discard low abundance ZOTUs #
negatives$total <- rowSums(negatives[,12:690])

drop <- c(rep(TRUE, 11), colSums(negatives[12:690]) >=10)
negatives <- negatives[,drop]

aggregate(negatives[, 22], list(negatives$group), mean)
