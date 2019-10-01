library(reshape)
nemdat <- read.csv("nem-zotus-merged-with-meta-minabd-10.csv")
pasdat <- read.csv("pas-zotus-merged-with-meta.csv")

esfn_nem <- nemdat[grep("ESFN", nemdat$group), ]
# Only samples which amplified.
esfn_nem_amp <- subset(esfn_nem, band.strength >=2)
# Leavs 264 samples 
esfn_pas <- pasdat[grep("ESFN", pasdat$group), ]

# Join Pasteuria and Metazoan datasets.
esfn.data <- merge(esfn_nem_amp, esfn_pas, by = "pcrid", all.x = TRUE)
esfn.data$band.pas <- esfn.data$band
esfn.data$band.nem <- esfn.data$band.strength
esfn.data$vol.pas <- esfn.data$vol
esfn.data$vol.nem <- esfn.data$volume
colnames(esfn.data)
esfn.data <- esfn.data[,-c(2:11,691:701,826,827)]
# Set N/A (samples which amplified from nematodes but not from Pasteuria) to 0)
esfn.data[is.na(esfn.data)] <- 0

colnames(esfn.data)

write.csv(esfn.data, "all-esfn-zotus.csv")

esfn.data$pas.reads <- rowSums(esfn.data[,681:804])
esfn.data$nem.reads <- rowSums(esfn.data[,2:680])

# discard suspicious samples (which contain multiple control sequences)
discard_list <- c(522,
                  524,
                  525,
                  530,
                  531,
                  532,
                  533,
                  538,
                  539,
                  540,
                  541,
                  546,
                  547,
                  548,
                  549
)

esfn.data <- esfn.data[!esfn.data$pcrid %in% discard_list,]
# leaves 249

# convert to reads per microlitre
esfn.data[,681:804] <- esfn.data[,681:804]/esfn.data$vol.pas
esfn.data[,2:680] <- esfn.data[,2:680]/esfn.data$vol.nem
colnames(esfn.data)

colSums(esfn.data[,681:804])

# Make a drop list
drop <- c(TRUE, colSums(esfn.data[,2:804]) >= 10, rep(TRUE, 6))
esfn.data <- esfn.data[, drop]

# write filtered table
write.csv(esfn.data, "all-esfn-zotus-post-filter.csv")
