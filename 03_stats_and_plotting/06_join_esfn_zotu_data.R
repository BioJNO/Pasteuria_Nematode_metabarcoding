library(reshape)
metdat <- read.csv("met-zotus-merged-with-meta-minabd-10.csv")
nemdat <- read.csv("nem-zotus-merged-with-meta-minabd-10.csv")
othdat <- read.csv("oth-zotus-merged-with-meta-minabd-10.csv")
pasdat <- read.csv("pas-zotus-merged-with-meta-minabd-10.csv")

# subset ESFN samples from Metazoan datasets.
esfn_met <- metdat[grep("ESFN", metdat$group), ]
esfn_nem <- nemdat[grep("ESFN", nemdat$group), ]
esfn_oth <- othdat[grep("ESFN", othdat$group), ]

# subset only ESFN samples which amplified.
esfn_met_amp <- subset(esfn_met, band.strength >=2)
esfn_nem_amp <- subset(esfn_nem, band.strength >=2)
esfn_oth_amp <- subset(esfn_oth, band.strength >=2)
# Leaves data from 264 samples 

# Subset ESFN samples from the Pasteuria dataset.
esfn_pas <- pasdat[grep("ESFN", pasdat$group), ]

# Join Pasteuria and Metazoan ESFN datasets.
esfn.data <- merge(esfn_nem_amp, esfn_pas, by = "pcrid", all.x = TRUE)
esfn.data$band.pas <- esfn.data$band
esfn.data$band.nem <- esfn.data$band.strength
esfn.data$vol.pas <- esfn.data$vol
esfn.data$vol.nem <- esfn.data$volume
colnames(esfn.data)
esfn.data <- esfn.data[,-c(2:11,588:598,674,675)]
colnames(esfn.data)
esfn.data <- esfn.data[,c(1:576,578:656,577)]
# Set N/A (samples which amplified from nematodes but not from Pasteuria) to 0.
esfn.data[is.na(esfn.data)] <- 0

colnames(esfn.data)

esfn.data$pas.reads <- rowSums(esfn.data[,577:651])
esfn.data$nem.reads <- rowSums(esfn.data[,2:576])

# discard suspicious samples (which contain multiple control sequences included
# in the controls).
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
# leaves data from 249 samples.

write.csv(esfn.data, "ESFN_NMDS/nem-pas-esfn-zotus.csv")
write.csv(esfn.data, "ESFN_NMDS/nem-pas-esfn-zotus.csv")
write.csv(esfn_oth_amp, "ESFN_NMDS/esfn_other_zotus.csv")

