# Load reshape ---------------------------------------------------------------
library(reshape)

# Prepare data ---------------------------------------------------------------
nemdat <- read.csv("nem-zotus-merged-with-meta-minabd-10.csv")
pasdat <- read.csv("pas-zotus-merged-with-meta.csv")

# Pull out ESFN samples 
esfn_nem <- nemdat[grep("ESFN", nemdat$group), ]
esfn_pas <- pasdat[grep("ESFN", pasdat$group), ]

# Merge Pasteuria and nematode ZOTU tables.
esfn.data <- merge(esfn_pas, esfn_nem, by = "pcrid")
colnames(esfn.data)

# Give common PCR metadata variables unique names
esfn.data$pas.sample <- esfn.data$sample.x
esfn.data$nem.sample <- esfn.data$sample.y
esfn.data$band.pas <- esfn.data$band
esfn.data$band.nem <- esfn.data$band.strength
esfn.data$vol.pas <- esfn.data$vol
esfn.data$vol.nem <- esfn.data$volume
colnames(esfn.data)

esfn.data <- esfn.data[, -c(2:12, 166:175)]
colnames(esfn.data)

# Discard ESFN samples which are uninformative (no metzoan amplification).
esfn.data.subset <- subset(esfn.data, band.nem >=2)
# leaves 103 samples total

# write.csv(esfn.data.subset, "subset_esfn_zotus_2.csv")

# Discard suspicious samples (which contain multiple plasmid control sequences, 9 total)
rownames(esfn.data.subset)
esfn.data.subset <- esfn.data.subset[-c(35,95:103),]

# convert to reads per microlitre #
colnames(esfn.data.subset)
esfn.data.subset[,2:154] <- esfn.data.subset[,2:154]/esfn.data.subset$vol.pas
esfn.data.subset[,155:833] <- esfn.data.subset[,155:833]/esfn.data.subset$vol.nem

# Discard low abundance ZOTUs #
esfn.data$pas.reads <- rowSums(esfn.data[,2:154])
esfn.data$nem.reads <- rowSums(esfn.data[,155:833])

drop <- c(TRUE, colSums(esfn.data.subset[2:833]) >=10, rep(TRUE, 6))
esfn.data.subset <- esfn.data.subset[,drop]

# Leaves 200 ZOTUs and and 93 samples

colnames(esfn.data.subset)

# write filtered table #
write.csv(esfn.data.subset, "all-zotus-post-filter.csv")

# remove metadata #
esfn.data.subset <- esfn.data.subset[,-c(1,195:200)]

colnames(esfn.data.subset)

# Calculate pairwise ZOTU correlation
results<-matrix(nrow=0,ncol=6)
options(warnings=-1)

for(b in 1:(dim(esfn.data.subset)[2]-1)){
  for(c in (b+1):(dim(esfn.data.subset[2]-1))){
    species1.ab<-sum(esfn.data.subset[,b])
    species2.ab<-sum(esfn.data.subset[,c])
    if(species1.ab >1 & species2.ab >1){
      test<-cor.test(esfn.data.subset[,b],esfn.data.subset[,c],method="spearman",na.action=na.rm, exact = TRUE)
      rho<-test$estimate
      p.value<-test$p.value
    }
    if(species1.ab <=1 | species2.ab <= 1){
      rho<-0
      p.value<-1
    }
    new.row<-c(names(esfn.data.subset)[b], names(esfn.data.subset)[c],rho,p.value,species1.ab,species2.ab)
    results<-rbind(results,new.row)			
  }
}

results<-data.frame(data.matrix(results))

names(results)<-c("taxa1", "taxa2", "rho","p.value","ab1","ab2")

lapply(results, class)

head(results)
colnames(results)

# save all spearmans rank correlations
write.csv(results, "zotu-all-correlations.csv")

# store  values as numeric
results$p.value <- as.character(results$p.value)
results$p.value <- as.numeric(results$p.value)
results$p.value

sig_filter <- subset(results, p.value <= 0.05)

sig_filter$p.adj <- p.adjust(sig_filter$p.value, method = "BH", n = length(results$p.value))

sig_filter$rho <- as.character(sig_filter$rho)
sig_filter$rho <- as.numeric(sig_filter$rho)
sig_filter$rho

neg_filter <- subset(sig_filter, rho < 0)
pos_filter <- subset(sig_filter, rho > 0)

sig_filter2 <- subset(sig_filter, p.adj <=0.05)

neg_filter2 <- subset(sig_filter2, rho < 0)
pos_filter2 <- subset(sig_filter2, rho > 0)

write.csv(sig_filter, "sig_correlations_pre_correction.csv")
write.csv(neg_filter2, "zotu-sig_neg_assoc_corrected_perul.csv")
write.csv(pos_filter2, "zotu-sig_pos_assoc_corrected_perul.csv")
write.csv(sig_filter2, "zotu-sig_assoc_post-correction_perul.csv")
