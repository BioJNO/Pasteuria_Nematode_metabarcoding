# Data -----------------------------------------------------------------------
metdat <- read.csv("metazoan_samples_as_rows-min10.csv")
nemdat <- read.csv("nematodes_samples_as_rows-min10.csv")
othdat <- read.csv("non_nems_samples_as_rows-min10.csv")

PCR_meta <- read.csv("PCR_data/nem-meta.csv")

# Rename sample ids so we can merge them.
colnames(metdat)[1] <- "pcrid"
colnames(nemdat)[1] <- "pcrid"
colnames(othdat)[1] <- "pcrid"

metdat$pcrid <- gsub('X', '', metdat$pcrid)
nemdat$pcrid <- gsub('X', '', nemdat$pcrid)
othdat$pcrid <- gsub('X', '', othdat$pcrid)

met_merged <- merge(PCR_meta, metdat, by = "pcrid")
nem_merged <- merge(PCR_meta, nemdat, by = "pcrid")
oth_merged <- merge(PCR_meta, othdat, by = "pcrid")

colnames(met_merged)
met_merged$met_reads <- rowSums(met_merged[,11:833])
nem_merged$met_reads <- rowSums(met_merged[,11:833])
oth_merged$met_reads <- rowSums(met_merged[,11:833])

# Write out merged and filtered file. 
write.csv(met_merged, "met-zotus-merged-with-meta-minabd-10.csv")
write.csv(nem_merged, "nem-zotus-merged-with-meta-minabd-10.csv")
write.csv(oth_merged, "oth-zotus-merged-with-meta-minabd-10.csv")
