# Data -----------------------------------------------------------------------
# Load ZOTU table, combined taxonomy table, and ZOTU metadata table.
zotutab <- read.csv("Linux_data/minee2-nematode_alpha1_zotutab.csv")
taxtab <- read.csv("Linux_data/minee2-nematode_taxonomy_combined.csv")
metatab <- read.csv("Linux_data/minee2-nematode_zotu_metadata.csv")

# Merge all ZOTU metadata.
meta <- merge(taxtab, metatab, by = "Zotu_ID")
setdiff(metatab$Zotu_ID, zotutab$Zotu_ID)
# missing metadat for a green algae (not sure why)
# Merge ZOTU and metadata tables.
merged <- merge(meta, zotutab, by="Zotu_ID")
# Write out merged Zotu table.
write.csv(merged, "Filtering/metazoan_ZOTU_tab_no_filter.csv")
colnames(merged)

esfn.only <- merged[,1:511]
colnames(esfn.only)
esfn.only$abundance <- rowSums(esfn.only[,8:511])
esfn.only <- esfn.only[esfn.only$abundance > 0, ]
# 834 - 803 = 31 unique to controls
write.csv(esfn.only, "Filtering/ESFN_ZOTU_tab_no_filter.csv")


# Data filtering -------------------------------------------------------------
# Filter low abundance ZOTUS.
colnames(merged)
merged$abundance <- rowSums(merged[,8:910])
merged_min10 <- merged[merged$abundance >= 10, ]
# leaves 823 ZOTUs

# Split nematode and non-nematode ZOTUs 
nem.only <- merged_min10[grepl("Nematoda", merged_min10$taxid),]
not_nem <- merged_min10[!grepl("Nematoda", merged_min10$taxid),]

# Write out filtered taxa table.
write.csv(merged_min10, "Filtering/all_metazoa_min10.csv")
write.csv(nem.only, "Filtering/nem_only_min10.csv")
write.csv(not_nem, "Filtering/non_nematode_min10.csv")

# Change data structure ------------------------------------------------------
# Flip samples to rows and ZOTUs to columns.
metazoan_samples_as_rows <- t(merged_min10[, c(5,8:910)])
nematodes_samples_as_rows <- t(nem.only[, c(5, 8:910)])
non_nems_samples_as_rows <- t(not_nem[, c(5, 8:910)])

# Set row one as column headers.
colnames(metazoan_samples_as_rows) <- metazoan_samples_as_rows[1,]
metazoan_samples_as_rows <- metazoan_samples_as_rows[-1,]

colnames(nematodes_samples_as_rows) <- nematodes_samples_as_rows[1,]
nematodes_samples_as_rows <- nematodes_samples_as_rows[-1,]

colnames(non_nems_samples_as_rows) <- non_nems_samples_as_rows[1,]
non_nems_samples_as_rows <- non_nems_samples_as_rows[-1,]

# Write out flipped otu tables as csv.
write.csv(metazoan_samples_as_rows,
          "metazoan_samples_as_rows-min10.csv")
write.csv(nematodes_samples_as_rows,
          "nematodes_samples_as_rows-min10.csv")
write.csv(non_nems_samples_as_rows,
          "non_nems_samples_as_rows-min10.csv")

