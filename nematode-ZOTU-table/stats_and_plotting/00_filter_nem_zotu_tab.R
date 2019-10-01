# Load ZOTU, taxonomy, and ZOTU metadata tables ------------------------------
zotutab <- read.csv("nematode_alpha1_zotutab.csv")
taxtab <- read.csv("nematode_taxonomy_combined.csv")
metatab <- read.csv("nematode_zotu_metadata.csv")

# Merge ZOTU table and all metadata ------------------------------------------
meta <- merge(taxtab, metatab, by = "Zotu_ID")
merged <- merge(meta, zotutab, by="Zotu_ID")

# Write out merged Zotu table.
write.csv(merged, "metazoan_ZOTU_tab_no_filter.csv")

# Filter ZOTUs ---------------------------------------------------------------
# Remove unwanted taxa.
merged <- merged[!grepl("Eukaryota;Alveolata", merged$taxid),]
merged <- merged[!grepl("Eukaryota;Metazoa;Annelida", merged$taxid),]
merged <- merged[!grepl("Eukaryota;Rhizaria", merged$taxid),]
merged <- merged[!grepl("Eukaryota;Metazoa;Arthropoda", merged$taxid),]
merged <- merged[!grepl("Eukaryota;Viridiplantae", merged$taxid),]
# Then grab anything at least more specific than Eukaryota.
merged <- merged[grepl("Eukaryota;", merged$taxid),]

# Filter low abundance ZOTUs
colnames(merged)
merged$abundance <- rowSums(merged[,8:910])
merged <- merged[merged$abundance >= 10, ]

# Filter sequence length.
merged <- merged[merged$seqlen >= 314, ]

# Write out filtered taxa table.
write.csv(merged, "metazoan_ZOTU_tab_post_filter.csv")

# Convert ZOTU table format. -------------------------------------------------
# Flip samples to rows and ZOTUs to columns.
metazoa_samples_as_rows <- t(merged[, c(5,8:911)])
# Set row one as column headers
colnames(metazoa_samples_as_rows) <- metazoa_samples_as_rows[1,]
metazoa_samples_as_rows <- metazoa_samples_as_rows[-1,]
# Write out flipped ZOTU table as csv.
write.csv(metazoa_samples_as_rows, "metazoa-zotu-tab-samples-as-rows-min10.csv")
