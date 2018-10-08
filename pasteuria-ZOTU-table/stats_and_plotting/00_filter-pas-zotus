# Load ZOTU, taxonomy, and ZOTU metadata tables ------------------------------
zotutab <- read.csv("zotutab.csv")
taxtab <- read.csv("taxonomy_combined.csv")
metatab <- read.csv("zotu_add_data.csv")

# Merge ZOTU table and all metadata ------------------------------------------
meta <- merge(taxtab, metatab, by = "Zotu_ID")
merged <- merge(meta, zotutab, by="Zotu_ID")

# Write out merged Zotu table.
write.csv(merged, "pasteuria_ZOTU_tab_no_filter.csv")

# Filter ZOTUs ---------------------------------------------------------------
# Exclude low abundance ZOTUs
pasdat <- merged[merged$abundance >= 10, ]
# Write out abundance filtered table.
write.csv(pasdat, "Filtered-pas-tab.csv")
# Cut filtering metadata.
pasdat <- pasdat[,-c(1:4,6:8)]

# Convert ZOTU table format. -------------------------------------------------
# Flip samples to rows and taxa to columns so that each sample is an
# observation and each ZOTU is a variable. 
pasdat_samples_as_rows <- t(pasdat)
# Set row one as column headers
colnames(pasdat_samples_as_rows) <- pasdat_samples_as_rows[1,]
pasdat_samples_as_rows <- pasdat_samples_as_rows[-1,]
# Write flipped ZOTU table to csv.
write.csv(pasdat_samples_as_rows, "Pas-zotu-tab-samples-as-rows.csv")
