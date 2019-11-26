# Load ZOTU, taxonomy, and ZOTU metadata tables ------------------------------
zotutab <- read.csv("Linux_data/pasteuria_alpha1_zotutab.csv")
taxtab <- read.csv("Linux_data/pasteuria_taxonomy_combined.csv")
metatab <- read.csv("Linux_data/pasteuria_zotu_metadata.csv")

# Merge ZOTU table and all metadata ------------------------------------------
# meta <- merge(taxtab, metatab, by = "Zotu_ID", all.y = TRUE)
meta <- merge(taxtab, metatab, by = "Zotu_ID")
merged <- merge(meta, zotutab, by="Zotu_ID")

# generate a column giving the total abundance in all sampels
colnames(merged)
merged$abundance <- rowSums(merged[8:682])

# esfn_only <- merged[,1:291]
# colnames(esfn_only)
# esfn_only$abundance <- rowSums(esfn_only[,8:291])
# 
# nsis_only <- merged[,c(1:7,292:389)]
# colnames(nsis_only)
# nsis_only$abundance <- rowSums(nsis_only[,8:105])


# Write out merged Zotu table.
write.csv(merged, "Filtering/pasteuria_ZOTU_tab_no_filter.csv")
write.csv(esfn_only, "Filtering/pasteuria_ZOTU_tab_ESFN_only_no_filter.csv")
write.csv(nsis_only, "Filtering/pasteuria_ZOTU_tab_NSIS_only_no_filter.csv")

# Filter ZOTUs ---------------------------------------------------------------
# Exclude low abundance ZOTUs
colnames(merged)
pasdat <- merged[merged$abundance >= 10, ]
# leaves 75

# Cut filtering metadata.
colnames(pasdat)
pasdat <- pasdat[,-c(1:4,6:7,683)]
colnames(pasdat)

# Convert ZOTU table format. -------------------------------------------------
# Flip samples to rows and taxa to columns so that each sample is an
# observation and each ZOTU is a variable. 
pasdat_samples_as_rows <- t(pasdat)
# Set row one as column headers
colnames(pasdat_samples_as_rows) <- pasdat_samples_as_rows[1,]
pasdat_samples_as_rows <- pasdat_samples_as_rows[-1,]

# Write flipped ZOTU table to csv.
write.csv(pasdat_samples_as_rows, "Filtering/Pas-zotu-tab-samples-as-rows.csv")
