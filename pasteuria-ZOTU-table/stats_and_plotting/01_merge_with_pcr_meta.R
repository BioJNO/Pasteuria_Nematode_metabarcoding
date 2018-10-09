# Merge ZOTU tab with PCR metadata

# Load data ------------------------------------------------------------------
# Load in flipped ZOTU table and PCR metadata tables
pasdat <- read.csv("Pas-zotu-tab-samples-as-rows.csv", header = TRUE)
pasmeta <- read.csv("pasteuria-meta.csv")


# Merge data -----------------------------------------------------------------
# Rename sample ids.
colnames(pasdat)[1] <- "pcrid"
# Remove 'X' from sample ID numeric.
pasdat$pcrid <- gsub('X', '', pasdat$pcrid)
# merge ZOTU data with PCR metdata
pas_merged <- merge(pasmeta, pasdat, by = "pcrid")
# Check that worked.
colnames(pas_merged)
# Write out pre-filtering table.
write.csv(pas_merged, "pre-filter-pas-zotu.csv")

# Filter ZOTU table based on PCR metadata ------------------------------------
# Remove 'dead' barcodes.
pas_merged <- pas_merged[!grepl("GATGGT", pas_merged$ftag),]
pas_merged <- pas_merged[!grepl("GTTCAG", pas_merged$rtag),]
pas_merged <- pas_merged[!grepl("AGCAAG", pas_merged$rtag),]
pas_merged <- pas_merged[!grepl("ACACGT", pas_merged$rtag),]

# Remove plates with PCR band to read mismatch.
pas_merged <- pas_merged[!grepl("2", pas_merged$Edit.tag.plate),]
pas_merged <- pas_merged[!grepl("8", pas_merged$Edit.tag.plate),]

# Remove contaminated barcode pairs.
compromised <- read.csv("pattern-otus.csv")
compromised <- compromised[grep("3", compromised$rtag.grp), ]
comp_list <- compromised$pcrid
pas_merged <- pas_merged[!grepl(paste(comp_list, collapse='|'),
                         pas_merged$pcrid),]

# Write out PCR filtered ZOTU table.
write.csv(pas_merged, "pas-zotus-merged-with-meta.csv")
