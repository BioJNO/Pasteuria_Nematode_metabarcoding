# Merge ZOTU tab with PCR metadata

# Load data ------------------------------------------------------------------
# Load in flipped ZOTU table and PCR metadata tables
pasdat <- read.csv("Filtering/Pas-zotu-tab-samples-as-rows.csv", header = TRUE)
pasmeta <- read.csv("PCR_data/pasteuria_PCR_metadata.csv")

# Merge data -----------------------------------------------------------------
# Rename sample ids.
colnames(pasdat)[1] <- "pcrid"
# Remove 'X' from sample ID numeric.
pasdat$pcrid <- sub('X', '', pasdat$pcrid)
# merge ZOTU data with PCR metdata
pas_merged <- merge(pasmeta, pasdat, by = "pcrid")
# Check that worked.
colnames(pas_merged)
# make a total abundance column 
pas_merged$abundance <- rowSums(pas_merged[12:86])
# make an abundance per ul column
pas_merged$perul <- pas_merged$abundance/pas_merged$vol
# Write out pre-filtering table.
write.csv(pas_merged, "Filtering/pre-filter-pas-zotu.csv")

# Filter ZOTU table based on PCR metadata ------------------------------------
# Remove barcodes which returned no product.
# This is largely uneccecary with updated pipeline as samples with no reads
# return no record
pas_merged <- pas_merged[!grepl("GATGGT", pas_merged$ftag),]
pas_merged <- pas_merged[!grepl("GTTCAG", pas_merged$rtag),]
pas_merged <- pas_merged[!grepl("AGCAAG", pas_merged$rtag),]
pas_merged <- pas_merged[!grepl("ACACGT", pas_merged$rtag),]

# Remove plates with PCR band to read mismatch.
pas_merged <- pas_merged[!grepl("2", pas_merged$Edit.tag.plate),]
pas_merged <- pas_merged[!grepl("8", pas_merged$Edit.tag.plate),]

# Write out PCR filtered ZOTU table.
write.csv(pas_merged, "pas-zotus-merged-with-meta-minabd-10.csv")
