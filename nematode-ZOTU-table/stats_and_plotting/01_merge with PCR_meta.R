# Merge ZOTU table with PCR metadata

# Load data ------------------------------------------------------------------
# Load in flipped ZOTU table and PCR metadata tables
nemdat <- read.csv("metazoa-zotu-tab-samples-as-rows-min10.csv")
nemmeta <- read.csv("nem-meta.csv")

# Rename sample ids so we can merge them
colnames(nemdat)[1] <- "pcrid"
# Remove 'X' from sample ID number.
nemdat$pcrid <- gsub('X', '', nemdat$pcrid)
nem_merged <- merge(nemmeta, nemdat, by = "pcrid")
colnames(nem_merged)
# Write out unfiltered ZOTU table.
write.csv(nem_merged, "nem-zotus-merged-with-meta-minabd-10.csv")

# Filter ZOTU table based on PCR metadata ------------------------------------
# Remove dead barcodes.
nem_merged <- merged[!grepl("CTGAAC", merged$ftag),]
nem_merged <- merged[!grepl("GATGGT", merged$ftag),]
nem_merged <- merged[!grepl("GACCTT", merged$rtag),]
nem_merged <- merged[!grepl("TGGACT", merged$rtag),]

# Write out PCR filtered ZOTU table.
write.csv(nem_merged, "nem-zotus-merged-with-meta-minabd-10.csv")
