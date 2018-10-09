# Load required packages -----------------------------------------------------
library(vegan)
library(dummies)
library(plyr)
library(ggplot2)
library(ggrepel)
library(grid)

# Prepare data ---------------------------------------------------------------
data <- read.csv("pas-zotus-merged-with-meta.csv",
                 header =T,
                 na.strings = c("", " ", "NA"))
# Pull out NSIS sample group
nsis_data <- data[grep("NSIS", data$group),]
colnames(nsis_data)

# Drop zotus with a total abundance below 10
drop <- c(rep(TRUE, 12), colSums(nsis_data[,13:165]) > 10)
nsis_data <- nsis_data[,drop]

# Load in nsis meta data master sheet
meta <- read.csv("NSIS-meta-master-sheet.csv")
colnames(meta)
meta <- meta[,1:270]
colnames(meta)

# All samples (including blanks) merged with NSIS soil metadata.
merged_zero                 <- merge(nsis_data, meta, by = "sample")
colnames(merged_zero)
write.csv(merged_zero, "zero_nsis_master_merged.csv")

# Some manual database sorting /trimming of uniformative/incomplete columns.
cleaned_merged <- read.csv("zero_nsis_master_merged.csv")

# Set up so that if upstream filtering is changed, manual database cleanup
# does not need to be re-done.
colnames(cleaned_merged)
clean_merge_meta <- cleaned_merged[,c(1,47:139)]
colnames(merged_zero)
raw_merge_zotu <- merged_zero[, c(1,3:101)]
final_merge <- merge(raw_merge_zotu, clean_merge_meta, by = "sample")
colnames(final_merge)
write.csv(final_merge, "final_merged_NSIS_ZOTU_table.csv")

# Check the class of all merged variables is correct.
lapply(final_merge, class)
colnames(final_merge)

# Exclude samples with no Pasteuria.
final_merge$total <- rowSums(final_merge[,13:100])
final_merge$perul <- final_merge$total/final_merge$vol
non_zero <- final_merge[final_merge$perul > 10, ]

# Subset main pit samples (with the most complete metadata).
main_pit_non_zero <- non_zero[grep("V0", non_zero$OBSERVATION.NO), ]

# To look at main pit variables only substitute all instances of the
# "non-zero" data frame in the following lines with the
# "main_pit_non_zero" data subset generated above.

# Remove columns with missing values.
complete_cols_non_zero <- non_zero[,!sapply(non_zero, function(x) any(is.na(x)))]

# Convert all ZOTUs to per ul.
colnames(complete_cols_non_zero)
perul_non_zero <- complete_cols_non_zero
perul_non_zero[,11:99] <- complete_cols_non_zero[,11:99]/complete_cols_non_zero$vol

# Split zotus and metadata.
colnames(perul_non_zero)    

non_zero_zotu <- perul_non_zero[,11:99]
non_zero_zotu$total <- rowSums(non_zero_zotu)
non_zero_meta <- complete_cols_non_zero[,100:148]
non_zero_meta$total <- non_zero_zotu$total
non_zero_zotu <- non_zero_zotu[,1:89]

# Assign factors dummy numerical values.
non_zero_meta_dummy <- dummy.data.frame(non_zero_meta,
                                        names = NULL,
                                        omit.constants = TRUE,
                                        dummy.classes = getOption("dummy.classes"))

# Run NMDS and store as data frame -------------------------------------------
NMDS.log <- log1p(non_zero_zotu)
sol <- metaMDS(NMDS.log, k=3, trymax = 100)
scrs <- as.data.frame(scores(sol, display = "sites"))
# Make sure to note down 
sol
# Quck check of NMDS.
plot(sol)

# Run envfit on metadata 
non_zero_ef <- envfit(sol, non_zero_meta_dummy, perm=999)
non_zero_ef

# Convert envfit object to data frame.
non_zero_ef_df <- data.frame((non_zero_ef$vectors)$arrows*sqrt(non_zero_ef$vectors$r),
                             (non_zero_ef$vectors)$r,
                             (non_zero_ef$vectors)$pvals)

# Rename envfit data frame column headers.
non_zero_ef_df <- rename(non_zero_ef_df,
                         c("X.non_zero_ef.vectors..r" = "r",
                           "X.non_zero_ef.vectors..pvals" = "p"))
                           
non_zero_ef_df <- cbind(non_zero_ef_df, Species = rownames(non_zero_ef_df))

# Adjust p. values for number of variables
non_zero_ef_df$p.adj <- p.adjust(non_zero_ef_df$p,
                                 method = "BH",
                                 n = length(non_zero_ef_df$p))

# Write envfit dataframes to file.
write.csv(non_zero_ef_df, "non_zero_ef.csv")

# For point sizing.
scrs <- cbind(scrs, total.amplicons = non_zero_meta$total)
scrs <- cbind(scrs, volume = non_zero_meta$vol)

# Significant variables after correction for false discovery.
scrs <- cbind(scrs, record.type = non_zero_meta$RECORD.TYPE)
scrs <- cbind(scrs, horizon = non_zero_meta$Horizon)

# Fit ZOTUs to NMDS plot.
vf <- envfit(sol, NMDS.log, perm=999)
vf

# Store as data frame.
spp.scrs <- data.frame((vf$vectors)$arrows*sqrt(vf$vectors$r),
                       (vf$vectors)$r,
                       (vf$vectors)$pvals)
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs <- rename(spp.scrs,
                   c("X.vf.vectors..r" = "r",
                     "X.vf.vectors..pvals" = "p"))
spp.scrs$p.adj <- p.adjust(spp.scrs$p,
                           method = "BH",
                           n = length(spp.scrs$p))

# Write out species fits.
write.csv(spp.scrs, "non_zero_species_fit.csv")

all_fits <- rbind(spp.scrs, non_zero_ef_df)
 
# Create subgroup of pval cuttoffs.
adjsig <- subset(all_fits, p.adj<0.05)

# Plots ----------------------------------------------------------------------

# Colour by record type (organic/mineral).
p <- ggplot(scrs) + 
     geom_point(mapping = aes(x=NMDS1,
                              y=NMDS2,
                              color = record.type,
                              size = total.amplicons)) +
     coord_fixed() +
     geom_segment(data = adjsig,
                  aes(x = 0, xend=NMDS1,
                      y = 0, yend = NMDS2),
                  arrow = arrow(length = unit(0.25, "cm")),
                  colour = "grey") + 
     geom_text_repel(data = adjsig,
                     aes(x=NMDS1, 
                         y = NMDS2,
                         label = Species),
                     size = 3)
p

ggsave("all_significant_variables_NMDS-by-size.svg",
       dpi=600,
       width = 14,
       height=9)

# Colour by horizon.
print(levels(scrs$horizon))
scrs$horizon <- factor(scrs$horizon, levels(scrs$horizon)[c(3, 4, 2, 5, 1)])

p <- ggplot(scrs) + 
     geom_point(mapping = aes(x=NMDS1,
                              y=NMDS2,
                              color = horizon,
                              size = total.amplicons)) + 
     coord_fixed() + 
     geom_segment(data = adjsig,
                  aes(x=0,
                      xend=NMDS1,
                      y=0,
                      yend = NMDS2),
                      arrow = arrow(length = unit(0.25, "cm")), colour = "grey") + 
     geom_text_repel(data = adjsig,
                     aes(x=NMDS1,
                         y = NMDS2,
                         label = Species),
                     size = 3)
p

ggsave("Horizon-by-size.svg",
        dpi=600,
        width = 14,
        height=9)
