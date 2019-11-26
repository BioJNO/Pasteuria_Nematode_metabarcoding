# Load required packages -----------------------------------------------------
library(vegan)
library(dummies)
library(plyr)
library(ggplot2)
library(ggrepel)
library(grid)
library(caret)
library(tidyr)

# Prepare data ---------------------------------------------------------------
data <- read.csv("pas-zotus-merged-with-meta-minabd-10.csv",
                 header =T,
                 na.strings = c("", " ", "NA"))
# Pull out NSIS sample group
nsis_data <- data[grep("NSIS", data$group),]
colnames(nsis_data)

# Drop zotus with a total abundance below 10 in this sample set
drop <- c(rep(TRUE, 12), colSums(nsis_data[,13:87]) > 10, rep(TRUE, 2))
nsis_data <- nsis_data[,drop]
colnames(nsis_data)

# Load in nsis meta data master sheet
meta <- read.csv("Soil_data/NSIS-meta-master-sheet.csv")
colnames(meta)
meta <- meta[,1:270]
colnames(meta)

# All samples (including blanks) merged with NSIS soil metadata.
merged_zero <- merge(nsis_data, meta, by = "sample")
colnames(merged_zero)
write.csv(merged_zero, "Mapping/NSIS_metadata.csv")

# # Some manual database sorting /trimming of uniformative/incomplete columns.
# colnames(merged_zero)
# clean_merged <- merged_zero[,-c(112:115,118:120,122:128,136,141,144:146,148,
#                                  150,151,154:159,163:166,168:173,177,180:245,
#                                  247:256,258:263,265:269, 271:276, 320:353,355:379)]
# colnames(clean_merged)
# write.csv(clean_merged, "NSIS_NMDS/clean_merged_nsis.csv")
# clean_merged <- read.csv("NSIS_NMDS/clean_merged_nsis.csv")
# colnames(clean_merged)

clean_merged <- read.csv("NSIS_NMDS/clean_merged_nsis.csv")
colnames(clean_merged)
clean_merged_meta <- clean_merged[,c(2:11,111:188)]
colnames(merged_zero)
clean_merged_zotus <- merged_zero[,1:80]

clean_merged <- merge(clean_merged_zotus, clean_merged_meta, by = "pcrid")

# Check the class of all merged variables is correct.
lapply(clean_merged, class)
colnames(clean_merged)

# remove empty samples.
# To acccount for background noise (we don't want to plot pasteuria where it wasn't really present)
# set the perul abundance > 0
non_zero <- subset(clean_merged, perul > 10)
write.csv(non_zero, "Mapping/NSIS_clean_merged_non_zero.csv")
write.csv(non_zero, "NSIS_NMDS//NSIS_clean_merged_non_zero.csv")

# Some samples have incomplete metadata which envfit does not like.
# Restrict soil metadata to variables which are complete for all samples.
complete_cols_non_zero <- non_zero[,!sapply(non_zero, function(x) any(is.na(x)))]

# Convert all ZOTUs to per ul.
colnames(complete_cols_non_zero)
perul_non_zero <- complete_cols_non_zero
perul_non_zero[,12:77] <- complete_cols_non_zero[,12:77]/complete_cols_non_zero$vol.x

# Split zotus and metadata.
colnames(perul_non_zero)
non_zero_zotu <- perul_non_zero[,12:77]
non_zero_meta <- complete_cols_non_zero[,c(79,92:105)]
non_zero_meta$total <- non_zero$abundance
colnames(non_zero_meta)
# remove names and NGRs, mineral/organic, and total reads
non_zero_meta_reduced <- non_zero_meta[,-c(7:8,13,1,16)]
# double check class
lapply(non_zero_meta_reduced, class)

# Assign factors dummy numerical values.
non_zero_meta_dummy <- dummy.data.frame(non_zero_meta_reduced,
                                        names = NULL,
                                        omit.constants = TRUE,
                                        dummy.classes = getOption("dummy.classes"))
colnames(non_zero_meta_dummy)

# check for independence of metadata variables.
# using caret generate a correlation matrix of variables
meta_cor <- cor(non_zero_meta_dummy, method = "spearman")
highCorr <- sum(abs(meta_cor[upper.tri(meta_cor)]) > .75)
highCorr
correlations <- as.data.frame(meta_cor)
write.csv(correlations, "NSIS_NMDS/NSIS_complete_col_metadata_spearman_correlations.csv")

# Run NMDS and store as data frame -------------------------------------------
set.seed(8)
colnames(non_zero_zotu)
NMDS.log <- log1p(non_zero_zotu)
sol <- metaMDS(NMDS.log, k=2, trymax = 100)
scrs <- as.data.frame(scores(sol, display = "sites"))
# Make sure to note down 
sol
# stress = 0.245235
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
write.csv(non_zero_ef_df, "NSIS_NMDS/NMDS_soil_envfit.csv")

# For point sizing.
scrs <- cbind(scrs, total.amplicons = non_zero_meta$total)

# Significant variables after correction for false discovery.
sig_meta <- subset(non_zero_ef_df, p.adj<0.05)

scrs <- cbind(scrs, horizon = non_zero_meta$NSIS2.FIELD.HORZ)
scrs <- cbind(scrs, perul = non_zero_meta$perul)

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
write.csv(spp.scrs, "NSIS_NMDS/NSIS_species_ordination_fit.csv")
all_fits <- rbind(spp.scrs, non_zero_ef_df)
 
# Create subgroup of pval cuttoffs.
adjsig <- subset(all_fits, p.adj<0.05)
adjsig$Species
adjsig$SPecies.sensible <- c("X1_Pasteuria.penetrans_1.0",
                             "X25_Pasteuria.goettingianae_0.96",
                             "X3_Pasteuria.hartismeri_0.99",
                             "Horizon A",
                             "Dry Bulk Density",
                             "Carbon",
                             "pH",
                             "Moisture")

# store a colourblind pallete
cb = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")

# Plots ----------------------------------------------------------------------

p <- ggplot(scrs) + 
     geom_point(mapping = aes(x=NMDS1,
                              y=NMDS2,
                              color = horizon,
                              size = perul)) +
     scale_color_manual(values = cb) +
     coord_fixed() +
     geom_segment(data = adjsig,
                  aes(x = 0, xend=NMDS1,
                      y = 0, yend = NMDS2),
                  arrow = arrow(length = unit(0.25, "cm")),
                  colour = "grey") + 
     geom_text(data = adjsig,
                     aes(x=NMDS1, y = NMDS2, label = SPecies.sensible,
                         angle = (180/pi) * atan(NMDS2/NMDS1)),
                     size = 3) 
p


ggsave("Figures/Fig4_NSIS_NDMS_all_significant_variables_post_BH.svg",
       dpi=600)
