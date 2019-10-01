# Packages -------------------------------------------------------------------
library(vegan)
library(dummies)
library(plyr)
library(ggplot2)
library(ggrepel)
library(caret)

# Load data -------------------------------------------------------------------
esfn_data <- read.csv("all-zotus-post-filter.csv")
esfn_meta <- read.csv("esfn_metadata_summary.csv")
# NB: farm coordinates stripped from public version of above file.
pcr_meta <- read.csv("ESFN_sample_PCR.csv")

merged <- merge(esfn_data, pcr_meta, by="pcrid")
write.csv(merged, "esfn_zotus_merged_with_pcr_metadata.csv")
merged <- merge(merged, esfn_meta, by="field")
# metadata at 91 sites leaving 133 samples

write.csv(merged, "esfn_metadata_zotus_merged.csv")

setdiff(pcr_meta$field, merged$field)

# Split ZOTUs and PCR metadata. 
colnames(merged)

# Some samples have incomplete metadata which envfit does not like.
# Restrict samples to those with complete metadata.
merged <- merged[,-337]
complete_rows<- na.omit(merged)
# leaves 91 samples

colnames(complete_rows)

rowSums(complete_rows[,4:329])

complete_rows_min100 <- complete_rows[rowSums(complete_rows[,4:329]) > 100, ]
colnames(complete_rows_min100)

all_zotu <- complete_rows_min100[,c(4:329)]
all_meta <- complete_rows_min100[,c(330:362)]
colnames(all_zotu)

nem.otu <- all_zotu[,1:305]
pas.otu <- all_zotu[,306:326]

# Ordination -----------------------------------------------------------------
# Run NMDS and store as data frame.
set.seed(17)
NMDS_log <- log1p(nem.otu)
zotu_nmds <- metaMDS(NMDS_log, k=2, trymax = 1000)
scrs <- as.data.frame(scores(zotu_nmds, display = "sites"))
plot(zotu_nmds)
zotu_nmds
# stress=0.21647

# Fit ZOTUs and PCR metadata to ordination -----------------------------------

# Assign dummy values to environmental variable factors
colnames(all_meta)
all_meta_sub <- all_meta[,-c(3:7,10,12:16)]

meta_dummy <- dummy.data.frame(all_meta_sub,
                               names = NULL,
                               omit.constants = TRUE,
                               dummy.classes = getOption("dummy.classes"))
colnames(meta_dummy)

# Fit all variables to nematode community ordination.
soil_ef <- envfit(zotu_nmds, meta_dummy, perm=999)
soil_ef

# Store fits as data frame. 
ef_df <- data.frame((soil_ef$vectors)$arrows*sqrt(soil_ef$vectors$r),
                    (soil_ef$vectors)$r,
                    (soil_ef$vectors)$pvals)
ef_df <- cbind(ef_df,
               Species = rownames(ef_df))
ef_df <- rename(ef_df,
                c("X.soil_ef.vectors..r" = "r",
                  "X.soil_ef.vectors..pvals" = "p"))

# Adjust the p value to account for the false discovery rate using
# the Benjamini-Hochberg porcedure.                 
ef_df$p.adj <- p.adjust(ef_df$p,
                        method = "BH",
                        n = length(ef_df$p))


write.csv(ef_df, "ESFN_nem_zotu_NMDS_env_fit_79_fit.csv")


# species fit all zOTUS
# Fit ZOTUs to NMDS plot.
vf <- envfit(zotu_nmds, NMDS_log, perm=999)
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
write.csv(spp.scrs, "ESFN_nem_ordination_nem_fit.csv")

# Fit Pasteuria ZOTUs to ordination
pf <- envfit(zotu_nmds, pas.otu, perm=999)
pf

pas.scrs <- data.frame((pf$vectors)$arrows*sqrt(pf$vectors$r),
                       (pf$vectors)$r,
                       (pf$vectors)$pvals)
pas.scrs <- cbind(pas.scrs, Species = rownames(pas.scrs))
pas.scrs <- rename(pas.scrs,
                   c("X.pf.vectors..r" = "r",
                     "X.pf.vectors..pvals" = "p"))
pas.scrs$p.adj <- p.adjust(pas.scrs$p,
                           method = "BH",
                           n = length(pas.scrs$p))

# combine fits 
all_fits <- rbind(spp.scrs, ef_df)

# restrict to significant after BH correction
adjsig <- subset(all_fits, p.adj<=0.03)
adjsig$Species
adjsig <- adjsig[!grepl("Acrobeloides", adjsig$Species),]
adjsig <- adjsig[!grepl("Cephalob", adjsig$Species),]
adjsig <- adjsig[!grepl("uncultured", adjsig$Species),]
adjsig <- adjsig[-c(40:42),]
adjsig <- subset(adjsig, r > 0.14)

# BH correction cuts off all Pasteuria ZOTUs so set threshold at uncorrected
# probability
pas.sig <- subset(pas.scrs, p<0.05)

scrs <- cbind(scrs, pasband = all_meta$band.pas)
scrs <- cbind(scrs, farm.type = all_meta$Farm_Type)
scrs <- cbind(scrs, soil.type = all_meta$soil_type)
scrs <- cbind(scrs, nem.abd = complete_rows_min100$nem.reads)

complete_rows_min100$per.ul <- complete_rows_min100$nem.reads/complete_rows_min100$vol.nem

scrs <- cbind(scrs, per.ul = complete_rows_min100$per.ul)
scrs <- cbind(scrs, Mg = all_meta$Mg_mg_per_l)
scrs <- cbind(scrs, delta_PDB = all_meta$delta_PDB)

cb = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")

# Plots ----------------------------------------------------------------------
p <- ggplot(scrs) +
     geom_point(mapping = aes(x=NMDS1, y=NMDS2, color=farm.type, size = per.ul)) +
     theme_dark() +
     scale_color_manual(values = cb) +
     coord_fixed() +
     expand_limits(x=-1, y=-1)
p

p <- p + geom_segment(data = TData, aes(xend = TData[ ,i], yend=TData[ ,j]),
                      x=0, y=0, colour="black",
                      arrow=arrow(angle=25, length=unit(0.25, "cm")))

p2 <- p + 
      geom_segment(data = pas.sig,
               aes(x=0, xend=NMDS1, y=0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey") +
      geom_text(data=pas.sig,
                    aes(x=NMDS1, y=NMDS2, label=Species,
                        angle = (180/pi) * atan(NMDS2/NMDS1)),
                    size=3,
                    colour="white")
p2

p3 <- p2 +
      geom_segment(data = adjsig,
                   aes(x=0, xend=NMDS1, y=0, yend = NMDS2),
                   arrow = arrow(length = unit(0.25, "cm")),
                   colour = "grey") +
      geom_text(data = adjsig,
                aes(x=NMDS1, y = NMDS2, label = Species,
                    angle = (180/pi) * atan(NMDS2/NMDS1)),
                size = 3,
                color = "white")
p3

ggsave("ESFN_nem_zotu_ordination_farm_type_nem_env_79_fit.svg",
      dpi=600,
      width = 14,
      height = 9)
