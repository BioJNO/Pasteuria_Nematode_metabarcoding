# Packages -------------------------------------------------------------------
library(vegan)
library(dummies)
library(plyr)
library(ggplot2)
library(ggrepel)
library(caret)

# Load data -------------------------------------------------------------------
esfn_nem_pas <- read.csv("ESFN_NMDS/nem-pas-esfn-zotus.csv")
esfn_oth <- read.csv("ESFN_NMDS/esfn_other_zotus.csv")
esfn_meta <- read.csv("Soil_data/esfn_metadata_summary.csv")
# NB: farm coordinates stripped from public version of above file.
pcr_meta <- read.csv("PCR_data/ESFN_sample_PCR.csv")

# Remove sequences which are not abundant at the same filter level as before in 
# just the ESFN dataset
# Make a drop list
colnames(esfn_nem_pas)
drop <- c(FALSE, TRUE, colSums(esfn_nem_pas[,3:652]) >= 10, rep(TRUE, 7))
esfn_nem_pas_min10 <- esfn_nem_pas[,drop]
colnames(esfn_nem_pas_min10)
# leaves data for 553 ZOTUs in 249 samples

semiquant <- esfn_nem_pas_min10
colnames(semiquant)
semiquant$vol.pas[semiquant$vol.pas == 0] <- 1
semiquant[,2:524] <- semiquant[,2:524]/semiquant$vol.nem
semiquant[,525:546] <- semiquant[,525:546]/semiquant$vol.pas

# next I need to subset the remaining samples with associated environmental 
# metadata. I'll do this by merging the metadata dataframe to the zotu 
# dataframe by sample field. 
merged <- merge(semiquant, pcr_meta, by="pcrid")
colnames(merged)
merged <- merge(merged, esfn_meta, by="field")
# metadata at 91 sites leaving 133 samples
colnames(merged)

write.csv(merged, "Mapping/ESFN_merged_with_metadata.csv")

# Some samples have incomplete metadata which envfit does not like.
# Restrict samples to those with complete metadata.
#
# Remove the subgroup column (all NA)
merged <- merged[,-556]

# Avoid plotting total nematode abundance.

# Lastly I want to include the non-nematode ZOTU data so I can fit it to the
# nematode community ordination.
all_data <- merge(merged, esfn_oth, by="pcrid")
colnames(all_data)

all_data[,593:840] <- all_data[,593:840]/all_data$vol.nem
write.csv(all_data, "Pairwise_correlations/all_esfn_zotus_perul.csv")

drop <- c(rep(TRUE, 2),
          colSums(all_data[,3:547]) >= 10,
          rep(TRUE, 45),
          colSums(all_data[,593:840]) >=10,
          FALSE)
all_data <- all_data[,drop]
colnames(all_data)

# get rid of that pesky subgroup column again
all_data <- all_data[,-270]

complete_rows<- na.omit(all_data)
# leaves 91 samples

colnames(complete_rows)
complete_rows_min100 <- complete_rows[rowSums(complete_rows[,3:232]) > 100, ]
# leaves 78 samples

write.csv(all_data, "ESFN_NMDS/esfn_zotus_merged_with_metadata.csv")


nem_zotu <- complete_rows_min100[,3:216]
pas_zotu <- complete_rows_min100[,217:232]
env_meta <- complete_rows_min100[,c(250,254:260)]
oth_zotu <- complete_rows_min100[,277:317]

# Ordination -----------------------------------------------------------------
# Run NMDS and store as data frame.
set.seed(17)
lognem <- log1p(nem_zotu)
zotu_nmds <- metaMDS(lognem, k=2, trymax = 1000)
scrs <- as.data.frame(scores(zotu_nmds, display = "sites"))
plot(zotu_nmds)
zotu_nmds
# stress = 0.2482174

# Fit ZOTUs and PCR metadata to ordination -----------------------------------
# Assign dummy values to environmental variable factors
meta_dummy <- dummy.data.frame(env_meta,
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

write.csv(ef_df, "ESFN_NMDS/ESFN_nem_zotu_NMDS_env_fit_perul.csv")


# Fit nematode ZOTUS back to NMDS
vf <- envfit(zotu_nmds, lognem, perm=999)
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
write.csv(spp.scrs, "ESFN_NMDS/ESFN_nem_ordination_nem_fit_perul.csv")

# Fit Pasteuria ZOTUs to NMDS
pf <- envfit(zotu_nmds, pas_zotu, perm=999)
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

write.csv(pas.scrs, "ESFN_NMDS/ESFN_nem_ordination_pas_fit_perul.csv")

# Fit non-nem ZOTUs to NMDS
oth_fit <- envfit(zotu_nmds, log1p(oth_zotu), perm=999)
oth_fit

oth_fit_df <- data.frame((oth_fit$vectors)$arrows*sqrt(oth_fit$vectors$r),
                         (oth_fit$vectors)$r,
                         (oth_fit$vectors)$pvals)
oth_fit_df <- cbind(oth_fit_df, Species = rownames(oth_fit_df))
oth_fit_df <- rename(oth_fit_df,
                     c("X.oth_fit.vectors..r" = "r",
                       "X.oth_fit.vectors..pvals" = "p"))
oth_fit_df$p.adj <- p.adjust(oth_fit_df$p,
                             method = "BH",
                             n = length(oth_fit_df$p))

write.csv(oth_fit_df, "ESFN_NMDS/ESFN_nem_ordination_non_nem_metazoan_fit_perul.csv")

# restrict to significant fits
env.sig <- subset(ef_df, p <=0.05)
nem.sig <- subset(spp.scrs, p.adj <=0.05)
nem.sig$Species
nem.sig <- nem.sig[c(1,4,19,24,28,35,46,60,71,72,77,86),]
pas.sig <- subset(pas.scrs, p.adj <=0.05)
oth.sig <- subset(oth_fit_df, p <=0.05)

scrs <- cbind(scrs, pasband = complete_rows_min100$band.pas)
scrs <- cbind(scrs, farm.type = complete_rows_min100$Farm_Type)
scrs <- cbind(scrs, sand = complete_rows_min100$Weight_Sand)
scrs <- cbind(scrs, silt = complete_rows_min100$Weight_Silt)
scrs <- cbind(scrs, nem.abd = complete_rows_min100$nem.reads)

complete_rows_min100$per.ul <- complete_rows_min100$met_reads.x/complete_rows_min100$vol.nem

scrs <- cbind(scrs, per.ul = complete_rows_min100$per.ul)

# store a colour blind friendly pallette
cb = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")

# Plots ----------------------------------------------------------------------
p <- ggplot(scrs) +
  geom_point(mapping = aes(x=NMDS1, y=NMDS2, color=(farm.type), size = per.ul)) +
  scale_color_manual(values = cb) +
  coord_fixed()

p

p2 <- p + 
  geom_segment(data = pas.sig,
               aes(x=0, xend=NMDS1, y=0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey") +
  geom_text(data=pas.sig,
            aes(x=NMDS1, y=NMDS2, label=Species,
                angle = (180/pi) * atan(NMDS2/NMDS1)),
            size=3,
            colour="black")
p2

p3 <- p2 +
  geom_segment(data = env.sig,
               aes(x=0, xend=NMDS1, y=0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey") +
  geom_text(data = env.sig,
            aes(x=NMDS1, y = NMDS2, label = Species,
                angle = (180/pi) * atan(NMDS2/NMDS1)),
            size = 3,
            color = "black")
p3

p4 <- p3 +
  geom_segment(data = oth.sig,
               aes(x=0, xend=NMDS1, y=0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey") +
  geom_text(data = oth.sig,
            aes(x=NMDS1, y = NMDS2, label = Species,
                angle = (180/pi) * atan(NMDS2/NMDS1)),
            size = 3,
            color = "black")
p4

ggsave("Figures/ESFN_nem_zotu_ordination_farm_type_pas_env_fit_semiquant.svg",
       dpi=600,
       width = 14,
       height = 9)

all_sig <- rbind(nem.sig, pas.sig, oth.sig, env.sig)


p5 <- p +
  geom_segment(data = all_sig,
               aes(x=0, xend=NMDS1, y=0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey") +
  geom_text(data = all_sig,
            aes(x=NMDS1, y = NMDS2, label = Species,
                angle = (180/pi) * atan(NMDS2/NMDS1)),
            size = 3,
            color = "black")
p5

p5 + expand_limits(x= -.75, y= -.75)

ggsave("Figures/ESFN_nem_zotu_ordination_all_fits_proportional.svg",
       dpi=600,
       width = 14,
       height = 9)
