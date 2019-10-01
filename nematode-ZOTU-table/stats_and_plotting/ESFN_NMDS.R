
library(vegan)
library(dummies)
library(plyr)
library(ggplot2)
library(ggrepel)

esfn.data <- read.csv("all-zotus-post-filter.csv")

colnames(esfn.data)

all.otu <- esfn.data[,c(3:195)]
all.meta <- esfn.data[,c(1,2,196:201)]

# run NMDS and store as data frame

NMDS.log <- log1p(all.otu)
sol <- metaMDS(NMDS.log, k=3, trymax = 100)
scrs <- as.data.frame(scores(sol, display = "sites"))

plot(sol)
sol

scrs <- cbind(scrs, pasband = all.meta$band.pas)
scrs <- cbind(scrs, sample = all.meta$nem.sample)

nem_meta <- read.csv("nem-meta.csv")
all_meta <- merge(all.meta, nem_meta, by="pcrid")


all.vf <- envfit(sol, all.otu, perm=999)
all.vf

meta_dummy         <- dummy.data.frame(all_meta, names = NULL, omit.constants = TRUE, dummy.classes = getOption("dummy.classes"))

meta.vf <- envfit(sol, meta_dummy, perm=999)
meta.vf

spp.scrs <- data.frame((meta.vf$vectors)$arrows*sqrt(meta.vf$vectors$r), (meta.vf$vectors)$r, (meta.vf$vectors)$pvals)
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs <- rename(spp.scrs, c("X.meta.vf.vectors..r" = "r", "X.meta.vf.vectors..pvals" = "p"))
spp.scrs$p.adj <- p.adjust(spp.scrs$p, method = "BH", n = length(spp.scrs$p))

adjsig <- subset(spp.scrs, p.adj<=0.05)

write.csv(adjsig, "esfn_zotu_nmds_PCRmetadata.csv")

p <- ggplot(scrs) + geom_point(mapping = aes(x=NMDS1, y=NMDS2, color=pasband)) + coord_fixed() + geom_segment(data = adjsig, aes(x=0, xend=NMDS1, y=0, yend = NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey") + geom_text(data = adjsig, aes(x=NMDS1, y = NMDS2, label = Species), size = 3)
p


ggsave("ESFN_NMDS_pasteuria_bandscore.svg", dpi=600, width = 14, height = 9)


# store as data frame #
spp.scrs <- data.frame((all.vf$vectors)$arrows*sqrt(all.vf$vectors$r), (all.vf$vectors)$r, (all.vf$vectors)$pvals)
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs <- rename(spp.scrs, c("X.all.vf.vectors..r" = "r", "X.all.vf.vectors..pvals" = "p"))
spp.scrs$p.adj <- p.adjust(spp.scrs$p, method = "BH", n = length(spp.scrs$p))

adjsig <- subset(spp.scrs, p.adj<=0.05)

write.csv(adjsig, "esfn_zotu_nmds.csv")

p <- ggplot(scrs) + geom_point(mapping = aes(x=NMDS1, y=NMDS2, color=pasband)) + coord_fixed() + geom_segment(data = adjsig, aes(x=0, xend=NMDS1, y=0, yend = NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey") + geom_text(data = adjsig, aes(x=NMDS1, y = NMDS2, label = Species), size = 3)
p


ggsave("ESFN_NMDS_pasteuria_bandscore.svg", dpi=600, width = 14, height = 9)
