# Split nematode control groups.

nemdat <- read.csv("nem-zotus-merged-with-meta-minabd-10.csv")

detlims <- nemdat[grep("detection", nemdat$subgroup), ]
juv_ratios <- nemdat[grep("juvenile ratios", nemdat$subgroup),]
plas_ratios <- nemdat[grep("plasmid ratios", nemdat$subgroup),]
replicates <- nemdat[grep("replicates", nemdat$subgroup),]

write.csv(replicates, "nematode-replicates.csv")

colnames(detlims)
#convert to perul 
detlims[12:690] <- detlims[12:690]/detlims$volume

drop <- c(rep(TRUE, 11), colSums(detlims[12:690]) > 500)
detlims <- detlims[,drop]

colSums(detlims[12:18])
colnames(detlims)

hapla_det <- detlims[grep("hapla.detection.limits", detlims$subgroup), ]
globo_det <- detlims[grep("palida.detection.limits", detlims$subgroup), ]
chitw_det <- detlims[grep("chitwoodi.detection.limits", detlims$subgroup), ]
javan_det <- detlims[grep("javanica.detection.limits", detlims$subgroup), ]

write.csv(hapla_det, "hapla_detection_limits.csv")
write.csv(globo_det, "globo_detection_limits.csv")
write.csv(chitw_det, "chitw_detection_limits.csv")
write.csv(javan_det, "javan_detection_limits.csv")

juv_ratios$X8_X13_M.javanica_0.995 <- juv_ratios$X8_Meloidogyne.incognita..southern.root.knot.nematode._0.995 + juv_ratios$X13_Meloidogyne.incognita..southern.root.knot.nematode._1.0

plas_ratios$X8_X13_M.javanica_0.995 <- plas_ratios$X8_Meloidogyne.incognita..southern.root.knot.nematode._0.995 + plas_ratios$X13_Meloidogyne.incognita..southern.root.knot.nematode._1.0

colnames(juv_ratios)

juv_ratios <- juv_ratios[, c(1:11,13,14,20,26,255)]

plas_ratios <- plas_ratios[, c(1:11,13,14,20,26,255)]

write.csv(juv_ratios, "nematode_juvenile_ratios.csv")
write.csv(plas_ratios, "nematode_plasmid_ratios.csv")





