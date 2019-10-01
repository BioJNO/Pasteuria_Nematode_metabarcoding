unfiltered <- read.csv("metazoan_ZOTU_tab_no_filter.csv", colClasses = c(rep("factor", 8), c(rep("numeric", 903))))

colnames(unfiltered)

# flip samples to rows and taxa to columns
samples_as_rows <- t(unfiltered[, c(3,9:911)])

# set row one as column headers
colnames(samples_as_rows) <- samples_as_rows[1,]
samples_as_rows <- samples_as_rows[-1,]

nemdat <- as.data.frame(samples_as_rows)

rownames(nemdat)

nemdat$pcrid <- rownames(nemdat)
nemdat$pcrid

nemdat$pcrid <- gsub('X', '', nemdat$pcrid)

nemmeta <- read.csv("nem-meta.csv")

merged <- merge(nemmeta, nemdat, by = "pcrid")

nem_merged <- merged[!grepl("CTGAAC", merged$ftag),]
nem_merged <- merged[!grepl("GATGGT", merged$ftag),]
nem_merged <- merged[!grepl("GACCTT", merged$rtag),]
nem_merged <- merged[!grepl("TGGACT", merged$rtag),]


nem_merged$pcrid
#discardlist <- ["79", "212", "213", "216", "217", "220", "221", "223", "224", "225"]
nem_merged <- nem_merged[-c(63,196:208)]

colnames(merged)

ESFN <- merged[grep("ESFN", merged$group), ]

colnames(ESFN)

write.csv(ESFN, "unfiltered_ESFN_ZOTUs.csv")
ESFN <- read.csv("unfiltered_ESFN_ZOTUs.csv")

colnames(ESFN)

drop <- c(rep(TRUE, 11), colSums(ESFN[12:851]) > 0)
ESFN <- ESFN[,drop]

colnames(ESFN)

tdat <- t(ESFN[12:819])

tdat <- as.data.frame(tdat)

lapply(tdat, class)

tdat$taxid <- rownames(tdat)
tdat$taxid
colnames(tdat)
tdat$total <- rowSums(tdat[,1:504])
# 808 ZOTUs
totalsum <-colSums(tdat[506])
totalsum
# 888907 reads

nem_count <- tdat[grep("Nematoda", tdat$taxid), ]
# 548 ZOTUs

nemsum <- colSums(nem_count[506])
nemsum
# 610009 reads

percent_nem <- (nemsum/totalsum)*100
percent_nem
# 68.6%

fungi_count <- tdat[grep("Fungi", tdat$taxid), ]
# 73 ZOTUs
fungi_sum <- colSums(fungi_count[505])
fungi_sum
# 10148 reads

percent_fungi <- (fungi_sum/totalsum)*100
percent_fungi
# 1.14%

arth_count <- tdat[grep("Arthropoda", tdat$taxid), ]
# 9 ZOTUs
arth_sum <- colSums(arth_count[505])
arth_sum
# 2050 reads

percent_arth <- (arth_sum/totalsum)*100
percent_arth
# 0.2%

alveo_count <- tdat[grep("Alveolata", tdat$taxid), ]
# 31 ZOTUs
alveo_sum <- colSums(alveo_count[505])
alveo_sum
# 8772 reads

percent_alveo <- (alveo_sum/totalsum)*100
percent_alveo
# 0.99%

not_nem <- tdat[!grepl("Nematoda", tdat$taxid), ]
not_nem$taxid
