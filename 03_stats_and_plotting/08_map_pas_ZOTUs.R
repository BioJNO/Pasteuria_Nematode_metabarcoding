# Load required packages -----------------------------------------------------
library(ggplot2)
library(rgdal)
library(maptools)
library(ggmap)
library(dplyr)
library(reshape2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(stringr)

# Plot a base R map of Scotland ----------------------------------------------
# Read shapefile data.
scotland = readOGR(dsn="Mapping/Scotland_caspcs_2001",
                   layer="scotland_caspcs_2001")
scotland@data$id = rownames(scotland@data)
scotland.points = fortify(scotland, region="id")
scotland.df = left_join(scotland.points, scotland@data, by="id")

# Plot the base map
basemap = ggplot(scotland.df, aes(long, lat, group=group))
basemap = basemap + geom_polygon(fill="gray75") +
          theme(panel.background = element_rect(fill = "gray95")) +
          coord_map("mercator")
basemap = basemap + coord_equal()
basemap

# Load in ZOTU data ----------------------------------------------------------
samples <- read.csv("Mapping/NSIS_metadata.csv")
nsis_zotus <- read.csv("Mapping/NSIS_clean_merged_non_zero.csv")
esfn.data <- read.csv("Mapping/ESFN_merged_with_metadata.csv")
colnames(samples)
colnames(esfn.data)
colnames(nsis_zotus)

# Filter the NSIS data frame to leave only the most abundant 
# for ease of plotting.
#
# First make a list where TRUE means a column will be kept,
# and FALSE means a column will be discarded. 
#
# Set TRUE explicitly for values you want to keep (metadata),
# and set a condition for the ZOTU columns (>1000).
sums.df <- as.data.frame(colSums(nsis_zotus[,14:79]))
drop <- c(rep(TRUE, 13), colSums(nsis_zotus[,14:79]) > 2000, rep(TRUE, 89))
nsis_zotus <- nsis_zotus[,drop]
colnames(nsis_zotus)
colSums(nsis_zotus[14:19])

# Change the structure of the NSIS data frame -------------------------------------
# Select only the ZOTU columns from the data frame.
mnsis <- nsis_zotus[,c(14:19)]
# Transpose the ZOTUs.
mnsis <- t(mnsis)
# Take the long format ZOTU table and stack the columns into one.
dat.m <- melt(mnsis, id.vars = "sample")
# Give the columns in the stacked dataframe sensible names.
dat.m <- dat.m %>% rename(colid=Var2, zotu=Var1, reads=value)
# Set nsis_zotus as rownames
nsis_zotus$colid <- rownames(nsis_zotus)
# Merge stacked dataframe with sample dataframe and metadata.
merged <- merge(dat.m, nsis_zotus, by = "colid")
# Filter the merged dataframe to exclude samples with no merged read pairs.
merged <- merged[merged$reads > 0, ]
colnames(merged)
# Convert total assembled reads to amplicons per microliter.
merged$amplicons_per_ul <- merged$reads/merged$vol.x

# Add ZOTU data to base map ---------------------------------------------------
# Extract the columns required for plotting: 
# zotuid, sample, easting, northing, and amplicons per ul.
colnames(merged)
nsis2.plotpoints <- merged[,c(2,6,37,38,112)]
colnames(nsis2.plotpoints)

# Filter the ESFN data frame --------------------------------------------------
colnames(esfn.data)
colSums(esfn.data[,527:548])
drop <- c(rep(TRUE, 3),
          rep(FALSE, 523),
		  colSums(esfn.data[,527:548]) > 1000,
		  rep(TRUE, 35))
esfn.data <- esfn.data[,drop]
colnames(esfn.data)

# Change the structure of the ESFN data frame ---------------------------------
# Select only the ZOTU columns from the data frame.
mesfn <- esfn.data[,c(4:6)]
# Transpose the ZOTUs.
mesfn <- t(mesfn)
# Take the long format ZOTU table and stack the columns into one.
dat.m <- melt(mesfn, id.vars = "sample")
# Give the columns in the stacked dataframe sensible names.
dat.m <- dat.m %>% rename(colid=Var2, zotu=Var1, reads=value)
# Set nsis_zotus as rownames
esfn.data$colid <- rownames(esfn.data)
# Merge stacked dataframe with sample dataframe and metadata.
merged <- merge(dat.m, esfn.data, by = "colid")
# Filter the merged dataframe to exclude samples with no merged read pairs.
merged <- merged[merged$reads > 0, ]
colnames(merged)
# Convert total assembled reads to amplicons per microliter.
merged$amplicons_per_ul <- merged$pas.reads/merged$vol.pas

# Join NSIS and ESFN data -----------------------------------------------------
colnames(merged)
# Extract the columns required for plotting: 
# zotuid, sample, easting, northing, and amplicons per ul.
esfn.plotpoints <- merged[,c(2,21,30,31,45)]
colnames(esfn.plotpoints)

# Rename ESFN coordinate columns.
esfn.plotpoints <- rename(esfn.plotpoints,
                          EASTING = GIS_X_V2,
                          NORTHING = GIS_Y_V2)
colnames(esfn.plotpoints)
# Add a sample pool column.
esfn.plotpoints$sample_pool <- "ESFN"

# rename NSIS2 coordinate columns.
nsis2.plotpoints <- rename(nsis2.plotpoints,
                         EASTING = NSIS2.EASTING,
                         NORTHING = NSIS2.NORTHING,
                         sample = sample.x)
colnames(nsis2.plotpoints)
nsis2.plotpoints$sample_pool <- "NSIS2"
nsis2.plotpoints$sample <- as.factor(nsis2.plotpoints$sample)

# Join NSIS and ESFN data tables.
all_plot_points <- rbind(esfn.plotpoints, nsis2.plotpoints)
all_plot_points$sample_pool <- as.factor(all_plot_points$sample_pool)
print(levels(all_plot_points$sample_pool))
# Re-order the levels of the sample pool column factors.
all_plot_points$sample_pool <- factor(all_plot_points$sample_pool,
                                      levels(all_plot_points$sample_pool)[c(2, 1)])

print(levels(all_plot_points$zotu))
all_plot_points$zotu <- factor(all_plot_points$zotu,
                               levels(all_plot_points$zotu)[c(2,3,1,4,5,6,7)])


# Plot all samples included in the dataset after filtering from NSIS2 samples.
p1 = basemap + geom_point(mapping = aes(NSIS2.EASTING, NSIS2.NORTHING),
                          color = "black",
                          fill = "white",
                          size = 3,
                          shape = 21,
                          data = samples)
p1

# Plot all samples included in the dataset after filtering from ESFN samples.
# NB: the above file providing exact coordinates of farms is not provided.
p2 =  p1 + geom_point(shape = 24,
                      size = 3,
                      mapping = aes(GIS_X_V2, GIS_Y_V2), 
                      color = "black",
                      fill = "white",
                      data = esfn.data)

p2
							   
# Plot ESFN and NSIS samples simultaneously.
p4 = p2 + geom_point(mapping = aes(EASTING, NORTHING,
                                   color = zotu,
                                   fill = zotu,
                                   shape = sample_pool,
                                   group = 1,
                                   size = (100*amplicons_per_ul)),
                                   alpha = 0.6,
                                   data=all_plot_points) + 
          scale_fill_brewer(type = "div", palette = "Paired") +
          scale_color_brewer(type = "div", palette = "Paired")
p4

# Save as a scalable vector graphic.
ggsave("Figures/NSIS_and_ESFN_map.svg", dpi=600)
