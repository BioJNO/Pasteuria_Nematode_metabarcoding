library(ggplot2)
library(reshape2)
library(tidyr)

id.mat <- read.csv("pas-sh-distance-matrix.csv")

long.plot.data <- melt(id.mat)

print(levels(long.plot.data$variable))
long.plot.data$variable <- factor(long.plot.data$variable,
                                  levels(long.plot.data$variable)[c(1,16,17,18,
                                                                    9,10,11,2,13,12,14,
                                                                    3,5,4,8,15,
                                                                    6,7)])

long.plot.data$variable <- factor(long.plot.data$variable,
                                  rev(levels(long.plot.data$variable)))


print(levels(long.plot.data$X))
long.plot.data$X <- factor(long.plot.data$X,
                                  levels(long.plot.data$X)[c(14,11,16,17,
                                                             5,12,3,18,9,10,8,
                                                             6,15,7,13,2,4,1)])

# long.plot.data <- long.plot.data[!grepl("ramosa", long.plot.data$X),]
# long.plot.data <- long.plot.data[!grepl("ramosa", long.plot.data$variable),]
# long.plot.data <- long.plot.data[!grepl("daqus", long.plot.data$X),]
# long.plot.data <- long.plot.data[!grepl("daqus", long.plot.data$variable),]

# Generate heatmaps -----------------------------------------------------------
# Using RColorBrewer pallette. Diverging pallettes give an appreciable scale
# from a "cold" to a "hot" colour.
#
# Some diverging colour blind friendly pallettes are: BrBG, RdYlBu, RdBu, PuOr,
# PRGn, and PiYG. Gone with RdYlBu because red and blue are "heat-y" colours.
identity.heat <- ggplot(long.plot.data, aes(x=X, y=variable, fill=value)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
identity.heat

# Save it.
ggsave("Pas_ref_identity_heat.svg",
       height = 9, width =15, dpi=600)
