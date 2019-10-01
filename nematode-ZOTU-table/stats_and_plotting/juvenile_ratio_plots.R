library(ggplot2)
library(reshape2)
library(plyr)



juv_ratios <- read.csv("nematode_juvenile_ratios.csv", header =T, na.strings = c("", " ", "NA"), colClasses = c(rep("factor",11), rep("numeric", 6)))

juv_ratios$total <- rowSums(juv_ratios[,12:17])

juv_ratios <- juv_ratios[juv_ratios$total > 500, ]

colnames(juv_ratios)

tjuv_ratio 				<- t(juv_ratios[12:17])
tjuv_ratio 				<- scale(tjuv_ratio, center = F, scale= colSums(tjuv_ratio))
dat.m 							<- melt(tjuv_ratio, id.vars = "sample")
dat.m 							<- rename(dat.m, c("Var2"="colid", "Var1"="ZOTU", "value"="proportion_of_assembled_reads"))
juv_ratios$colid 				<- rownames(juv_ratios)
merged 							<- merge(dat.m, juv_ratios, by = "colid")
dat.m 							<- merged[,1:14]
dat.m$pcrid 					<- as.numeric(dat.m$pcrid)

dat.m$X <- as.numeric(dat.m$X)

glob <- dat.m[grep("Glob", dat.m$sample), ]
java <- dat.m[grep("Jav", dat.m$sample), ]
chit <- dat.m[grep("Chit", dat.m$sample), ]
hapl <- dat.m[grep("Hap", dat.m$sample), ]
comu <- dat.m[grep("Com", dat.m$sample), ]

write.csv(glob, "globodera_juvenile_prop.csv")

glob <- read.csv("globodera_juvenile_prop.csv")

glob_juv_ratio_plot 						<- ggplot(glob, aes(x = X, y = proportion_of_assembled_reads, fill=ZOTU, order=ZOTU)) + geom_bar(stat='identity')
glob_juv_ratio_plot

java_juv_ratio_plot 						<- ggplot(java, aes(x = X, y = proportion_of_assembled_reads, fill=ZOTU, order=ZOTU)) + geom_bar(stat='identity')
java_juv_ratio_plot

chit_juv_ratio_plot 						<- ggplot(chit, aes(x = X, y = proportion_of_assembled_reads, fill=ZOTU, order=ZOTU)) + geom_bar(stat='identity')
chit_juv_ratio_plot

hapl_juv_ratio_plot 						<- ggplot(hapl, aes(x = X, y = proportion_of_assembled_reads, fill=ZOTU, order=ZOTU)) + geom_bar(stat='identity')
hapl_juv_ratio_plot

comu_juv_ratio_plot 						<- ggplot(comu, aes(x = sample, y = proportion_of_assembled_reads, fill=ZOTU, order=ZOTU)) + geom_bar(stat='identity')
comu_juv_ratio_plot


comu_juv_ratio_plot <- ggplot(comu, aes(x = sample, y = proportion_of_assembled_reads, fill=ZOTU, order=ZOTU)) + geom_boxplot(aes(fill = ZOTU)) + xlab("Sample") + ylab("Proportion of assembled reads")
comu_juv_ratio_plot 

ggsave("community_plots.svg", dpi=600, width = 8, height = 5)

