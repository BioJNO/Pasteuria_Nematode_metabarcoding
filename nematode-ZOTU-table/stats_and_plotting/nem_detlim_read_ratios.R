library(ggplot2)

# M. hapla #

hapla_det <- read.csv("hapla_detection_limits.csv")

# re-order sample factors (in order of number of juveniles as opposed to alphabetical)
print(levels(hapla_det$sample))
hapla_det$sample 					<- factor(hapla_det$sample, levels(hapla_det$sample)[c(2, 3, 1, 4)])

hapla_det_box					<- ggplot(data = hapla_det, aes(sample, (X4_Meloidogyne.hapla_1.0 + X16_Meloidogyne_0.995))) + geom_boxplot(aes(fill = sample)) + xlab("Number of juveniles") + ylab("Assembled reads per ul") + geom_point(aes(colour = sample))

hapla_det_box

hapla_det$total <- hapla_det$X4_Meloidogyne.hapla_1.0 + hapla_det$X16_Meloidogyne_0.995

# remove string from sample ID to give numeric individual nem value #
hapla_det$juveniles <- gsub("hapla ", "", hapla_det$sample)
hapla_det$juveniles <- as.numeric(hapla_det$juveniles)

cor.test(hapla_det$juveniles, hapla_det$total, method = "spearman")

ggsave("hapla_detection_per_ul.svg", dpi= 600, width = 6, height = 4)

# Globodera # 

globo_det <- read.csv("globo_detection_limits.csv")

# re-order sample factors (in order of number of juveniles as opposed to alphabetical)
print(levels(globo_det$sample))
globo_det$sample 					<- factor(globo_det$sample, levels(globo_det$sample)[c(2, 3, 1, 4)])
colnames(globo_det)

globo_det_box					<- ggplot(data = globo_det, aes(sample, X1_Heterodera_0.995)) + geom_boxplot(aes(fill = sample)) + xlab("Number of juveniles") + ylab("Amplicons per ul") + geom_point(aes(colour = sample))

globo_det_box

globo_det$juveniles <- gsub("pal ", "", globo_det$sample)
globo_det$juveniles <- as.numeric(globo_det$juveniles)

cor.test(globo_det$juveniles, globo_det$X1_Heterodera_0.995, method = "spearman")

ggsave("globo_detection_per_ul.svg", dpi= 600, width = 6, height = 4)

# javanica # 

javan_det <- read.csv("javan_detection_limits.csv")

# re-order sample factors (in order of number of juveniles as opposed to alphabetical)
print(levels(javan_det$sample))
javan_det$sample 					<- factor(javan_det$sample, levels(javan_det$sample)[c(2, 3, 1, 4)])

javan_det$total <- javan_det$X8_Meloidogyne.javanica_1.0 + javan_det$X13_Meloidogyne.incognita..southern.root.knot.nematode._1.0

javan_det_box					<- ggplot(data = javan_det, aes(sample, total)) + geom_boxplot(aes(fill = sample)) + xlab("Number of juveniles") + ylab("Amplicons per ul") + geom_point(aes(colour = sample))

javan_det_box

javan_det$juveniles <- gsub("jav ", "", javan_det$sample)
javan_det$juveniles <- as.numeric(javan_det$juveniles)

cor.test(javan_det$juveniles, javan_det$total, method = "spearman")

ggsave("javan_detection_per_ul.svg", dpi= 600, width = 6, height = 4)

# chitwoodi # 

chitw_det <- read.csv("chitw_detection_limits.csv")

# re-order sample factors (in order of number of juveniles as opposed to alphabetical)
print(levels(chitw_det$sample))
chitw_det$sample 					<- factor(chitw_det$sample, levels(chitw_det$sample)[c(3, 4, 1, 5, 2)])

chitw_det <- chitw_det[-9,]

chitw_det_box					<- ggplot(data = chitw_det, aes(sample, X5_Meloidogyne.chitwoodi_1.0)) + geom_boxplot(aes(fill = sample)) + xlab("Number of juveniles") + ylab("Assembled reads per ul") + geom_point(aes(colour = sample))

chitw_det_box

chitw_det$juveniles <- gsub("chit ", "", chitw_det$sample)
chitw_det$juveniles <- as.numeric(chitw_det$juveniles)

cor.test(chitw_det$juveniles, chitw_det$X5_Meloidogyne.chitwoodi_1.0, method = "spearman")

ggsave("chitw_detection_per_ul.svg", dpi= 600, width = 6, height = 4)
