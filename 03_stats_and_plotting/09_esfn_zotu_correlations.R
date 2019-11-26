# Packages -------------------------------------------------------------------
library(reshape)

# Data -----------------------------------------------------------------------
esfn_data <- read.csv("Pairwise_correlations/all_esfn_ZOTUS.csv")
colnames(esfn_data)

# Remove metadata.
zotus <- esfn_data[,-c(1:2)]
colnames(zotus)

# Pairwise ZOTU correlations -------------------------------------------------
results<-matrix(nrow=0,ncol=6)
options(warnings=-1)

# Function to determine the pairwise spearmans rank correlation
# between ZOTUs 
for(b in 1:(dim(zotus)[2]-1)){
  for(c in (b+1):(dim(zotus[2]-1))){
    species1.ab<-sum(zotus[,b])
    species2.ab<-sum(zotus[,c])
    if(species1.ab >1 & species2.ab >1){
      test<-cor.test(zotus[,b],
                     zotus[,c],
                     method="spearman",
                     na.action=na.rm,
                     exact = TRUE)
      rho<-test$estimate
      p.value<-test$p.value
    }
    if(species1.ab <=1 | species2.ab <= 1){
      rho<-0
      p.value<-1
    }
    new.row<-c(names(zotus)[b],
               names(zotus)[c],
               rho,
               p.value,
               species1.ab,
               species2.ab)
    results<-rbind(results,new.row)			
  }
}

results <- data.frame(data.matrix(results))

names(results)<-c("taxa1", "taxa2", "rho","p.value","ab1","ab2")

lapply(results, class)

head(results)
colnames(results)

# Save all spearmans rank correlations.
write.csv(results, "Pairwise_correlations/zotu-all-spearman-pairwise-correlations.csv")

# Store  values as numeric.
results$p.value <- as.character(results$p.value)
results$p.value <- as.numeric(results$p.value)
results$p.value

sig_filter <- subset(results, p.value <= 0.05)

# Correct for the false discovery rate. 
sig_filter$p.adj <- p.adjust(sig_filter$p.value,
                             method = "BH",
                             n = length(results$p.value))

sig_filter$rho <- as.character(sig_filter$rho)
sig_filter$rho <- as.numeric(sig_filter$rho)
sig_filter$rho

neg_filter <- subset(sig_filter, rho < 0)
pos_filter <- subset(sig_filter, rho > 0)

sig_filter2 <- subset(sig_filter, p.adj <=0.05)

neg_filter2 <- subset(sig_filter2, rho < 0)
pos_filter2 <- subset(sig_filter2, rho > 0)

write.csv(sig_filter, "Pairwise_correlations/zotu_sig_spearman_pre_correction.csv")
write.csv(sig_filter2, "Pairwise_correlations/zotu-sig_spearman_assoc_post-correction.csv")

pas_sig <- sig_filter2[grep("Pasteuria", sig_filter2$taxa1),]
pas_sig <- pas_sig[!grepl("Pasteuria", pas_sig$taxa2),]

write.csv(pas_sig, "Pairwise_correlations//all_significant_pas_correlations.csv")

esfn_data_X5 <- esfn_data[esfn_data$X5_Meloidogyne_1.0 > 1, ]

esfn_data_X5$X5_Meloidogyne_1.0
