# Heatmap and Sankey Diagram tutorial
# By: Sheng Dong
# Edited: 6/27/2023

#############################  Prepare input file  #############################

# Set working directory
setwd("D:/My Dropbox/Dropbox/Bioinformatics and visualization tutorial/PICRUSt2")

# Read in the stratified output file from PICRUSt2 (big file, failed to upload)
strat <- read.table(file = 'pred_metagenome_strat.tsv', sep = '\t', header = TRUE)
# Separate the huge table into 5 sub table and subset manually
# Tested code running automatically but failed (may think about this again later)
s1 <- as.data.frame(strat[1:500000,])
s2 <- as.data.frame(strat[500001:1000000,])
s3 <- as.data.frame(strat[1000001:1500000,])
s4 <- as.data.frame(strat[1500001:2000000,])
s5 <- as.data.frame(strat[2000001:2385023,])

write.csv(s1,file = "s1.csv",row.names = F)
write.csv(s2,file = "s2.csv",row.names = F)
write.csv(s3,file = "s3.csv",row.names = F)
write.csv(s4,file = "s4.csv",row.names = F)
write.csv(s5,file = "s5.csv",row.names = F)

# Manually subset the data (s1-s5) to get one data set with long format of KO link to ASVs
data <- read.csv('prepare data for alluvial diagram.csv')
count<- data$HLNR_ED+data$HLNR_NA
data1 <- cbind(data[,1:2], count)

# Add taxa info
taxa <- read.csv("taxa info.csv", header = TRUE)
taxa1 <- data.frame(taxa$Taxonomy)
taxonomy <- data.frame()

for (i in 1:nrow(data)) {
  
  taxonomy <- rbind(taxonomy, taxa1[as.numeric(which(taxa$X == data1[i,2])),])
}
names(taxonomy) <- c("tax")

library(tidyr)
df.taxonomy <-separate(taxonomy,tax,into = c("kingdom","phylum","class","order","family","genus"),sep = ";")
data2 <- cbind.data.frame(data1, df.taxonomy)
# remove rows with count = 0
data3 <- data2[!(data2$count == 0),]
write.csv(data3,file='prepare data for alluvial diagram_taxa.csv', row.names=FALSE)

######################### Plot Sankey/Alluvial Diagram #########################
# Call R package
library(ggalluvial)

# Prepare input file
alluv <- read.csv('prepare data for alluvial diagram_taxa.csv')
colnames(alluv)[1:2] <-  c("ko", "asv")
# Subset targeted genes
alluv1 <- alluv[(alluv$ko == "K00217" | alluv$ko == "K04098"),]

# Set output format
pdf(file="alluv.pdf")
# Set colors
library(RColorBrewer)
colourCount = length(unique(alluv1$genus))
getPalette = colorRampPalette(brewer.pal(12, "Set3"))

ggplot(data = alluv1,
       aes(axis1 = ko,   
           axis2 = genus,
           y = 0.2)) +
  geom_alluvium(aes(fill = genus), width = 1/30) +
  geom_stratum(width = 1/30, fill = "grey", colour = "white") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)), fontface = "bold", size = 2) +
  scale_x_discrete(limits = c("ko", "genus"),
                   expand = c(0.05, 0.05, 0.05, 0.2)) +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme_void() +
  theme(legend.position = "none")
dev.off()


################################### Plot Heatmap ###############################

# Read in a list of targeted genes and results with abundance
test <- read.csv("heatmap.csv", header = TRUE)

all_ko <- read.table(file = 'pred_metagenome_unstrat.tsv', sep = '\t', header = TRUE)

# Find out the abundance for those targeted genes
find_ko <- test[,1]
index <- all_ko$function.%in%find_ko
index <- data.frame(index)
test_ko <- cbind.data.frame(all_ko, index)
test_ko <- subset(test_ko, test_ko$index == "TRUE")
# Output for manual calculation or additional calculation in R
write.csv(test_ko,file='heatmap_input.csv')


# Call R library
library(pheatmap)
dataset <- read.csv("heatmap_input_ManualCal.csv",header = TRUE, row.names = 1)
# Generate heatmap figure for certain groups
aero <- dataset[,c(10:17)]

bk <- c(seq(-3,2, by = 0.01)) # adjust heatmap color

pdf(file="aero.pdf")
aero <- pheatmap(aero, 
                 cluster_rows = F,
                 cluster_cols = F,
                 show_rownames = T,
                 show_colnames = T,
                 scale = "row", 
                 color = colorRampPalette(c("navy", "white","firebrick3"))(length(bk))
)
dev.off()














