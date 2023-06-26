# Network analysis tutorial #1:prepare input files for Gephi
# By: Sheng Dong
# Edited: 6/20/2023

#################### STEP 1 Prepare input file for Gephi #####################
 
# Prepare an input file by saving the absolute abundance data as a .txt file

# Set working directory
setwd("D:/My Dropbox/Dropbox/Bioinformatics and visualization tutorial/Network analysis")

# Call R packages (Install packages)
library(psych)

# Import the absolute abundance data and calculate the relative abundance
asv <- read.delim('Network analysis tutorial 1_input1.txt',sep='\t',row.names=1)
asv <- t(asv)
asv <- asv/rowSums(asv)
asv <- t(asv)

# Preliminary data screening (optional)
# Filter out ASV with very low abundance
asv <- asv[which(rowSums(asv) >= 0.01), ]
# Add other filter conditions (optional)
#asv1 <- asv
#asv1[asvs1>0] <- 1
#asv <- asv[which(rowSums(asv1) >= 12), ]
asv <- t(asv)
# Output the filtered relative abundance data and input again (optional)
write.table(asv,file='Network analysis tutorial 1_output1.txt',sep='\t')
asv <- read.delim('Network analysis tutorial 1_output1.txt', row.names = 1)

# Calculate correlation coefficient matrix
occor <-  corr.test(asv,use="pairwise",method="spearman",adjust="fdr",alpha=.05)
occor.r <- occor$r #Get R value
occor.p <- occor$p #Get p value

# Set the threshold of the existence of interactions between ASVs and make all the other data in the R matrix as zero
occor.r[occor.p>0.01|occor.r<0.8] <- 0
# Make the values on diagonal all zeros (self-related)
diag(occor.r) <- 0


# Output the matrix of network file
write.csv(occor.r,file='Network analysis tutorial 1_output2.csv')

#################### STEP 2 Modify Edge and Node file in Gephi #####################

# Use the output2 file as an input in Gephi
# File -> Input spreadsheet -> Next -> Finish
# Graph type -> undirected -> OK
# Data table -> Node/edge -> Export

# Import the Edge file into R
edge <-read.csv('Network analysis tutorial 1_input2_Edge.csv')
edge[which(edge$Weight > 0),'pn'] <- 'p'
edge[which(edge$Weight < 0),'pn'] <- 'n'
write.csv(edge,file='Network analysis tutorial 1_output3_Edge.csv')
# Delete the first column manually in Excel

# Modify Node file
# Add column of other label/genus name/etc.
taxa <- read.csv("taxa info.csv", header = TRUE)
taxa1 <- data.frame(taxa$Taxonomy)

node <-read.csv('Network analysis tutorial 1_input2_Node.csv')

taxonomy <- data.frame()

for (i in 1:nrow(node)) {
  
  taxonomy <- rbind(taxonomy, taxa1[as.numeric(which(taxa$X == node[i,1])),])
}
names(taxonomy) <- c("tax")


library(tidyr)
df.taxonomy <-separate(taxonomy,tax,into = c("kingdom","phylum","class","order","family","genus"),sep = ";")
node <- cbind.data.frame(node, df.taxonomy)
write.csv(node,file='Network analysis tutorial 1_output3_Node.csv')
# Delete the first column manually


# Go back to Gephi
# Import spreadsheet: the node file 
# Import spreadsheet: the edge file (Append to existing workspace)

# Run statistics in Gephi
# Adjust Appearance in Gephi (Layout, nodes and Edges color, size, etc.)

