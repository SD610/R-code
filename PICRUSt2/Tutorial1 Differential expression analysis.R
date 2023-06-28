# Differential expression analysis tutorial using edgeR for PICRUSt2 results
# By: Sheng Dong
# Edited: 6/25/2023

##########  STEP 1 Prepare input file for certain set of experiments ###########

# Set working directory
setwd("D:/My Dropbox/Dropbox/Bioinformatics and visualization tutorial/PICRUSt2")

# Use the unstratified outputfile from PICRUSt2 as the input file
unstrat <- read.table(file = 'pred_metagenome_unstrat.tsv', sep = '\t', header = TRUE, row.names=1, stringsAsFactors = FALSE)

# Call R packages (Install packages)
library(edgeR)

# Take Nitrate-reducing set as an example
# Extract Nitrate-reducing experiment data from the raw data

# Control and ED treatment
ED <- unstrat[,c(15,16,11,12)]
# Control and NA treatment
NaAt <- unstrat[,c(15,16,13,14)]
# Create group variables
group_list <- factor(c(rep("control",2),rep("treatment",2)))
# Find gene expressed in at least two samples
df_ED <- ED[rowSums(cpm(ED) > 1) >= 2,]
df_NA <- NaAt[rowSums(cpm(NaAt) > 1) >= 2,]
# Normalization
df_ED <- DGEList(counts = df_ED, group = group_list)
df_ED <- calcNormFactors(df_ED)

df_NA <- DGEList(counts = df_NA, group = group_list)
df_NA <- calcNormFactors(df_NA)
# Estimate a common dispersion value across all genes and tagwise dispersion values
df_ED <- estimateCommonDisp(df_ED)
df_ED <- estimateTagwiseDisp(df_ED)

df_NA <- estimateCommonDisp(df_NA)
df_NA <- estimateTagwiseDisp(df_NA)
# Find genes with significantly different abundance
et_ED <- exactTest(df_ED)
tTag_ED <- topTags(et_ED, n=nrow(df_ED)) # may adjust the parameter to set the items to be shown and display in a certain order

et_NA <- exactTest(df_NA)
tTag_NA <- topTags(et_NA, n=nrow(df_NA))

DEG_ED=as.data.frame(tTag_ED)
DEG_NA=as.data.frame(tTag_NA)
# Set threshold
logFC_cutoff <- 2
k1_ED = (DEG_ED$FDR < 0.05)&(DEG_ED$logFC < -logFC_cutoff)
k2_ED = (DEG_ED$FDR < 0.05)&(DEG_ED$logFC > logFC_cutoff)

k1_NA = (DEG_NA$FDR < 0.05)&(DEG_NA$logFC < -logFC_cutoff)
k2_NA = (DEG_NA$FDR < 0.05)&(DEG_NA$logFC > logFC_cutoff)

DEG_ED$change = ifelse(k1_ED,"DOWN",ifelse(k2_ED,"UP","NOT"))
DEG_NA$change = ifelse(k1_NA,"DOWN",ifelse(k2_NA,"UP","NOT"))

table(DEG_ED$change)
table(DEG_NA$change)

NR_ED <- DEG_ED
NR_ED$ID <- rownames(DEG_ED)
NR_ED

NR_NA <- DEG_NA
NR_NA$ID <- rownames(DEG_NA)
NR_NA

####### STEP 2 Extract the genes increased or decreased significantly  #########

# Extract up sub dataset, merge ED and NA, and remove duplicate genes
NR_ED_up <- subset(NR_ED, change == "UP")
NR_NA_up <- subset(NR_NA, change == "UP")

NR_up <- as.data.frame(c(NR_ED_up$ID,NR_NA_up$ID))
NR_up <- NR_up[!duplicated(NR_up),]

# Add enzyme description and pathway to the UP genes
int_gene <- as.data.frame(NR_up)

info <- read.delim(file.choose()) # Input the description and pathway file (KO1-4.txt)

PWYL1 <- data.frame()
PWYL2 <- data.frame()
PWY <- data.frame()
kodescription <-data.frame()

ko <- info[,1]

for (i in 1:nrow(int_gene)){
  rownum = match(int_gene$NR_up[i], ko)
  PWYL1 <- rbind.data.frame(PWYL1, info[rownum,2])
  PWYL2 <- rbind.data.frame(PWYL2, info[rownum,3])
  PWY <- rbind.data.frame(PWY, info[rownum,4])
  kodescription <- rbind.data.frame(kodescription, info[rownum,5])
}
names(PWYL1) <- c("Pathways level 1")
names(PWYL2) <- c("Pathways level 2")
names(PWY) <- c("Pathway")
names(kodescription) <- c("kO Description")
int_gene2 <- cbind.data.frame(int_gene,PWYL1,PWYL2,PWY,kodescription)
write.csv(int_gene2,file = "gene and pathway.csv",row.names = F)

