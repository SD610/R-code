# Screening and Visualization of NTA data
# By: Sheng Dong
# Edited: 3/25/2024

# Using 82FTOH aerobic robins data as an example

# Input the spreadsheet

OriTable <- read.csv("82FTOH_NonTargetedCompounds_aerobic_Robins_arranged.csv", header = TRUE)

# Set criteria to lower the total number of candidates

# Keep the rows with molecular weight less than 500 (8:2FTOH MW = 464.12)
Table_MW <- subset(OriTable, (OriTable[,4] <= 500))

# Rule out possible artifacts (find big deviation in peak area of replicates)
# Only in live group
select_label_1 <- data.frame() 

for(i in 1:nrow(Table_MW)){
  
  if(max(Table_MW[i,14:16])/min(Table_MW[i,14:16]) <= 2 & max(Table_MW[i,17:19])/min(Table_MW[i,17:19]) <= 2 & max(Table_MW[i,20:22])/min(Table_MW[i,20:22]) <= 2) {
    select_label_1 <- rbind(select_label_1, 1)
  } else {
    select_label_1 <- rbind(select_label_1, 0)  
  }
  
}
Table_rep <- cbind.data.frame(Table_MW, select_label_1)
# Keep the rows with small variation in peak areas in live treatment
Table_rep <- subset(Table_rep, (select_label_1 == 1))



########################################################################
# Visualization
# Plot peak area by time

library("ggpubr")
library("ggplot2")

for(i in 1:500){
  live_cpd <- cbind.data.frame(as.factor(c(0, 0, 0, 28, 28, 28, 90, 90, 90)), as.factor(rep("live", times=9)),data.frame(t(Table_rep[i:i, 14:22])))
  names(live_cpd) <- c("time","treatment", "area")
  
  abiotic_cpd <- cbind.data.frame(as.factor(c(0, 0, 0, 28, 28, 28, 90, 90, 90)), as.factor(rep("abiotic", times=9)),data.frame(t(Table_rep[i:i, 23:31])))
  names(abiotic_cpd) <- c("time","treatment", "area")
  
  positive_cpd <- cbind.data.frame(as.factor(c(0, 0, 0, 28, 28, 28, 90, 90, 90)), as.factor(rep("positive", times=9)),data.frame(t(Table_rep[i:i, 32:40])))
  names(positive_cpd) <- c("time","treatment", "area")
  
  
  cpd <- rbind.data.frame(live_cpd, abiotic_cpd, positive_cpd)
  
  meanplot <- ggline(cpd, x = "time", y = "area", color = "treatment", add = c("mean_se", "jitter"), merge = TRUE,  ylab = "Peak area count", xlab = "Time (day)", main = Table_rep[i:i,2])
  
  print(i) 
  print(meanplot)
  
}




