# To generate boxplot for the PFAS components in FT-based AFFF
# By: Sheng Dong
# Edited: 7/25/2024

#Import the data: use FT_AFFF.csv
FTAFFF<-read.csv(file.choose(), header=T)
#Remove blank rows
FTAFFF<-na.omit(FTAFFF)

###Start using ggplot###
library(ggplot2)
library(plyr)
library(scales)

level_order <- c("PFCA (3)","PFCA (4)","PFCA (5)","PFCA (6)","PFCA (7)","PFCA (8)","PFCA (9)",         
                 "PFCA (11)","PFSA (4)","PFSA (6)","PFSA (8)","PFSA (9)","PFSA (10)","PFPAs (8/10)","n:2 FTOH (6)",     
                 "n:2 FTOH (8)","n:2 FTOH (10)","n:2 FTCA (8)","n:2 FTS (4)","n:2 FTS (6)","n:2 FTS (8)",     
                 "FOSA-based","n:2 FtTAoS (4)","n:2 FtTAoS (6)","n:2 FtTAoS (8)","n:2 FTAB (6)","n:2 FTAB (8)","n:2 FTAB (10)",  
                 "n:2 FTAB (12)","n:2 FtSaAm (6)","n:2 FtSaAm (8)","n:2 FtTHN (6)","n:2 H-FTB (5/7/9)","n:3 FTB (5/7/9)")


#reorder the boxplot in ascendind conc., based on median concentration
ggplot(FTAFFF, aes(x= factor(PFAS, level = level_order), y=conc))+
  ggtitle("PFAS in FT-based AFFF formulations") +
  xlab("PFAS") + ylab("Concentration (mg/L)") +
  #set outlier 
  geom_boxplot(outlier.colour = "red", outlier.shape = 1)+ 
  # axis label rotate
  theme (axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  #log scale y-scale, etc.
  scale_y_continuous(trans = "log10", breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  #Adding jitter
  geom_jitter(color="blue", size=0.7, alpha=0.5)+
#Add Mean value for each group
  stat_summary(fun.y=mean, geom="point", shape=2, size=2, color="darkred", fill="red")
#Save the picture
ggsave('PFAS in FT-based AFFF formulations.tiff', units="in", width=10, height=5, dpi=300, compression = 'lzw')
