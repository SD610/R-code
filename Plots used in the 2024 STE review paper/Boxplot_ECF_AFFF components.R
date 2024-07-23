#Import the data
ECFAFFF<-read.csv(file.choose(), header=T)
library(ggplot2)
library(plyr)
library(scales)

ECFAFFF <- na.omit(ECFAFFF)

level_order <- c("PFCA (3)","PFCA (4)","PFCA (5)","PFCA (6)","PFCA (7)","PFCA (8)","PFCA (9)","PFCA (10)",
                 "PFCA (11)","PFCA (12)","PFSA (2)","PFSA (3)","PFSA (4)", "PFSA (5)","PFSA (6)","PFSA (7)",
                 "PFSA (8)","PFSA (9)","PFSA (10)","PFPA (8/10)","n:2 FTS (4/6/8)","FASA (4/6/8)","Me/Et-FOSA(A)",
                 "Me/Et-FOSE", "FtTAoS (6)","AmPr-FASA (4)","AmPr-FASA (5)","AmPr-FASA (6)","AmPr-FASA (7)",
                 "AmPr-FASA (8)","AmPr-FASA-PrA (4)","AmPr-FASA-PrA (5)","AmPr-FASA-PrA (6)","AmPr-FASA-PrA (7)",
                 "AmPr-FASA-PrA (8)")


#ggplot(ECFAFFF, aes(x=reorder(PFAS,conc,FUN=median), y=conc))+ #reorder the boxplot in ascendind conc., based on median concentration
ggplot(ECFAFFF, aes(x= factor(PFAS, level = level_order), y=conc))+
    ggtitle("PFAS in ECF-based AFFF formulations") +
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
ggsave('PFAS in ECF-based AFFF formulations.tiff', units="in", width=10, height=5, dpi=300, compression = 'lzw')
