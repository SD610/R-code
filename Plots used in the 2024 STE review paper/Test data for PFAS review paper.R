setwd("/Users/shengdong/Desktop")
df <- read.csv("Test data for PFAS review paper.csv", header = F)
df <- df[3:10,]
colnames(df) <- df[1,]
df <- df[-1,]

# Keep the columns with at least one value >= 0.01
df_sub <- df[,colSums(df >= 0.01) > 0]

library(ggplot2)
library(reshape2)
library(dplyr)

#convert data frame from a "wide" format to a "long" format
df_long = melt(df_sub, id = c("Name"))
df_long$value <- as.numeric(df_long$value)
df_long$value <- df_long$value * 100
df_long$value[df_long$value < 1] <- NA
df_long$value <- round(df_long$value,1) 

#colnames(df2) <- c("Compound", "Group","value")

# May be useful https://r-graph-gallery.com/320-the-basis-of-bubble-plot.html; https://jkzorz.github.io/2019/06/05/Bubble-plots.html

# Bubble plot test 2
bp2 = ggplot(df_long,aes(x=Name,y=variable, size=value, color=Name)) + 
  geom_point(alpha=0.5) + 
  theme(legend.key=element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color="grey"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        axis.text.x = element_text(color = "black", size = 9, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(color = "black", size = 9), 
        legend.position = "right") + 
  scale_y_discrete(limits = rev(sort(df_long$variable)) ) +
  theme(plot.margin=unit(c(.1,.1,.1,.1),"cm")) + 
  labs(title= "Relative abundance at and above 1%") +
  geom_text(data = df_long, aes(x=Name,y=variable, label = value), color = I(alpha("black", 0.85)), size = 3 ) +
  scale_size(range = c(3, 20), name="Relative abundance (%)")
bp2

# Bubble plot test 1 (doesn't work, need to assign colors)
# bp1 = ggplot(df_long, aes(x = Name, y = variable)) + 
#  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) + 
#  scale_size_continuous(limits = c(0.000001, 100), range = c(1,17), breaks = c(1,10,50,75)) + 
#  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
#  theme(legend.key=element_blank(), 
#        axis.text.x = element_text(color = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
#        axis.text.y = element_text(color = "black", face = "bold", size = 11), 
#        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
#        legend.title = element_text(size = 12, face = "bold"), 
#        panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = 1.2), 
#        legend.position = "right") +  
#  scale_fill_manual(values = colors, guide = FALSE) + 
#  scale_y_discrete(limits = rev(levels(df_long$variable))) 
#bp1

# Pie chart
# May be useful https://www.statology.org/ggplot-pie-chart/; https://r-graph-gallery.com/piechart-ggplot2.html; https://www.statmethods.net/graphs/pie.html

# Pie chart test using ECF-based AFFF data

EA <- df[2,]
rownames(EA) <- EA[,1]
EA <- EA[,-1]
EA <- as.data.frame(t(EA))

EA$`ECF-based AFFF ` <- as.numeric(EA$`ECF-based AFFF `)
EA$`ECF-based AFFF ` <-  EA$`ECF-based AFFF ` * 100
#EA$`ECF-based AFFF ` <- round(EA$`ECF-based AFFF `,1) 

cutoff <- which(EA$`ECF-based AFFF `< 1)

# create a new row with the sum of the values less than 0.01
other_row <- data.frame(value = sum(EA$`ECF-based AFFF`[cutoff]))
rownames(other_row) <- "Other"

# remove rows with values less than 0.01
EA <- EA[-cutoff, ]

# add the "Other" row to the dataframe
EA <- rbind(EA, other_row)

pc <- ggplot(data, aes(x="", y=amount, fill=category)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) 


