# To generate bubble plot and pie chart for PFAS composition
# By: Sheng Dong
# Edited: 7/25/2024

#Import data
setwd("/Users/shengdong/Desktop")
df <- read.csv("Test data for PFAS review paper.csv", header = T)
df[df == 0] <- NA
df[df == ""] <- NA

library(ggplot2)
library(reshape2)
library(dplyr)

water <- df[1:2,]
solid<- df[3:5,]

# Remove columns with all values as NA
water <- water[, colSums(is.na(water)) < nrow(water)]
solid <- solid[, colSums(is.na(solid)) < nrow(solid)]

# Convert data frame from a "wide" format to a "long" format
wl = melt(water, id = c("Name"))
sl = melt(solid, id = c("Name"))

wl$value <- as.numeric(wl$value)
sl$value <- as.numeric(sl$value)

wl$value <- round(wl$value,2) 
sl$value <- round(sl$value,2) 

colnames(wl) <- c("Matrix","Compound","Concentration")
colnames(sl) <- c("Matrix","Compound","Concentration")

# Output the long format file and add detection frequency manually (the manually modifed files are already prepared in the folder)
write.csv(wl, file = "Matrix in water.csv", row.names = FALSE)
write.csv(sl, file = "Matrix in solid.csv", row.names = FALSE)

# Read the new files in
wl = read.csv("Matrix in water-1.csv", header = T)
sl = read.csv("Matrix in solid-1.csv", header = T)

# Define the desired order of the Compound variable
cmp_wl <- c("PFPrA","PFBA","PFPeA","PFHxA","PFHpA","PFOA","PFNA","PFDA", "PFUdA","PFDoA","PFEtS","PFPrS","PFBS","PFPeS","PFHxS","PFHpS","PFOS","PFNS","PFDS","X4.2.FTS","X6.2.FTS","X8.2.FTS","FHxSA","FOSA","X6.2.FTAB","AmPr.FBSA","AmPr.FPeSA","AmPr.FHxSA","AmPr.FBSA.PrA","AmPr.FPeSA.PrA","AmPr.FHxSA.PrA")
cmp_sl <- c("PFBA","PFPeA","PFHxA","PFHpA","PFOA","PFNA","PFDA","PFUdA","PFDoA","PFTrDA","PFTeDA","PFPrS","PFBS","PFPeS","PFHxS","PFHpS","PFOS","PFNS","PFDS","X6.2.FTUA","X8.2.FTUA","X4.2.FTS","X6.2.FTS","X8.2.FTS","X10.2.FTS","X5.3.acid","X7.3.acid","FHxSA","FOSA","FOSAA","N.EtFOSA","N.MeFOSAA","X6.2.FTAB","X8.2.FTAB","X10.2.FTAB","X12.2.FTAB","X6.2.FtSaAm","X8.2.FtSaAm","AmPr.FHxSA","AmPr.FOSA","TAmPr.FOSA","TAmPr.FHxSA")

# Convert Compound to a factor with desired levels
wl$Compound <- factor(wl$Compound, levels = cmp_wl)
sl$Compound <- factor(sl$Compound, levels = cmp_sl)

# Reorder the data frame by the Occupation variable
wl <- wl[order(wl$Compound),]
sl <- sl[order(sl$Compound),]

# Bubble plot
# Show Cm on the plot
bp_water = ggplot(wl,aes(x=Matrix,y=Compound, size=Concentration, color=Matrix)) + 
  geom_point(alpha=0.5) + 
  theme(legend.key=element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color="grey"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        axis.text.x = element_text(color = "black", size = 9, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(color = "black", size = 9), 
        legend.position = "right") + 
  scale_y_discrete(limits = rev(sort(wl$Compound)) ) +
  theme(plot.margin=unit(c(.1,.1,.1,.1),"cm")) + 
  labs(title= "Cm") +
  geom_text(data = wl, aes(x=Matrix,y=Compound, label = Concentration), color = I(alpha("black", 0.85)), size = 3 ) +
  scale_size(range = c(2, 15), name="Cm (ng/L)")
bp_water

#Show DF on the plot
bp_water1 = ggplot(wl,aes(x=Matrix,y=Compound, size=Concentration, color=Matrix)) + 
  geom_point(alpha=0.5) + 
  theme(legend.key=element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color="grey"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        axis.text.x = element_text(color = "black", size = 9, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(color = "black", size = 9), 
        legend.position = "right") + 
  scale_y_discrete(limits = rev(sort(wl$Compound)) ) +
  theme(plot.margin=unit(c(.1,.1,.1,.1),"cm")) + 
  labs(title= "Cm") +
  geom_text(data = wl, aes(x=Matrix,y=Compound, label = DF), color = I(alpha("black", 0.85)), size = 3 ) +
  scale_size(range = c(2, 15), name="Cm (ng/L)")
bp_water1


# Show both Cm and DF on the plot
bp_water2 = ggplot(wl,aes(x=Matrix,y=Compound, size=Concentration, color=Matrix)) + 
  geom_point(alpha=0.5) + 
  theme(legend.key=element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color="grey"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        axis.text.x = element_text(color = "black", size = 9, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(color = "black", size = 9), 
        legend.position = "right") + 
  scale_y_discrete(limits = rev(sort(wl$Compound)) ) +
  theme(plot.margin=unit(c(.1,.1,.1,.1),"cm")) + 
  labs(title= "Cm") +
  geom_text(data = wl, aes(x=Matrix,y=Compound, label = paste(Concentration, DF, sep = "\n")), color = I(alpha("black", 0.85)), size = 3 ) +
  scale_size(range = c(2, 15), name="Cm (ng/L)")
bp_water2

# Show DF on the plot
bp_solid = ggplot(sl,aes(x=Matrix,y=Compound, size=Concentration, color=Matrix)) + 
  geom_point(alpha=0.5) + 
  theme(legend.key=element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color="grey"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        axis.text.x = element_text(color = "black", size = 9, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(color = "black", size = 9), 
        legend.position = "right") + 
  scale_y_discrete(limits = rev(sort(sl$Compound)) ) +
  theme(plot.margin=unit(c(.1,.1,.1,.1),"cm")) + 
  labs(title= "Cm") +
  geom_text(data = sl, aes(x=Matrix,y=Compound, label = DF), color = I(alpha("black", 0.85)), size = 3 ) +
  scale_size(range = c(2, 15), name="Cm (µg/kg)")
bp_solid

# Change the type from character to factor
sl$Matrix <- factor(sl$Matrix)
# Define a desired order of x variable
x_order <- c("Surface soil", "Subsurface soil", "Sediment")
# Reorder x variable
sl$Matrix <- factor(sl$Matrix, levels = x_order)
# Check the order
# levels(sl$Matrix)
# Replot
bp_solid_1 = ggplot(sl,aes(x=Matrix,y=Compound, size=Concentration, color=Matrix)) + 
  geom_point(alpha=0.5) + 
  theme(legend.key=element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color="grey"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        axis.text.x = element_text(color = "black", size = 9, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(color = "black", size = 9), 
        legend.position = "right") + 
  scale_y_discrete(limits = rev(sort(sl$Compound)) ) +
  theme(plot.margin=unit(c(.1,.1,.1,.1),"cm")) + 
  labs(title= "Cm") +
  geom_text(data = sl, aes(x=Matrix,y=Compound, label = DF), color = I(alpha("black", 0.85)), size = 3 ) +
  scale_size(range = c(2, 15), name="Cm (µg/kg)")
bp_solid_1

# Create a bubble plot with more bubbles and scales in the legend
bp_solid_2 = ggplot(sl,aes(x=Matrix,y=Compound, size=Concentration, color=Matrix)) + 
  geom_point(alpha=0.5) + 
  theme(legend.key=element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color="grey"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        axis.text.x = element_text(color = "black", size = 9, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(color = "black", size = 9), 
        legend.position = "right") + 
  scale_y_discrete(limits = rev(sort(sl$Compound)) ) +
  theme(plot.margin=unit(c(.1,.1,.1,.1),"cm")) + 
  labs(title= "Cm") +
  geom_text(data = sl, aes(x=Matrix,y=Compound, label = DF), color = I(alpha("black", 0.85)), size = 6 ) +
  scale_size_continuous(range = c(2, 15), 
                        name="Cm (µg/kg)",
                        breaks = c(10, 50, 100, 300, 600),
                        limits = c(min(sl$Concentration), max(sl$Concentration)))
bp_solid_2


bp_water_3 = ggplot(wl,aes(x=Matrix,y=Compound, size=Concentration, color=Matrix)) + 
  geom_point(alpha=0.5) + 
  theme(legend.key=element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color="grey"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        axis.text.x = element_text(color = "black", size = 9, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(color = "black", size = 9), 
        legend.position = "right") + 
  scale_y_discrete(limits = rev(sort(wl$Compound)) ) +
  theme(plot.margin=unit(c(.1,.1,.1,.1),"cm")) + 
  labs(title= "Cm") +
  geom_text(data = wl, aes(x=Matrix,y=Compound, label = DF), color = I(alpha("black", 0.85)), size = 6 ) +
  scale_size_continuous(range = c(2, 15), 
                        name="Cm (ng/L)",
                        breaks = c(10, 100, 1000, 5000, 10000),
                        limits = c(min(wl$Concentration), max(wl$Concentration)))
bp_water_3


# set resolution (dpi) and dimensions (inch) for the output file
res <- 300 
width <- 5
height <- 10 
# save the plot as a high-resolution JPG file
ggsave("Matrix in water.jpg", plot = bp_water_3, width = width, height = height, dpi = res)
ggsave("Matrix in solid.jpg", plot = bp_solid_2, width = width, height = height, dpi = res)

# Combine plots in one file
library(gridExtra)
cp <- grid.arrange(bp_water_3, bp_solid_2, nrow = 1)
ggsave("Matrix.jpg", plot = cp, width = 10, height = 10, dpi = 300)


# Pie chart
# May be helpful https://developer.aliyun.com/article/928960

df1 <- read.csv("Data for pie chart.csv", header = T)

# Separate each matrix
gw <- df1[df1$Matrix == "Groundwater",]
sw <- df1[df1$Matrix == "Surface water",]
ss <- df1[df1$Matrix == "Surface soil",]
subss <- df1[df1$Matrix == "Subsurface soil",]
sd <- df1[df1$Matrix == "Sediment",]
ecf <- df1[df1$Matrix == "ECF-based AFFF",]
ft <- df1[df1$Matrix == "FT-based AFFF",]

# Calculate sum for each class
gwd <- aggregate(Median.Conc ~ Class, gw, sum)
swd <- aggregate(Median.Conc ~ Class, sw, sum) 
ssd <- aggregate(Median.Conc ~ Class, ss, sum)
subssd <- aggregate(Median.Conc ~ Class, subss, sum)
sdd <- aggregate(Median.Conc ~ Class, sd, sum)
ecfd  <- aggregate(Median.Conc ~ Class, ecf, sum)
ftd <- aggregate(Median.Conc ~ Class, ft, sum)

# Pie chart
library(RColorBrewer)
library(dplyr)
library(graphics)
library(ggplot2)

labs_gw <- paste0(round(gwd$Median.Conc/sum(gwd$Median.Conc)*100,1), "%", sep="") 
labs_sw <- paste0(round(swd$Median.Conc/sum(swd$Median.Conc)*100,1), "%", sep="") 
labs_ss <- paste0(round(ssd$Median.Conc/sum(ssd$Median.Conc)*100,1), "%", sep="") 
labs_subss <- paste0(round(subssd$Median.Conc/sum(subssd$Median.Conc)*100,1), "%", sep="") 
labs_sd <- paste0(round(sdd$Median.Conc/sum(sdd$Median.Conc)*100,1), "%", sep="")
labs_ecf <- paste0(round(ecfd$Median.Conc/sum(ecfd$Median.Conc)*100,1), "%", sep="")
labs_ft <- paste0(round(ftd$Median.Conc/sum(ftd$Median.Conc)*100,1), "%", sep="")

# Set color for different fraction
#FFB6C1: PFCA
#FFDAB9: PFSA
#90EE90: PFPA
#FFFFE0: primary ECF
#ADD8E6: primary FT
#D3D3D3: secondary ECF
#E6E6FA: secondary FT
    
# Generate plot

my_colors <- c("#FFB6C1", "#FFDAB9", "#FFFFE0", "#ADD8E6","#D3D3D3","#E6E6FA")
names(my_colors) <- labs_gw
jpeg("Groundwater.jpeg", res = 300, width = 1000, height = 600)
par(mar = c(1, 1, 1, 3))
pie(gwd$Median.Conc, labels = labs_gw, init.angle = 90, col = my_colors[match(labs_gw, names(my_colors))], border = NA, main = "Groundwater",cex = 0.4, cex.main = 0.5)
dev.off()

my_colors <- c("#FFB6C1", "#FFDAB9", "#ADD8E6","#D3D3D3","#E6E6FA")
names(my_colors) <- labs_sw
jpeg("Surface water.jpeg", res = 300, width = 1000, height = 600)
par(mar = c(1, 1, 1, 3))
pie(swd$Median.Conc,labels=labs_sw, init.angle=90,col = my_colors[match(labs_sw, names(my_colors))], border = NA, main = "Surface water",cex = 0.4, cex.main = 0.5)
dev.off()

my_colors <- c("#FFB6C1", "#FFDAB9", "#FFFFE0", "#ADD8E6","#D3D3D3","#E6E6FA")
names(my_colors) <- labs_ss
jpeg("Surface soil.jpeg", res = 300, width = 1000, height = 600)
par(mar = c(1, 1, 1, 3))
pie(ssd$Median.Conc,labels=labs_ss, init.angle=90, col =  my_colors[match(labs_ss, names(my_colors))], border = NA, main = "Surface soil", cex = 0.4, cex.main = 0.5)
dev.off()

my_colors <- c("#FFB6C1", "#FFDAB9", "#FFFFE0", "#ADD8E6","#D3D3D3","#E6E6FA")
names(my_colors) <- labs_subss
jpeg("Subsurface soil.jpeg", res = 300, width = 1000, height = 600)
par(mar = c(1, 1, 1, 3))
pie(subssd$Median.Conc,labels=labs_subss, init.angle=90, col =  my_colors[match(labs_subss, names(my_colors))], border = NA, main = "Subsurface soil", cex = 0.4, cex.main = 0.5)
dev.off()

my_colors <- c("#FFB6C1", "#FFDAB9", "#ADD8E6","#D3D3D3","#E6E6FA")
names(my_colors) <- labs_sd
jpeg("Sediment.jpeg", res = 300, width = 1000, height = 600)
par(mar = c(1, 1, 1, 3))
pie(sdd$Median.Conc,labels=labs_sd, init.angle=90, col =  my_colors[match(labs_sd, names(my_colors))], border = NA, main = "Sediment", cex = 0.4, cex.main = 0.5)
dev.off()

my_colors <- c("#FFB6C1", "#90EE90", "#FFDAB9", "#FFFFE0", "#ADD8E6","#D3D3D3","#E6E6FA")
names(my_colors) <- labs_ecf
jpeg("ECF-based AFFF.jpeg", res = 300, width = 1000, height = 600)
par(mar = c(1, 1, 1, 3))
pie(ecfd$Median.Conc,labels=labs_ecf, init.angle=90, col =  my_colors[match(labs_ecf, names(my_colors))], border = NA, main = "ECF-based AFFF", cex = 0.4, cex.main = 0.5)
dev.off()

my_colors <- c("#FFB6C1", "#90EE90", "#FFDAB9", "#ADD8E6","#D3D3D3","#E6E6FA")
names(my_colors) <- labs_ft
jpeg("FT-based AFFF.jpeg", res = 300, width = 1000, height = 600)
par(mar = c(1, 1, 1, 3))
pie(ftd$Median.Conc,labels=labs_ft, init.angle=90, col =  my_colors[match(labs_ft, names(my_colors))], border = NA, main = "FT-based AFFF", cex = 0.4, cex.main = 0.5)
dev.off()

my_colors <- c("#FFB6C1", "#90EE90", "#FFDAB9", "#FFFFE0", "#ADD8E6","#D3D3D3","#E6E6FA")
names(my_colors) <- labs_ecf
jpeg("ECF-based AFFF with legend.jpeg", res = 300, width = 1000, height = 600)
par(mar = c(1, 1, 1, 3))
pie(ecfd$Median.Conc,labels=labs_ecf, init.angle=90, col =  my_colors[match(labs_ecf, names(my_colors))], border = NA, main = "ECF-based AFFF", cex = 0.4, cex.main = 0.5)
# Add a legend to the chart
par(lwd = 0.5) # set line width
legend("right", legend = ecfd$Class, fill = my_colors, cex = 0.6, lty = "blank", border = "black", x.intersp = 0)
dev.off()









