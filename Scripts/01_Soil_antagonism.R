################################################
# This script allows to reproduce the Figure 1 #

setwd(C:/[[YOUR PATH]]/Cuticular_bacteria_isolated_from_Meloidogyne_graminicola)

### Libraries
library(ggsignif) 
library(tidyverse)

Data <- read.csv("Data/Soil_antagonism.csv", sep = ",", dec = ".", header = TRUE)

Data[,"Samples"] <- gsub("^CT_F", "CT F", Data[,"Samples"])
Data[,"Samples"] <- gsub("^CA_F", "CA F", Data[,"Samples"])

### Layers for boxplot
Colors <- c("orangered", "orangered", "lightskyblue", "lightskyblue")
Fills <- c("orangered", "white", "lightskyblue", "white")

### Boxplot
p1 <- ggplot(Data)+
  geom_boxplot(aes(x = Samples, y = Percentage_of_mobile_larvae, fill = Samples), alpha = 0.65, color = "black")+
  scale_fill_manual("Samples", values = Fills, guide = "none")+
  geom_signif(aes(x = Samples, y = Percentage_of_mobile_larvae), test = wilcox.test, comparisons = list(c("CT", "CA"), c("CT Filtered", "CT"), c("CA Filtered", "CA")), map_signif_level = TRUE, y_position = c(165, 145, 145), textsize = 5)+
  geom_point(aes(x = Samples, y = Percentage_of_mobile_larvae), color = "black")+
  stat_summary(aes(x = Samples, y = Percentage_of_mobile_larvae), fun = mean, colour = "red",  geom = "point", shape = 3, size = 2)+
  ylab("Mobile J2")+
  scale_y_continuous(breaks = seq(0, 200, 50))+
  ylim(0,175)+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 12, color = Colors), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 13), legend.title = element_blank())
p1

ggsave(plot = p1, dpi = 1000, device = "jpeg", width = 7, height = 6, filename = "Results/Figure 1.jpeg")
