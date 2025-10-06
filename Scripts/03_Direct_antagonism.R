########################################################
# This script allows to reproduce the Figure 3A and 3B #


setwd(C:/[[YOUR PATH]]/Cuticular_bacteria_isolated_from_Meloidogyne_graminicola)


######## Libraries
library(ggtext)
library(multcompView)
library(patchwork)
library(tidyverse)

#### Create a function to determine the statistical groups based on the p-values
tri.to.squ <- function(x)
{
  rn <- row.names(x)
  cn <- colnames(x)
  an <- unique(c(cn,rn))
  myval <-  x[!is.na(x)]
  mymat <-  matrix(1, nrow = length(an), ncol = length(an), dimnames = list(an,an))
  for(ext in 1:length(cn))
  {
    for(int in 1:length(rn))
    {
      if(is.na(x[row.names(x) == rn[int], colnames(x) == cn[ext]])) next
      mymat[row.names(mymat) == rn[int], colnames(mymat) == cn[ext]] <- x[row.names(x) == rn[int], colnames(x) == cn[ext]]
      mymat[row.names(mymat) == cn[ext], colnames(mymat) == rn[int]] <- x[row.names(x) == rn[int], colnames(x) == cn[ext]]
    }
  }
  return(mymat)
}

### All genera
Data <- read.csv("Data/Nematicidal effect - All genus.csv", dec = ".", sep = ";", header = TRUE)

Data[grep("Paenibacillus", Data[,"Bacteria"]), "Strains"] <- "<i>Paenibacillus</i> CT02"
Data[grep("Methylobact", Data[,"Bacteria"]), "Strains"] <- "<i>Methylobacterium</i> CA03"
Data[grep("Bosea", Data[,"Bacteria"]), "Strains"] <- "<i>Bosea</i> CT05"
Data[grep("Rothia", Data[,"Bacteria"]), "Strains"] <- "<i>Rothia</i> CTRL01"
Data[grep("Janibact", Data[,"Bacteria"]), "Strains"] <- "<i>Janibacter</i> CA04"
Data[grep("Microbacterium", Data[,"Bacteria"]), "Strains"] <- "<i>Microbacterium</i> CTRL04"
Data[grep("Priestia", Data[,"Bacteria"]), "Strains"] <- "<i>Priestia</i> CTRL03"
Data[grep("Microco", Data[,"Bacteria"]), "Strains"] <- "<i>Micrococcus</i> CA05"
Data[grep("Sphingom", Data[,"Bacteria"]), "Strains"] <- "<i>Sphingomonas</i> CTRL06"
Data[grep("Panto", Data[,"Bacteria"]), "Strains"] <- "<i>Pantoea</i> CA02"
Data[grep("Burkho", Data[,"Bacteria"]), "Strains"] <- "<i>Burkholderia</i> CT03"
Data[grep("Acineto", Data[,"Bacteria"]), "Strains"] <- "<i>Acinetobacter</i> CA01"
Data[grep("Agrobact", Data[,"Bacteria"]), "Strains"] <- "<i>Agrobacterium</i> CT06"

#### Crossed larvae
# Compares the percentage of mobility between samples with a Wilcoxon-Mann-Whitney test with an adjustment of p-value using a FDR method
Stat <- pairwise.wilcox.test(Data$Percent_of_crossed_J2s, Data$Strains, p.adjust.method = "fdr")

# Determines the statistical groups
Mymat <- tri.to.squ(Stat$p.value)

Letters <- multcompLetters(Mymat, compare = "<=", threshold = 0.05, Letters = letters)
Letters <- data.frame(group = names(Letters$Letters), letter = Letters$Letters)

## Colors
Order <- c("Saline solution", "<i>Paenibacillus</i> CT02", "<i>Methylobacterium</i> CA03", "<i>Bosea</i> CT05", "<i>Rothia</i> CTRL01", "<i>Janibacter</i> CA04", "<i>Microbacterium</i> CTRL04", "<i>Priestia</i> CTRL03", "<i>Micrococcus</i> CA05", "<i>Sphingomonas</i> CTRL06", "<i>Pantoea</i> CA02", "<i>Burkholderia</i> CT03", "<i>Acinetobacter</i> CA01", "<i>Agrobacterium</i> CT06")

Fills <- c("white", "grey70", "grey70", "grey70", "grey70", "grey70", "grey70", "grey70", "grey70", "grey70", "grey70", "grey70", "grey70", "grey70")

## Plots the boxplots
p1 <- ggplot(Data)+
  geom_boxplot(aes(x = fct_relevel(Strains, Order), y = Percent_of_crossed_J2s, fill = fct_relevel(Strains, Order)), alpha = 0.65, color = "black")+
  scale_fill_manual("Bacteria", values = Fills, guide = "none")+
  geom_point(aes(x = Strains, y = Percent_of_crossed_J2s, color = reorder(Strains, - Percent_of_crossed_J2s, mean)), color = "black")+
  stat_summary(aes(x = Strains, y = Percent_of_crossed_J2s), fun = mean, colour = "red", geom = "point", shape = 3, size = 2)+
  geom_text(data = Letters, aes(label = letter, y = 125, x = group), colour = "black", size = 4)+
  ylab("J2 that passed through the sieve")+
  scale_y_continuous(breaks = seq(0, 150, 25))+
  ggtitle(" A")+
  xlab("")+
  guides(alpha = FALSE)+
  theme_bw()+
  theme(axis.text.x = element_markdown(size = 10, angle = 30, hjust = 1), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12), plot.title = element_text(face = "bold", hjust = -0.05))
p1

### Agrobacterium
Data <- read.csv("Data/Nematicidal effect - Agrobacterium.csv", dec = ".", sep = ";", header = TRUE)

Data[grep("Agrobacterium CT06", Data[,"Bacteria"]), "Strains"] <- "<i>Agrobacterium</i> CT06"
Data[grep("Agrobacterium CT07", Data[,"Bacteria"]), "Strains"] <- "<i>Agrobacterium</i> CT07"
Data[grep("Agrobacterium CT08", Data[,"Bacteria"]), "Strains"] <- "<i>Agrobacterium</i> CT08"
Data[grep("Agrobacterium CT09", Data[,"Bacteria"]), "Strains"] <- "<i>Agrobacterium</i> CT09"
Data[grep("Agrobacterium CT10", Data[,"Bacteria"]), "Strains"] <- "<i>Agrobacterium</i> CT10"

Stat <- pairwise.wilcox.test(Data$Percent_of_crossed_J2s, Data$Strains, p.adjust.method = "fdr")

# Determines the statistical groups
Mymat <- tri.to.squ(Stat$p.value)

Letters <- multcompLetters(Mymat, compare = "<=", threshold = 0.05, Letters = letters)
Letters <- data.frame(group = names(Letters$Letters), letter = Letters$Letters)

## Colors
Fills <- c("white", "grey70", "grey70", "grey70", "grey70", "grey70")

## Plots the boxplots
p2 <- ggplot(Data)+
  geom_boxplot(aes(reorder(Strains, - Percent_of_crossed_J2s, mean), y = Percent_of_crossed_J2s, fill = reorder(Strains, - Percent_of_crossed_J2s, mean)), alpha = 0.65, color = "black")+
  scale_fill_manual("Strains", values = Fills, guide = "none")+
  geom_point(aes(x = Strains, y = Percent_of_crossed_J2s, color = reorder(Strains, - Percent_of_crossed_J2s, mean)), color = "black")+
  stat_summary(aes(x = Strains, y = Percent_of_crossed_J2s), fun = mean, colour = "red", geom = "point", shape = 3, size = 2)+
  geom_text(data = Letters, aes(label = letter, y = 125, x = group), colour = "black", size = 4)+
  ylab("J2 that passed through the sieve")+
  scale_y_continuous(breaks = seq(0, 150, 25))+
  ggtitle(" B")+
  xlab("")+
  guides(alpha = FALSE)+
  theme_bw()+
  theme(axis.text.x = element_markdown(size = 10, angle = 30, hjust = 1), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(face = "bold", hjust = -0.05))
p2

### Figure 3
Model <- "AAB"

p10 <- p1 + p2 + plot_layout(design = Model)
p10

ggsave(plot = p10, dpi = 1000, device = "jpeg", width = 12, height = 6, filename = "Results/Figure 3.jpeg")
