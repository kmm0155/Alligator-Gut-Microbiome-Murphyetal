####Beta_Diversity####

setwd("/Users/kamurphy/Desktop")
getwd()

#install these packages
install.packages('ggplot2')
install.packages('vegan')
install.packages('devtools')
install.packages('tidyverse')
install.packages('data.table')

#call on these libraries
library(ggplot2)
library(vegan)
library(devtools)
library(tidyverse)
library(data.table)

#####NMDS Plot########

#https://chrischizinski.github.io/rstats/vegan-ggplot2/
#https://ourcodingclub.github.io/tutorials/ordination/

#Plot for Figure 2C

#Call on your taxa file
datum_rc=read.csv(file.choose())
head(datum_rc)

#Call on metadata file that was also used for QIIME analyses
#Make sure your metadata file is in the same order as your taxa file
metadata=read.csv(file.choose())
head(metadata)

#Remove control rows from metadata file
metadata <- metadata[-c(14,17), ]

# transpose all but the first column (name)
datum_rc <- as.data.frame(t(datum_rc[,-1]))

#Change NAs in data set to 0
datum_rc[is.na(datum_rc)] <- 0

#Check the file
ncol(datum_rc)
nrow(datum_rc)
head(datum_rc)

#This is the NMDS analysis
#Default for Vegan is Bray-Curtis
datum_rc.dis <- metaMDS(datum_rc, k = 3, try = 50, trace = 1)
head(datum_rc.dis)
summary(datum_rc.dis)
stressplot(datum_rc.dis)
datum_rc.dis$stress

##Let's put these points together with the metadata sheet 
datum_rc.sub=as.data.frame(datum_rc.dis$points)
head(datum_rc.sub)
datum_rc.final=cbind(datum_rc.sub, Treatment = metadata$subject)
head(datum_rc.final)

#Plot for Figure 2B
datum_rc.final$TREAT_2 = factor(datum_rc.final$Treatment,c("Control", "Low", "High"))

tiff("Figure2.tiff", width = 8, height = 4, units = 'in', res = 300)

ggplot(datum_rc.final, aes(x=MDS1, y=MDS2, fill = TREAT_2)) +
  #stat_ellipse(size = 1) + 
  geom_polygon(aes(color = TREAT_2), alpha=0.30, show.legend = NA) +
  geom_point(aes(shape= TREAT_2, color = TREAT_2), size = 3.5) +
  scale_shape_manual(values=c("High" = 17, "Low" = 15, "Control" = 19), name = "Legend") +
  scale_fill_manual(values=c("High" = "#CC6677", "Low" = "#DDCC77", "Control" = "#44AA99"), name= "Legend") +
  scale_color_manual(values=c("High" = "#CC6677", "Low" = "#DDCC77", "Control" = "#44AA99"),name= "Legend") +
  theme_bw() +
  labs(title = "nMDS(Bray-Curtis) Alligator OTUs", legend = "Legend")

dev.off()

#####ANOSIM#####
#https://jkzorz.github.io/2019/06/11/ANOSIM-test.html

ano = anosim(datum_rc, metadata$subject, distance = "bray", permutations = 9999)
ano

#####PERMANOVA#####
#https://rpubs.com/hafezahmad/948799

#Preview test of similarity
datum_rc.veg <- vegdist(datum_rc, method="bray")
datum_rc.veg

#Summarize values
datum_rc.div<-adonis2(datum_rc.veg~as.factor(metadata$subject), data=datum_rc, permutations=9999, method="bray")
datum_rc.div