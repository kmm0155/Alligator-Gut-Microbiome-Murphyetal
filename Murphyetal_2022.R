####AlligatorExp####

setwd("/Users/kaitlynmurphy/Desktop")
getwd()

#install these packages
install.packages('lme4')
install.packages('nlme')
install.packages('ggplot2')
install.packages('vegan')
install.packages('devtools')
install.packages('tidyverse')
install.packages('lsmeans')

#call on these libraries
library(lme4)
library(nlme)
library(ggplot2)
library(vegan)
library(devtools)
library(tidyverse)
library(lsmeans)

#####Rarefaction Curve####

#https://fromthebottomoftheheap.net/2015/04/16/drawing-rarefaction-curves-with-custom-colours/

#Call on your data set for the rarefaction curve. This data frame can be found in "Alligator_Sequencing_Data.xlsx"
datum_rc=read.csv(file.choose())
head(datum_rc)

BCI2 <- datum_rc[1:23, ]
raremax <- min(rowSums(BCI2))
raremax

#raremax is the minimum sample count achieved over the 23 samples. 
#We will rarefy the sample counts to this value.

pars <- expand.grid(col = "black", lty = "solid", stringsAsFactors = FALSE)
head(pars)

#Note that I save the output from rarecurve() in object out. This object contains everything we need to draw our own version of the plot if we wish. 
#For example, we could use fewer colors and alter the line thickness1 instead to make up the required number of combinations.

#Plot for Supplementary Figure 1

tiff("rarefaction.tiff", width = 8, height = 4, units = 'in', res = 300)

out <- with(pars[1:23, ],
            rarecurve(BCI2, step = 20, col = "black",
                      lty = "solid", label = TRUE, ylab = "Rarefied No. of bacterial groups", xlab = "Observed No. of bacterial groups"))
dev.off()

#####################Shannons############

#http://rstudio-pubs-static.s3.amazonaws.com/11473_29c2de401c55441ca5350862ebd04456.html

#This is the main dataframe from Alligator_Phenotypic_Data
datum=read.csv(file.choose())
head(datum)

resultssh=lm(SHANNONS~TREAT,data=datum,na.action=na.omit)
summary(resultssh)

resultssh1=lm(SHANNONS~relevel(TREAT,ref = "High"),data=datum,na.action=na.omit)
summary(resultssh1)

#Plot for Figure 2A

datum$TREAT1 = factor(datum$TREAT,c("Control", "Low", "High"))

tiff("Shannons.tiff", width = 4, height = 4, units = 'in', res = 300)

boxplot(SHANNONS~TREAT1,data=datum, notch=FALSE, 
        col=(c("#44AA99", "#DDCC77", "#CC6677")),
        xlab="Treatment", ylab="Shannon's Diversity Index",
        outline=FALSE)
stripchart(SHANNONS ~ TREAT1, vertical = TRUE, data = datum, 
           method = "jitter", add=TRUE, pch = 20)

dev.off()

##########Pielou's evenness##########

#https://cran.r-project.org/web/packages/vegan/vignettes/diversity-vegan.pdf
#You will use the data set from the rarefaction curve you've already called on here

H<-diversity(datum_rc)

J <- H/log(specnumber(datum_rc))
J

#added output to my original datafile
res.aovp <- aov(PIE ~ TREAT, data = datum, na.action=na.omit)
# Summary of the analysis
summary(res.aovp)

resultssh=lm(PIE~TREAT,data=datum,na.action=na.omit)
summary(resultssh)

resultssh1=lm(PIE~relevel(TREAT,ref = "High"),data=datum,na.action=na.omit)
summary(resultssh1)

#Plot for Supplementary Figure 2B

datum$TREAT1 = factor(datum$TREAT,c("Control", "Low", "High"))

tiff("Pielous.tiff", width = 4, height = 4, units = 'in', res = 300)

boxplot(PIE~TREAT1,data=datum, notch=FALSE, 
        col=(c("#44AA99", "#DDCC77", "#CC6677")),
        xlab="Treatment", ylab="Pielou's Evenness Index",
        outline=FALSE)
stripchart(PIE ~ TREAT1, vertical = TRUE, data = datum, 
           method = "jitter", add=TRUE, pch = 20)

dev.off()

#########OTU#######

resultsotu=lm(OT~TREAT,data=datum,na.action=na.omit)
summary(resultsotu)
confint(resultsotu)
  
resultsotu1=lm(OT~relevel(TREAT,ref = "High"),data=datum,na.action=na.omit)
summary(resultsotu1)

#Plot for Figure 2A
  
datum$TREAT1 = factor(datum$TREAT,c("Control", "Low", "High"))

tiff("OTU.tiff", width = 4, height = 4, units = 'in', res = 300)

boxplot(OT~TREAT1,data=datum, notch=FALSE, 
         col=(c("#44AA99", "#DDCC77", "#CC6677")),
          xlab="Treatment", ylab="Operational taxonomic units (OTUs)",
          outline=FALSE)
stripchart(OT ~ TREAT1, vertical = TRUE, data = datum, 
           method = "jitter", add=TRUE, pch = 20)

dev.off()

##############Abundances###########

#Call on the dataset from Alligator_Sequence_Data, tab "Main_LongFormat"
taxa=read.csv(file.choose())
head(taxa)

resultstaxa=lme(ABUND~TREAT,random=~1|ID,data=phyla,na.action=na.omit)
summary(resultstaxa)

###Phyla###

#Call on the dataset from Alligator_Taxonomy_Data, tab "Phylum"
phyla=read.csv(file.choose())
head(phyla)

#Firmicutes
# Compute the analysis of variance
res.aov1 <- aov(Firmicutes ~ TREAT, data = phyla, na.action=na.omit)
# Summary of the analysis
summary(res.aov1)

firm=lm(Firmicutes ~ TREAT, data = phyla, na.action=na.omit)
summary(firm)
firm1=lm(Firmicutes ~ relevel(TREAT,ref = "High"), data = phyla, na.action=na.omit)
summary(firm1)

plot(Firmicutes~TREAT, data=phyla)

#Fusobacteria
# Compute the analysis of variance
res.aov2 <- aov(Fusobacteria ~ TREAT, data = phyla, na.action=na.omit)
# Summary of the analysis
summary(res.aov2)

fuso=lm(Fusobacteria ~ TREAT, data = phyla, na.action=na.omit)
summary(fuso)
fuso1=lm(Fusobacteria ~ relevel(TREAT,ref = "High"), data = phyla, na.action=na.omit)
summary(fuso1)

plot(Fusobacteria~TREAT, data=phyla)

#Proteobacteria
# Compute the analysis of variance
res.aov3 <- aov(Protobacteria ~ TREAT, data = phyla, na.action=na.omit)
# Summary of the analysis
summary(res.aov3)

TukeyHSD(res.aov3)

proteo=lm(Protobacteria ~ TREAT, data = phyla, na.action=na.omit)
summary(proteo)
proteo1=lm(Protobacteria ~ relevel(TREAT,ref = "High"), data = phyla, na.action=na.omit)
summary(proteo1)

plot(Proteobacteria~TREAT, data=phyla)

TukeyHSD(res.aov3)

#Plot
#https://stackoverflow.com/questions/22305023/how-to-get-a-barplot-with-several-variables-side-by-side-grouped-by-a-factor
#http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization

#Call on dataset in long format from "Alligator_Taxonomy_Data", tab "Phyla_Long_Format
phyla1=read.csv(file.choose())
head(phyla1)

phyla1$TREAT_1 = factor(phyla1$TREAT,c("Control", "Low", "High"))

tiff("phyla.tiff", width = 4, height = 5, units = 'in', res = 300)

ggplot(data=phyla1,aes(x=TREAT_1,y=ABUND,fill=factor(PHYLA))) +
  geom_bar(stat="identity",position="dodge") +
  theme_bw() +
  theme(legend.position="bottom") +
  scale_fill_manual(name="Phylum",
                    breaks=c(1,2,3),
                    values=c("grey", "orange","blue"),
                    labels=c("Firmicutes", "Fusobacteria", "Proteobacteria")) +
  xlab("Treatment")+ylab("OTU Abundance")

dev.off()

###Genus###

#Call on the dataset from Alligator_Taxonomy_Data, tab "Genera"
genus=read.csv(file.choose())
head(genus)

#Clostridium
# Compute the analysis of variance
res.aov4 <- aov(Clostridium ~ TREAT, data = genus, na.action=na.omit)
# Summary of the analysis
summary(res.aov4)

clos=lm(Clostridium ~ TREAT, data = genus, na.action=na.omit)
summary(clos)
clos1=lm(Clostridium ~ relevel(TREAT,ref = "High"), data = genus, na.action=na.omit)
summary(clos1)

#Cetobacterium
# Compute the analysis of variance
res.aov5 <- aov(Cetobacterium ~ TREAT, data = genus, na.action=na.omit)
# Summary of the analysis
summary(res.aov5)

ceto=lm(Cetobacterium ~ TREAT, data = genus, na.action=na.omit)
summary(clos)
ceto1=lm(Cetobacterium ~ relevel(TREAT,ref = "High"), data = genus, na.action=na.omit)
summary(ceto1)

#Defluviitalea
# Compute the analysis of variance
res.aov6 <- aov(Defluviitalea ~ TREAT, data = genus, na.action=na.omit)
# Summary of the analysis
summary(res.aov6)

defl=lm(Defluviitalea ~ TREAT, data = genus, na.action=na.omit)
summary(defl)
defl1=lm(Defluviitalea ~ relevel(TREAT,ref = "High"), data = genus, na.action=na.omit)
summary(defl1)

#Helicobacter
# Compute the analysis of variance
res.aov6 <- aov(Helicobacter ~ TREAT, data = genus, na.action=na.omit)
# Summary of the analysis
summary(res.aov6)

helico=lm(Helicobacter ~ TREAT, data = genus, na.action=na.omit)
summary(helico)
helico1=lm(Helicobacter ~ relevel(TREAT,ref = "High"), data = genus, na.action=na.omit)
summary(helico1)

#Call on the dataset from Alligator_Taxonomy_Data, tab "Genera_Long_Format"
genus1=read.csv(file.choose())
head(genus1)

genus1$TREAT_1 = factor(genus1$TREAT,c("Control", "Low", "High"))

tiff("genus.tiff", width = 4, height = 5, units = 'in', res = 300)

ggplot(data=genus1,aes(x=TREAT_1,y=ABUND,fill=factor(GENERA))) +
  geom_bar(stat="identity",position="dodge") +
  theme_bw() +
  theme(legend.position="right") +
  scale_fill_manual(name="Genus",
                    breaks=c(1,2,3,4),
                    values=c("grey", "orange", "blue", "black"),
                    labels=c("Clostridium", "Cetobacterium", "Defluviitalea", "Helicobacter")) +
  xlab("Treatment")+ylab("OTU Abundance")

dev.off()

############Piechart##################

#Plot for Figure 2C

#https://www.displayr.com/how-to-make-a-pie-chart-in-r/

#Control

#Call on the dataset from "Alligator_Taxonomy_Data", tab "PieChart_Control
Control=read.csv(file.choose())
head(Control)

safe_colorblind_palette <- c("#44AA99", "#999933", "#882255", "#661100","#52854C", "#6699CC", "#888888",
                             "#88CCEE", "#CC6677", "#332288", "#DDCC77", "#117733", "#D55E00")

tiff("Control.tiff", width = 8 , height = 3, units = 'in', res = 300)

ggplot(Control, aes(x= "", y=ABUND, fill=GROUP)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    scale_fill_manual(values=safe_colorblind_palette) +
    labs(x = NULL, y = NULL, fill = NULL, title = "Control OTUs") +
    theme_classic() +
    theme(axis.line = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    plot.title = element_text(hjust = 0.5, color = "black")) +
     theme(legend.text = element_text(size= 5.5), 
            legend.key.size = unit(0.5, 'cm'), 
            legend.key.height = unit(0.5, 'cm'), 
            legend.key.width = unit(0.5, 'cm'))

dev.off()

#Low Treatment

#Call on the dataset from "Alligator_Taxonomy_Data", tab "PieChart_Low

Low=read.csv(file.choose())
head(Low)

tiff("Low.tiff", width = 8, height = 3, units = 'in', res = 300)

ggplot(Low, aes(x= "", y=ABUND, fill=GROUP)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values=safe_colorblind_palette) +
  labs(x = NULL, y = NULL, fill = NULL, title = "Low OTUs") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "black")) +
  theme(legend.text = element_text(size= 5.5), 
        legend.key.size = unit(0.5, 'cm'), 
        legend.key.height = unit(0.5, 'cm'), 
        legend.key.width = unit(0.5, 'cm'))

dev.off()

#High Treatment

#Call on the dataset from "Alligator_Taxonomy_Data", tab "PieChart_High
High=read.csv(file.choose())
head(High)

tiff("High.tiff", width = 8, height = 3, units = 'in', res = 300)

ggplot(High, aes(x= "", y=ABUND, fill=GROUP)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values=safe_colorblind_palette) +
  labs(x = NULL, y = NULL, fill = NULL, title = "High OTUs") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "black")) +
  theme(legend.text = element_text(size= 5.5), 
        legend.key.size = unit(0.5, 'cm'), 
        legend.key.height = unit(0.5, 'cm'), 
        legend.key.width = unit(0.5, 'cm'))

dev.off()

#####################NMDS Plot################################

#https://chrischizinski.github.io/rstats/vegan-ggplot2/
#https://ourcodingclub.github.io/tutorials/ordination/

#Plot for Figure 2C

#Call on metadata file that was also used for QIIME analyses
metadata=read.csv(file.choose())

#Use dataset used for rarefaction curve

names(datum_rc)
#for the number of species
ncol(datum_rc)
#for the number of samples
nrow(datum_rc)
nrow(metadata)

#Change NA to 0
datum_rc[is.na(datum_rc)] <- 0

#make sure samples are in same order as metadata
head(metadata)
head(datum_rc)

#Transpose it because Vegan is gay and likes it setup this way (Sample = rows, taxnomic=columns)
#vt_MDS= as.data.frame(t(sheet2))
#otherwise run:
Vt_MDS = datum_rc
head (Vt_MDS)

##Actually run the nMDS, default is Bray, but you can change that (see manpage)
###########IF THERE IS A PROBLEM AT THIS STEP ABOUT NUMERIC DATA, go back to your original file and delete the first column of the number series###
#V.nMDSAG5 <- metaMDS(Vt_MDS, k = 2, trymax = 100, trace = F)- the default (bray-curtis) is the same as this
V.nMDSAG5 <- metaMDS(Vt_MDS)
head(V.nMDSAG5)
summary(V.nMDSAG5)
stressplot(V.nMDSAG5)

##Let's put these points together with the metadata sheet so we can graph some shit 
V.nMDSnumbersAG5=as.data.frame(V.nMDSAG5$points)
head(V.nMDSnumbersAG5)
V.nMDSnummetaAG5=cbind(V.nMDSnumbersAG5, Treatment = metadata$subject)
#V.nMDSnummetaAG5=cbind(V.nMDSnumbersAG5,metadata$subject)
head(V.nMDSnummetaAG5)

#################################plot by Treatment
tiff("Figure2.tiff", width = 8, height = 4, units = 'in', res = 300)

ggplot(V.nMDSnummetaAG5, aes(x=MDS1, y=MDS2, color= factor(metadata$subject))) +
  #stat_ellipse(size = 1) + 
  geom_polygon(data=V.nMDSnummetaAG5,aes(x=MDS1,y=MDS2,fill= factor(metadata$subject),group= factor(metadata$subject)),alpha=0.30,show.legend = NA) + # add the convex hulls
  scale_fill_manual(values=c("High" = "#CC6677", "Low" = "#DDCC77", "Control" = "#44AA99")) +
  geom_point(aes(shape = factor(metadata$subject)), size = 3.5) +
  scale_colour_manual(values=c("High" = "#CC6677", "Low" = "#DDCC77", "Control" = "#44AA99")) +
  theme_bw() +
  labs(title = "nMDS(Bray-Curtis) Alligator OTUs")

dev.off()

#ANOSIM
#https://jkzorz.github.io/2019/06/11/ANOSIM-test.html

m_com = as.matrix(datum_diversity)
pc=cbind(datum_diversity, Treatment = metadata$subject)
head(pc)

ano = anosim(m_com, pc$Treatment, distance = "bray", permutations = 9999)
ano

####Estradiol####

results=lm(ESTRADIOL~TREAT+SEX+WEIGHT,data=datum,na.action=na.omit)
summary(results)
results=lm(ESTRADIOL~relevel(TREAT,ref = "High")+SEX+WEIGHT,data=datum)
summary(results)

resultssex=lm(ESTRADIOL~TREAT+relevel(SEX, ref= "M")+WEIGHT,data=datum)
summary(resultssex)

datum$TREAT_2 = factor(datum$TREAT,c("Control", "Low", "High"))

tiff("Figure1.tiff", width = 4, height = 4, units = 'in', res = 300)

boxplot(ESTRADIOL~TREAT_2, data=datum,
        col=(c("#44AA99", "#DDCC77", "#CC6677")),
        xlab="Treatment", ylab="Estradiol concentration (pg/mL)",
        outline=FALSE)
stripchart(ESTRADIOL~TREAT_2, vertical = TRUE, data = datum, 
           method = "jitter", add=TRUE, pch = 20)

dev.off()
