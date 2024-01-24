# Data Analysis
#Project: p451 run180802
#User(s): Stefanie Stadelmann (stadelms@student.ethz.ch)
#Data type: AmpSeq PE300 NodA (acyltransferase)
#Date   : 23.08.2018

rm(list=ls()) # clean/reset environment 
#list.files(getwd()) # Show the content of my working directory.

## Working directory
# e.g. setwd("your/path")

setwd("T:/ETH/Data/Pot_experiment/Rhizobia/nodA_work/test")

list.files(getwd())
wd <- getwd()

## Load homebrew-R functions
library(parallel)
source("Rfunctions.R")

## Load libraries 
#source("https://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
library(phyloseq)
#biocLite("ggplot2")
library(ggplot2)
library(vegan)
library(permute)
library(lattice)
#biocLite("microbiome")
library(microbiome)
#biocLite("ampvis2")
#source("https://bioconductor.org/biocLite.R")
#biocLite("ampvis2")
library(ampvis2)
library(knitr)
library(magrittr)
library(ggpubr)
library("RColorBrewer")


otufile      <- "p451_run180802_NodA_ZOTU_c99_Count_NCBI_nodA.txt"
mapfile      <- "p451_runrun180802_nodA_MapFile.txt"
treefile     <- "p451_run180802_NodA_ZOTU_c99_MSA.tre"
#refseqfile  <- "p451_run180802_NodA_ZOTU_c99.fa"

## Import data into Phyloseq
d <- import_qiime(otufilename = otufile, mapfilename = mapfile, treefilename = treefile)
d

## Summary
summarize_phyloseq(d)
# => Sparsity = 0.69 

## Number of counts per Taxa
sort(taxa_sums(d))

## Number of counts per Sample
sort(sample_sums(d))

## Read-Count Distribution - Careful this are still raw counts!

#pdf("Raw_Counts_per_sample.pdf", paper="a4")
#png("Raw_Counts_per_sample.png", width = 700, height = 350, units = "px")

plot(sample_sums(d), xaxt= "n", xlab="Sample", ylab="Number of Reads", pch=19, cex=c(sample_sums(d)/25000), col=rgb(0,0,1,alpha=0.5), main="Read-Counts per Sample", ylim=range(0,max(sample_sums(d))*1.1))
ac <- rep("blue",length(sample_sums(d)))
ac[get_variable(d, "SID") == "m"] <- "green"
points(sample_sums(d), pch=3, col=ac, cex=2)
v <- sample_names(d)
axis(side = 1, at = seq(1,length(sample_sums(d))), labels = v, tck=-.02, las=2, col.axis = 1)

#dev.off()

### Samples with low counts
head(sort(sample_sums(d)),20)
low.counts  <- head(sort(sample_sums(d)),4)
low.names   <- attr(low.counts, "names")
low.names
# => Samples with (very) low counts: r209  r074  r402  r265.   

# Remove samples with low counts:
all.names   <- sample_names(d)
high.names  <- all.names[!(all.names %in% low.names)]
d.high.S  <- prune_samples(high.names, d)

### Samples with missing meta data

### Rarefaction curve

#We can run all samples together but this might be a little bit messy. I would recommend to split the data into groups. 

D <- d.high.S 

## All sampels (this is a bit messy)
ggrare(D, step = 1000, color = "Farm", label = "X.SampleID", se = FALSE)

## Samples with low counts
D.low <- prune_samples(sample_sums(D) < 5000, D)
ggrare(D.low, step = 1000, color = "Farm", label = "X.SampleID", se = FALSE)

# > Samples below 2000 counts are problematic because diversity is underestimated. 

## Remove samples with counts below 3000
D.b2k <- prune_samples(sample_sums(D) > 3000, D)

#Note: It is important that a sample reaches the plateau otherwise we underestimate the diversity of the sample.

## Remove problematic OTUs

## Show tree
plot_tree(D.b2k, label.tips="taxa_names", "treeonly", nodeplotblank, ladderize=TRUE)

## The bad
badZOTU <- c("ZOTU152","ZOTU143","ZOTU133","ZOTU119","ZOTU124","ZOTU156","ZOTU160","ZOTU134","ZOTU44","ZOTU153","ZOTU128","ZOTU104","ZOTU127","ZOTU127","ZOTU103","ZOTU151","ZOTU157","ZOTU169","ZOTU146","ZOTU144")

## All
allZOTU  <- taxa_names(D.b2k)

## The good
goodZOUT <- allZOTU[!(allZOTU %in% badZOTU)]

## Keep only the good ones
D.b2k.goodZOTU <- prune_taxa(goodZOUT, D.b2k)
plot_tree(D.b2k.goodZOTU , label.tips="taxa_names", "treeonly", nodeplotblank, ladderize=TRUE)

D <- D.b2k.goodZOTU

##----------------------------------
#relative OTU abundance
##-------------------------------------
Dt <- transform_sample_counts(D, function(OTU) OTU/sum(OTU))

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

##-------------------------------------------------------------------------------------
### Ordination analysis
##------------------------------------------------------------------------------------
library(reshape2)
library(dplyr)
library(plyr)
library(tidyverse)

#Only by fertilization
fert = subset_samples(Dt, Fertilization == "fertilized")
unfert = subset_samples(Dt, Fertilization == "unfertilized")

fert.matrix = distance(fert, method="bray")
unfert.matrix = distance(unfert, method="bray")

fert.dbb = melt(as.matrix(fert.matrix))
fertilized = rep("fertilized", nrow(fert.dbb))
fert.anova = cbind(fert.dbb, fertilized)
fert.anova$Fertilization = fert.anova$fertilized

unfert.dbb = melt(as.matrix(unfert.matrix))
unfertilized = rep("unfertilized", nrow(unfert.dbb))
unfert.anova = cbind(unfert.dbb, unfertilized)
unfert.anova$Fertilization = unfert.anova$unfertilized

fertilization.dbb = rbind.fill(fert.anova, unfert.anova)
fertilization.dbb$BrayCurtis = fertilization.dbb$value

lmu = lm(fertilization.dbb$BrayCurtis ~ fertilization.dbb$Fertilization)
anova(lmu)
summary(lmu)
plot(fertilization.dbb$BrayCurtis ~ fertilization.dbb$Fertilization)

#Plot with 2xSE
library(Rmisc)
tgc <- summarySE(fertilization.dbb, measurevar="BrayCurtis", groupvars=c("Fertilization"))
tgc
ggplot(tgc, aes(x=Fertilization, y=BrayCurtis, colour = Fertilization)) + 
  geom_errorbar(aes(ymin=BrayCurtis-2*se, ymax=BrayCurtis+2*se), width=.1) +
  geom_point(size = 3) + theme_bw() + theme_classic()


#Only by soil
cult = subset_samples(Dt, Soil == "cultivated")
wild = subset_samples(Dt, Soil == "wild")
mix = subset_samples(Dt, Soil == "mixed")

cult.matrix = distance(cult, method="bray")
wild.matrix = distance(wild, method="bray")
mix.matrix = distance(mix, method="bray")

cult.dbb = melt(as.matrix(cult.matrix))
cultivated = rep("cultivated", nrow(cult.dbb))
cult.anova = cbind(cult.dbb, cultivated)
cult.anova$Soil = cult.anova$cultivated
mix.dbb = melt(as.matrix(mix.matrix))
mixed = rep("mixed", nrow(mix.dbb))
mix.anova = cbind(mix.dbb, mixed)
mix.anova$Soil = mix.anova$mixed
wild.dbb = melt(as.matrix(wild.matrix))
wild = rep("wild", nrow(wild.dbb))
wild.anova = cbind(wild.dbb, wild)
wild.anova$Soil = wild.anova$wild

soil.dbb = rbind.fill(cult.anova, mix.anova, wild.anova)
soil.dbb$BrayCurtis = soil.dbb$value

lmu = lm(soil.dbb$BrayCurtis ~ soil.dbb$Soil)
anova(lmu)
summary(lmu)
plot(soil.dbb$BrayCurtis ~ soil.dbb$Soil)
#Tukey's tests
a1 = aov(soil.dbb$BrayCurtis ~ soil.dbb$Soil)
tuk.soil = TukeyHSD(x=a1, "soil.dbb$Soil", conf.level=0.95)

#Plot with 2xSE
tgc <- summarySE(soil.dbb, measurevar="BrayCurtis", groupvars=c("Soil"))
tgc
ggplot(tgc, aes(x=Soil, y=BrayCurtis, colour = Soil)) + 
  geom_errorbar(aes(ymin=BrayCurtis-2*se, ymax=BrayCurtis+2*se), width=.1) +
  geom_point(size = 3)+ theme_bw() + theme_classic()



#Fertilization x Soil
c.fert = subset_samples(cult, Fertilization == "fertilized") 
w.fert = subset_samples(wild, Fertilization == "fertilized")
m.fert = subset_samples(mix, Fertilization == "fertilized")
c.unfert = subset_samples(cult, Fertilization == "unfertilized")
w.unfert = subset_samples(wild, Fertilization == "unfertilized")
m.unfert = subset_samples(mix, Fertilization == "unfertilized")

cult.fert.matrix = distance(c.fert, method="bray")
wild.fert.matrix = distance(w.fert, method="bray")
mix.fert.matrix = distance(m.fert, method="bray")

cultfert.dbb = melt(as.matrix(cult.fert.matrix))
cultivated.f = rep("cult.fert", nrow(cultfert.dbb))
cultfert.anova = cbind(cultfert.dbb, cultivated.f)
cultfert.anova$Soil = cultfert.anova$cultivated.f

mixfert.dbb = melt(as.matrix(mix.fert.matrix))
mixed.f = rep("mix.fert", nrow(mixfert.dbb))
mixfert.anova = cbind(mixfert.dbb, mixed.f)
mixfert.anova$Soil = mixfert.anova$mixed.f

wildfert.dbb = melt(as.matrix(wild.fert.matrix))
wild.f = rep("wild.fert", nrow(wildfert.dbb))
wildfert.anova = cbind(wildfert.dbb, wild.f)
wildfert.anova$Soil = wildfert.anova$wild.f

soil.dbb = rbind.fill(cultfert.anova, mixfert.anova, wildfert.anova)
soil.dbb$BrayCurtis = soil.dbb$value

lmu = lm(soil.dbb$BrayCurtis ~ soil.dbb$Soil)
anova(lmu)
summary(lmu)
plot(soil.dbb$BrayCurtis ~ soil.dbb$Soil)
#Tukey's tests
a1 = aov(soil.dbb$BrayCurtis ~ soil.dbb$Soil)
tuk.soil = TukeyHSD(x=a1, "soil.dbb$Soil", conf.level=0.95)

#Plot with 2xSE
tgc <- summarySE(soil.dbb, measurevar="BrayCurtis", groupvars=c("Soil"))
tgc
ggplot(tgc, aes(x=Soil, y=BrayCurtis, colour = Soil)) + 
  geom_errorbar(aes(ymin=BrayCurtis-2*se, ymax=BrayCurtis+2*se), width=.1) +
  geom_point(size = 3)

#Unfertilized
cult.unfert.matrix = distance(c.unfert, method="bray")
wild.unfert.matrix = distance(w.unfert, method="bray")
mix.unfert.matrix = distance(m.unfert, method="bray")

cultunfert.dbb = melt(as.matrix(cult.unfert.matrix))
cultivated.u = rep("cult.unfert", nrow(cultunfert.dbb))
cultunfert.anova = cbind(cultunfert.dbb, cultivated.u)
cultunfert.anova$Soil = cultunfert.anova$cultivated.u

mixunfert.dbb = melt(as.matrix(mix.unfert.matrix))
mixed.u = rep("mix.unfert", nrow(mixunfert.dbb))
mixunfert.anova = cbind(mixunfert.dbb, mixed.u)
mixunfert.anova$Soil = mixunfert.anova$mixed.u

wildunfert.dbb = melt(as.matrix(wild.unfert.matrix))
wild.u = rep("wild.unfert", nrow(wildunfert.dbb))
wildunfert.anova = cbind(wildunfert.dbb, wild.u)
wildunfert.anova$Soil = wildunfert.anova$wild.u

soil.dbb.u = rbind.fill(cultunfert.anova, mixunfert.anova, wildunfert.anova)
soil.dbb.u$BrayCurtis = soil.dbb.u$value

lmu = lm(soil.dbb.u$BrayCurtis ~ soil.dbb.u$Soil)
anova(lmu)
summary(lmu)
plot(soil.dbb.u$BrayCurtis ~ soil.dbb.u$Soil)
#Tukey's tests
a1 = aov(soil.dbb.u$BrayCurtis ~ soil.dbb.u$Soil)
tuk.soil = TukeyHSD(x=a1, "soil.dbb.u$Soil", conf.level=0.95)

#Plot with 2xSE
tgc <- summarySE(soil.dbb.u, measurevar="BrayCurtis", groupvars=c("Soil"))
tgc
ggplot(tgc, aes(x=Soil, y=BrayCurtis, colour = Soil)) + 
  geom_errorbar(aes(ymin=BrayCurtis-2*se, ymax=BrayCurtis+2*se), width=.1) +
  geom_point(size = 3)

#Fert and unfert together with soils
soil.dbb.all = rbind.fill(soil.dbb, soil.dbb.u)
soil.dbb.all = cbind(soil.dbb.all, c(rep("fertilized", 5638), rep("unfertilized", 5555)))
soil.dbb.all$Fertilization = soil.dbb.all$`c(rep("fertilized", 5638), rep("unfertilized", 5555))`
#Plot with 2xSE
tgc <- summarySE(soil.dbb.all, measurevar="BrayCurtis", groupvars=c("Soil", "Fertilization"))
tgc
ggplot(tgc, aes(x=Soil, y=BrayCurtis, colour = Fertilization)) + 
  geom_errorbar(aes(ymin=BrayCurtis-2*se, ymax=BrayCurtis+2*se), width=.1) +
  geom_point(size = 3)



##MOST ABUNDANT TAXA IN CULT-MIXED-WILD TO EXPLAIN BRAY-CURTIS PATTERNS
Dtc = subset_samples(Dt, Soil == "cultivated")
Dtm = subset_samples(Dt, Soil == "mixed")
Dtw = subset_samples(Dt, Soil == "wild")
#Choose only taxa with cumulative relabund >0.5%
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtc, f1, A=2)
prune_taxa(wh1, Dtc)
sort(taxa_sums(Dtc)) #ZOTU66, ZOTU138, ZOTU472, ZOTU17, ZOTU287, ZOTU7, ZOTU25, ZOTU54, ZOTU26, ZOTU68, ZOTU20, ZOTU47, ZOTU30, ZOTU59
topZOTU <- c("ZOTU66", "ZOTU138", "ZOTU472", "ZOTU17", "ZOTU287", "ZOTU7", "ZOTU25", "ZOTU54", "ZOTU26", "ZOTU68", "ZOTU20", "ZOTU47", "ZOTU30", "ZOTU59")
gyrB.c.top <- prune_taxa(topZOTU, Dtc)

f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtm, f1, A=2)
prune_taxa(wh1, Dtm)
sort(taxa_sums(Dtm)) #ZOTU59, ZOTU472, ZOTU7, ZOTU138, ZOTU54, ZOTU25, ZOTU96, ZOTU287, ZOTU68, ZOTU26, ZOTU66, ZOTU41, ZOTU30, ZOTU17  
topZOTU <- c("ZOTU59", "ZOTU472", "ZOTU7", "ZOTU138", "ZOTU54", "ZOTU25", "ZOTU96", "ZOTU287", "ZOTU68", "ZOTU26", "ZOTU66", "ZOTU41", "ZOTU30", "ZOTU17")
gyrB.m.top <- prune_taxa(topZOTU, Dtm)

#gyrB.w <- otu_table(Dtw)
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtw, f1, A=2)
prune_taxa(wh1, Dtw)
sort(taxa_sums(Dtw)) #ZOTU22, ZOTU7, ZOTU20, ZOTU68, ZOTU17, ZOTU287, ZOTU54, ZOTU472, ZOTU50, ZOTU138, ZOTU26 
topZOTU <- c("ZOTU22", "ZOTU7", "ZOTU20", "ZOTU68", "ZOTU17", "ZOTU287", "ZOTU54", "ZOTU472", "ZOTU50", "ZOTU138", "ZOTU26")
gyrB.w.top <- prune_taxa(topZOTU, Dtw)

#gyrB.c.top = merge_samples(gyrB.c.top, "Soil")
#plot_bar(gyrB.c.top, fill="OTU")
#gyrB.m.top = merge_samples(gyrB.m.top, "Soil")
#plot_bar(gyrB.m.top, fill="OTU")
#gyrB.w.top = merge_samples(gyrB.w.top, "Soil")
#plot_bar(gyrB.w.top, fill="OTU")

#Merge into single phyloseq object
all.soil = merge_phyloseq(otu_table(gyrB.c.top), otu_table(gyrB.m.top), otu_table(gyrB.w.top), sample_data(Dt))
all.soil.merged = merge_samples(all.soil, "Soil")
all.soil.merged = transform_sample_counts(all.soil.merged, function(OTU) OTU/sum(OTU))
plot_bar(all.soil.merged, fill="OTU")




##BY FARM
#BLOMFONTEIN
Dt.blo = subset_samples(Dt, Farm == "blo")
Dtc = subset_samples(Dt.blo, Soil == "cultivated")
Dtm = subset_samples(Dt.blo, Soil == "mixed")
Dtw = subset_samples(Dt.blo, Soil == "wild")
#Choose only taxa with cumulative relabund >0.5%
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtc, f1, A=2)
prune_taxa(wh1, Dtc)
sort(taxa_sums(Dtc)) #ZOTU138, ZOTU47, ZOTU59, ZOTU7,  ZOTU25, ZOTU54
topZOTU <- c("ZOTU138", "ZOTU47", "ZOTU59", "ZOTU7","ZOTU25", "ZOTU54")
gyrB.c.top <- prune_taxa(topZOTU, Dtc)
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtm, f1, A=2)
prune_taxa(wh1, Dtm)
sort(taxa_sums(Dtm)) #ZOTU26, ZOTU41, ZOTU59, ZOTU138, ZOTU54, ZOTU25
topZOTU <- c("ZOTU26", "ZOTU41", "ZOTU59", "ZOTU138", "ZOTU54", "ZOTU25")
gyrB.m.top <- prune_taxa(topZOTU, Dtm)
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtw, f1, A=2)
prune_taxa(wh1, Dtw)
sort(taxa_sums(Dtw)) #ZOTU50, ZOTU20, ZOTU287, ZOTU54, ZOTU25
topZOTU <- c("ZOTU50", "ZOTU20", "ZOTU287", "ZOTU54", "ZOTU25")
gyrB.w.top <- prune_taxa(topZOTU, Dtw)
#Merge into single phyloseq object
all.soil = merge_phyloseq(otu_table(gyrB.c.top), otu_table(gyrB.m.top), otu_table(gyrB.w.top), sample_data(Dt))
all.soil.merged = merge_samples(all.soil, "Soil")
all.soil.merged = transform_sample_counts(all.soil.merged, function(OTU) OTU/sum(OTU))
plot_bar(all.soil.merged, fill="OTU")

#DOBBELARSKOP
Dt.dob = subset_samples(Dt, Farm == "dob")
Dtc = subset_samples(Dt.dob, Soil == "cultivated")
Dtm = subset_samples(Dt.dob, Soil == "mixed")
Dtw = subset_samples(Dt.dob, Soil == "wild")
#Choose only taxa with cumulative relabund >0.5%
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtc, f1, A=2)
prune_taxa(wh1, Dtc)
sort(taxa_sums(Dtc)) 
topZOTU <- c("ZOTU17", "ZOTU287", "ZOTU7", "ZOTU25", "ZOTU54")
gyrB.c.top <- prune_taxa(topZOTU, Dtc)
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtm, f1, A=2)
prune_taxa(wh1, Dtm)
sort(taxa_sums(Dtm)) 
topZOTU <- c("ZOTU25", "ZOTU54")
gyrB.m.top <- prune_taxa(topZOTU, Dtm)
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtw, f1, A=2)
prune_taxa(wh1, Dtw)
sort(taxa_sums(Dtw)) 
topZOTU <- c("ZOTU17", "ZOTU54", "ZOTU25")
gyrB.w.top <- prune_taxa(topZOTU, Dtw)
#Merge into single phyloseq object
all.soil = merge_phyloseq(otu_table(gyrB.c.top), otu_table(gyrB.m.top), otu_table(gyrB.w.top), sample_data(Dt))
all.soil.merged = merge_samples(all.soil, "Soil")
all.soil.merged = transform_sample_counts(all.soil.merged, function(OTU) OTU/sum(OTU))
plot_bar(all.soil.merged, fill="OTU")

#LANDSKLOF
Dt.lan = subset_samples(Dt, Farm == "lan")
Dtc = subset_samples(Dt.lan, Soil == "cultivated")
Dtm = subset_samples(Dt.lan, Soil == "mixed")
Dtw = subset_samples(Dt.lan, Soil == "wild")
#Choose only taxa with cumulative relabund >0.5%
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtc, f1, A=2)
prune_taxa(wh1, Dtc)
sort(taxa_sums(Dtc)) 
topZOTU <- c("ZOTU138", "ZOTU66", "ZOTU7", "ZOTU472", "ZOTU25", "ZOTU54")
gyrB.c.top <- prune_taxa(topZOTU, Dtc)
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtm, f1, A=2)
prune_taxa(wh1, Dtm)
sort(taxa_sums(Dtm)) 
topZOTU <- c("ZOTU472", "ZOTU138","ZOTU25", "ZOTU54", "ZOTU7")
gyrB.m.top <- prune_taxa(topZOTU, Dtm)
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtw, f1, A=2)
prune_taxa(wh1, Dtw)
sort(taxa_sums(Dtw)) 
topZOTU <- c("ZOTU287", "ZOTU26", "ZOTU20","ZOTU17", "ZOTU54", "ZOTU25")
gyrB.w.top <- prune_taxa(topZOTU, Dtw)
#Merge into single phyloseq object
all.soil = merge_phyloseq(otu_table(gyrB.c.top), otu_table(gyrB.m.top), otu_table(gyrB.w.top), sample_data(Dt))
all.soil.merged = merge_samples(all.soil, "Soil")
all.soil.merged = transform_sample_counts(all.soil.merged, function(OTU) OTU/sum(OTU))
plot_bar(all.soil.merged, fill="OTU")

#MATARAKOPPIES
Dt.mat = subset_samples(Dt, Farm == "mat")
Dtc = subset_samples(Dt.mat, Soil == "cultivated")
Dtm = subset_samples(Dt.mat, Soil == "mixed")
Dtw = subset_samples(Dt.mat, Soil == "wild")
#Choose only taxa with cumulative relabund >0.5%
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtc, f1, A=2)
prune_taxa(wh1, Dtc)
sort(taxa_sums(Dtc)) 
topZOTU <- c("ZOTU20", "ZOTU7", "ZOTU287", "ZOTU25", "ZOTU54")
gyrB.c.top <- prune_taxa(topZOTU, Dtc)
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtm, f1, A=2)
prune_taxa(wh1, Dtm)
sort(taxa_sums(Dtm)) 
topZOTU <- c("ZOTU59","ZOTU25", "ZOTU54")
gyrB.m.top <- prune_taxa(topZOTU, Dtm)
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtw, f1, A=2)
prune_taxa(wh1, Dtw)
sort(taxa_sums(Dtw)) 
topZOTU <- c("ZOTU287", "ZOTU22", "ZOTU54", "ZOTU25")
gyrB.w.top <- prune_taxa(topZOTU, Dtw)
#Merge into single phyloseq object
all.soil = merge_phyloseq(otu_table(gyrB.c.top), otu_table(gyrB.m.top), otu_table(gyrB.w.top), sample_data(Dt))
all.soil.merged = merge_samples(all.soil, "Soil")
all.soil.merged = transform_sample_counts(all.soil.merged, function(OTU) OTU/sum(OTU))
plot_bar(all.soil.merged, fill="OTU")

#MELKKRAAL
Dt.mel = subset_samples(Dt, Farm == "mel")
Dtc = subset_samples(Dt.mel, Soil == "cultivated")
Dtm = subset_samples(Dt.mel, Soil == "mixed")
Dtw = subset_samples(Dt.mel, Soil == "wild")
#Choose only taxa with cumulative relabund >0.5%
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtc, f1, A=2)
prune_taxa(wh1, Dtc)
sort(taxa_sums(Dtc)) 
topZOTU <- c("ZOTU472", "ZOTU17", "ZOTU25", "ZOTU54")
gyrB.c.top <- prune_taxa(topZOTU, Dtc)
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtm, f1, A=2)
prune_taxa(wh1, Dtm)
sort(taxa_sums(Dtm)) 
topZOTU <- c("ZOTU138","ZOTU25", "ZOTU54")
gyrB.m.top <- prune_taxa(topZOTU, Dtm)
f1  <- filterfun_sample(topk(5))
wh1 <- genefilter_sample(Dtw, f1, A=2)
prune_taxa(wh1, Dtw)
sort(taxa_sums(Dtw)) 
topZOTU <- c("ZOTU17","ZOTU7", "ZOTU68", "ZOTU54", "ZOTU25")
gyrB.w.top <- prune_taxa(topZOTU, Dtw)
#Merge into single phyloseq object
all.soil = merge_phyloseq(otu_table(gyrB.c.top), otu_table(gyrB.m.top), otu_table(gyrB.w.top), sample_data(Dt))
all.soil.merged = merge_samples(all.soil, "Soil")
all.soil.merged = transform_sample_counts(all.soil.merged, function(OTU) OTU/sum(OTU))
plot_bar(all.soil.merged, fill="OTU")
