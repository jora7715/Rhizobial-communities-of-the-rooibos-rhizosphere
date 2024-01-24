# Data Analysis
#Project: p451 run180802
#User(s): Stefanie Stadelmann (stadelms@student.ethz.ch)
#Data type: AmpSeq PE300 NodA (acyltransferase)
#Date   : 23.08.2018

rm(list=ls()) # clean/reset environment 
setwd("T:/ETH/Data/Pot_experiment/Rhizobia/nodA_work/R")

list.files(getwd())
wd <- getwd()

## Load homebrew-R functions
library(parallel)
source("Rfunctions.R")

## Load libraries 
#source("https://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
library(phyloseq)
library(vegan)
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

## Import data into Phyloseq
otufile      <- "p451_run180802_NodA_ZOTU_c99_Count_NCBI_nodA.txt"
mapfile      <- "p451_runrun180802_nodA_MapFile.txt"
treefile     <- "p451_run180802_NodA_ZOTU_c99_MSA.tre"
#refseqfile  <- "p451_run180802_NodA_ZOTU_c99.fa"

d <- import_qiime(otufilename = otufile, mapfilename = mapfile, treefilename = treefile)
d

##---------------------------------------------------------
## Samples with low counts and removal of bad OTUs
##----------------------------------------------------------
head(sort(sample_sums(d)),20)
low.counts  <- head(sort(sample_sums(d)),4)
low.names   <- attr(low.counts, "names")
low.names
# => Samples with (very) low counts: r209  r074  r402  r265.   

# Remove samples with low counts:
all.names   <- sample_names(d)
high.names  <- all.names[!(all.names %in% low.names)]
d.high  <- prune_samples(high.names, d)

### Samples with missing meta data

# Remove samples without meta data
d.high.S <- subset_samples(d.high, SID != "m")

### Rarefaction curve

#We can run all samples together but this might be a little bit messy. I would recommend to split the data into groups. 

D <- d.high.S 

## All samples 
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
#Transform sample counts
##-------------------------------------
Dt <- transform_sample_counts(D, function(OTU) OTU/sum(OTU))

# square root transformation (Hellinger)
sqrt.Dt <- transform_sample_counts(Dt, function(OTU) sqrt(OTU))



#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

####Find the most abundant taxa per plant####
#Goes through a phyloseq object, picks out the most abundant taxa and gives the abundance for each
#and identifies which taxa is most abundant for which sample
# Function to find the top most abundant ASV/OTU per sample, number of ASV/OTU per sample depends on user input. This is particularly relevant for finding strain level community structure when a few genera dominates the communities, for  example "is it a single variant of Pseudomonas dominating all the samples?"
find.top.asv <- function(x,num){
  # x <- Dt
  # num <- 3
  otu <- t(otu_table(x))
  tax <- t(tax_table(x))
  j1 <- apply(otu,1,sort,index.return=T, decreasing=T) # modifying which.max to return a list of sorted index
  j2 <- lapply(j1,'[[',"ix") # select for index
  
  # l <- data.frame(unique(tax@.Data[unlist(j2),]))
  m <- data.frame(otu@.Data[,unique(unlist(j2))])
  n <- apply(m,1,sort,index.return=T, decreasing=T) %>%
    lapply('[[',"ix") %>%  # Extract index
    lapply(head,n=num) # This to returns the top x tax
  
  p <- list()
  for(i in 1:length(n)){
    p[[i]]<- colnames(m)[n[[i]]]
  }
  m$taxa <- p
  return(m)
}
test = find.top.asv(Dt,3)
write.table(as.matrix(test), "top_abundant_taxa.txt")
toptax = read.csv("top_abundant_taxa.csv", header = T, sep=",", row.names = 1) 
toptax = as.data.frame(toptax)
tax123 = toptax[,94:96]
sample_data(Dt) = cbind(sample_data(Dt), tax123)

#Subset phyloseq object to keep only most abundant ZOTUs
dominant.OTU <- subset(otu_table(Dt), rownames(otu_table(Dt)) %in% c("ZOTU43", "ZOTU6", "ZOTU28", "ZOTU23","ZOTU99", "ZOTU7","ZOTU72","ZOTU18", "ZOTU27","ZOTU22","ZOTU155","ZOTU12","ZOTU42","ZOTU21","ZOTU97","ZOTU86","ZOTU112","ZOTU31","ZOTU26","ZOTU19","ZOTU63","ZOTU50","ZOTU98","ZOTU76","ZOTU40","ZOTU14","ZOTU33","ZOTU166","ZOTU45"))
dom.physeq <- merge_phyloseq(dominant.OTU, tax_table(Dt), sample_data(Dt))

###INDICATOR SPECIES ANALYSIS###
#install.packages("indicspecies")
library(indicspecies)
#Indicator species by farm
indicspecies.farm = multipatt(as.data.frame(t(otu_table(Dt))), cluster = sample_data(Dt)$Farm, control = how(nperm = 9999))
summary(indicspecies.farm)
write.csv(summary(indicspecies.farm), "indic_spp_farm.csv")
#Indicator species by soil origin
indicspecies.soil = multipatt(as.data.frame(t(otu_table(Dt))), cluster = sample_data(Dt)$Soil, control = how(nperm = 9999))
indsoil = summary(indicspecies.soil)
#Indicator species by Farm and Soil 
#write.csv(sample_data(Dt), "data_Dt_concat.csv")
newDt = read.csv("data_Dt_concat.csv", header = T, sep = ",")
indicspecies.fs = multipatt(as.data.frame(t(otu_table(Dt))), cluster = newDt$Farm.Soil, control = how(nperm = 9999))
indsoil = summary(indicspecies.fs)
#Indicator species by fertilization
indicspecies.fert = multipatt(as.data.frame(t(otu_table(Dt))), cluster = sample_data(Dt)$Fertilization, control = how(nperm = 9999))
summary(indicspecies.fert)
#Indicator species by Farm and Soil within unfertilized
Dtu = subset_samples(Dt, Fertilization == "unfertilized")
Dtu = subset_samples(Dtu, Farm == "mel")
#newDtu = read.csv("data_Dtu_concat.csv", header = T, sep = ",")
newDtu = read.csv("data_Dtu_Mel_concat.csv", header = T, sep = ",")
indicspecies.unfert = multipatt(as.data.frame(t(otu_table(Dtu))), cluster = newDtu$Soil, control = how(nperm = 9999))
indunfert = summary(indicspecies.unfert)


#Indicator species by plant response category
indicspecies.DM = multipatt(as.data.frame(t(otu_table(Dt))), cluster = sample_data(Dt)$plant.DM.class, control = how(nperm = 9999))
summary(indicspecies.DM)

##TREE OF INDICATOR SPECIES AND IMPORT INTO ITOL##
#prune OTUs excluded in indicator species analysis
#FARM
indicfarm = read.table("list_farm.txt", sep = "", header = T)
species_indicator = indicfarm[,1]
big_tree = phy_tree(Dt)
pruned.tree.indic = drop.tip(big_tree,big_tree$tip.label[-match(species_indicator, big_tree$tip.label)])
plot(pruned.tree.indic)
#SOIL
indicsoil = read.table("list_soil.txt", sep = "", header = T)
species_indicator = indicsoil[,1]
big_tree = phy_tree(Dt)
pruned.tree.indicsoil = drop.tip(big_tree,big_tree$tip.label[-match(species_indicator, big_tree$tip.label)])
plot(pruned.tree.indicsoil)
#DM
indicdm = read.table("indicspp_DM.txt", sep = "", header = T)
species_indicator = indicdm[,1]
big_tree = phy_tree(Dt)
pruned.tree.indicdm = drop.tip(big_tree,big_tree$tip.label[-match(species_indicator, big_tree$tip.label)])
plot(pruned.tree.indicsoil)
#save as newick
#write.tree(pruned.tree.indic, "indicspp_tree_farm.tre")
#write.tree(pruned.tree.indicsoil, "indicspp_tree_soil.tre")

indicmelkunfert = read.table("list_melk_unfert.txt", sep = "", header = T)
indsp.melkunfert = indicmelkunfert[,1]

##-------------------------------------------------------------------------------------
#correlation of ZOTU relative abundance to plant response using indicator species
##-------------------------------------------------------------------------------------
OTU_id = as.matrix(rownames(tax_table(Dt)))
indicfarm = read.table("list_farm.txt", sep = "", header = T)
species_indicator = indicfarm[,1]
#OTU_id = as.matrix(rownames(tax_table(Dtu)))
ind.spp.subset = subset_taxa(Dt, OTU_id%in%as.matrix(species_indicator))
#ind.spp.subset = subset_taxa(Dtu, OTU_id%in%as.matrix(indsp.melkunfert))
#Farm
df1 <- data.frame(sample_data(ind.spp.subset))
df2 <- data.frame(otu_table(ind.spp.subset))
df2.t <- t(df2) #transpose df2
df_corr <- cbind(df1, df2.t) #combine by column (r001...)
#Soil
df1 <- data.frame(sample_data(ind.spp.subset))
df2 <- data.frame(otu_table(ind.spp.subset))
df2.t <- t(df2) #transpose df2
df_corr <- cbind(df1, df2.t) #combine by column (r001...)
#DM
df1 <- data.frame(sample_data(ind.spp.subset))
df2 <- data.frame(otu_table(ind.spp.subset))
df2.t <- t(df2) #transpose df2
df_corr <- cbind(df1, df2.t) #combine by column (r001...)
###Scatter plot
##Farm
qplot(x=ZOTU45, y=plant.DM, color=Farm, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
##Soil
#!!!OTU166 indicator for mixed from Blo cultivations
qplot(x=ZOTU166, y=plant.DM, color=Soil, shape = Farm, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#!!!OTU19 indicator for mixed from Dob cultivations
qplot(x=ZOTU19, y=Leaf.d15N, color=Soil, shape = Farm, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#!!!OTU17 indicator for mixed from Lan wild only
qplot(x=ZOTU17, y=plant.DM, color=Soil, shape = Farm, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 

#Plots using indicator ZOTUs for plant DM: ZOTU61, ZOTU62, ZOTU161, ZOTU165
qplot(x=ZOTU61, y=plant.DM, color=Soil, shape = Farm, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
qplot(x=ZOTU62, y=plant.DM, color=Soil, shape = Farm, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
qplot(x=ZOTU161, y=plant.DM, color=Soil, shape = Farm, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
qplot(x=ZOTU165, y=plant.DM, color=Soil, shape = Farm, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 

qplot(x=ZOTU149, y=plant.DM, color=Soil, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 



