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
  x <- Dt
  num <- 3
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
#write.table(as.matrix(test), "top_abundant_taxa.txt")
toptax = read.csv("top_abundant_taxa.csv", header = T, sep=",", row.names = 1) 
# toptax = as.data.frame(t(toptax))
tax123 = toptax[,94]
tax123 = toptax[,94:96]
sample_data(Dt) = cbind(sample_data(Dt), tax123)
filtering = tax123
filtering = c(as.character(tax123[,1]), as.character(tax123[,2]), as.character(tax123[,3]))
dominant.OTU <- subset(otu_table(Dt), rownames(otu_table(Dt)) %in% filtering)
dom.physeq <- merge_phyloseq(dominant.OTU, tax_table(Dt), sample_data(Dt), phy_tree(Dt))

#Proportion of reads top3 dominant OTUs by treatment
#Fertilized
fert = subset_samples(dom.physeq, Fertilization =="fertilized")
fert = as.data.frame(cbind(sample_data(fert), as.data.frame(sample_sums(fert))))
fert$RelativeAbundance = fert$`sample_sums(fert)`
#unfertilized
unfert = subset_samples(dom.physeq, Fertilization =="unfertilized")
unfert = as.data.frame(cbind(sample_data(unfert), as.data.frame(sample_sums(unfert))))
unfert$RelativeAbundance = unfert$`sample_sums(unfert)`
#Cultivated
cult = subset_samples(dom.physeq, Soil =="cultivated")
cult = as.data.frame(cbind(sample_data(cult), as.data.frame(sample_sums(cult))))
cult$RelativeAbundance = cult$`sample_sums(cult)`
#Mixed
mix = subset_samples(dom.physeq, Soil =="mixed")
mix = as.data.frame(cbind(sample_data(mix), as.data.frame(sample_sums(mix))))
mix$RelativeAbundance = mix$`sample_sums(mix)`
#Wild
wild = subset_samples(dom.physeq, Soil =="wild")
wild = as.data.frame(cbind(sample_data(wild), as.data.frame(sample_sums(wild))))
wild$RelativeAbundance = wild$`sample_sums(wild)`
#All
dom.physeq.all = as.data.frame(cbind(sample_data(dom.physeq), as.data.frame(sample_sums(dom.physeq))))
dom.physeq.all$RelativeAbundance = dom.physeq.all$`sample_sums(dom.physeq)`
##Plots
#All
plot(dom.physeq.all$RelativeAbundance ~ dom.physeq.all$Fertilization)
plot(dom.physeq.all$RelativeAbundance ~ dom.physeq.all$Soil)
plot(dom.physeq.all$RelativeAbundance ~ dom.physeq.all$Farm)
plot(dom.physeq.all$RelativeAbundance ~ dom.physeq.all$Farm.Soil)
#Fertilization
plot(fert$RelativeAbundance ~ fert$Soil)
plot(unfert$RelativeAbundance ~ unfert$Soil)
#Cult
plot(cult$RelativeAbundance ~ cult$Fertilization)
plot(cult$RelativeAbundance ~ cult$Farm)
#Mixed
plot(mix$RelativeAbundance ~ mix$Fertilization)
plot(mix$RelativeAbundance ~ mix$Farm)
#Wild
plot(wild$RelativeAbundance ~ wild$Fertilization)
plot(wild$RelativeAbundance ~ wild$Farm)

##BY FARM
#Blomfontein
blo = subset_samples(dom.physeq, Farm =="blo")
blo.c = subset_samples(dom.physeq, Soil =="cultivated")
blo.m = subset_samples(dom.physeq, Soil =="mixed")
blo.w = subset_samples(dom.physeq, Soil =="wild")
#OTU tables Blomfontein
otus.blo = as.data.frame(otu_table(blo))
otus.blo.c = as.data.frame(otu_table(blo.c))
otus.blo.m = as.data.frame(otu_table(blo.m))
otus.blo.w = as.data.frame(otu_table(blo.w))

plot_heatmap(blo.c, method = "NMDS", distance = "bray", sample.order = "Fertilization")
plot_heatmap(blo.m, method = "NMDS", distance = "bray", sample.order = "Fertilization")
plot_heatmap(blo.w, method = "NMDS", distance = "bray", sample.order = "Fertilization")

otus.blo.c = sort()




#Subset phyloseq object to keep only most abundant ZOTUs
a = subset_samples(Dt, tax1 == "ZOTU43") #Dobbelarskop wild ZOTU
b = subset_samples(Dt, tax1 == "ZOTU6") #Blomfontein cultivated ZOTU
c = subset_samples(Dt, tax1 == "ZOTU28") #Dobbelarskop wild ZOTU
d = subset_samples(Dt, tax1 == "ZOTU23") #Dobbelarskop cultivated ZOTU
e = subset_samples(Dt, tax1 == "ZOTU99") #Blomfontein wild ZOTU
f = subset_samples(Dt, tax1 == "ZOTU7") #Melkkraal all ZOTU
g = subset_samples(Dt, tax1 == "ZOTU72") #Melkkraal cultivated ZOTU in 1 sample
h = subset_samples(Dt, tax1 == "ZOTU18") #Melkkraal cultivated ZOTU in 2 samples
i = subset_samples(Dt, tax1 == "ZOTU27") #Melkkraal mostly cultivated ZOTU 
j = subset_samples(Dt, tax1 == "ZOTU22") #Widespread mostly wild ZOTU 
k = subset_samples(Dt, tax1 == "ZOTU155") #Landskloof and Matarakopies cultivation ZOTU
l = subset_samples(Dt, tax1 == "ZOTU12") #!!!!Matarakopies mixed fertilized ZOTU
m = subset_samples(Dt, tax1 == "ZOTU42") #Widespread ZOTU 
n = subset_samples(Dt, tax1 == "ZOTU21") #Melkkraal wild ZOTU 
o = subset_samples(Dt, tax1 == "ZOTU97") #Landskloof cultivated ZOTU in 1 sample 
p = subset_samples(Dt, tax1 == "ZOTU86") #Landskloof cultivated ZOTU 
q = subset_samples(Dt, tax1 == "ZOTU112") #Landskloof cultivated ZOTU
r = subset_samples(Dt, tax1 == "ZOTU31") #Blomfontein wild ZOTU in 2 samples
s = subset_samples(Dt, tax1 == "ZOTU26") #Blomfontein cultivated fertilized ZOTU
t = subset_samples(Dt, tax1 == "ZOTU19") #Dobbelarskop mixed fertilized ZOTU
u = subset_samples(Dt, tax1 == "ZOTU63") #Matarakopies wild unfertilized ZOTU in 1 sample
v = subset_samples(Dt, tax1 == "ZOTU50") #Blomfontein mixed fertilized ZOTU in 1 sample
w = subset_samples(Dt, tax1 == "ZOTU98") #Blomfontein mixed unfertilized ZOTU in 1 sample
x = subset_samples(Dt, tax1 == "ZOTU76") #Blomfontein wild unfertilized ZOTU in 1 sample
y = subset_samples(Dt, tax1 == "ZOTU40") #Dobbelarskop cultivated ZOTU in 2 samples
z = subset_samples(Dt, tax1 == "ZOTU14") #Matarakoppies wild ZOTU 
ba = subset_samples(Dt, tax1 == "ZOTU33") #Melkkraal mixed unfertilized ZOTU in 1 sample
bb = subset_samples(Dt, tax1 == "ZOTU166") #Blomfontein mixed unfertilized ZOTU in 1 sample
bc = subset_samples(Dt, tax1 == "ZOTU45") #Melkkraal cultivated fertilized ZOTU in 1 sample

dominant.OTU <- subset(otu_table(Dt), rownames(otu_table(Dt)) %in% c("ZOTU43", "ZOTU6", "ZOTU28", "ZOTU23","ZOTU99", "ZOTU7","ZOTU72","ZOTU18", "ZOTU27","ZOTU22","ZOTU155","ZOTU12","ZOTU42","ZOTU21","ZOTU97","ZOTU86","ZOTU112","ZOTU31","ZOTU26","ZOTU19","ZOTU63","ZOTU50","ZOTU98","ZOTU76","ZOTU40","ZOTU14","ZOTU33","ZOTU166","ZOTU45"))
dom.physeq <- merge_phyloseq(dominant.OTU, tax_table(Dt), sample_data(Dt), phy_tree(Dt))

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

#Recalculate evenness for only dominant strains
even <- evenness(dom.physeq, "all")
pielou = even$pielou
sample_data(dom.physeq) = as.data.frame(cbind(sample_data(dom.physeq), pielou))

#Ordination
#NMDS farm
ord <- ordinate(dom.physeq, "NMDS", "bray")
ord <- ordinate(dom.physeq, "MDS", "unifrac", weighted = T)
plot_ordination(dom.physeq, ord, color = "Farm") +
  geom_point(size = 5)
#NMDS soil
plot_ordination(dom.physeq, ord, color = "Soil") +
  geom_point(size = 5)
#NMDS Farm-Soil
p1 <- plot_ordination(physeq=dom.physeq, ordination=ord, shape="Soil", color="Farm")
p1 + geom_point(size = 5) + ggtitle("Ordination Plots (NMDS / Bray-Curtis)")
#NMDS Fert-Soil
p1 <- plot_ordination(physeq=dom.physeq, ordination=ord, shape="Soil", color="Fertilization")
p1 + geom_point(size = 5) + ggtitle("Ordination Plots (NMDS / Bray-Curtis)")

#Facet wrap
#By farm
p2 = plot_ordination(dom.physeq, ord, color="Soil", shape="Fertilization")
print(p2)
p2 + facet_wrap(~Farm) + geom_point(size=4) + ggtitle("NMDS for the different Farms")
#By soil origin
p3 = plot_ordination(dom.physeq, ord, color="Farm", shape="Fertilization")
p3
p3 + facet_wrap(~Soil) + geom_point(size=4) + ggtitle("NMDS for the different soil origins")

#PERMANOVA
meta = meta(dom.physeq)
adonis(distance(Dt, method="bray")~Soil*Farm, data=meta)
adonis(distance(Dt, method="bray")~Farm*Fertilization, data=meta)
adonis(distance(Dt, method="bray")~Soil*Fertilization, data=meta)

###INDICATOR SPECIES ANALYSIS###
#install.packages("indicspecies")
library(indicspecies)
#Indicator species by farm
indicspecies.farm = multipatt(as.data.frame(t(otu_table(dom.physeq))), cluster = sample_data(dom.physeq)$Farm, control = how(nperm = 9999))
summary(indicspecies.farm)
write.csv(summary(indicspecies.farm), "indic_spp_farm_dominant.csv")
#Indicator species by soil origin
indicspecies.soil = multipatt(as.data.frame(t(otu_table(dom.physeq))), cluster = sample_data(dom.physeq)$Soil, control = how(nperm = 9999))
indsoil = summary(indicspecies.soil)
#write.csv(summary(indicspecies.soil), "indic_spp_soil_dominant.csv")
# #Indicator species by Farm and Soil 
# #write.csv(sample_data(Dt), "data_Dt_concat.csv")
# newDt = read.csv("data_Dt_concat.csv", header = T, sep = ",")
# indicspecies.fs = multipatt(as.data.frame(t(otu_table(Dt))), cluster = newDt$Farm.Soil, control = how(nperm = 9999))
# indsoil = summary(indicspecies.fs)
#Indicator species by fertilization
indicspecies.fert = multipatt(as.data.frame(t(otu_table(dom.physeq))), cluster = sample_data(dom.physeq)$Fertilization, control = how(nperm = 9999))
summary(indicspecies.fert)
#Indicator species by Farm and Soil within unfertilized
# Dtu = subset_samples(dom.physeq, Fertilization == "unfertilized")
# Dtu = subset_samples(Dtu, Farm == "mel")
# #newDtu = read.csv("data_Dtu_concat.csv", header = T, sep = ",")
# newDtu = read.csv("data_Dtu_Mel_concat.csv", header = T, sep = ",")
# indicspecies.unfert = multipatt(as.data.frame(t(otu_table(Dtu))), cluster = Dtu$Soil, control = how(nperm = 9999))
# indunfert = summary(indicspecies.unfert)


##-------------------------------------------------------------------------------------
#correlation of ZOTU relative abundance to plant response using indicator species
##-------------------------------------------------------------------------------------
OTU_id = as.matrix(rownames(tax_table(dom.physeq)))
indicfarm = read.table("list_farm_dom.txt", sep = "", header = T)
species_indicator = indicfarm[,1]
#OTU_id = as.matrix(rownames(tax_table(Dtu)))
ind.spp.subset = subset_taxa(dom.physeq, OTU_id%in%as.matrix(species_indicator))
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

###Scatter plots
##Evenness
#Pielou-DM!!!
qplot(x=pielou, y=plant.DM, color=Soil, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#Pielou-tiller length!!!
qplot(x=pielou, y=tiller.total.length, color=Soil, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#Pielou-nodule number!!!
qplot(x=pielou, y=nodule.number, color=Soil, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#Pielou-root diameter!!!
qplot(x=pielou, y=root.diameter, color=Soil, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#Pielou-delta15N!!!No effects
qplot(x=pielou, y=Leaf.d15N, color=Soil, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#Pielou-Nodule 15N!!!
qplot(x=pielou, y=Nodule.15N, color=Soil, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#Pielou-Nodule total N!!!
qplot(x=pielou, y=Nodule.totalN, color=Soil, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 


##Farm
qplot(x=pielou, y=plant.DM, color=Farm, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 


##Soil
#!!!OTU166 indicator for mixed from Blo cultivations
qplot(x=ZOTU166, y=plant.DM, color=Soil, shape = Farm, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#!!!OTU19 indicator for mixed from Dob cultivations
qplot(x=ZOTU19, y=Leaf.d15N, color=Soil, shape = Farm, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#!!!OTU17 indicator for mixed from Lan wild only
qplot(x=ZOTU17, y=plant.DM, color=Soil, shape = Farm, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 


###TESTING TREATMENT EFFECTS ON EVENNESS##
##Testing (ANOVA + TUKEY'S HSD)

#Pielou-Soil
dom.even = sample_data(dom.physeq)
lm = lm(dom.even$pielou ~ dom.even$Soil)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "dom.even$Soil", conf.level=0.95)
#Boxplot
boxplot(dom.even$pielou ~ dom.even$Soil, na.rm = T, ylab = "Evenness (Pielou's index)", xlab = "Soil")   

#Pielou-Farm
dom.even = sample_data(dom.physeq)
lm = lm(dom.even$pielou ~ dom.even$Farm)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "dom.even$Farm", conf.level=0.95)
#Boxplot
boxplot(dom.even$pielou ~ dom.even$Farm, na.rm = T, ylab = "Evenness (Pielou's index)", xlab = "Farm")   

#Pielou-Fertilization
dom.even = sample_data(dom.physeq)
lm = lm(dom.even$pielou ~ dom.even$Fertilization)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "dom.even$Fertilization", conf.level=0.95)
#Boxplot
boxplot(dom.even$pielou ~ dom.even$Fertilization, na.rm = T, ylab = "Evenness (Pielou's index)", xlab = "Fertilization")   

#Pielou-DM
dom.even = sample_data(dom.physeq)
lm = lm(dom.even$plant.DM ~ dom.even$pielou)
lmu = aov(lm)
summary(lmu)

#Pielou-leaf15N
dom.even = sample_data(dom.physeq)
lm = lm(dom.even$Leaf.d15N ~ dom.even$pielou)
lmu = aov(lm)
summary(lmu)

#Pielou-nodule15N
dom.even = sample_data(dom.physeq)
lm = lm(dom.even$Nodule.totalN ~ dom.even$pielou)
lmu = aov(lm)
summary(lmu)

#Pielou-leaf N content
dom.even = sample_data(dom.physeq)
lm = lm(dom.even$Leaf.N.content ~ dom.even$pielou)
lmu = aov(lm)
summary(lmu)


##TREE OF INDICATOR SPECIES AND IMPORT INTO ITOL##
#prune OTUs excluded in indicator species analysis
#FARM
OTU_id = as.matrix(rownames(tax_table(dom.physeq)))
indicfarm = read.table("list_farm_dom.txt", sep = "", header = T)
species_indicator = indicfarm[,1]
big_tree = phy_tree(dom.physeq)
pruned.tree.indic = drop.tip(big_tree,big_tree$tip.label[-match(species_indicator, big_tree$tip.label)])
plot(pruned.tree.indic)

#SOIL
indicsoil = read.table("indic_spp_soil_dominant.txt", sep = "", header = T)
species_indicator = indicsoil[,1]
big_tree = phy_tree(dom.physeq)
pruned.tree.indicsoil = drop.tip(big_tree,big_tree$tip.label[-match(c("ZOTU166", "ZOTU26", "ZOTU63", "ZOTU155","ZOTU19","ZOTU98","ZOTU33"), big_tree$tip.label)])
plot(pruned.tree.indicsoil)

#DM
indicdm = read.table("indicspp_DM.txt", sep = "", header = T)
species_indicator = indicdm[,1]
big_tree = phy_tree(Dt)
pruned.tree.indicdm = drop.tip(big_tree,big_tree$tip.label[-match(as.data.frame(species_indicator), big_tree$tip.label)])
plot(pruned.tree.indicdm)
#save as newick
write.tree(pruned.tree.indic, "indicspp_tree_farm_dominant.tre")
write.tree(pruned.tree.indicsoil, "indicspp_tree_soil_dominant.tre")

indicmelkunfert = read.table("list_melk_unfert.txt", sep = "", header = T)
indsp.melkunfert = indicmelkunfert[,1]
