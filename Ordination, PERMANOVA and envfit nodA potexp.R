# Data Analysis
#Project: p451 run180802
#User(s): Stefanie Stadelmann (stadelms@student.ethz.ch)
#Data type: AmpSeq PE300 NodA (acyltransferase)
#Date   : 23.08.2018

rm(list=ls()) # clean/reset environment 
setwd("E:/ETH/Data/Pot_experiment/Rhizobia/nodA_work/test")

list.files(getwd())
wd <- getwd()

## Load homebrew-R functions
library(parallel)
source("Rfunctions.r")

## Load libraries 
#source("https://bioconductor.org/biocLite.R")
#biocLite("vegan")
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
#d.high.S <- subset_samples(d.high, SID != "m")

### Rarefaction curve

#We can run all samples together but this might be a little bit messy. I would recommend to split the data into groups. 

d <- d.high 

## All samples 
ggrare(d, step = 200, color = "Fertilization", label = "X.SampleID", se = FALSE)

## Samples with low counts
d <- prune_samples(sample_sums(d) < 200, d)
ggrare(d, step = 1000, color = "Farm", label = "X.SampleID", se = FALSE)

## Remove problematic OTUs

## Show tree
#plot_tree(D.b2k, label.tips="taxa_names", "treeonly", nodeplotblank, ladderize=TRUE)

## The bad
badZOTU <- c("ZOTU152","ZOTU143","ZOTU133","ZOTU119","
             4","ZOTU156","ZOTU160","ZOTU134","ZOTU44","ZOTU153","ZOTU128","ZOTU104","ZOTU127","ZOTU127","ZOTU103","ZOTU151","ZOTU157","ZOTU169","ZOTU146","ZOTU144")

## All
allZOTU  <- taxa_names(d)

## The good
goodZOUT <- allZOTU[!(allZOTU %in% badZOTU)]

## Keep only the good ones
D.b2k.goodZOTU <- prune_taxa(goodZOUT, d)
plot_tree(D.b2k.goodZOTU , label.tips="taxa_names", "treeonly", nodeplotblank, ladderize=TRUE)
D <- D.b2k.goodZOTU
D = subset_samples(D, X.SampleID != "r498")

##----------------------------------
#Transform sample counts
##-------------------------------------
Dt <- transform_sample_counts(D, function(OTU) OTU/sum(OTU))

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


#Keep only OTUs with >1% relative abundance of the total samples
Dt = filter_taxa(Dt, function(OTU) sum(OTU) > .001, TRUE)

# square root transformation (Hellinger)
sqrt.Dt <- transform_sample_counts(Dt, function(OTU) sqrt(OTU))

#scale to even sampling depth (data transformation)
min_depth <- min(sample_sums(D))
D.scale = transform_sample_counts(D, function(x) min_depth * x/sum(x))

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

Dt = subset_samples(Dt, Soil=="mixed")
# Ordinate data with PCoA
ord <- ordinate(Dt, "NMDS", "bray")
ord <- ordinate(D, "NMDS", "bray")


plot_ordination(D, ord,  color = "Farm") +
  geom_point(size = 5)+ theme_bw() + theme_classic()

##NMDS has >0.25 stress!!!!
#NMDS farm
library(wesanderson)
ord <- ordinate(Dt, "NMDS", "bray")
ord <- ordinate(Dt, "MDS")
# Multidimensional scaling (MDS / PCoA)
plot_ordination(Dt, ord, color = "Farm") + scale_color_manual(values=wes_palette(n=5, name="Cavalcanti1")) + geom_point(size = 5) + theme_bw() 

p + stat_ellipse(type="norm")
facet_wrap(~Farm)

#NMDS soil
ord <- ordinate(Dt, "NMDS", "bray")
ord <- ordinate(Dt, "MDS", "unifrac", weighted = T)
# Multidimensional scaling (MDS / PCoA)
plot_ordination(Dt, ord, color = "Soil") +
  geom_point(size = 5) + theme_bw()

#NMDS fertilization
ord <- ordinate(Dt, "NMDS", "bray")
ord <- ordinate(Dt, "MDS", "unifrac", weighted = T)
# Multidimensional scaling (MDS / PCoA)
plot_ordination(Dt, ord, color = "Fertilization") + geom_point(size = 5)

#NMDS Soil and Farm
ord <- ordinate(Dt, method="NMDS", distance="bray")
ord <- ordinate(Dt, "MDS", "unifrac", weighted = T)
p1 <- plot_ordination(physeq=Dt, ordination=ord, shape="Soil", color="Farm")
p1 + geom_point(size = 5) + ggtitle("Ordination Plots (NMDS / Bray-Curtis)")

#NMDS Farm and Fertilization
ord <- ordinate(Dt, method="NMDS", distance="bray")
p1 <- plot_ordination(physeq=Dt, ordination=ord, shape="Fertilization", color="Farm")
p1 + geom_point(size = 5) + ggtitle("Ordination Plots (NMDS / Bray-Curtis)")

#NMDS Fertilization and Soil
ord <- ordinate(Dt, method="NMDS", distance="bray")
p1 <- plot_ordination(physeq=Dt, ordination=ord, shape="Soil", color="Fertilization")
p1 + geom_point(size = 5) + ggtitle("Ordination Plots (NMDS / Bray-Curtis)")

#Facet wrap
#By farm
p2 = plot_ordination(Dt, ord, color="Soil")
p2 + facet_wrap(~Farm) + geom_point(size=4) + ggtitle("NMDS for the different Farms")
#By soil origin
p3 = plot_ordination(Dt, ord, color="Farm", shape="Fertilization")
p3
p3 + facet_wrap(~Soil) + geom_point(size=4) + ggtitle("NMDS for the different soil origins")

#add the ellipses to the NMDS
ordn=ordinate(Dt, "NMDS", "bray")
p6=plot_ordination(Dt, ordn, color="Soil", shape="Farm")
p6 + geom_point(size=4) + ggtitle("All samples") + facet_wrap(~Farm)
p6 + stat_ellipse(type="norm") + theme_bw()



#--------------------------------------------------------------------------------------
#PERMANOVA
#--------------------------------------------------------------------------------------
library(pairwiseAdonis)
#Interactions
meta = meta(Dt)
meta = na.omit(meta)

bcmatrix = distance(Dt, method = "bray")
bcmatrix = na.omit(bcmatrix)

adonis <- adonis(bcmatrix ~ Block / Farm / Soil * Fertilization, strata = meta$Block, data = meta)

adonis(distance(Dt, method = "bray") ~ Soil, data=meta, permutations = 9999)
adonis(vegdist(otu_table(Dt), "bray")~Farm*Fertilization, data=meta, permutations = 9999)
adonis(vegdist(t(as.data.frame(otu_table(Dt))), method="bray")~Soil*Fertilization, data=meta, permutations = 9999)
adonis(distance(Dt, method="bray")~Soil*Fertilization, data=meta, permutations = 9999)
adonis(distance(Dt, method="bray")~Altitude, data=meta, permutations = 9999)
adonis(distance(Dt, method="bray")~Latitude, data=meta, permutations = 9999)
adonis(distance(Dt, method="bray")~Longitude, data=meta, permutations = 9999)
adonis(distance(Dt, method="bray")~cluster.presence, data=meta, permutations = 9999)

##Pairwise adonis!!!!!
data.pa = as.data.frame(sample_data(Dt))
dist.pa = distance(Dt, method="bray")
pairwise.adonis(dist.pa, factors = data.pa$Soil, sim.method = "bray", p.adjust.m = "none")
pairwise.adonis(dist.pa, factors = data.pa$Farm, sim.method = "bray", p.adjust.m = "none")
pairwise.adonis(dist.pa, factors = data.pa$Fertilization, sim.method = "bray", p.adjust.m = "none")
pairwise.adonis(dist.pa, factors = data.pa$Farm.Soil, sim.method = "bray", p.adjust.m = "none")


#new Permanova from Jean-Claude

#scale to even sampling depth (data transformation)
min_depth <- min(sample_sums(D))
D.scale = transform_sample_counts(D, function(x) min_depth * x/sum(x))

bdist <- distance(D.scale, "bray")
dataframe <- data.frame(sample_data(D))

fert <- as(sample_data(D.scale), "data.frame")[,"Fertilization"]
soil <- as(sample_data(D.scale), "data.frame")[,"Soil"]
farm <- as(sample_data(D.scale), "data.frame")[,"Farm"]

ordn=ordinate(D.scale, "NMDS", "bray")
p6=plot_ordination(D.scale, ordn, color="Soil", shape="Soil")
p6 + geom_point(size=4) + ggtitle("All samples") + facet_wrap(~Farm)
p6 + stat_ellipse(type="norm") + theme_bw()
#Adonis Test
adonis.bdist.fert <- adonis(bdist~fert)
adonis.bdist.soil <- adonis(bdist~soil)
adonis.bdist.farm <- adonis(bdist~farm)

adonis.bdist.run_fert_soil_farm <- adonis(bdist ~ fert + soil + farm)

print(adonis.bdist.run_fert_soil_farm$aov.tab)
capture.output(adonis.bdist.run_fert_soil_farm$aov.tab, file="permanova.doc")

#Homogeneity of dispersion test
betatax.fert <- betadisper(bdist, fert)
betatax.soil <- betadisper(bdist, soil)
betatax.farm <- betadisper(bdist, farm)


#--------------------------------------------------------------------------------------
#REGRESSION OF NMDS SCORES AND PLANT RESPONSES 
#--------------------------------------------------------------------------------------

#Obtain NMDS axis scores
library(vegan)
ord.dt <- ordinate(Dt, method="NMDS", distance="bray")
scores.ord = scores(ord.dt, choices =c(1,2), display = "sites")
Dtu = subset_samples(Dt, Fertilization == "unfertilized")

#Get data frame of plant responses and unify with NMDS scores
Plant.response = as.data.frame(sample_data(Dt))
Plant.response.scores = cbind(scores.ord, Plant.response)
write.csv(Plant.response.scores, "NMDS_plant_response_dbb.csv")

#Plot interesting regressions 
#Plant.DM
reg = ggplot(Plant.response.scores, aes(x=NMDS1, y=plant.DM, color = Soil)) +
  geom_point(size=5) + geom_smooth(method="lm", se = T) 
reg
#Above DM
reg = ggplot(Plant.response.scores, aes(x=NMDS1, y=above.DM, color = Soil)) +
  geom_point(size=5) + geom_smooth(method="lm", se = T) 
reg
#Leaf DM
reg = ggplot(Plant.response.scores, aes(x=NMDS1, y=leaf.DM, color = Soil)) +
  geom_point(size=5) + geom_smooth(method="lm", se = T) 
reg
#Stem DM
reg = ggplot(Plant.response.scores, aes(x=NMDS1, y=stem.DM, color = Soil)) +
  geom_point(size=5) + geom_smooth(method="lm", se = T) 
reg
#Root DM
reg = ggplot(Plant.response.scores, aes(x=NMDS1, y=root.DM, color = Soil)) +
  geom_point(size=5) + geom_smooth(method="lm", se = T) 
reg
#Ratio.above.root
reg = ggplot(Plant.response.scores, aes(x=NMDS2, y=Ratio.above.root, color = Soil)) +
  geom_point(size=5) + geom_smooth(method="lm", se = T) 
reg
#Leaf.Pconc
reg = ggplot(Plant.response.scores, aes(x=NMDS1, y=Leaf.Pconc..mg.g., color = Soil)) +
  geom_point(size=5) + geom_smooth(method="lm", se = T) 
reg
#Leaf.Pcont
reg = ggplot(Plant.response.scores, aes(x=NMDS2, y=Leaf.Pcont, color = Soil)) +
  geom_point(size=5) + geom_smooth(method="lm", se = T) 
reg
#Leaf.Nconc
reg = ggplot(Plant.response.scores, aes(x=NMDS1, y=Leaf.Nconc.mg.g., color = Soil)) +
  geom_point(size=5) + geom_smooth(method="lm", se = T) 
reg
#Leaf.NPratio
reg = ggplot(Plant.response.scores, aes(x=NMDS2, y=N.P.ratio, color = Soil)) +
  geom_point(size=5) + geom_smooth(method="lm", se = T) 
reg
#Leafd13C
reg = ggplot(Plant.response.scores, aes(x=NMDS2, y=Leaf.d13C, color = Soil)) +
  geom_point(size=5) + geom_smooth(method="lm", se = T) 
reg
#Leaf.Cu.conc
reg = ggplot(Plant.response.scores, aes(x=NMDS1, y=Leaf.Cu..ug.g., color = Soil)) +
  geom_point(size=5) + geom_smooth(method="lm", se = T) 
reg
#Leaf.15N
reg = ggplot(Plant.response.scores, aes(x=NMDS1, y=Leaf.15N.standardized, color = Soil)) +
  geom_point(size=5) + geom_smooth(method="lm", se = T) 
reg
#Root.15N
reg = ggplot(Plant.response.scores, aes(x=NMDS1, y=Root.15N, color = Soil)) +
  geom_point(size=5) + geom_smooth(method="lm", se = T) 
reg
#!!!Nodule15N
reg = ggplot(Plant.response.scores, aes(x=NMDS2, y=Nodule.15N, color = Soil)) +
  geom_point(size=5) + geom_smooth(method="lm", se = T) 
reg


##ENVIRONMENTAL VARIABLE FITTING TO NMDS##
Dtu = subset_samples(Dt, Fertilization == "unfertilized")
ord = ordinate(Dt, "NMDS", "bray")
ord = ordinate(Dt, "MDS", "unifrac", weighted = T)
#All variables together with PCoA!!!!!!!!
envfit(ord ~ sample_data(Dt)$ln.tiller.total.length + sample_data(Dt)$ln.leaf.DM + sample_data(Dt)$ln.above.DM + sample_data(Dt)$ln.root.DM + sample_data(Dt)$ln.plant.DM + sample_data(Dt)$ln.nodule.number + sample_data(Dt)$ln.nodules.rootDM + sample_data(Dt)$ln.Leaf.Pcont + sample_data(Dt)$ln.Leaf.N.content + sample_data(Dt)$DeltaDelta.leaf.15N + sample_data(Dt)$ln.Leaf.Mn..ug.g., na.rm = T, permutations = 9999)
envfit(ord$vectors ~ sample_data(Dt)$Altitude + sample_data(Dt)$Latitude + sample_data(Dt)$Longitude, na.rm = T, permutations = 9999)

#Pot soil N concentration
soilN = envfit(ord, sample_data(Dt)$Pot.soil.N.conc, permutations = 9999)
#Pot soil 15N
envfit(ord, sample_data(Dt)$Pot.soil.15N, permutations = 9999)
#Pot soil P concentration
envfit(ord, sample_data(Dt)$Pot.soil.P.conc, permutations = 9999)
#Pot soil NP ratio
envfit(ord, sample_data(Dt)$Pot.NPratio, permutations = 9999)
#pH
envfit(ord, sample_data(Dt)$pH, permutations = 9999, na.rm = T)
#Altitude
envfit(ord, sample_data(Dt)$Altitude, permutations = 9999, na.rm = T)
#Latitude
envfit(ord, sample_data(Dt)$Latitude, permutations = 9999, na.rm = T)
#Longitude
envfit(ord, sample_data(Dt)$Longitude, permutations = 9999, na.rm = T)
#CEC
envfit(ord, sample_data(Dt)$CEC.cmol.kg., permutations = 9999, na.rm = T)
#WHC
envfit(ord, sample_data(Dt)$WHC...10kPa., permutations = 9999, na.rm = T)
#Coarse sand
envfit(ord, sample_data(Dt)$Coarse.sand..., permutations = 9999, na.rm = T)




##PAIRWISE NMDS DIFFERENCES AND PLANT RESPONSES
Dtu = subset_samples(Dt, Fertilization == "unfertilized")
#Plant responses
response = sample_data(Dt)
#Euclidean distance matrices (for correlation plots)
dist.ordination.NMDS1 = dist(scores.ord[,1], method = "euclidean")
dist.ordination.NMDS2 = dist(scores.ord[,2], method = "euclidean")


#Distance matrices for environmental factors
#Geography
dist.latitude = dist(response$Latitude, method = "euclidean", na.action(na.omit)) #latitude
dist.longitude = dist(response$Longitude, method = "euclidean", na.action(na.omit)) #longitude
dist.km = read.csv("geodist_original.csv", header = T, row.names = 1)
dist.altitude = dist(response$Altitude, method = "euclidean", na.action(na.omit)) #longitude

#Soil factors
dist.CEC = dist(response$CEC.cmol.kg., method = "euclidean", na.action(na.omit)) #CEC
dist.WHC = dist(response$WHC...10kPa., method = "euclidean", na.action(na.omit)) #WHC
dist.sand = dist(response$Coarse.sand..., method = "euclidean", na.action(na.omit)) #sand coarseness
dist.pH = dist(response$pH, method = "euclidean", na.action(na.omit)) #pH
dist.soilN = dist(response$Pot.soil.N.conc, method = "euclidean", na.action(na.omit)) #soil N conc
dist.soilP = dist(response$Pot.soil.P.conc, method = "euclidean", na.action(na.omit)) #soil P conc
dist.soilCa = dist(response$Soil.Ca.mg.g., method = "euclidean", na.action(na.omit)) #soil Ca conc
dist.soilK = dist(response$Soil.K.mg.g., method = "euclidean", na.action(na.omit)) #soil K conc
dist.soilMg = dist(response$Soil.Mg.mg.g., method = "euclidean", na.action(na.omit)) #soil Mg conc
dist.soilMn = dist(response$Soil.Mn.mg.g., method = "euclidean", na.action(na.omit)) #soil Mn conc
dist.orgC = dist(response$OrgC.soil.g, method = "euclidean", na.action(na.omit)) #soil organic C
#Plant response
#Biomass accumulation
dist.tillerlength = dist(response$tiller.total.length, method = "euclidean", na.action(na.omit)) #Tiller total length
dist.leafDM = dist(response$leaf.DM, method = "euclidean", na.action(na.omit)) #Leaf dry matter
dist.stemDM = dist(response$stem.DM, method = "euclidean", na.action(na.omit)) #Stem dry matter
dist.aboveDM = dist(response$above.DM, method = "euclidean", na.action(na.omit)) #Above dry matter
dist.rootDM = dist(response$root.DM, method = "euclidean", na.action(na.omit)) #Root dry matter
dist.plantDM = dist(response$plant.DM, method = "euclidean", na.action(na.omit)) #Total plant dry matter
dist.ratioAB = dist(response$Ratio.above.root, method = "euclidean", na.action(na.omit)) #Ratio above-belowground
dist.rootdiam = dist(response$root.diameter, method = "euclidean", na.action(na.omit)) #Root diameter
#Nodulation
dist.nodulenumber = dist(response$nodule.number, method = "euclidean", na.action(na.omit)) #Nodule number
dist.nodules.rootDM = dist(response$nodules.rootDM, method = "euclidean", na.action(na.omit)) #Nodules per root DM
dist.noduleweight = dist(response$total.nodule.weight.mg, method = "euclidean", na.action(na.omit)) #Total nodule weight
#Stoichiometry
#Leaf
dist.leafPconc = dist(response$Leaf.Pconc..mg.g., method = "euclidean", na.action(na.omit)) #Foliar P concentration
dist.leafPcont = dist(response$Leaf.Pcont, method = "euclidean", na.action(na.omit)) #Foliar P content
dist.leafNconc = dist(response$Leaf.Nconc.mg.g., method = "euclidean", na.action(na.omit)) #Foliar N concentration
dist.leafNcont = dist(response$Leaf.N.content, method = "euclidean", na.action(na.omit)) #Foliar N content
dist.leafNP = dist(response$N.P.ratio, method = "euclidean", na.action(na.omit)) #Foliar NP ratio
dist.leaf15N = dist(response$Leaf.d15N, method = "euclidean", na.action(na.omit)) #Foliar 15N signature
dist.leaf13C = dist(response$Leaf.d13C, method = "euclidean", na.action(na.omit)) #Foliar 13C signature
dist.leafMn = dist(response$Leaf.Mn..ug.g., method = "euclidean", na.action(na.omit)) #Foliar Mn concentration
dist.leafFe = dist(response$Leaf.Fe..ug.g., method = "euclidean", na.action(na.omit)) #Foliar Fe concentration
dist.leafCu = dist(response$Leaf.Cu..ug.g., method = "euclidean", na.action(na.omit)) #Foliar Cu concentration
dist.leafZn = dist(response$Leaf.Zn..ug.g., method = "euclidean", na.action(na.omit)) #Foliar Zn concentration
#Root
dist.rootN = dist(response$Root.N.mg.g., method = "euclidean", na.action(na.omit)) #Root N concentration
dist.root15N = dist(response$Root.15N, method = "euclidean", na.action(na.omit)) #Root 15N signature
dist.root13C = dist(response$Root.13C, method = "euclidean", na.action(na.omit)) #Root 13C signature
dist.rootP = dist(response$rootP.mg.g., method = "euclidean", na.action(na.omit)) #Root P concentration
dist.rootCa = dist(response$rootCa.mg.g., method = "euclidean", na.action(na.omit)) #Root Ca concentration
dist.rootFe = dist(response$rootFe.mg.g., method = "euclidean", na.action(na.omit)) #Root Fe concentration
dist.rootK = dist(response$rootK.mg.g., method = "euclidean", na.action(na.omit)) #Root K concentration
dist.rootMg = dist(response$rootMg.mg.g., method = "euclidean", na.action(na.omit)) #Root Mg concentration
dist.rootMn = dist(response$rootMn.mg.g., method = "euclidean", na.action(na.omit)) #Root Mn concentration
#Nodule
dist.nodule15N = dist(response$Nodule.15N, method = "euclidean", na.action(na.omit)) #Nodule 15N signature
dist.nodule13C = dist(response$Nodule.13C, method = "euclidean", na.action(na.omit)) #Nodule 13C signature
dist.noduleN = dist(response$Nodule.totalN, method = "euclidean", na.action(na.omit)) #Nodule N content
dist.noduleCN = dist(response$Nodule.CNratio, method = "euclidean", na.action(na.omit)) #Nodule CN ratio


###MANTEL TEST### !!!!!Best correlations with Geographical distance and soil N,P!!!!!
aa = t(as.data.frame(otu_table(Dt)))
bdist <- vegdist(aa, "bray")
wunif.dist <- distance(Dt, "unifrac", weighted = T)
#Geography
#library(ape)
test.wunif.lat = mantel(dist.latitude, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.lat = mantel(dist.latitude, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.lon = mantel(dist.longitude, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.lon = mantel(dist.longitude, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.km = mantel(dist.km, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.km = mantel(dist.km, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.alt = mantel(dist.altitude, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.alt = mantel(dist.altitude, bdist, method = "pearson", permutations = 9999, na.rm = T)
##PLOT MANTEL CORREL
#correlog.km = mantel.correlog(D.geo = dist.km, D.eco = bdist, r.type = "spearman", nperm = 999)
#plot(correlog.km, alpha = 0.15)
plot(x=as.matrix(dist.km), y=as.matrix(bdist), geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 

#Soil factors
test.wunif.cec = mantel(dist.CEC, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.cec = mantel(dist.CEC, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.whc = mantel(dist.WHC, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.whc = mantel(dist.WHC, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.sand = mantel(dist.sand, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.sand = mantel(dist.sand, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.ph = mantel(dist.pH, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.ph = mantel(dist.pH, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.N = mantel(dist.soilN, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.N = mantel(dist.soilN, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.P = mantel(dist.soilP, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.P = mantel(dist.soilP, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.Ca = mantel(dist.soilCa, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.Ca = mantel(dist.soilCa, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.K = mantel(dist.soilK, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.K = mantel(dist.soilK, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.Mg = mantel(dist.soilMg, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.Mg = mantel(dist.soilMg, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.Mn = mantel(dist.soilMn, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.Mn = mantel(dist.soilMn, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.orgC = mantel(dist.orgC, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.orgC = mantel(dist.orgC, bdist, method = "pearson", permutations = 9999, na.rm = T)

#Plant response
test.wunif.tiller = mantel(dist.tillerlength, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.tiller = mantel(dist.tillerlength, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.leafDM = mantel(dist.leafDM, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.leafDM = mantel(dist.leafDM, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.stemDM = mantel(dist.stemDM, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.stemDM = mantel(dist.stemDM, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.aboveDM = mantel(dist.aboveDM, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.aboveDM = mantel(dist.aboveDM, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.rootDM = mantel(dist.rootDM, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.rootDM = mantel(dist.rootDM, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.plantDM = mantel(dist.plantDM, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.plantDM = mantel(dist.plantDM, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.ratioAB = mantel(dist.ratioAB, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.ratioAB = mantel(dist.ratioAB, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.rootdiam = mantel(dist.rootdiam, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.rootdiam = mantel(dist.rootdiam, bdist, method = "pearson", permutations = 9999, na.rm = T)
#Nodulation
test.wunif.nodulenumber = mantel(dist.nodulenumber, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.nodulenumber = mantel(dist.nodulenumber, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.nodules.rootDM = mantel(dist.nodules.rootDM, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.nodules.rootDM = mantel(dist.nodules.rootDM, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.noduleweight = mantel(dist.noduleweight, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.noduleweight = mantel(dist.noduleweight, bdist, method = "pearson", permutations = 9999, na.rm = T)
#Leaf stoichiometry
test.wunif.leafPconc = mantel(dist.leafPconc, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.leafPconc = mantel(dist.leafPconc, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.leafPcont = mantel(dist.leafPcont, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.leafPcont = mantel(dist.leafPcont, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.leafNconc = mantel(dist.leafNconc, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.leafNconc = mantel(dist.leafNconc, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.leafNcont = mantel(dist.leafNcont, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.leafNcont = mantel(dist.leafNcont, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.leafNP = mantel(dist.leafNP, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.leafNP = mantel(dist.leafNP, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.leaf15N = mantel(dist.leaf15N, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.leaf15N = mantel(dist.leaf15N, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.leaf13C = mantel(dist.leaf13C, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.leaf13C = mantel(dist.leaf13C, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.leafMn = mantel(dist.leafMn, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.leafMn = mantel(dist.leafMn, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.leafFe = mantel(dist.leafFe, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.leafFe = mantel(dist.leafFe, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.leafCu = mantel(dist.leafCu, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.leafCu = mantel(dist.leafCu, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.leafZn = mantel(dist.leafZn, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.leafZn = mantel(dist.leafZn, bdist, method = "pearson", permutations = 9999, na.rm = T)
#Root stoichiometry
test.wunif.rootN = mantel(dist.rootN, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.rootN = mantel(dist.rootN, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.root15N = mantel(dist.root15N, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.root15N = mantel(dist.root15N, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.root13C = mantel(dist.root13C, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.root13C = mantel(dist.root13C, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.rootP = mantel(dist.rootP, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.rootP = mantel(dist.rootP, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.rootCa = mantel(dist.rootCa, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.rootCa = mantel(dist.rootCa, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.rootFe = mantel(dist.rootFe, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.rootFe = mantel(dist.rootFe, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.rootK = mantel(dist.rootK, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.rootK = mantel(dist.rootK, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.rootMg = mantel(dist.rootMg, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.rootMg = mantel(dist.rootMg, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.rootMn = mantel(dist.rootMn, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.rootMn = mantel(dist.rootMn, bdist, method = "pearson", permutations = 9999, na.rm = T)
#Nodule stoichiometry
test.wunif.nodule15N = mantel(dist.nodule15N, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.nodule15N = mantel(dist.nodule15N, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.nodule13C = mantel(dist.nodule13C, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.nodule13C = mantel(dist.nodule13C, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.noduleN = mantel(dist.noduleN, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.noduleN = mantel(dist.noduleN, bdist, method = "pearson", permutations = 9999, na.rm = T)
test.wunif.noduleCN = mantel(dist.noduleCN, wunif.dist, method = "pearson", permutations = 9999, na.rm = T)
test.bray.noduleCN = mantel(dist.noduleCN, bdist, method = "pearson", permutations = 9999, na.rm = T)



##EVENNESS ANALYSIS##
Dtu = subset_samples(Dt, Fertilization == "unfertilized")
Dtf = subset_samples(Dt, Fertilization == "fertilized")

#Evenness considering all taxa
even <- evenness(Dt, "all")
pielou = even$pielou
simpson = even$simpson
sample_data(Dt) = as.data.frame(cbind(sample_data(Dt), pielou, simpson))

###TESTING TREATMENT EFFECTS ON EVENNESS##
##Testing (ANOVA + TUKEY'S HSD)
even = sample_data(Dt)
#Pielou-Soil
lm = lm(even$pielou ~ even$Soil)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "even$Soil", conf.level=0.95)
#Boxplot
boxplot(even$pielou ~ even$Soil, na.rm = T, ylab = "Evenness (Pielou's index)", xlab = "Soil")   

#Pielou-Farm
lm = lm(even$pielou ~ even$Soil*even$Farm)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "even$Soil:even$Farm", conf.level=0.95)
#Boxplot
boxplot(even$pielou ~ even$Farm.Soil, na.rm = T, ylab = "Evenness (Pielou's index)", xlab = "Farm")   

#Pielou-Fertilization
lm = lm(even$pielou ~ even$Fertilization)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "even$Fertilization", conf.level=0.95)
#Boxplot
boxplot(even$pielou ~ even$Fertilization, na.rm = T, ylab = "Evenness (Pielou's index)", xlab = "Fertilization")   

#Pielou-DM !!!!!!!!
even <- evenness(Dtu, "all")
pielou = even$pielou
simpson = even$simpson
sample_data(Dtu) = as.data.frame(cbind(sample_data(Dtu), pielou, simpson))
even = sample_data(Dtu)

lm = lm(even$plant.DM ~ even$pielou)
lmu = aov(lm)
summary(lmu)

#Pielou-leaf15N
lm = lm(even$Leaf.d15N ~ even$pielou)
lmu = aov(lm)
summary(lmu)

#Pielou-nodule15N
lm = lm(even$Nodule.totalN ~ even$pielou)
lmu = aov(lm)
summary(lmu)

#Pielou-leaf N content!!!!!!!
lm = lm(even$Leaf.N.content ~ even$pielou)
lmu = aov(lm)
summary(lmu)

#Pielou-leaf P content
lm = lm(even$Leaf.Pcont ~ even$pielou)
lmu = aov(lm)
summary(lmu)

###Scatter plots
##Evenness
#Pielou-DM!!!
p=qplot(x=pielou, y=plant.DM, color=Soil, data=as.data.frame(even), geom="point") + geom_point(size=5) + geom_smooth(method = "lm") 
p + facet_wrap(~Soil)
rcorr = cor.test(even$pielou, even$plant.DM, method = "pearson")
#Pielou-tiller length!!!
qplot(x=pielou, y=tiller.total.length, color=Soil, data=as.data.frame(even), geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
rcorr = cor.test(even$pielou, even$tiller.total.length, method = "pearson")
#Pielou-nodule number!!!
qplot(x=pielou, y=nodule.number, color=Soil, data=as.data.frame(even), geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
rcorr = cor.test(even$pielou, even$nodule.number, method = "pearson")
#Pielou-root diameter!!!
qplot(x=pielou, y=root.diameter, color=Soil, data=as.data.frame(even), geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
rcorr = cor.test(even$pielou, even$root.diameter, method = "pearson")
#Pielou-P content!!!!
p = qplot(x=pielou, y=Leaf.Pcont, color=Soil, data=as.data.frame(even), geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
p + facet_wrap(~Soil)
rcorr = cor.test(even$pielou, even$Leaf.Pcont, method = "pearson")

##Richness-Diversity
even$richness = richness(Dtu)[,1]
even$shannon = diversities(otu_table(Dtu), index = "shannon")[,]
#Richness NOT SIGNIFICANT
lm = lm(even$plant.DM ~ even$richness)
lmu = aov(lm)
summary(lmu)
#Shannon index !!!!!!!!
lm = lm(even$plant.DM ~ even$shannon)
lmu = aov(lm)
summary(lmu)
#Boxplot
boxplot(even$shannon ~ even$Soil, na.rm = T, ylab = "Diversity (Shannon's H index)", xlab = "Soil")   
#Shannon-Soil
lm = lm(even$shannon ~ even$Soil)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "even$Soil", conf.level=0.95)
#Shannon-Farm
lm = lm(even$shannon ~ even$Farm)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "even$Farm", conf.level=0.95)
#Boxplot
boxplot(even$shannon ~ even$Farm.Soil, na.rm = T, ylab = "Evenness (Pielou's index)", xlab = "Farm")   
#Shannon-Fertilization
lm = lm(even$shannon ~ even$Fertilization)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "even$Fertilization", conf.level=0.95)
#Boxplot
boxplot(even$shannon ~ even$Fertilization, na.rm = T, ylab = "Evenness (Pielou's index)", xlab = "Fertilization")   




#Evenness for only dominant strains
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
write.table(as.matrix(test), "top_abundant_taxa.txt")
toptax = read.csv("top_abundant_taxa.csv", header = T, sep=",", row.names = 1) 
# toptax = as.data.frame(t(toptax))
tax123 = toptax[,94:96]
sample_data(Dt) = cbind(sample_data(Dt), tax123)
dominant.OTU <- subset(otu_table(Dt), rownames(otu_table(Dt)) %in% c("ZOTU43", "ZOTU6", "ZOTU28", "ZOTU23","ZOTU99", "ZOTU7","ZOTU72","ZOTU18", "ZOTU27","ZOTU22","ZOTU155","ZOTU12","ZOTU42","ZOTU21","ZOTU97","ZOTU86","ZOTU112","ZOTU31","ZOTU26","ZOTU19","ZOTU63","ZOTU50","ZOTU98","ZOTU76","ZOTU40","ZOTU14","ZOTU33","ZOTU166","ZOTU45"))
dom.physeq <- merge_phyloseq(dominant.OTU, tax_table(Dt), sample_data(Dt), phy_tree(Dt))







##EVENNESS CALCULATION
even <- evenness(Dt, "all")
pielou = even$pielou
sample_data(Dt) = as.data.frame(cbind(sample_data(Dt), pielou))

###TESTING TREATMENT EFFECTS ON EVENNESS##
##Testing (ANOVA + TUKEY'S HSD)
#Pielou-Soil
dom.even = sample_data(Dtu)
lm = lm(dom.even$pielou ~ dom.even$Soil)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "dom.even$Soil", conf.level=0.95)
#Boxplot
boxplot(dom.even$pielou ~ dom.even$Soil, na.rm = T, ylab = "Evenness (Pielou's index)", xlab = "Soil")   

#Pielou-Farm
dom.even = sample_data(Dt)
lm = lm(dom.even$pielou ~ dom.even$Farm)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "dom.even$Farm", conf.level=0.95)
#Boxplot
boxplot(dom.even$pielou ~ dom.even$Farm, na.rm = T, ylab = "Evenness (Pielou's index)", xlab = "Farm")   

#Pielou-Fertilization
dom.even = sample_data(Dt)
lm = lm(dom.even$pielou ~ dom.even$Fertilization)
lmu = aov(lm)
anova(lmu)
tukey = TukeyHSD(lmu, "dom.even$Fertilization", conf.level=0.95)
#Boxplot
boxplot(dom.even$pielou ~ dom.even$Fertilization, na.rm = T, ylab = "Evenness (Pielou's index)", xlab = "Fertilization")   

#Pielou-DM
Dtu = subset_samples(Dt, Fertilization == "unfertilized")
dom.even = sample_data(Dtu)
lm = lm(dom.even$plant.DM ~ dom.even$pielou)
lmu = aov(lm)
anova(lmu)
qplot(x=pielou, y=plant.DM, color=Soil, data=dom.even, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
cor.test(x=dom.even$pielou, y=dom.even$plant.DM, method = "pearson")
#Pielou-leaf15N
dom.even = sample_data(Dtu)
lm = lm(dom.even$Leaf.d15N ~ dom.even$pielou)
lmu = aov(lm)
summary(lmu)

#Pielou-nodule15N
dom.even = sample_data(Dt)
lm = lm(dom.even$Nodule.totalN ~ dom.even$pielou)
lmu = aov(lm)
summary(lmu)

#Pielou-leaf N content
dom.even = sample_data(Dtu)
lm = lm(dom.even$Leaf.N.content ~ dom.even$pielou)
lmu = aov(lm)
summary(lmu)
qplot(x=pielou, y=Leaf.N.content, color=Soil, data=dom.even, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
cor.test(x=dom.even$pielou, y=dom.even$Leaf.N.content, method = "pearson")
#Pielou-leaf P content
dom.even = sample_data(Dtu)
lm = lm(dom.even$Leaf.Pcont ~ dom.even$pielou)
lmu = aov(lm)
summary(lmu)
qplot(x=pielou, y=Leaf.Pcont, color=Soil, data=dom.even, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
cor.test(x=dom.even$pielou, y=dom.even$Leaf.Pcont, method = "pearson")

###Scatter plots
##Evenness
#Pielou-DM!!!
qplot(x=pielou, y=plant.DM, color=Soil, data=dom.even, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#Pielou-tiller length!!!
qplot(x=pielou, y=tiller.total.length, color=Soil, data=dom.even, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#Pielou-nodule number!!!
qplot(x=pielou, y=nodule.number, color=Soil, data=dom.even, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#Pielou-root diameter!!!
qplot(x=pielou, y=root.diameter, color=Soil, data=dom.even, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#Pielou-delta15N!!!No effects
qplot(x=pielou, y=Leaf.d15N, color=Soil, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#Pielou-Nodule 15N!!!
qplot(x=pielou, y=Nodule.15N, color=Soil, data=dom.even, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#Pielou-Nodule total N!!!
qplot(x=pielou, y=Nodule.totalN, color=Soil, data=dom.even, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 




##DIST TO CENTROID ANALYSIS TO RELATE TO SOIL FACTORS WITHOUT VIOLATING STATISTICAL ASSUMPTIONS
bray = vegdist(t(otu_table(Dt)), method="bray", na.omit=T)
aa = as.data.frame(sample_data(Dt))
blo = subset(aa, Farm == "blo")
blo = rownames(blo)
dob = subset(aa, Farm == "dob")
dob = rownames(dob)
lan = subset(aa, Farm == "lan")
lan = rownames(lan)
mat = subset(aa, Farm == "mat")
mat = rownames(mat)
mel = subset(aa, Farm == "mel")
mel = rownames(mel)

library(usedist)
library(cba)
library(funfuns)

centroids.blo.dob = dist_between_centroids(bray, blo, dob)
centroids.blo.lan = dist_between_centroids(bray, blo, lan)
centroids.blo.mat = dist_between_centroids(bray, blo, mat)
centroids.blo.mel = dist_between_centroids(bray, blo, mel)
centroids.dob.lan = dist_between_centroids(bray, dob, lan)
centroids.dob.mat = dist_between_centroids(bray, dob, mat)
centroids.dob.mel = dist_between_centroids(bray, dob, mel)
centroids.lan.mat = dist_between_centroids(bray, lan, mat)
centroids.lan.mel = dist_between_centroids(bray, lan, mel)
centroids.mat.mel = dist_between_centroids(bray, mat, mel)

cent.distances = as.data.frame(rbind(centroids.blo.dob, centroids.blo.lan, centroids.blo.mat, centroids.blo.mel,
                       centroids.dob.lan, centroids.dob.mat, centroids.dob.mel, centroids.lan.mat,
                       centroids.lan.mel, centroids.mat.mel))
cent.distances$Distance = cent.distances$V1 
km.distances = as.data.frame(rbind(9.06, 14.25, 17.48, 14.39, 7.67, 12.60, 22.55, 5.23, 28.89, 31.04))
km.distances$Km = km.distances$V1
aa = cbind(cent.distances, km.distances)

#Centroids + soil parameters => Farm x Soil interaction
bloc = subset(aa, Farm.Soil == "blocultivated")
bloc = rownames(bloc)
dobc = subset(aa, Farm.Soil == "dobcultivated")
dobc = rownames(dobc)
lanc = subset(aa, Farm.Soil == "lancultivated")
lanc = rownames(lanc)
matc = subset(aa, Farm.Soil == "matcultivated")
matc = rownames(matc)
melc = subset(aa, Farm.Soil == "melcultivated")
melc = rownames(melc)
blom = subset(aa, Farm.Soil == "blomixed")
blom = rownames(blom)
dobm = subset(aa, Farm.Soil == "dobmixed")
dobm = rownames(dobm)
lanm = subset(aa, Farm.Soil == "lanmixed")
lanm = rownames(lanm)
matm = subset(aa, Farm.Soil == "matmixed")
matm = rownames(matm)
melm = subset(aa, Farm.Soil == "melmixed")
melm = rownames(melm)
blow = subset(aa, Farm.Soil == "blowild")
blow = rownames(blow)
dobw = subset(aa, Farm.Soil == "dobwild")
dobw = rownames(dobw)
lanw = subset(aa, Farm.Soil == "lanwild")
lanw = rownames(lanw)
matw = subset(aa, Farm.Soil == "matwild")
matw = rownames(matw)
melw = subset(aa, Farm.Soil == "melwild")
melw = rownames(melw)

library(usedist)
library(cba)
library(funfuns)

centroids.bloc.dobc = dist_between_centroids(bray, bloc, dobc)
centroids.bloc.lanc = dist_between_centroids(bray, bloc, lanc)
centroids.bloc.matc = dist_between_centroids(bray, bloc, matc)
centroids.bloc.melc = dist_between_centroids(bray, bloc, melc)
centroids.dobc.lanc = dist_between_centroids(bray, dobc, lanc)
centroids.dobc.matc = dist_between_centroids(bray, dobc, matc)
centroids.dobc.melc = dist_between_centroids(bray, dobc, melc)
centroids.lanc.matc = dist_between_centroids(bray, lanc, matc)
centroids.lanc.melc = dist_between_centroids(bray, lanc, melc)
centroids.matc.melc = dist_between_centroids(bray, matc, melc)
centroids.blow.dobw = dist_between_centroids(bray, blow, dobw)
centroids.blow.lanw = dist_between_centroids(bray, blow, lanw)
centroids.blow.matw = dist_between_centroids(bray, blow, matw)
centroids.blow.melw = dist_between_centroids(bray, blow, melw)
centroids.dobw.lanw = dist_between_centroids(bray, dobw, lanw)
centroids.dobw.matw = dist_between_centroids(bray, dobw, matw)
centroids.dobw.melw = dist_between_centroids(bray, dobw, melw)
centroids.lanw.matw = dist_between_centroids(bray, lanw, matw)
centroids.lanw.melw = dist_between_centroids(bray, lanw, melw)
centroids.matw.melw = dist_between_centroids(bray, matw, melw)
centroids.blom.dobm = dist_between_centroids(bray, blom, dobm)
centroids.blom.lanm = dist_between_centroids(bray, blom, lanm)
centroids.blom.matm = dist_between_centroids(bray, blom, matm)
centroids.blom.melm = dist_between_centroids(bray, blom, melm)
centroids.dobm.lanm = dist_between_centroids(bray, dobm, lanm)
centroids.dobm.matm = dist_between_centroids(bray, dobm, matm)
centroids.dobm.melm = dist_between_centroids(bray, dobm, melm)
centroids.lanm.matm = dist_between_centroids(bray, lanm, matm)
centroids.lanm.melm = dist_between_centroids(bray, lanm, melm)
centroids.matm.melm = dist_between_centroids(bray, matm, melm)

centroids.bloc.blow = dist_between_centroids(bray, bloc, blow)
centroids.bloc.blom = dist_between_centroids(bray, bloc, blom)
centroids.blow.blom = dist_between_centroids(bray, blow, blom)
centroids.dobc.dobw = dist_between_centroids(bray, dobc, dobw)
centroids.dobc.dobm = dist_between_centroids(bray, dobc, dobm)
centroids.dobm.dobw = dist_between_centroids(bray, dobm, dobw)
centroids.lanc.lanw = dist_between_centroids(bray, lanc, lanw)
centroids.lanc.lanm = dist_between_centroids(bray, lanc, lanm)
centroids.lanm.lanw = dist_between_centroids(bray, lanm, lanw)
centroids.matc.matw = dist_between_centroids(bray, matc, matw)
centroids.matc.matm = dist_between_centroids(bray, matc, matm)
centroids.matm.matw = dist_between_centroids(bray, matm, matw)
centroids.melc.melw = dist_between_centroids(bray, melc, melw)
centroids.melc.melm = dist_between_centroids(bray, melc, melm)
centroids.melm.melw = dist_between_centroids(bray, melm, melw)

centroids.bloc.dobm = dist_between_centroids(bray, bloc, dobm)
centroids.bloc.lanm = dist_between_centroids(bray, bloc, lanm)
centroids.bloc.matm = dist_between_centroids(bray, bloc, matm)
centroids.bloc.melm = dist_between_centroids(bray, bloc, melm)
centroids.dobc.lanm = dist_between_centroids(bray, dobc, lanm)
centroids.dobc.matm = dist_between_centroids(bray, dobc, matm)
centroids.dobc.melm = dist_between_centroids(bray, dobc, melm)
centroids.lanc.matm = dist_between_centroids(bray, lanc, matm)
centroids.lanc.melm = dist_between_centroids(bray, lanc, melm)
centroids.matc.melm = dist_between_centroids(bray, matc, melm)

centroids.blom.dobc = dist_between_centroids(bray, blom, dobc)
centroids.blom.lanc = dist_between_centroids(bray, blom, lanc)
centroids.blom.matc = dist_between_centroids(bray, blom, matc)
centroids.blom.melc = dist_between_centroids(bray, blom, melc)
centroids.dobm.lanc = dist_between_centroids(bray, dobm, lanc)
centroids.dobm.matc = dist_between_centroids(bray, dobm, matc)
centroids.dobm.melc = dist_between_centroids(bray, dobm, melc)
centroids.lanm.matc = dist_between_centroids(bray, lanm, matc)
centroids.lanm.melc = dist_between_centroids(bray, lanm, melc)
centroids.matm.melc = dist_between_centroids(bray, matm, melc)

centroids.bloc.dobw = dist_between_centroids(bray, bloc, dobw)
centroids.bloc.lanw = dist_between_centroids(bray, bloc, lanw)
centroids.bloc.matw = dist_between_centroids(bray, bloc, matw)
centroids.bloc.melw = dist_between_centroids(bray, bloc, melw)
centroids.dobc.lanw = dist_between_centroids(bray, dobc, lanw)
centroids.dobc.matw = dist_between_centroids(bray, dobc, matw)
centroids.dobc.melw = dist_between_centroids(bray, dobc, melw)
centroids.lanc.matw = dist_between_centroids(bray, lanc, matw)
centroids.lanc.melw = dist_between_centroids(bray, lanc, melw)
centroids.matc.melw = dist_between_centroids(bray, matc, melw)

centroids.blow.dobc = dist_between_centroids(bray, blow, dobc)
centroids.blow.lanc = dist_between_centroids(bray, blow, lanc)
centroids.blow.matc = dist_between_centroids(bray, blow, matc)
centroids.blow.melc = dist_between_centroids(bray, blow, melc)
centroids.dobw.lanc = dist_between_centroids(bray, dobw, lanc)
centroids.dobw.matc = dist_between_centroids(bray, dobw, matc)
centroids.dobw.melc = dist_between_centroids(bray, dobw, melc)
centroids.lanw.matc = dist_between_centroids(bray, lanw, matc)
centroids.lanw.melc = dist_between_centroids(bray, lanw, melc)
centroids.matw.melc = dist_between_centroids(bray, matw, melc)

centroids.blow.dobm = dist_between_centroids(bray, blow, dobm)
centroids.blow.lanm = dist_between_centroids(bray, blow, lanm)
centroids.blow.matm = dist_between_centroids(bray, blow, matm)
centroids.blow.melm = dist_between_centroids(bray, blow, melm)
centroids.dobw.lanm = dist_between_centroids(bray, dobw, lanm)
centroids.dobw.matm = dist_between_centroids(bray, dobw, matm)
centroids.dobw.melm = dist_between_centroids(bray, dobw, melm)
centroids.lanw.matm = dist_between_centroids(bray, lanw, matm)
centroids.lanw.melm = dist_between_centroids(bray, lanw, melm)
centroids.matw.melm = dist_between_centroids(bray, matw, melm)

centroids.blom.dobw = dist_between_centroids(bray, blom, dobw)
centroids.blom.lanw = dist_between_centroids(bray, blom, lanw)
centroids.blom.matw = dist_between_centroids(bray, blom, matw)
centroids.blom.melw = dist_between_centroids(bray, blom, melw)
centroids.dobm.lanw = dist_between_centroids(bray, dobm, lanw)
centroids.dobm.matw = dist_between_centroids(bray, dobm, matw)
centroids.dobm.melw = dist_between_centroids(bray, dobm, melw)
centroids.lanm.matw = dist_between_centroids(bray, lanm, matw)
centroids.lanm.melw = dist_between_centroids(bray, lanm, melw)
centroids.matm.melw = dist_between_centroids(bray, matm, melw)


cent.distances = as.data.frame(rbind(centroids.bloc.dobc, centroids.bloc.lanc, centroids.bloc.matc, centroids.bloc.melc,centroids.dobc.lanc,centroids.dobc.matc,centroids.dobc.melc,
                                     centroids.lanc.matc, centroids.lanc.melc, centroids.matc.melc, centroids.blow.dobw, centroids.blow.lanw, centroids.blow.matw, centroids.blow.melw,centroids.dobw.lanw,centroids.dobw.matw,centroids.dobw.melw,
                                     centroids.lanw.matw, centroids.lanw.melw, centroids.matw.melw, centroids.blom.dobm, centroids.blom.lanm, centroids.blom.matm, centroids.blom.melm,centroids.dobm.lanm,centroids.dobm.matm,centroids.dobm.melm,
                                     centroids.lanm.matm, centroids.lanm.melm, centroids.matm.melm, centroids.bloc.blow, centroids.bloc.blom, centroids.blow.blom, centroids.dobc.dobw, centroids.dobc.dobm, centroids.dobm.dobw,
                                     centroids.lanc.lanw, centroids.lanc.lanm, centroids.lanm.lanw,centroids.matc.matw, centroids.matc.matm, centroids.matm.matw, centroids.melc.melw, centroids.melc.melm, centroids.melm.melw,
                                     centroids.bloc.dobm, centroids.bloc.lanm, centroids.bloc.matm, centroids.bloc.melm, centroids.dobc.lanm, centroids.dobc.matm, centroids.dobc.melm,
                                     centroids.lanc.matm, centroids.lanc.melm, centroids.matc.melm, centroids.blom.dobc, centroids.blom.lanc,centroids.blom.matc,centroids.blom.melc,
                                     centroids.dobm.lanc, centroids.dobm.matc, centroids.dobm.melc, centroids.lanm.matc,centroids.lanm.melc,
                                     centroids.matm.melc, centroids.bloc.dobw, centroids.bloc.lanw, centroids.bloc.matw, centroids.bloc.melw, centroids.dobc.lanw,centroids.dobc.matw,
                                     centroids.dobc.melw, centroids.lanc.matw,centroids.lanc.melw,centroids.matc.melw,centroids.blow.dobc, centroids.blow.lanc, centroids.blow.matc,
                                     centroids.blow.melc, centroids.dobw.lanc, centroids.dobw.matc,centroids.dobw.melc,centroids.lanw.matc,centroids.lanw.melc,centroids.matw.melc,
                                     centroids.blow.dobm, centroids.blow.lanm, centroids.blow.matm, centroids.blow.melm, centroids.dobw.lanm, centroids.dobw.matm,centroids.dobw.melm,centroids.lanw.matm,centroids.lanw.melm,centroids.matw.melm,
                                     centroids.blom.dobw, centroids.blom.lanw, centroids.blom.matw, centroids.blom.melw, centroids.dobm.lanw,centroids.dobm.matw,centroids.dobm.melw,
 
                                                                         centroids.lanm.matw, centroids.lanm.melw, centroids.matm.melw))
p = qplot(x=aa$Km, y=aa$Distance, geom="point") + geom_point(size=5) 
+ geom_smooth(method="lm", se = T) 
p + theme_bw() + theme_classic()
rcorr = cor.test(aa$Distance, aa$Km, method = "pearson")

p = qplot(x=aa$Km, y=aa$Distance, geom="point") + geom_point(size=5) 
p + theme_bw() + theme_classic()
rcorr = cor.test(aa$Distance, aa$Km, method = "pearson")

##CENTROID AND SOIL PARAMETERS DISTANCES IN SINGLE DATAFRAME 
soil.params = read.csv("SoilParamsFarmSoil.csv", header = T)

#CEC
a = as.matrix(dist(soil.params$CEC.cmol.kg.))
a = data.frame(as.table(a))[lower.tri(a, diag = TRUE), ]

levels(a$Var1) <- c(levels(a$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
 "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
 "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 
levels(a$Var2) <- c(levels(a$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 

a$Var1[a$Var1==1]  <- "blocultivated" 
a$Var2[a$Var2==1]  <- "blocultivated" 
a$Var1[a$Var1==2]  <- "blomixed" 
a$Var2[a$Var2==2]  <- "blomixed" 
a$Var1[a$Var1==3]  <- "blowild" 
a$Var2[a$Var2==3]  <- "blowild" 
a$Var1[a$Var1==4]  <- "dobcultivated" 
a$Var2[a$Var2==4]  <- "dobcultivated" 
a$Var1[a$Var1==5]  <- "dobmixed" 
a$Var2[a$Var2==5]  <- "dobmixed" 
a$Var1[a$Var1==6]  <- "dobwild" 
a$Var2[a$Var2==6]  <- "dobwild" 
a$Var1[a$Var1==7]  <- "lancultivated" 
a$Var2[a$Var2==7]  <- "lancultivated" 
a$Var1[a$Var1==8]  <- "lanmixed" 
a$Var2[a$Var2==8]  <- "lanmixed" 
a$Var1[a$Var1==9]  <- "lanwild" 
a$Var2[a$Var2==9]  <- "lanwild" 
a$Var1[a$Var1==10]  <- "matcultivated" 
a$Var2[a$Var2==10]  <- "matcultivated" 
a$Var1[a$Var1==11]  <- "matmixed" 
a$Var2[a$Var2==11]  <- "matmixed" 
a$Var1[a$Var1==12]  <- "matwild" 
a$Var2[a$Var2==12]  <- "matwild" 
a$Var1[a$Var1==13]  <- "melcultivated" 
a$Var2[a$Var2==13]  <- "melcultivated" 
a$Var1[a$Var1==14]  <- "melmixed" 
a$Var2[a$Var2==14]  <- "melmixed" 
a$Var1[a$Var1==15]  <- "melwild" 
a$Var2[a$Var2==15]  <- "melwild" 

#WHC
b = as.matrix(dist(soil.params$WHC...10kPa.))
b = data.frame(as.table(b))[lower.tri(b, diag = TRUE), ]

levels(b$Var1) <- c(levels(b$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 
levels(b$Var2) <- c(levels(b$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 

b$Var1[b$Var1==1]  <- "blocultivated" 
b$Var2[b$Var2==1]  <- "blocultivated" 
b$Var1[b$Var1==2]  <- "blomixed" 
b$Var2[b$Var2==2]  <- "blomixed" 
b$Var1[b$Var1==3]  <- "blowild" 
b$Var2[b$Var2==3]  <- "blowild" 
b$Var1[b$Var1==4]  <- "dobcultivated" 
b$Var2[b$Var2==4]  <- "dobcultivated" 
b$Var1[b$Var1==5]  <- "dobmixed" 
b$Var2[b$Var2==5]  <- "dobmixed" 
b$Var1[b$Var1==6]  <- "dobwild" 
b$Var2[b$Var2==6]  <- "dobwild" 
b$Var1[b$Var1==7]  <- "lancultivated" 
b$Var2[b$Var2==7]  <- "lancultivated" 
b$Var1[b$Var1==8]  <- "lanmixed" 
b$Var2[b$Var2==8]  <- "lanmixed" 
b$Var1[b$Var1==9]  <- "lanwild" 
b$Var2[b$Var2==9]  <- "lanwild" 
b$Var1[b$Var1==10]  <- "matcultivated" 
b$Var2[b$Var2==10]  <- "matcultivated" 
b$Var1[b$Var1==11]  <- "matmixed" 
b$Var2[b$Var2==11]  <- "matmixed" 
b$Var1[b$Var1==12]  <- "matwild" 
b$Var2[b$Var2==12]  <- "matwild" 
b$Var1[b$Var1==13]  <- "melcultivated" 
b$Var2[b$Var2==13]  <- "melcultivated" 
b$Var1[b$Var1==14]  <- "melmixed" 
b$Var2[b$Var2==14]  <- "melmixed" 
b$Var1[b$Var1==15]  <- "melwild" 
b$Var2[b$Var2==15]  <- "melwild" 


#Sandiness
c = as.matrix(dist(soil.params$Coarse.sand...))
c = data.frame(as.table(c))[lower.tri(c, diag = TRUE), ]

levels(c$Var1) <- c(levels(c$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 
levels(c$Var2) <- c(levels(c$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 

c$Var1[c$Var1==1]  <- "blocultivated" 
c$Var2[c$Var2==1]  <- "blocultivated" 
c$Var1[c$Var1==2]  <- "blomixed" 
c$Var2[c$Var2==2]  <- "blomixed" 
c$Var1[c$Var1==3]  <- "blowild" 
c$Var2[c$Var2==3]  <- "blowild" 
c$Var1[c$Var1==4]  <- "dobcultivated" 
c$Var2[c$Var2==4]  <- "dobcultivated" 
c$Var1[c$Var1==5]  <- "dobmixed" 
c$Var2[c$Var2==5]  <- "dobmixed" 
c$Var1[c$Var1==6]  <- "dobwild" 
c$Var2[c$Var2==6]  <- "dobwild" 
c$Var1[c$Var1==7]  <- "lancultivated" 
c$Var2[c$Var2==7]  <- "lancultivated" 
c$Var1[c$Var1==8]  <- "lanmixed" 
c$Var2[c$Var2==8]  <- "lanmixed" 
c$Var1[c$Var1==9]  <- "lanwild" 
c$Var2[c$Var2==9]  <- "lanwild" 
c$Var1[c$Var1==10]  <- "matcultivated" 
c$Var2[c$Var2==10]  <- "matcultivated" 
c$Var1[c$Var1==11]  <- "matmixed" 
c$Var2[c$Var2==11]  <- "matmixed" 
c$Var1[c$Var1==12]  <- "matwild" 
c$Var2[c$Var2==12]  <- "matwild" 
c$Var1[c$Var1==13]  <- "melcultivated" 
c$Var2[c$Var2==13]  <- "melcultivated" 
c$Var1[c$Var1==14]  <- "melmixed" 
c$Var2[c$Var2==14]  <- "melmixed" 
c$Var1[c$Var1==15]  <- "melwild" 
c$Var2[c$Var2==15]  <- "melwild" 


#Organic carbon
d = as.matrix(dist(soil.params$OrgC.soil.g))
d = data.frame(as.table(d))[lower.tri(d, diag = TRUE), ]

levels(d$Var1) <- c(levels(d$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 
levels(d$Var2) <- c(levels(d$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 

d$Var1[d$Var1==1]  <- "blocultivated" 
d$Var2[d$Var2==1]  <- "blocultivated" 
d$Var1[d$Var1==2]  <- "blomixed" 
d$Var2[d$Var2==2]  <- "blomixed" 
d$Var1[d$Var1==3]  <- "blowild" 
d$Var2[d$Var2==3]  <- "blowild" 
d$Var1[d$Var1==4]  <- "dobcultivated" 
d$Var2[d$Var2==4]  <- "dobcultivated" 
d$Var1[d$Var1==5]  <- "dobmixed" 
d$Var2[d$Var2==5]  <- "dobmixed" 
d$Var1[d$Var1==6]  <- "dobwild" 
d$Var2[d$Var2==6]  <- "dobwild" 
d$Var1[d$Var1==7]  <- "lancultivated" 
d$Var2[d$Var2==7]  <- "lancultivated" 
d$Var1[d$Var1==8]  <- "lanmixed" 
d$Var2[d$Var2==8]  <- "lanmixed" 
d$Var1[d$Var1==9]  <- "lanwild" 
d$Var2[d$Var2==9]  <- "lanwild" 
d$Var1[d$Var1==10]  <- "matcultivated" 
d$Var2[d$Var2==10]  <- "matcultivated" 
d$Var1[d$Var1==11]  <- "matmixed" 
d$Var2[d$Var2==11]  <- "matmixed" 
d$Var1[d$Var1==12]  <- "matwild" 
d$Var2[d$Var2==12]  <- "matwild" 
d$Var1[d$Var1==13]  <- "melcultivated" 
d$Var2[d$Var2==13]  <- "melcultivated" 
d$Var1[d$Var1==14]  <- "melmixed" 
d$Var2[d$Var2==14]  <- "melmixed" 
d$Var1[d$Var1==15]  <- "melwild" 
d$Var2[d$Var2==15]  <- "melwild" 


#pH
E = as.matrix(dist(soil.params$pH))
E = data.frame(as.table(E))[lower.tri(E, diag = TRUE), ]

levels(E$Var1) <- c(levels(E$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 
levels(E$Var2) <- c(levels(E$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 

E$Var1[E$Var1==1]  <- "blocultivated" 
E$Var2[E$Var2==1]  <- "blocultivated" 
E$Var1[E$Var1==2]  <- "blomixed" 
E$Var2[E$Var2==2]  <- "blomixed" 
E$Var1[E$Var1==3]  <- "blowild" 
E$Var2[E$Var2==3]  <- "blowild" 
E$Var1[E$Var1==4]  <- "dobcultivated" 
E$Var2[E$Var2==4]  <- "dobcultivated" 
E$Var1[E$Var1==5]  <- "dobmixed" 
E$Var2[E$Var2==5]  <- "dobmixed" 
E$Var1[E$Var1==6]  <- "dobwild" 
E$Var2[E$Var2==6]  <- "dobwild" 
E$Var1[E$Var1==7]  <- "lancultivated" 
E$Var2[E$Var2==7]  <- "lancultivated" 
E$Var1[E$Var1==8]  <- "lanmixed" 
E$Var2[E$Var2==8]  <- "lanmixed" 
E$Var1[E$Var1==9]  <- "lanwild" 
E$Var2[E$Var2==9]  <- "lanwild" 
E$Var1[E$Var1==10]  <- "matcultivated" 
E$Var2[E$Var2==10]  <- "matcultivated" 
E$Var1[E$Var1==11]  <- "matmixed" 
E$Var2[E$Var2==11]  <- "matmixed" 
E$Var1[E$Var1==12]  <- "matwild" 
E$Var2[E$Var2==12]  <- "matwild" 
E$Var1[E$Var1==13]  <- "melcultivated" 
E$Var2[E$Var2==13]  <- "melcultivated" 
E$Var1[E$Var1==14]  <- "melmixed" 
E$Var2[E$Var2==14]  <- "melmixed" 
E$Var1[E$Var1==15]  <- "melwild" 
E$Var2[E$Var2==15]  <- "melwild" 

#N content
f = as.matrix(dist(soil.params$Pot.soil.N.conc))
f = data.frame(as.table(f))[lower.tri(f, diag = TRUE), ]

levels(f$Var1) <- c(levels(f$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 
levels(f$Var2) <- c(levels(f$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 

f$Var1[f$Var1==1]  <- "blocultivated" 
f$Var2[f$Var2==1]  <- "blocultivated" 
f$Var1[f$Var1==2]  <- "blomixed" 
f$Var2[f$Var2==2]  <- "blomixed" 
f$Var1[f$Var1==3]  <- "blowild" 
f$Var2[f$Var2==3]  <- "blowild" 
f$Var1[f$Var1==4]  <- "dobcultivated" 
f$Var2[f$Var2==4]  <- "dobcultivated" 
f$Var1[f$Var1==5]  <- "dobmixed" 
f$Var2[f$Var2==5]  <- "dobmixed" 
f$Var1[f$Var1==6]  <- "dobwild" 
f$Var2[f$Var2==6]  <- "dobwild" 
f$Var1[f$Var1==7]  <- "lancultivated" 
f$Var2[f$Var2==7]  <- "lancultivated" 
f$Var1[f$Var1==8]  <- "lanmixed" 
f$Var2[f$Var2==8]  <- "lanmixed" 
f$Var1[f$Var1==9]  <- "lanwild" 
f$Var2[f$Var2==9]  <- "lanwild" 
f$Var1[f$Var1==10]  <- "matcultivated" 
f$Var2[f$Var2==10]  <- "matcultivated" 
f$Var1[f$Var1==11]  <- "matmixed" 
f$Var2[f$Var2==11]  <- "matmixed" 
f$Var1[f$Var1==12]  <- "matwild" 
f$Var2[f$Var2==12]  <- "matwild" 
f$Var1[f$Var1==13]  <- "melcultivated" 
f$Var2[f$Var2==13]  <- "melcultivated" 
f$Var1[f$Var1==14]  <- "melmixed" 
f$Var2[f$Var2==14]  <- "melmixed" 
f$Var1[f$Var1==15]  <- "melwild" 
f$Var2[f$Var2==15]  <- "melwild" 


#P content
g = as.matrix(dist(soil.params$Pot.soil.P.conc))
g = data.frame(as.table(g))[lower.tri(g, diag = TRUE), ]

levels(g$Var1) <- c(levels(g$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 
levels(g$Var2) <- c(levels(g$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 

g$Var1[g$Var1==1]  <- "blocultivated" 
g$Var2[g$Var2==1]  <- "blocultivated" 
g$Var1[g$Var1==2]  <- "blomixed" 
g$Var2[g$Var2==2]  <- "blomixed" 
g$Var1[g$Var1==3]  <- "blowild" 
g$Var2[g$Var2==3]  <- "blowild" 
g$Var1[g$Var1==4]  <- "dobcultivated" 
g$Var2[g$Var2==4]  <- "dobcultivated" 
g$Var1[g$Var1==5]  <- "dobmixed" 
g$Var2[g$Var2==5]  <- "dobmixed" 
g$Var1[g$Var1==6]  <- "dobwild" 
g$Var2[g$Var2==6]  <- "dobwild" 
g$Var1[g$Var1==7]  <- "lancultivated" 
g$Var2[g$Var2==7]  <- "lancultivated" 
g$Var1[g$Var1==8]  <- "lanmixed" 
g$Var2[g$Var2==8]  <- "lanmixed" 
g$Var1[g$Var1==9]  <- "lanwild" 
g$Var2[g$Var2==9]  <- "lanwild" 
g$Var1[g$Var1==10]  <- "matcultivated" 
g$Var2[g$Var2==10]  <- "matcultivated" 
g$Var1[g$Var1==11]  <- "matmixed" 
g$Var2[g$Var2==11]  <- "matmixed" 
g$Var1[g$Var1==12]  <- "matwild" 
g$Var2[g$Var2==12]  <- "matwild" 
g$Var1[g$Var1==13]  <- "melcultivated" 
g$Var2[g$Var2==13]  <- "melcultivated" 
g$Var1[g$Var1==14]  <- "melmixed" 
g$Var2[g$Var2==14]  <- "melmixed" 
g$Var1[g$Var1==15]  <- "melwild" 
g$Var2[g$Var2==15]  <- "melwild" 


#NP ratio
h = as.matrix(dist(soil.params$Pot.NPratio))
h = data.frame(as.table(h))[lower.tri(h, diag = TRUE), ]

levels(h$Var1) <- c(levels(h$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 
levels(h$Var2) <- c(levels(h$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 

h$Var1[h$Var1==1]  <- "blocultivated" 
h$Var2[h$Var2==1]  <- "blocultivated" 
h$Var1[h$Var1==2]  <- "blomixed" 
h$Var2[h$Var2==2]  <- "blomixed" 
h$Var1[h$Var1==3]  <- "blowild" 
h$Var2[h$Var2==3]  <- "blowild" 
h$Var1[h$Var1==4]  <- "dobcultivated" 
h$Var2[h$Var2==4]  <- "dobcultivated" 
h$Var1[h$Var1==5]  <- "dobmixed" 
h$Var2[h$Var2==5]  <- "dobmixed" 
h$Var1[h$Var1==6]  <- "dobwild" 
h$Var2[h$Var2==6]  <- "dobwild" 
h$Var1[h$Var1==7]  <- "lancultivated" 
h$Var2[h$Var2==7]  <- "lancultivated" 
h$Var1[h$Var1==8]  <- "lanmixed" 
h$Var2[h$Var2==8]  <- "lanmixed" 
h$Var1[h$Var1==9]  <- "lanwild" 
h$Var2[h$Var2==9]  <- "lanwild" 
h$Var1[h$Var1==10]  <- "matcultivated" 
h$Var2[h$Var2==10]  <- "matcultivated" 
h$Var1[h$Var1==11]  <- "matmixed" 
h$Var2[h$Var2==11]  <- "matmixed" 
h$Var1[h$Var1==12]  <- "matwild" 
h$Var2[h$Var2==12]  <- "matwild" 
h$Var1[h$Var1==13]  <- "melcultivated" 
h$Var2[h$Var2==13]  <- "melcultivated" 
h$Var1[h$Var1==14]  <- "melmixed" 
h$Var2[h$Var2==14]  <- "melmixed" 
h$Var1[h$Var1==15]  <- "melwild" 
h$Var2[h$Var2==15]  <- "melwild" 


#Ca
i = as.matrix(dist(soil.params$Soil.Ca.mg.g.))
i = data.frame(as.table(i))[lower.tri(i, diag = TRUE), ]

levels(i$Var1) <- c(levels(i$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 
levels(i$Var2) <- c(levels(i$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 

i$Var1[i$Var1==1]  <- "blocultivated" 
i$Var2[i$Var2==1]  <- "blocultivated" 
i$Var1[i$Var1==2]  <- "blomixed" 
i$Var2[i$Var2==2]  <- "blomixed" 
i$Var1[i$Var1==3]  <- "blowild" 
i$Var2[i$Var2==3]  <- "blowild" 
i$Var1[i$Var1==4]  <- "dobcultivated" 
i$Var2[i$Var2==4]  <- "dobcultivated" 
i$Var1[i$Var1==5]  <- "dobmixed" 
i$Var2[i$Var2==5]  <- "dobmixed" 
i$Var1[i$Var1==6]  <- "dobwild" 
i$Var2[i$Var2==6]  <- "dobwild" 
i$Var1[i$Var1==7]  <- "lancultivated" 
i$Var2[i$Var2==7]  <- "lancultivated" 
i$Var1[i$Var1==8]  <- "lanmixed" 
i$Var2[i$Var2==8]  <- "lanmixed" 
i$Var1[i$Var1==9]  <- "lanwild" 
i$Var2[i$Var2==9]  <- "lanwild" 
i$Var1[i$Var1==10]  <- "matcultivated" 
i$Var2[i$Var2==10]  <- "matcultivated" 
i$Var1[i$Var1==11]  <- "matmixed" 
i$Var2[i$Var2==11]  <- "matmixed" 
i$Var1[i$Var1==12]  <- "matwild" 
i$Var2[i$Var2==12]  <- "matwild" 
i$Var1[i$Var1==13]  <- "melcultivated" 
i$Var2[i$Var2==13]  <- "melcultivated" 
i$Var1[i$Var1==14]  <- "melmixed" 
i$Var2[i$Var2==14]  <- "melmixed" 
i$Var1[i$Var1==15]  <- "melwild" 
i$Var2[i$Var2==15]  <- "melwild" 


#K
j = as.matrix(dist(soil.params$Soil.K.mg.g.))
j = data.frame(as.table(j))[lower.tri(j, diag = TRUE), ]

levels(j$Var1) <- c(levels(j$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 
levels(j$Var2) <- c(levels(j$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 

j$Var1[j$Var1==1]  <- "blocultivated" 
j$Var2[j$Var2==1]  <- "blocultivated" 
j$Var1[j$Var1==2]  <- "blomixed" 
j$Var2[j$Var2==2]  <- "blomixed" 
j$Var1[j$Var1==3]  <- "blowild" 
j$Var2[j$Var2==3]  <- "blowild" 
j$Var1[j$Var1==4]  <- "dobcultivated" 
j$Var2[j$Var2==4]  <- "dobcultivated" 
j$Var1[j$Var1==5]  <- "dobmixed" 
j$Var2[j$Var2==5]  <- "dobmixed" 
j$Var1[j$Var1==6]  <- "dobwild" 
j$Var2[j$Var2==6]  <- "dobwild" 
j$Var1[j$Var1==7]  <- "lancultivated" 
j$Var2[j$Var2==7]  <- "lancultivated" 
j$Var1[j$Var1==8]  <- "lanmixed" 
j$Var2[j$Var2==8]  <- "lanmixed" 
j$Var1[j$Var1==9]  <- "lanwild" 
j$Var2[j$Var2==9]  <- "lanwild" 
j$Var1[j$Var1==10]  <- "matcultivated" 
j$Var2[j$Var2==10]  <- "matcultivated" 
j$Var1[j$Var1==11]  <- "matmixed" 
j$Var2[j$Var2==11]  <- "matmixed" 
j$Var1[j$Var1==12]  <- "matwild" 
j$Var2[j$Var2==12]  <- "matwild" 
j$Var1[j$Var1==13]  <- "melcultivated" 
j$Var2[j$Var2==13]  <- "melcultivated" 
j$Var1[j$Var1==14]  <- "melmixed" 
j$Var2[j$Var2==14]  <- "melmixed" 
j$Var1[j$Var1==15]  <- "melwild" 
j$Var2[j$Var2==15]  <- "melwild" 


#Mg
k = as.matrix(dist(soil.params$Soil.Mg.mg.g.))
k = data.frame(as.table(k))[lower.tri(k, diag = TRUE), ]

levels(k$Var1) <- c(levels(k$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 
levels(k$Var2) <- c(levels(k$Var1), "blocultivated", "blomixed", "blowild", "dobcultivated",
                    "dobmixed", "dobwild", "lancultivated", "lanmixed", "lanwild", "matcultivated",
                    "matmixed", "matwild", "melcultivated", "melmixed", "melwild") 

k$Var1[k$Var1==1]  <- "blocultivated" 
k$Var2[k$Var2==1]  <- "blocultivated" 
k$Var1[k$Var1==2]  <- "blomixed" 
k$Var2[k$Var2==2]  <- "blomixed" 
k$Var1[k$Var1==3]  <- "blowild" 
k$Var2[k$Var2==3]  <- "blowild" 
k$Var1[k$Var1==4]  <- "dobcultivated" 
k$Var2[k$Var2==4]  <- "dobcultivated" 
k$Var1[k$Var1==5]  <- "dobmixed" 
k$Var2[k$Var2==5]  <- "dobmixed" 
k$Var1[k$Var1==6]  <- "dobwild" 
k$Var2[k$Var2==6]  <- "dobwild" 
k$Var1[k$Var1==7]  <- "lancultivated" 
k$Var2[k$Var2==7]  <- "lancultivated" 
k$Var1[k$Var1==8]  <- "lanmixed" 
k$Var2[k$Var2==8]  <- "lanmixed" 
k$Var1[k$Var1==9]  <- "lanwild" 
k$Var2[k$Var2==9]  <- "lanwild" 
k$Var1[k$Var1==10]  <- "matcultivated" 
k$Var2[k$Var2==10]  <- "matcultivated" 
k$Var1[k$Var1==11]  <- "matmixed" 
k$Var2[k$Var2==11]  <- "matmixed" 
k$Var1[k$Var1==12]  <- "matwild" 
k$Var2[k$Var2==12]  <- "matwild" 
k$Var1[k$Var1==13]  <- "melcultivated" 
k$Var2[k$Var2==13]  <- "melcultivated" 
k$Var1[k$Var1==14]  <- "melmixed" 
k$Var2[k$Var2==14]  <- "melmixed" 
k$Var1[k$Var1==15]  <- "melwild" 
k$Var2[k$Var2==15]  <- "melwild" 

parameters = cbind(a,b[,3],c[,3],d[,3],E[,3],f[,3],g[,3],h[,3],i[,3],j[,3],k[,3])

write.csv(cent.distances, "centroid_distances_gyrB.csv")
write.csv(parameters, "soilparams_distances_nodA.csv")


###IMPORT AGAIN DATABASE READY FOR ANALYSES
soils = read.csv("soilparams_distances_nodA_gyrB.csv", header = T)
+ geom_smooth(method="lm", se = T) 

p = qplot(x=soils$WHC...10kPa., y=soils$Centroid.nodA, geom="point") + geom_point(size=5) 
p + theme_bw() + theme_classic()
rcorr = cor.test(soils$WHC...10kPa., soils$Centroid.nodA, method = "pearson")

p = qplot(x=soils$Pot.soil.N.conc, y=soils$Centroid.gyrB, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
p + theme_bw() + theme_classic()
rcorr = cor.test(soils$WHC...10kPa., soils$Centroid.gyrB, method = "pearson")

#nodA correlations
#!!Mg
qplot(x=Soil.Mg.mg.g., y=Centroid.nodA, data=soils, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
cor(soils$Centroid.nodA, soils$Soil.Mg.mg.g.)
lm = lm(soils$Centroid.nodA ~ soils$Soil.Mg.mg.g.)
summary(lm)
#gyrB correlations
#!!WHC gyrB
qplot(x=WHC...10kPa., y=Centroid.gyrB, data=soils, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
cor(soils$Centroid.gyrB, soils$WHC...10kPa.)
lm = lm(soils$Centroid.gyrB ~ soils$WHC...10kPa.)
summary(lm)
#!!Soil N content
qplot(x=Pot.soil.N.conc, y=Centroid.gyrB, data=soils, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
cor(soils$Centroid.gyrB, soils$Pot.soil.N.conc)
lm = lm(soils$Centroid.gyrB ~ soils$Pot.soil.N.conc)
summary(lm)
#!!Soil K content
qplot(x=Soil.K.mg.g., y=Centroid.gyrB, data=soils, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
cor(soils$Centroid.gyrB, soils$Soil.K.mg.g.)
lm = lm(soils$Centroid.gyrB ~ soils$Soil.K.mg.g.)
summary(lm)
#!!Soil Mg content
qplot(x=Soil.Mg.mg.g., y=Centroid.gyrB, data=soils, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
cor(soils$Centroid.gyrB, soils$Soil.Mg.mg.g.)
lm = lm(soils$Centroid.gyrB ~ soils$Soil.Mg.mg.g.)
summary(lm)










Dt = subset_samples(Dt, Soil == "mixed")
bray = vegdist(t(otu_table(Dt)), method="bray", na.omit=T)
aa = as.data.frame(sample_data(Dt))
blo = subset(aa, Farm == "blo")
blo = rownames(blo)
dob = subset(aa, Farm == "dob")
dob = rownames(dob)
lan = subset(aa, Farm == "lan")
lan = rownames(lan)
mat = subset(aa, Farm == "mat")
mat = rownames(mat)
mel = subset(aa, Farm == "mel")
mel = rownames(mel)

library(usedist)
library(cba)
library(funfuns)

centroids.blo.dob = dist_between_centroids(bray, blo, dob)
centroids.blo.lan = dist_between_centroids(bray, blo, lan)
centroids.blo.mat = dist_between_centroids(bray, blo, mat)
centroids.blo.mel = dist_between_centroids(bray, blo, mel)
centroids.dob.lan = dist_between_centroids(bray, dob, lan)
centroids.dob.mat = dist_between_centroids(bray, dob, mat)
centroids.dob.mel = dist_between_centroids(bray, dob, mel)
centroids.lan.mat = dist_between_centroids(bray, lan, mat)
centroids.lan.mel = dist_between_centroids(bray, lan, mel)
centroids.mat.mel = dist_between_centroids(bray, mat, mel)

cent.distances = as.data.frame(rbind(centroids.blo.dob, centroids.blo.lan, centroids.blo.mat, centroids.blo.mel,
                                     centroids.dob.lan, centroids.dob.mat, centroids.dob.mel, centroids.lan.mat,
                                     centroids.lan.mel, centroids.mat.mel))
cent.distances$Distance = cent.distances$V1 
km.distances = as.data.frame(rbind(9.06, 14.25, 17.48, 14.39, 7.67, 12.60, 22.55, 5.23, 28.89, 31.04))
km.distances$Km = km.distances$V1
aa = cbind(cent.distances, km.distances)




lm = lm(aa$Distance ~ aa$Km)
summary(lm)








