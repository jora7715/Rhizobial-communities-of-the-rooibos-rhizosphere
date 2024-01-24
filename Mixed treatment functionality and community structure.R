
rm(list=ls()) # clean/reset environment 
setwd("T:/ETH/Data/Pot_experiment/Rhizobia/nodA_work/test")

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

D <- d.high 

## All samples 
ggrare(D, step = 1000, color = "Farm", label = "X.SampleID", se = FALSE)

## Samples with low counts
D.low <- prune_samples(sample_sums(D) < 3000, D)
ggrare(D.low, step = 1000, color = "Farm", label = "X.SampleID", se = FALSE)

# > Samples below 2000 counts are problematic because diversity is underestimated. 

## Remove samples with counts below 3000
D.b2k <- prune_samples(sample_sums(D) > 3000, D)

#Note: It is important that a sample reaches the plateau otherwise we underestimate the diversity of the sample.

## Remove problematic OTUs

## Show tree
#plot_tree(D.b2k, label.tips="taxa_names", "treeonly", nodeplotblank, ladderize=TRUE)

## The bad
badZOTU <- c("ZOTU152","ZOTU143","ZOTU133","ZOTU119","
             4","ZOTU156","ZOTU160","ZOTU134","ZOTU44","ZOTU153","ZOTU128","ZOTU104","ZOTU127","ZOTU127","ZOTU103","ZOTU151","ZOTU157","ZOTU169","ZOTU146","ZOTU144")

## All
allZOTU  <- taxa_names(D.b2k)

## The good
goodZOUT <- allZOTU[!(allZOTU %in% badZOTU)]

## Keep only the good ones
D.b2k.goodZOTU <- prune_taxa(goodZOUT, D.b2k)
plot_tree(D.b2k.goodZOTU , label.tips="taxa_names", "treeonly", nodeplotblank, ladderize=TRUE)
D <- D.b2k.goodZOTU
D = subset_samples(D, X.SampleID != "r498")

##----------------------------------
#Transform sample counts
##-------------------------------------
Dt <- transform_sample_counts(D, function(OTU) OTU/sum(OTU))


###WE HAVE 3 TYPES OF MIXED COMMUNITIES: 
##EQUAL TO C/W -> BLO
##MORE SIMILAR TO W -> DOB, LAN
##INTERMEDIATE -> MAT, MEL
#We will study the effects of mixing on functioning separately
setwd("T:/ETH/Data/Pot_experiment/Nutrients")
pr = read.csv("Final_database_ambient.csv", header = T)
library(lme4)
library(nlme)

#TEST IF THE MIXES FROM DIFFERENT FARMS LEAD TO DIFFERENT DM
prm = subset(pr, pr$Soil=="mixed")
fit = lme(plant.DM ~ Farm, random= ~1|Block, data = prm)
lmu = anova(fit)
lm = plot(prm$plant.DM~ prm$Farm)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "pr$Farm", conf.level=0.95)

##BLO
blo = subset(pr, Farm == "blo")
blo = blo[,c("Block", "Fertilization", "Soil", "plant.DM", "cluster.presence")]
blo = na.omit(blo)
fit = lme(plant.DM ~ Fertilization*Soil, random= ~1|Block, data = blo)
lmu = anova(fit)
p.adjust(lmu$`p-value`, method = "bonferroni")

library(wesanderson)
p10 <- ggplot(blo, aes(Fertilization, plant.DM, fill = Soil)) + geom_boxplot() + theme_classic()
p10 = p10 + scale_fill_manual(values=wes_palette(n=3, name="Cavalcanti1"))
p10 + geom_point(size = 3, position=position_dodge(width=0.75),aes(group=Soil, colour = cluster.presence)) + scale_color_grey()
#p = qplot(x=root.DM, y=nodule.number, data = pr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
#p + theme_bw() + theme_classic() + scale_color_grey()


##DOB
dob = subset(pr, Farm == "dob")
dob = dob[,c("Block", "Fertilization", "Soil", "plant.DM", "cluster.presence")]
dob = na.omit(dob)
fit = lme(plant.DM ~ Fertilization*Soil, random= ~1|Block, data = dob)
lmu = anova(fit)

p.adjust(lmu$`p-value`, method = "bonferroni")
library(wesanderson)
p10 <- ggplot(dob, aes(Fertilization, plant.DM, fill = Soil)) + geom_boxplot() + theme_classic()
p10 = p10 + scale_fill_manual(values=wes_palette(n=3, name="Cavalcanti1"))
p10 + geom_point(size = 3, position=position_dodge(width=0.75),aes(group=Soil, colour = cluster.presence)) + scale_color_grey()
#p = qplot(x=root.DM, y=nodule.number, data = pr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
#p + theme_bw() + theme_classic() + scale_color_grey()


##lan
lan = subset(pr, Farm == "lan")
lan = lan[,c("Block", "Fertilization", "Soil", "plant.DM", "cluster.presence")]
lan = na.omit(lan)
fit = lme(plant.DM ~ Fertilization*Soil, random= ~1|Block, data = lan)
lmu = anova(fit)

p.adjust(lmu$`p-value`, method = "bonferroni")
library(wesanderson)
p10 <- ggplot(lan, aes(Soil, plant.DM, fill = Soil)) + geom_boxplot() + theme_classic()
p10 = p10 + scale_fill_manual(values=wes_palette(n=3, name="Cavalcanti1"))
p10 + geom_point(size = 3, position=position_dodge(width=0.75),aes(group=Soil, colour = cluster.presence)) + scale_color_grey()
#p = qplot(x=root.DM, y=nodule.number, data = pr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
#p + theme_bw() + theme_classic() + scale_color_grey()



##mat
mat = subset(pr, Farm == "mat")
mat = mat[,c("Block", "Fertilization", "Soil", "plant.DM", "cluster.presence")]
mat = na.omit(mat)
fit = lme(plant.DM ~ Fertilization*Soil, random= ~1|Block, data = mat)
lmu = anova(fit)

p.adjust(lmu$`p-value`, method = "bonferroni")
library(wesanderson)
p10 <- ggplot(mat, aes(Fertilization, plant.DM, fill = Soil)) + geom_boxplot() + theme_classic()
p10 = p10 + scale_fill_manual(values=wes_palette(n=3, name="Cavalcanti1"))
p10 + geom_point(size = 3, position=position_dodge(width=0.75),aes(group=Soil, colour = cluster.presence)) + scale_color_grey()
#p = qplot(x=root.DM, y=nodule.number, data = pr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
#p + theme_bw() + theme_classic() + scale_color_grey()



##mel
mel = subset(pr, Farm == "mel")
mel = mel[,c("Block", "Fertilization", "Soil", "plant.DM", "cluster.presence")]
mel = na.omit(mel)
fit = lme(plant.DM ~ Soil, random= ~1|Block, data = mel)
lmu = anova(fit)

p.adjust(lmu$`p-value`, method = "bonferroni")
library(wesanderson)
p10 <- ggplot(mel, aes(Soil, plant.DM, fill = Soil)) + geom_boxplot() + theme_classic()
p10 = p10 + scale_fill_manual(values=wes_palette(n=3, name="Cavalcanti1"))
p10 + geom_point(size = 3, position=position_dodge(width=0.75),aes(group=Soil, colour = cluster.presence)) + scale_color_grey()
#p = qplot(x=root.DM, y=nodule.number, data = pr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
#p + theme_bw() + theme_classic() + scale_color_grey()


##TUKEY TESTS
#Soil
lm = lm(pr$plant.DM~ pr$Soil)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "pr$Soil", conf.level=0.95)



##NMDS showing differences fertilized vs unfertilized (we know fertilization equalizes functionality, so it should equalize microbial community as well)
##-------------------------------------------------------------------------------------
### Ordination analysis
##------------------------------------------------------------------------------------
Dt = subset_samples(Dt, Soil=="mixed")
# Ordinate data with PCoA
ord <- ordinate(Dt, "NMDS", "bray")
ord <- ordinate(D, "NMDS", "bray")


plot_ordination(Dt, ord,  color = "Farm", shape = "Soil") +
  geom_point(size = 5)+ theme_bw() + theme_classic() + scale_color_manual(values=wes_palette(n=5, name="Cavalcanti1")) + stat_ellipse(type="norm", level = 0.85) 

+ facet_wrap(~Farm)
##NMDS has >0.25 stress!!!!
#NMDS farm
library(wesanderson)
ord <- ordinate(Dt, "NMDS", "bray")
ord <- ordinate(Dt, "MDS")
# Multidimensional scaling (MDS / PCoA)
plot_ordination(Dt, ord, color = "Farm") + scale_color_manual(values=wes_palette(n=5, name="Cavalcanti1")) + geom_point(size = 5) + theme_bw() 

p + stat_ellipse(type="norm")
facet_wrap(~Farm)


#--------------------------------------------------------------------------------------
#PERMANOVA WITHIN FARM
#--------------------------------------------------------------------------------------
Dt <- transform_sample_counts(D, function(OTU) OTU/sum(OTU))
library(pairwiseAdonis)
library(vegan)
#Interactions
a = meta(Dt)
library(phyloseq)
bcmatrix = phyloseq::distance(Dt, method = "bray")
bcmatrix = na.omit(bcmatrix)
#adonis(t(as.data.frame(otu_table(Dt))) ~ Soil, data = a, permutations = 999, method = "bray")
adonis <- adonis(bcmatrix ~ Soil/Farm + Fertilization*Soil/Farm , strata = a$Farm, data = a, permutations = 9999)

##PERMANOVA ON EACH FARM TO COMPARE SHIFT DUE TO FERTILIZATION
blo = subset_samples(Dt, Farm == "blo")
a = meta(blo)
bcmatrix = phyloseq::distance(blo, method = "bray")
adonis(bcmatrix ~ Soil * Fertilization, strata =a$Block, data = a, permutations = 9999)

dob = subset_samples(Dt, Farm == "dob")
a = meta(dob)
bcmatrix = phyloseq::distance(dob, method = "bray")
adonis(bcmatrix ~ Soil * Fertilization, strata =a$Block, data = a, permutations = 9999)

lan = subset_samples(Dt, Farm == "lan")
a = meta(lan)
bcmatrix = phyloseq::distance(lan, method = "bray")
adonis(bcmatrix ~ Soil * Fertilization, strata =a$Block, data = a, permutations = 9999)

mat = subset_samples(Dt, Farm == "mat")
a = meta(mat)
bcmatrix = phyloseq::distance(mat, method = "bray")
adonis(bcmatrix ~ Soil * Fertilization, strata =a$Block, data = a, permutations = 9999)

mel = subset_samples(Dt, Farm == "mel")
a = meta(mel)
bcmatrix = phyloseq::distance(mel, method = "bray")
adonis(bcmatrix ~ Soil * Fertilization, strata =a$Block, data = a, permutations = 9999)


#PERMANOVA TESTING IF FERTILIZATION AFFECTS MORE STRONGLY THE MIXED SOILS THAN THE OTHERS
blo = subset_samples(Dt, Farm == "blo")
blo = subset_samples(blo, Soil == "wild")
a = meta(blo)
bcmatrix = phyloseq::distance(blo, method = "bray")
adonis(bcmatrix ~ Fertilization, strata =a$Block, data = a, permutations = 999)

dob = subset_samples(Dt, Farm == "dob")
dob = subset_samples(dob, Soil == "cultivated")
a = meta(dob)
bcmatrix = phyloseq::distance(dob, method = "bray")
adonis(bcmatrix ~ Fertilization, strata =a$Block, data = a, permutations = 999)

lan = subset_samples(Dt, Farm == "lan")
lan = subset_samples(lan, Soil == "cultivated")
a = meta(lan)
bcmatrix = phyloseq::distance(lan, method = "bray")
adonis(bcmatrix ~ Fertilization, strata =a$Block, data = a, permutations = 999)

mat = subset_samples(Dt, Farm == "mat")
mat = subset_samples(mat, Soil == "cultivated")
a = meta(mat)
bcmatrix = phyloseq::distance(mat, method = "bray")
adonis(bcmatrix ~ Fertilization, strata =a$Block, data = a, permutations = 999)

mel = subset_samples(Dt, Farm == "mel")
mel = subset_samples(mel, Soil == "cultivated")
a = meta(mel)
bcmatrix = phyloseq::distance(mel, method = "bray")
adonis(bcmatrix ~ Fertilization, strata =a$Block, data = a, permutations = 999)


##BETADISPER: DOES FERTILIZATION OR MIXING LEAD TO MORE ALTERNATIVE COMMUNITIES?
bdist = phyloseq::distance(Dt, method = "bray")
fert <- as(sample_data(Dt), "data.frame")[,"Fertilization"]
soil <- as(sample_data(Dt), "data.frame")[,"Soil"]
betatax.fert <- betadisper(bdist, fert)
betatax.soil <- betadisper(bdist, soil)

#BETADISPER BY FARM
blo = subset_samples(Dt, Farm == "blo")
bdist = phyloseq::distance(blo, method = "bray")
fert <- as(sample_data(blo), "data.frame")[,"Fertilization"]
soil <- as(sample_data(blo), "data.frame")[,"Soil"]
betatax.fert <- betadisper(bdist, fert)
betatax.soil <- betadisper(bdist, soil)

dob = subset_samples(Dt, Farm == "dob")
bdist = phyloseq::distance(dob, method = "bray")
fert <- as(sample_data(dob), "data.frame")[,"Fertilization"]
soil <- as(sample_data(dob), "data.frame")[,"Soil"]
betatax.fert <- betadisper(bdist, fert)
betatax.soil <- betadisper(bdist, soil)

lan = subset_samples(Dt, Farm == "lan")
bdist = phyloseq::distance(lan, method = "bray")
fert <- as(sample_data(lan), "data.frame")[,"Fertilization"]
soil <- as(sample_data(lan), "data.frame")[,"Soil"]
betatax.fert <- betadisper(bdist, fert)
betatax.soil <- betadisper(bdist, soil)

mat = subset_samples(Dt, Farm == "mat")
bdist = phyloseq::distance(mat, method = "bray")
fert <- as(sample_data(mat), "data.frame")[,"Fertilization"]
soil <- as(sample_data(mat), "data.frame")[,"Soil"]
betatax.fert <- betadisper(bdist, fert)
betatax.soil <- betadisper(bdist, soil)

mel = subset_samples(Dt, Farm == "mel")
bdist = phyloseq::distance(mel, method = "bray")
fert <- as(sample_data(mel), "data.frame")[,"Fertilization"]
soil <- as(sample_data(mel), "data.frame")[,"Soil"]
betatax.fert <- betadisper(bdist, fert)
betatax.soil <- betadisper(bdist, soil)


##BETADISPER: DOES FERTILIZATION MAKE MIXED COMMUNITIES CONVERGE WITHIN TREATMENT?
mix = subset_samples(Dt, Soil == "mixed")
bdist = phyloseq::distance(mix, method = "bray")
fert <- as(sample_data(mix), "data.frame")[,"Fertilization"]
betatax.fert <- betadisper(bdist, fert)

blo = subset_samples(mix, Farm == "blo")
bdist = phyloseq::distance(blo, method = "bray")
fert <- as(sample_data(blo), "data.frame")[,"Fertilization"]
betatax.fert <- betadisper(bdist, fert)

dob = subset_samples(mix, Farm == "dob")
bdist = phyloseq::distance(dob, method = "bray")
fert <- as(sample_data(dob), "data.frame")[,"Fertilization"]
betatax.fert <- betadisper(bdist, fert)

lan = subset_samples(mix, Farm == "lan")
bdist = phyloseq::distance(lan, method = "bray")
fert <- as(sample_data(lan), "data.frame")[,"Fertilization"]
betatax.fert <- betadisper(bdist, fert)

mat = subset_samples(mix, Farm == "mat")
bdist = phyloseq::distance(mat, method = "bray")
fert <- as(sample_data(mat), "data.frame")[,"Fertilization"]
betatax.fert <- betadisper(bdist, fert)

mel = subset_samples(mix, Farm == "mel")
bdist = phyloseq::distance(mel, method = "bray")
fert <- as(sample_data(mel), "data.frame")[,"Fertilization"]
betatax.fert <- betadisper(bdist, fert)

##Pairwise adonis to quantify differences between soil and fertilization treatments within farm
library(pairwiseAdonis)
blo = subset_samples(Dt, Farm == "blo" & Fertilization == "fertilized")
a = meta(blo)
dist.pa = distance(blo, method="bray")
pairwise.adonis(dist.pa, factors = a$Soil, sim.method = "bray", p.adjust.m = "none")

dob = subset_samples(Dt, Farm == "dob"& Fertilization == "fertilized")
a = meta(dob)
dist.pa = distance(dob, method="bray")
pairwise.adonis(dist.pa, factors = a$Soil, sim.method = "bray", p.adjust.m = "none")

lan = subset_samples(Dt, Farm == "lan"& Fertilization == "fertilized")
a = meta(lan)
dist.pa = distance(lan, method="bray")
pairwise.adonis(dist.pa, factors = a$Soil, sim.method = "bray", p.adjust.m = "none")

mat = subset_samples(Dt, Farm == "mat"& Fertilization == "fertilized")
a = meta(mat)
dist.pa = distance(mat, method="bray")
pairwise.adonis(dist.pa, factors = a$Soil, sim.method = "bray", p.adjust.m = "none")

mel = subset_samples(Dt, Farm == "mel"& Fertilization == "fertilized")
a = meta(mel)
dist.pa = distance(mel, method="bray")
pairwise.adonis(dist.pa, factors = a$Soil, sim.method = "bray", p.adjust.m = "none")


#FERTILIZATION EFFECTS ON MIXED SOILS WITHIN FARM
Dt = subset_samples(Dt, Soil == "wild")
blo = subset_samples(Dt, Farm == "blo")
a = meta(blo)
dist.pa = distance(blo, method="bray")
pairwise.adonis(dist.pa, factors = a$Fertilization, sim.method = "bray", p.adjust.m = "none")

dob = subset_samples(Dt, Farm == "dob")
a = meta(dob)
dist.pa = distance(dob, method="bray")
pairwise.adonis(dist.pa, factors = a$Fertilization, sim.method = "bray", p.adjust.m = "none")

lan = subset_samples(Dt, Farm == "lan")
a = meta(lan)
dist.pa = distance(lan, method="bray")
pairwise.adonis(dist.pa, factors = a$Fertilization, sim.method = "bray", p.adjust.m = "none")

mat = subset_samples(Dt, Farm == "mat")
a = meta(mat)
dist.pa = distance(mat, method="bray")
pairwise.adonis(dist.pa, factors = a$Fertilization, sim.method = "bray", p.adjust.m = "none")

mel = subset_samples(Dt, Farm == "mel")
a = meta(mel)
dist.pa = distance(mel, method="bray")
pairwise.adonis(dist.pa, factors = a$Fertilization, sim.method = "bray", p.adjust.m = "none")


###INDICATOR SPECIES ANALYSIS###
Dt <- transform_sample_counts(D, function(OTU) OTU/sum(OTU))
#install.packages("indicspecies")
library(indicspecies)
#Indicator species of soils by farm
blo = subset_samples(Dt, Farm == "blo")
a = meta(blo)
indicspecies.blo = multipatt(as.data.frame(t(otu_table(blo))), cluster = a$Soil, control = how(nperm = 9999))
summary(indicspecies.blo)

dob = subset_samples(Dt, Farm == "dob")
a = meta(dob)
indicspecies.dob = multipatt(as.data.frame(t(otu_table(dob))), cluster = a$Soil, control = how(nperm = 9999))
summary(indicspecies.dob)

lan = subset_samples(Dt, Farm == "lan")
a = meta(lan)
indicspecies.lan = multipatt(as.data.frame(t(otu_table(lan))), cluster = a$Soil, control = how(nperm = 9999))
summary(indicspecies.lan)

mat = subset_samples(Dt, Farm == "mat")
a = meta(mat)
indicspecies.mat = multipatt(as.data.frame(t(otu_table(mat))), cluster = a$Soil, control = how(nperm = 9999))
summary(indicspecies.mat)

mel = subset_samples(Dt, Farm == "mel")
a = meta(mel)
indicspecies.mel = multipatt(as.data.frame(t(otu_table(mel))), cluster = a$Soil, control = how(nperm = 9999))
summary(indicspecies.mel)


#Indicator species of fertilization by farm
blo = subset_samples(Dt, Farm == "blo")
a = meta(blo)
indicspecies.blo = multipatt(as.data.frame(t(otu_table(blo))), cluster = a$Fertilization, control = how(nperm = 9999))
summary(indicspecies.blo)

dob = subset_samples(Dt, Farm == "dob")
a = meta(dob)
indicspecies.dob = multipatt(as.data.frame(t(otu_table(dob))), cluster = a$Fertilization, control = how(nperm = 9999))
summary(indicspecies.dob)

lan = subset_samples(Dt, Farm == "lan")
a = meta(lan)
indicspecies.lan = multipatt(as.data.frame(t(otu_table(lan))), cluster = a$Fertilization, control = how(nperm = 9999))
summary(indicspecies.lan)

mat = subset_samples(Dt, Farm == "mat")
a = meta(mat)
indicspecies.mat = multipatt(as.data.frame(t(otu_table(mat))), cluster = a$Fertilization, control = how(nperm = 9999))
summary(indicspecies.mat)

mel = subset_samples(Dt, Farm == "mel")
a = meta(mel)
indicspecies.mel = multipatt(as.data.frame(t(otu_table(mel))), cluster = a$Fertilization, control = how(nperm = 9999))
summary(indicspecies.mel)

##PLOT CHANGES IN RELATIVE ABUNDANCE BY OTU TO COMPARE THE MIXED TO THE OTHER TREATMENTS
blom = subset_samples(blo, Soil == "mixed")
blow = subset_samples(blo, Soil == "wild")
bloc = subset_samples(blo, Soil == "cultivated")

##DESEQ DIFFERENTIAL ABUNDANCE ANALYSIS
library(DESeq2)
Dt = filter_taxa(D, function(OTU) sum(OTU) > 100, TRUE)
#SOILS
soil = phyloseq_to_deseq2(Dt, ~ 0 + Soil)
soil$Soil = relevel(soil$Soil, "mixed")
mdl = DESeq(soil, test="Wald", fitType="parametric")
res = results(mdl, contrast = c("Soil","cultivated","wild"))
ggplot(as.data.frame(res), aes(x=rownames(res), y=log2FoldChange)) + geom_point(size=4) +
  theme_classic() + geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))
write.csv(res, "cult_wild.csv")
#FERTILIZATION
fert = phyloseq_to_deseq2(Dt, ~ 0 + Fertilization)
fert$Fertilization = relevel(fert$Fertilization, "unfertilized")
mdl = DESeq(fert, test="Wald", fitType="parametric")
res = results(mdl, contrast = c("Fertilization","unfertilized","fertilized"))
ggplot(as.data.frame(res), aes(x=rownames(res), y=log2FoldChange)) + geom_point(size=4) +
  theme_classic() + geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))
write.csv(res, "unfert_fert.csv")
#SOILS WITHIN FERTILIZATION TREATMENT
fsub = subset_samples(Dt, Fertilization == "fertilized")
fsub = phyloseq_to_deseq2(fsub, ~ 0 + Soil)
fsub$Soil = relevel(fsub$Soil, "mixed")
mdl = DESeq(fsub, test="Wald", fitType="parametric")
res = results(mdl, contrast = c("Soil","cultivated","wild"))
ggplot(as.data.frame(res), aes(x=rownames(res), y=log2FoldChange)) + geom_point(size=4) +
  theme_classic() + geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))
write.csv(res, "fert_cultivated_wild.csv")
#FERTILIZATION WITHIN MIX
mix = subset_samples(Dt, Soil == "mixed")
fsub = phyloseq_to_deseq2(mix, ~ 0 + Fertilization)
fsub$Fertilization = relevel(fsub$Fertilization, "unfertilized")
mdl = DESeq(fsub, test="Wald", fitType="parametric")
res = results(mdl, contrast = c("Fertilization","unfertilized","fertilized"))
ggplot(as.data.frame(res), aes(x=rownames(res), y=log2FoldChange)) + geom_point(size=4) +
  theme_classic() + geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))
write.csv(res, "mix_unfert_fert.csv")
#SOILS WITHIN FARM
mel = subset_samples(Dt, Farm == "mel")
fsub = phyloseq_to_deseq2(mel, ~ 0 + Soil)
fsub$Soil = relevel(fsub$Soil, "mixed")
mdl = DESeq(fsub, test="Wald", fitType="parametric")
res = results(mdl, contrast = c("Soil","cultivated","wild"))
ggplot(as.data.frame(res), aes(x=rownames(res), y=log2FoldChange)) + geom_point(size=4) +
  theme_classic() + geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))
write.csv(res, "mel_cultivated_wild.csv")
#FERTILIZATION WITHIN FARM
mel = subset_samples(Dt, Farm == "mel")
fsub = phyloseq_to_deseq2(mel, ~ 0 + Fertilization)
fsub$Fertilization = relevel(fsub$Fertilization, "unfertilized")
mdl = DESeq(fsub, test="Wald", fitType="parametric")
res = results(mdl, contrast = c("Fertilization","unfertilized","fertilized"))
ggplot(as.data.frame(res), aes(x=rownames(res), y=log2FoldChange)) + geom_point(size=4) +
  theme_classic() + geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))
write.csv(res, "mel_unfert_fert.csv")
#FERTILIZED X SOIL WITHIN FARM
mel = subset_samples(Dt, Farm == "mel" & Fertilization == "unfertilized" & Soil == "cultivated")
fsub = phyloseq_to_deseq2(mel, ~ 0 + Block)
fsub$Soil = relevel(fsub$Soil, "mixed")
mdl = DESeq(fsub, test="Wald", fitType="parametric")
res = results(mdl, contrast = c("Soil","cultivated","wild"))
ggplot(as.data.frame(res), aes(x=rownames(res), y=log2FoldChange)) + geom_point(size=4) +
  theme_classic() + geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))
write.csv(res, "mel_fertilized_cultivated_wild.csv")
#UNFERTILIZED X SOIL WITHIN FARM
mel = subset_samples(Dt, Farm == "blo" & Soil == "cultivated")
fsub = phyloseq_to_deseq2(mel, ~ 0 + Fertilization)
fsub$Fertilization = relevel(fsub$Fertilization, "unfertilized")
mdl = DESeq(fsub, test="Wald", fitType="parametric")
res = results(mdl, contrast = c("Fertilization","unfertilized","fertilized"))
res$OTU = rownames(res)
res = res[order(res$log2FoldChange),]
res$OTU <- factor(res$OTU, levels = res$OTU[order(res$log2FoldChange)])
ggplot(as.data.frame(res), aes(x=OTU, y=log2FoldChange)) + geom_point(size=4) +
  theme_classic() + geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))
write.csv(res, "blo_cultivated_unfert_fert.csv")

#GET NORMALIZED COUNTS
counts(mdl, normalized=TRUE)


##DIVERGENCE INDEX PLOT BASED ON DIFFERENCES IN THE R2 OF UNFERT - FERT SOIL COMPARISONS
dgc = read.csv("Divergence_index_plot.csv", header = T)
ggplot(dgc, aes(x=Farm, y=Divergence_index, color = Soil_comparison)) + geom_point(size=4) + theme_classic() 

ddb = read.csv("Convergence_divergence_Pairwiseadonis.csv", header = T)  
ddb = subset(ddb, Fertilization == "fertilized")
ggplot(ddb, aes(x=Farm, y=R2, color = Soil_comparison)) + geom_point(size=4) + theme_classic() 

##NORMALIZED ABUNDANCE PLOT FOR OVERDOMINANT ZOTUS
ddb = read.csv("zotu12.csv", header = T)  
#ddb = subset(ddb, Fertilization == "unfertilized")
ggplot(ddb, aes(x=Soil, y=Abundance)) + geom_point(size=4) +
  theme_classic() + geom_errorbar(aes(ymin=Abundance-SE, ymax=Abundance+SE), width = .1)

##FERTILIZATION EFFECTS ON DIFFERENTIAL ABUNDANCES
ddb = read.csv("Differential_abundance_fertilization.csv", header = T)  
ggplot(ddb, aes(x=Fertilization, y=nOTUs, shape = Soil_comparison)) + geom_point(size=4) +
  theme_classic() 
 
ggplot(ddb, aes(x=Fertilization, y=log2fold)) + theme_classic() + geom_errorbar(aes(ymin=log2fold-sd, ymax=log2fold+sd), width = .1)
png("boxplot.png", width = 1000, height = 2500, res = 300)
boxplot(log2fold ~ Fertilization, data = ddb)
dev.off()
a = lm(log2fold ~ Fertilization, data = ddb)
summary(a)


setwd("E:/ETH/Data/Pot_experiment/Rhizobia/nodA_work/test")
ddb = read.csv("Differential_abundance_fertilization.csv", header = T)  
ddb = subset(ddb, Farm == "mel")
ggplot(ddb, aes(x=Fertilization, y=nOTUs, shape = Soil_comparison)) + geom_point(size=4) +
  theme_classic() 


##DIST TO CENTROID ANALYSIS TO RELATE TO SOIL FACTORS WITHOUT VIOLATING STATISTICAL ASSUMPTIONS
aa = as.data.frame(sample_data(Dt))

library(usedist)
library(cba)
library(funfuns)

#BLO
blo = subset_samples(Dt, Farm == "blo")
blo.UC = subset_samples(blo, Fertilization == "unfertilized" & Soil == "cultivated")
blo.UM = subset_samples(blo, Fertilization == "unfertilized" & Soil == "mixed")
blo.UW = subset_samples(blo, Fertilization == "unfertilized" & Soil == "wild")
blo.FC = subset_samples(blo, Fertilization == "fertilized" & Soil == "cultivated")
blo.FM = subset_samples(blo, Fertilization == "fertilized" & Soil == "mixed")
blo.FW = subset_samples(blo, Fertilization == "fertilized" & Soil == "wild")
blo.UC.names = rownames(as.data.frame(sample_data(blo.UC)))
blo.UM.names = rownames(as.data.frame(sample_data(blo.UM)))
blo.UW.names = rownames(as.data.frame(sample_data(blo.UW)))
blo.FC.names = rownames(as.data.frame(sample_data(blo.FC)))
blo.FM.names = rownames(as.data.frame(sample_data(blo.FM)))
blo.FW.names = rownames(as.data.frame(sample_data(blo.FW)))

bray.blo = vegdist(t(otu_table(blo)), method="bray", na.omit=T)

centroids.U.CW = dist_between_centroids(bray.blo, blo.FC.names, blo.FW.names)
centroids.U.MC = dist_between_centroids(bray.blo, blo.FC.names, blo.FM.names)
centroids.U.MW = dist_between_centroids(bray.blo, blo.FW.names, blo.FM.names)

#DOB
dob = subset_samples(Dt, Farm == "dob")
dob.UC = subset_samples(dob, Fertilization == "unfertilized" & Soil == "cultivated")
dob.UM = subset_samples(dob, Fertilization == "unfertilized" & Soil == "mixed")
dob.UW = subset_samples(dob, Fertilization == "unfertilized" & Soil == "wild")
dob.FC = subset_samples(dob, Fertilization == "fertilized" & Soil == "cultivated")
dob.FM = subset_samples(dob, Fertilization == "fertilized" & Soil == "mixed")
dob.FW = subset_samples(dob, Fertilization == "fertilized" & Soil == "wild")
dob.UC.names = rownames(as.data.frame(sample_data(dob.UC)))
dob.UM.names = rownames(as.data.frame(sample_data(dob.UM)))
dob.UW.names = rownames(as.data.frame(sample_data(dob.UW)))
dob.FC.names = rownames(as.data.frame(sample_data(dob.FC)))
dob.FM.names = rownames(as.data.frame(sample_data(dob.FM)))
dob.FW.names = rownames(as.data.frame(sample_data(dob.FW)))

bray.dob = vegdist(t(otu_table(dob)), method="bray", na.omit=T)

centroids.U.CW = dist_between_centroids(bray.dob, dob.FC.names, dob.FW.names)
centroids.U.MC = dist_between_centroids(bray.dob, dob.FC.names, dob.FM.names)
centroids.U.MW = dist_between_centroids(bray.dob, dob.FW.names, dob.FM.names)

#LAN
lan = subset_samples(Dt, Farm == "lan")
lan.UC = subset_samples(lan, Fertilization == "unfertilized" & Soil == "cultivated")
lan.UM = subset_samples(lan, Fertilization == "unfertilized" & Soil == "mixed")
lan.UW = subset_samples(lan, Fertilization == "unfertilized" & Soil == "wild")
lan.FC = subset_samples(lan, Fertilization == "fertilized" & Soil == "cultivated")
lan.FM = subset_samples(lan, Fertilization == "fertilized" & Soil == "mixed")
lan.FW = subset_samples(lan, Fertilization == "fertilized" & Soil == "wild")
lan.UC.names = rownames(as.data.frame(sample_data(lan.UC)))
lan.UM.names = rownames(as.data.frame(sample_data(lan.UM)))
lan.UW.names = rownames(as.data.frame(sample_data(lan.UW)))
lan.FC.names = rownames(as.data.frame(sample_data(lan.FC)))
lan.FM.names = rownames(as.data.frame(sample_data(lan.FM)))
lan.FW.names = rownames(as.data.frame(sample_data(lan.FW)))

bray.lan = vegdist(t(otu_table(lan)), method="bray", na.omit=T)

centroids.U.CW = dist_between_centroids(bray.lan, lan.UC.names, lan.UW.names)
centroids.U.MC = dist_between_centroids(bray.lan, lan.UC.names, lan.UM.names)
centroids.U.MW = dist_between_centroids(bray.lan, lan.UW.names, lan.UM.names)

#MAT
mat = subset_samples(Dt, Farm == "mat")
mat.UC = subset_samples(mat, Fertilization == "unfertilized" & Soil == "cultivated")
mat.UM = subset_samples(mat, Fertilization == "unfertilized" & Soil == "mixed")
mat.UW = subset_samples(mat, Fertilization == "unfertilized" & Soil == "wild")
mat.FC = subset_samples(mat, Fertilization == "fertilized" & Soil == "cultivated")
mat.FM = subset_samples(mat, Fertilization == "fertilized" & Soil == "mixed")
mat.FW = subset_samples(mat, Fertilization == "fertilized" & Soil == "wild")
mat.UC.names = rownames(as.data.frame(sample_data(mat.UC)))
mat.UM.names = rownames(as.data.frame(sample_data(mat.UM)))
mat.UW.names = rownames(as.data.frame(sample_data(mat.UW)))
mat.FC.names = rownames(as.data.frame(sample_data(mat.FC)))
mat.FM.names = rownames(as.data.frame(sample_data(mat.FM)))
mat.FW.names = rownames(as.data.frame(sample_data(mat.FW)))

bray.mat = vegdist(t(otu_table(mat)), method="bray", na.omit=T)

centroids.U.CW = dist_between_centroids(bray.mat, mat.FC.names, mat.FW.names)
centroids.U.MC = dist_between_centroids(bray.mat, mat.FC.names, mat.FM.names)
centroids.U.MW = dist_between_centroids(bray.mat, mat.FW.names, mat.FM.names)

#MEL
mel = subset_samples(Dt, Farm == "mel")
mel.UC = subset_samples(mel, Fertilization == "unfertilized" & Soil == "cultivated")
mel.UM = subset_samples(mel, Fertilization == "unfertilized" & Soil == "mixed")
mel.UW = subset_samples(mel, Fertilization == "unfertilized" & Soil == "wild")
mel.FC = subset_samples(mel, Fertilization == "fertilized" & Soil == "cultivated")
mel.FM = subset_samples(mel, Fertilization == "fertilized" & Soil == "mixed")
mel.FW = subset_samples(mel, Fertilization == "fertilized" & Soil == "wild")
mel.UC.names = rownames(as.data.frame(sample_data(mel.UC)))
mel.UM.names = rownames(as.data.frame(sample_data(mel.UM)))
mel.UW.names = rownames(as.data.frame(sample_data(mel.UW)))
mel.FC.names = rownames(as.data.frame(sample_data(mel.FC)))
mel.FM.names = rownames(as.data.frame(sample_data(mel.FM)))
mel.FW.names = rownames(as.data.frame(sample_data(mel.FW)))

bray.mel = vegdist(t(otu_table(mel)), method="bray", na.omit=T)

centroids.U.CW = dist_between_centroids(bray.mel, mel.UC.names, mel.UW.names)
centroids.U.MC = dist_between_centroids(bray.mel, mel.UC.names, mel.UM.names)
centroids.U.MW = dist_between_centroids(bray.mel, mel.UW.names, mel.UM.names)

##PLOT DISTANCES TO CENTROIDS TO VISUALIZE CONVERGENCE/DIVERGENCE
cent = read.csv("Centroids_soils_fertilization.csv", header = T)
#blo
cent.blo = subset(cent, Farm == "blo")
ggplot(cent.blo, aes(x=Fertilization, y=Centroid_distance, shape = Comparison)) + geom_point(size=10) +
  theme_classic() 
#dob
cent.dob = subset(cent, Farm == "dob")
ggplot(cent.dob, aes(x=Fertilization, y=Centroid_distance, shape = Comparison)) + geom_point(size=10) +
  theme_classic() 
#lan
cent.lan = subset(cent, Farm == "lan")
ggplot(cent.lan, aes(x=Fertilization, y=Centroid_distance, shape = Comparison)) + geom_point(size=10) +
  theme_classic() 
#mat
cent.mat = subset(cent, Farm == "mat")
ggplot(cent.mat, aes(x=Fertilization, y=Centroid_distance, shape = Comparison)) + geom_point(size=10) +
  theme_classic()
#mel
cent.mel = subset(cent, Farm == "mel")
ggplot(cent.mel, aes(x=Fertilization, y=Centroid_distance, shape = Comparison)) + geom_point(size=10) +
  theme_classic()



##CORRELATIONS BETWEEN PLANT RESPONSES AND REL ABUNDANCE OF ZOTUS OVER-EXPRESSED IN MIXED, FERTILIZED SOILS
fert = subset_samples(Dt, Fertilization == "fertilized")
response = as.data.frame(meta(fert))
otus = as.data.frame(t(get_sample(fert, c("ZOTU19", "ZOTU21", "ZOTU28", "ZOTU40", "ZOTU68", "ZOTU12", "ZOTU43", "ZOTU93", "ZOTU33", "ZOTU34", "ZOTU61", "ZOTU14", "ZOTU49", "ZOTU36", "ZOTU166", "ZOTU164", "ZOTU17", "ZOTU6", "ZOTU23", "ZOTU173", "ZOTU52", "ZOTU155", "ZOTU50", "ZOTU105", "ZOTU70", "ZOTU148"))))

##PLANT DM
r19 = cor.test(response$plant.DM, otus$ZOTU19, method = "pearson")                     
r21 = cor.test(response$plant.DM, otus$ZOTU21, method = "pearson")                     
r28 = cor.test(response$plant.DM, otus$ZOTU28, method = "pearson")                     
r40 = cor.test(response$plant.DM, otus$ZOTU40, method = "pearson")                     
r68 = cor.test(response$plant.DM, otus$ZOTU68, method = "pearson")                     
r12 = cor.test(response$plant.DM, otus$ZOTU12, method = "pearson")#r = 0.153, p = 0.014                     
r43 = cor.test(response$plant.DM, otus$ZOTU43, method = "pearson")#r = -0.185, p = 0.003
r93 = cor.test(response$plant.DM, otus$ZOTU93, method = "pearson")
r33 = cor.test(response$plant.DM, otus$ZOTU33, method = "pearson")
r34 = cor.test(response$plant.DM, otus$ZOTU34, method = "pearson")
r61 = cor.test(response$plant.DM, otus$ZOTU61, method = "pearson")
r14 = cor.test(response$plant.DM, otus$ZOTU14, method = "pearson")
r49 = cor.test(response$plant.DM, otus$ZOTU49, method = "pearson")
r36 = cor.test(response$plant.DM, otus$ZOTU36, method = "pearson")
r166 = cor.test(response$plant.DM, otus$ZOTU166, method = "pearson")
r164 = cor.test(response$plant.DM, otus$ZOTU164, method = "pearson")
r17 = cor.test(response$plant.DM, otus$ZOTU17, method = "pearson")
r6 = cor.test(response$plant.DM, otus$ZOTU6, method = "pearson")
r23 = cor.test(response$plant.DM, otus$ZOTU23, method = "pearson")
r173 = cor.test(response$plant.DM, otus$ZOTU173, method = "pearson")
r52 = cor.test(response$plant.DM, otus$ZOTU52, method = "pearson")
r155 = cor.test(response$plant.DM, otus$ZOTU155, method = "pearson")
r50 = cor.test(response$plant.DM, otus$ZOTU50, method = "pearson")
r105 = cor.test(response$plant.DM, otus$ZOTU105, method = "pearson")
r70 = cor.test(response$plant.DM, otus$ZOTU70, method = "pearson")
r148 = cor.test(response$plant.DM, otus$ZOTU148, method = "pearson")

##LEAF N
r19 = cor.test(response$Leaf.N.content, otus$ZOTU19, method = "pearson")                     
r21 = cor.test(response$Leaf.N.content, otus$ZOTU21, method = "pearson")                     
r28 = cor.test(response$Leaf.N.content, otus$ZOTU28, method = "pearson")                     
r40 = cor.test(response$Leaf.N.content, otus$ZOTU40, method = "pearson")                     
r68 = cor.test(response$Leaf.N.content, otus$ZOTU68, method = "pearson")                     
r12 = cor.test(response$Leaf.N.content, otus$ZOTU12, method = "pearson")#r = 0.153, p = 0.014                     
r43 = cor.test(response$Leaf.N.content, otus$ZOTU43, method = "pearson")#r = -0.185, p = 0.003
r93 = cor.test(response$Leaf.N.content, otus$ZOTU93, method = "pearson")
r33 = cor.test(response$Leaf.N.content, otus$ZOTU33, method = "pearson")
r34 = cor.test(response$Leaf.N.content, otus$ZOTU34, method = "pearson")
r61 = cor.test(response$Leaf.N.content, otus$ZOTU61, method = "pearson")
r14 = cor.test(response$Leaf.N.content, otus$ZOTU14, method = "pearson")
r49 = cor.test(response$Leaf.N.content, otus$ZOTU49, method = "pearson")
r36 = cor.test(response$Leaf.N.content, otus$ZOTU36, method = "pearson")
r166 = cor.test(response$Leaf.N.content, otus$ZOTU166, method = "pearson")
r164 = cor.test(response$Leaf.N.content, otus$ZOTU164, method = "pearson")
r17 = cor.test(response$Leaf.N.content, otus$ZOTU17, method = "pearson")
r6 = cor.test(response$Leaf.N.content, otus$ZOTU6, method = "pearson")
r23 = cor.test(response$Leaf.N.content, otus$ZOTU23, method = "pearson")
r173 = cor.test(response$Leaf.N.content, otus$ZOTU173, method = "pearson")
r52 = cor.test(response$Leaf.N.content, otus$ZOTU52, method = "pearson")
r155 = cor.test(response$Leaf.N.content, otus$ZOTU155, method = "pearson")
r50 = cor.test(response$Leaf.N.content, otus$ZOTU50, method = "pearson")
r105 = cor.test(response$Leaf.N.content, otus$ZOTU105, method = "pearson")
r70 = cor.test(response$Leaf.N.content, otus$ZOTU70, method = "pearson")
r148 = cor.test(response$Leaf.N.content, otus$ZOTU148, method = "pearson")

##NODULES PER ROOT DM
r19 = cor.test(response$nodules.rootDM, otus$ZOTU19, method = "pearson")                     
r21 = cor.test(response$nodules.rootDM, otus$ZOTU21, method = "pearson")                     
r28 = cor.test(response$nodules.rootDM, otus$ZOTU28, method = "pearson")                     
r40 = cor.test(response$nodules.rootDM, otus$ZOTU40, method = "pearson")                     
r68 = cor.test(response$nodules.rootDM, otus$ZOTU68, method = "pearson")                     
r12 = cor.test(response$nodules.rootDM, otus$ZOTU12, method = "pearson")#r = 0.153, p = 0.014                     
r43 = cor.test(response$nodules.rootDM, otus$ZOTU43, method = "pearson")#r = -0.185, p = 0.003
r93 = cor.test(response$nodules.rootDM, otus$ZOTU93, method = "pearson")
r33 = cor.test(response$nodules.rootDM, otus$ZOTU33, method = "pearson")
r34 = cor.test(response$nodules.rootDM, otus$ZOTU34, method = "pearson")
r61 = cor.test(response$nodules.rootDM, otus$ZOTU61, method = "pearson")
r14 = cor.test(response$nodules.rootDM, otus$ZOTU14, method = "pearson")
r49 = cor.test(response$nodules.rootDM, otus$ZOTU49, method = "pearson")
r36 = cor.test(response$nodules.rootDM, otus$ZOTU36, method = "pearson")
r166 = cor.test(response$nodules.rootDM, otus$ZOTU166, method = "pearson")
r164 = cor.test(response$nodules.rootDM, otus$ZOTU164, method = "pearson")
r17 = cor.test(response$nodules.rootDM, otus$ZOTU17, method = "pearson")
r6 = cor.test(response$nodules.rootDM, otus$ZOTU6, method = "pearson")
r23 = cor.test(response$nodules.rootDM, otus$ZOTU23, method = "pearson")
r173 = cor.test(response$nodules.rootDM, otus$ZOTU173, method = "pearson")
r52 = cor.test(response$nodules.rootDM, otus$ZOTU52, method = "pearson")
r155 = cor.test(response$nodules.rootDM, otus$ZOTU155, method = "pearson")
r50 = cor.test(response$nodules.rootDM, otus$ZOTU50, method = "pearson")
r105 = cor.test(response$nodules.rootDM, otus$ZOTU105, method = "pearson")
r70 = cor.test(response$nodules.rootDM, otus$ZOTU70, method = "pearson")
r148 = cor.test(response$nodules.rootDM, otus$ZOTU148, method = "pearson")

##NODULE NUMBER
r19 = cor.test(response$nodule.number, otus$ZOTU19, method = "pearson")                     
r21 = cor.test(response$nodule.number, otus$ZOTU21, method = "pearson")                     
r28 = cor.test(response$nodule.number, otus$ZOTU28, method = "pearson")                     
r40 = cor.test(response$nodule.number, otus$ZOTU40, method = "pearson")                     
r68 = cor.test(response$nodule.number, otus$ZOTU68, method = "pearson")                     
r12 = cor.test(response$nodule.number, otus$ZOTU12, method = "pearson")#r = 0.153, p = 0.014                     
r43 = cor.test(response$nodule.number, otus$ZOTU43, method = "pearson")#r = -0.185, p = 0.003
r93 = cor.test(response$nodule.number, otus$ZOTU93, method = "pearson")
r33 = cor.test(response$nodule.number, otus$ZOTU33, method = "pearson")
r34 = cor.test(response$nodule.number, otus$ZOTU34, method = "pearson")
r61 = cor.test(response$nodule.number, otus$ZOTU61, method = "pearson")
r14 = cor.test(response$nodule.number, otus$ZOTU14, method = "pearson")
r49 = cor.test(response$nodule.number, otus$ZOTU49, method = "pearson")
r36 = cor.test(response$nodule.number, otus$ZOTU36, method = "pearson")
r166 = cor.test(response$nodule.number, otus$ZOTU166, method = "pearson")
r164 = cor.test(response$nodule.number, otus$ZOTU164, method = "pearson")
r17 = cor.test(response$nodule.number, otus$ZOTU17, method = "pearson")
r6 = cor.test(response$nodule.number, otus$ZOTU6, method = "pearson")
r23 = cor.test(response$nodule.number, otus$ZOTU23, method = "pearson")
r173 = cor.test(response$nodule.number, otus$ZOTU173, method = "pearson")
r52 = cor.test(response$nodule.number, otus$ZOTU52, method = "pearson")
r155 = cor.test(response$nodule.number, otus$ZOTU155, method = "pearson")
r50 = cor.test(response$nodule.number, otus$ZOTU50, method = "pearson")
r105 = cor.test(response$nodule.number, otus$ZOTU105, method = "pearson")
r70 = cor.test(response$nodule.number, otus$ZOTU70, method = "pearson")
r148 = cor.test(response$nodule.number, otus$ZOTU148, method = "pearson")
