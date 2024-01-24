# Data Analysis
#Project: p451 run180802
#User(s): Stefanie Stadelmann (stadelms@student.ethz.ch)
#Data type: AmpSeq PE300 NodA (acyltransferase)
#Date   : 23.08.2018

rm(list=ls()) # clean/reset environment 
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

#################### Summary and exploration ####################

summarize_phyloseq(d)
# => Sparsity = 0.69 

# Number of counts per Taxa
sort(taxa_sums(d))

# Number of counts per Sample
sort(sample_sums(d))

# Read-Count Distribution - Careful this are still raw counts!
plot(sample_sums(d), xaxt= "n", xlab="Sample", ylab="Number of Reads", pch=19, cex=c(sample_sums(d)/25000), col=rgb(0,0,1,alpha=0.5), main="Read-Counts per Sample", ylim=range(0,max(sample_sums(d))*1.1))
ac <- rep("blue",length(sample_sums(d)))
ac[get_variable(d, "SID") == "m"] <- "green"
points(sample_sums(d), pch=3, col=ac, cex=2)
v <- sample_names(d)
axis(side = 1, at = seq(1,length(sample_sums(d))), labels = v, tck=-.02, las=2, col.axis = 1)

#pdf("Raw_Counts_per_sample.pdf", paper="a4")
#png("Raw_Counts_per_sample.png", width = 700, height = 350, units = "px")

#dev.off()

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
D <- d.high.S 

### Remove problematic OTUs
## Show tree
plot_tree(D, label.tips="taxa_names", "treeonly", nodeplotblank, ladderize=TRUE)
## The bad

badZOTU <- c("ZOTU152","ZOTU143","ZOTU133","ZOTU119","ZOTU124","ZOTU156","ZOTU160","ZOTU134","ZOTU44","ZOTU153","ZOTU128","ZOTU104","ZOTU127","ZOTU127","ZOTU103","ZOTU151","ZOTU157","ZOTU169","ZOTU146","ZOTU144")
## All
allZOTU  <- taxa_names(D)
## The good
goodZOUT <- allZOTU[!(allZOTU %in% badZOTU)]
## Keep only the good ones
D.b2k.goodZOTU <- prune_taxa(goodZOUT, D)
plot_tree(D.b2k.goodZOTU , label.tips="taxa_names", "treeonly", nodeplotblank, ladderize=TRUE)
D <- D.b2k.goodZOTU

##---------------------------------------------------------
## Intra-Extrapolation
##----------------------------------------------------------
library(metagMisc)
library(dada2)
library(ALDEx2)
library(metagenomeSeq)
library(DESeq2)
library(iNEXT)
## Phyloseq intra-extrapolation
#phyloseq_inext(D, Q = 0, curve_type = "diversity",
#              correct_singletons = FALSE, endpoint = NULL, knots = 10,
#              multithread = FALSE, show_CI = TRUE, show_sample_labels = FALSE,
#              show_plot = TRUE, justDF = FALSE, add_raw_data = TRUE)
#phyloseq_inext(D, Q = 0, curve_type = "coverage",
 #             correct_singletons = FALSE, endpoint = NULL, knots = 10,
  #            multithread = FALSE, show_CI = TRUE, show_sample_labels = FALSE,
   #           show_plot = TRUE, justDF = FALSE, add_raw_data = TRUE)
nodA_richness = phyloseq_inext(D, Q = 0, justDF = TRUE, conf = 0.95, endpoint = median(sample_sums(D)))
#nodA_shannon = phyloseq_inext(D, Q = 1, justDF = TRUE, conf = 0.95)
nodA_simpson = phyloseq_inext(D, Q = 2, justDF = TRUE, conf = 0.95, endpoint = median(sample_sums(D)))
write.csv(nodA_richness, "richness_iNEXT_nodA_medianendpoint.csv")
#write.csv(nodA_shannon, "shannon_iNEXT_nodA.csv")
write.csv(nodA_simpson, "simpson_iNEXT_nodA_medianendpoint.csv")


#Standardise the diversity estimates to a community coverage of 0.99:
#otus = as.data.frame(otu_table(D))
#AA = estimateD(otus, datatype="abundance", base="coverage", level=0.99, conf=0.95)
#AAA = phyloseq_inext(D, Q = 0, justDF = TRUE)
#write.csv(AAA, "inter-extrapol_nodA.csv")
##Intra Extrapolation iNext
#library(iNEXT)
#otus = as.data.frame(t(as.data.frame(otu_table(D))))
#samples = as.data.frame(sample_data(D))
#input = as.data.frame(t(as.data.frame(cbind(otus, samples))))
#input = as.data.frame(cbind(otus, samples))
#write.csv(input, "xxx.csv")
#m = c(200, 1200, 3000, 6000, 10000, 30000)
#out <- iNEXT(as.data.frame(t(input)), q=c(0), datatype="abundance", size=m, knots = 10, conf = 0.95, nboot = 50)
#ggiNEXT(out, se=TRUE, facet.var="none", color.var = "none")  

### Rarefaction curve
## All samples 
ggrare(D, step = 1000, color = "Farm", label = "X.SampleID", se = FALSE)

## Samples with low counts
D.low <- prune_samples(sample_sums(D) < 500, D)
ggrare(D.low, step = 500, color = "Farm", label = "X.SampleID", se = FALSE)

# > Samples below 2000 counts are problematic because diversity is underestimated. 

## Remove samples with counts below 3000
D.b2k <- prune_samples(sample_sums(D) > 500, D)
D = D.b2k






##----------------------------------------------------------
### Alpha diversity after inter- extrapolation
##---------------------------------------------------------
setwd("T:/ETH/Data/Pot_experiment/Rhizobia/nodA_work/R")
rich.nodA = read.csv("richness_iNEXT_nodA_medianendpoint.csv", sep = ",", header = T)
simp.nodA = read.csv("simpson_iNEXT_nodA_medianendpoint.csv", sep = ",", header = T)
#alpha.div = as.data.frame(cbind(rich.nodA, simp.nodA$Simpson))

library(lme4)
library(nlme)
library(car)

rich.nodA = subset(rich.nodA, Fertilization == "fertilized")
blo = subset(rich.nodA, Farm == "blo")
plot(blo$Richness ~ blo$Soil)
fit = lme(Simpson ~ Farm, random= ~1|Block, data=as.data.frame(simp.nodA))
anova(fit)
dob = subset(rich.nodA, Farm == "dob")
plot(dob$Richness ~ dob$Soil)
fit = lme(Richness ~ Soil, random= ~1|Block, data=dob)
anova(fit)
lan = subset(rich.nodA, Farm == "lan")
plot(lan$Richness ~ lan$Soil)
fit = lme(Richness ~ Soil, random= ~1|Block, data=lan)
anova(fit)
mat = subset(rich.nodA, Farm == "mat")
plot(mat$Richness ~ mat$Soil)
fit = lme(Richness ~ Soil, random= ~1|Block, data=mat)
anova(fit)
mel = subset(rich.nodA, Farm == "mel")
plot(mel$Richness ~ mel$Soil)
fit = lme(Richness ~ Soil, random= ~1|Block, data=mel)
anova(fit)


rich.nodA = rich.nodA[,c("Block", "Farm", "Richness")]
rich.nodA = na.omit(rich.nodA)
fit = lme(Richness ~ Farm, random= ~1|Block, data=rich.nodA)
anova(fit)

p = qplot(x=Farm, y=Richness, data = rich.nodA, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
p + theme_bw() + theme_classic() + scale_color_grey()
rcorr = cor.test(rich.nodA$root.DM, rich.nodA$Richness, method = "pearson")

p = qplot(x=nodule.number, y=Richness, data = rich.nodA, geom="point", colour = Fertilization) + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
p + theme_bw() + theme_classic() + scale_color_grey()
rcorr = cor.test(rich.nodA$nodule.number, rich.nodA$Richness, method = "pearson")

simp.nodA = simp.nodA[,c("Block", "Fertilization", "Soil", "Farm", "Simpson")]
simp.nodA = na.omit(simp.nodA)
fit = lme(Simpson ~ Farm*Fertilization*Soil, random= ~1|Block, data=simp.nodA)
anova(fit)

p = qplot(x=root.DM, y=Simpson, data = simp.nodA, geom="point", colour = Fertilization) + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
p + theme_bw() + theme_classic() + scale_color_grey()
rcorr = cor.test(simp.nodA$nodules.rootDM, simp.nodA$Simpson, method = "pearson")

p = qplot(x=nodule.number, y=Simpson, data = simp.nodA, geom="point", colour = Fertilization) + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
p + theme_bw() + theme_classic() + scale_color_grey()

rcorr = cor.test(simp.nodA$root.DM, simp.nodA$Simpson, method = "pearson")









##RICHNESS AND EVENNESS AGAINST PLANT NUTRITIONAL STATUS
rich.gyrB = read.csv("richness_iNEXT_gyrB_medianendpoint_OK.csv", sep = ",", header = T)
simp.gyrB = read.csv("simpson_iNEXT_gyrB_medianendpoint_new.csv", sep = ",", header = T)
#alpha.div = as.data.frame(cbind(rich.gyrB, simp.gyrB$Simpson))
rich.gyrB$Richness = rich.gyrB$qD
simp.gyrB$Simpson = simp.gyrB$qD

rich.gyrB = subset(rich.gyrB, Fertilization == "fertilized")
p = qplot(x=N.P.ratio, y=Richness, data = rich.gyrB, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
p + theme_bw() + theme_classic() + scale_color_grey()
rcorr = cor.test(rich.gyrB$root.DM, rich.gyrB$Richness, method = "pearson")


















## Alpha diversity
alpha <- diversities(D, index = "all")
kable(head(alpha))

## Dominance
domin <- dominance(D, index = "all")
kable(head(domin))

## Richness
rich <- richness(D)
kable(head(rich))

## Evenness
even <- evenness(D, "all")
kable(head(even))

## Meta data
meta  <- meta(D)

## add index to meta data 
meta$Shannon        <- alpha$shannon 
meta$InverseSimpson <- alpha$inverse_simpson
meta$DominanceSimpson <- domin$simpson
#$richness <- rich

#diversity measure
p = plot_richness(D, x="Soil", color="Farm", measures=c("Chao1", "Shannon"))
p + geom_point(size=5, alpha=0.7)
p


##----------------------------------
#relative OTU abundance
##-------------------------------------
Dt <- transform_sample_counts(D, function(OTU) OTU/sum(OTU))

# square root transformation (Hellinger)

sqrt.Dt <- transform_sample_counts(Dt, function(OTU) sqrt(OTU))

##----------------------------------------------------------------------------------
##Richness
##-------------------------------------------------------------------------------------
metaD <- meta(D)
alphaD= estimate_richness(D, split=TRUE, measures = NULL)

#put the two dataframes together
df.alphameta <- cbind(metaD, alphaD)
df.alphameta = cbind(df.alphameta, df.alphameta$Observed/df.alphameta$nodule.number)
df.alphameta$richness.nodule = df.alphameta$`df.alphameta$Observed/df.alphameta$nodule.number`

library(ggplot2)
#Boxplot
p10 <- ggplot(df.alphameta, aes(Soil, Observed)) + geom_boxplot()
p = p10 + labs(y= "Observed Species Richness", title = "OTU Richness for Soiltype") + theme_classic()
p

p10 <- ggplot(df.alphameta, aes(Fertilization, Observed)) + geom_boxplot()
p = p10 + labs(y= "Observed Species Richness", title = "OTU Richness for Soiltype") + theme_classic()
p

p10 <- ggplot(df.alphameta, aes(Fertilization, richness.nodule)) + geom_boxplot()
p = p10 + labs(y= "Observed Species Richness", title = "OTU Richness for Soiltype") + theme_classic()
p

p10 <- ggplot(df.alphameta, aes(Soil, richness.nodule)) + geom_boxplot()
p = p10 + labs(y= "Observed Species Richness", title = "OTU Richness for Soiltype") + theme_classic()
p



#Farm
lm = lm(df.alphameta$Observed ~ df.alphameta$Farm)
lmu = aov(lm)
anova(lmu)
Tuk = TukeyHSD(lmu, "df.alphameta$Farm", conf.level=0.95)
#Soil
lm = lm(df.alphameta$Observed ~ df.alphameta$Soil)
lmu = aov(lm)
anova(lmu)
Tuk = TukeyHSD(lmu, "df.alphameta$Soil", conf.level=0.95)
#Fertilization
lm = lm(df.alphameta$Observed ~ df.alphameta$Fertilization)
lmu = aov(lm)
anova(lmu)
Tuk = TukeyHSD(lmu, "df.alphameta$Fertilization", conf.level=0.95)
#Soil x Fertilization
lm = lm(df.alphameta$Observed ~ df.alphameta$Soil*df.alphameta$Fertilization)
lmu = aov(lm)
anova(lmu)
Tuk = TukeyHSD(lmu, "df.alphameta$Soil:df.alphameta$Fertilization", conf.level=0.95)

p20 <- ggplot(df.alphameta, aes(Soil, Observed, colour=Farm)) + geom_boxplot()
p = p20 + labs(y= "Observed Species Richness", title = "OTU Richness for Soiltype") + theme_classic()
p



#for each Farm
blo=subset_samples(D, Farm=="blo")
dob=subset_samples(D, Farm=="dob")
mel=subset_samples(D, Farm=="mel")
lan=subset_samples(D, Farm=="lan")
mat=subset_samples(D, Farm=="mat")

#estimated richness per Farm
alpha.blo= estimate_richness(blo, split=TRUE)
alpha.dob= estimate_richness(dob, split=TRUE)
alpha.mel= estimate_richness(mel, split=TRUE)
alpha.lan= estimate_richness(lan, split=TRUE)
alpha.mat= estimate_richness(mat, split=TRUE)

#plots
plot_richness(D, x = "Soil", color = "Fertilization", title = "Alpha diversity (Richness)", measures="Observed")

#blo
plot_richness(blo, x = "Soil", color = "Fertilization", title = "Richness Blo", measures="Observed")
#dob
plot_richness(dob, x = "Soil", color = "Fertilization", title = "Richness Dob", measures="Observed")
#mel
plot_richness(mel, x = "Soil", color = "Fertilization", title = "Richness Mel", measures="Observed")
#lan
plot_richness(lan, x = "Soil", color = "Fertilization", title = "Richness Lan", measures="Observed")
#mat
plot_richness(mat, x = "Soil", color = "Fertilization", title = "Richness Mat", measures="Observed")








#Fertilized
D.fert=subset_samples(D, Fertilization=="fertilized")
metaDfert <- meta(D.fert)
D.fert.alpha=estimate_richness(D.fert)
D.fert.alpha = cbind(metaDfert, D.fert.alpha)

#Unfertilized
D.unfert=subset_samples(D, Fertilization=="unfertilized")
metaDunfert <- meta(D.unfert)
D.unfert.alpha =estimate_richness(D.unfert)
D.unfert.alpha= cbind(metaDunfert, D.unfert.alpha)

#Fertilized vs. Unfertilized Boxplot
p11 <- ggplot(D.fert.alpha, aes(Soil, Observed, colour=Farm)) + geom_boxplot()
p = p11 + labs(y= "Observed Species Richness", title = "OTU Richness Fertilized") + theme_classic()
p

p12 <- ggplot(D.unfert.alpha, aes(Soil, Observed, colour=Farm)) + geom_boxplot()
p = p12 + labs(y= "Observed Species Richness", title = "OTU Richness Unfertilized") + theme_classic()
p

#table
library(dplyr)
df.alphameta %>% 
  group_by(Farm) %>% 
  summarise(Observed = mean(Observed))

df.alphameta %>% 
  group_by(Farm, Soil) %>%                            # multiple group columns
  summarise(mean_Observed = mean(Observed), sd(Observed))  # multiple summary columns

df.alphameta %>%
  group_by(Fertilization) %>%
  summarise(mean_Observed = mean(Observed), sd(Observed))

df.alphameta %>%
  group_by(Soil) %>%
  summarise(mean_Observed = mean(Observed), sd(Observed))

obs_soil <- subset.data.frame(df.alphameta, select=c(Observed, Soil))
write.table(obs_soil, file = "observed_soil.csv", sep = ",", col.names = NA,
            qmethod = "double")
cult <- subset.data.frame(obs_soil, Soil = cultivated)  


##----------------------------------------------------------
##Pairwise comparisons
##----------------------------------------------------------
###farms
farm <- levels(meta$Farm) # get the variables

## Make pairwise lists that we want to compare
farm.pairs <- combn(seq_along(farm), 2, simplify = FALSE, FUN = function(i)farm[i])
print(farm.pairs)

## Choose colors
n.color <- length(unique(farm))
mycolor <- brewer.pal(n = n.color, name = "RdBu")

## Violin plot with non-parametric test (Wilcoxon test)
#png("AlphaDiversity.png", width = 700, height = 350, units = "px")
p1 <- ggviolin(meta, x = "Farm", y = "Shannon", add = "boxplot",
               fill = "Farm", palette = mycolor)
p1 <- p1 + stat_compare_means(comparisons = farm.pairs) 
p1
#dev.off()
###Soil
soil <- levels(meta$Soil) # get the variables

## make a pairwise list that we want to compare
soil.pairs <- combn(seq_along(soil), 2, simplify = FALSE, FUN = function(i)soil[i])
print(soil.pairs)

## Choose colors
n.color <- length(unique(soil))
mycolor <-brewer.pal(n = n.color, name = "RdBu")

## Violin plot with non-parametric test (Wilcoxon test)
#png("AlphaDiversity.png", width = 700, height = 350, units = "px")
p1 <- ggviolin(meta, x = "Soil", y = "Shannon", add = "boxplot",
               fill = "Soil", palette = mycolor)
p1 <- p1 + stat_compare_means(comparisons = soil.pairs) 
p1
#dev.off()
## Violin plot with non-parametric test (Wilcoxon test)
#png("AlphaDiversity.png", width = 700, height = 350, units = "px")
p1 <- ggviolin(meta, x = "Fertilization", y = "Shannon", add = "boxplot",
               fill = "Fertilization", palette = mycolor)
p1 <- p1 + stat_compare_means(comparisons = fertilization.pairs)
print(p1)
#dev.off()
###Fertilization
fertilization <- levels(meta$Fertilization) # get the variables

## make a pairwise list that we want to compare
fertilization.pairs <- combn(seq_along(fertilization), 2, simplify = FALSE, FUN = function(i)fertilization[i])
print(fertilization.pairs)

## Choose colors
n.color <- length(unique(fertilization))
mycolor <- heat.colors(n=n.color, alpha = 1)


