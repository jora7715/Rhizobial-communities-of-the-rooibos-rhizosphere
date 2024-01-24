### ------------------------
### Setup                    
### ------------------------

## clean/reset environment 
rm(list=ls()) 

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

## Set working directory
setwd("T:/ETH/Data/Pot_experiment/Rhizobia/gyrB_work/R")
wd <- getwd()

list.files(wd)

## Functions
source("Rfunctions.r")

### ------------------------
### Data Import in R
### ------------------------

otufile     <- "p451_run190307_GryB_ZOTU_c99_Count_Sintax.txt"
mapfile     <- "p451_run190307_MapFile.txt"
treefile    <- "p451_run190307_GryB_ZOTU_c99_MSA.tre"

d <- import_qiime(otufilename = otufile, mapfilename = mapfile, treefilename = treefile)
summarize_phyloseq(d)

## Corect problem with taxa names import
rank_names(d)

## Counts per OTU
summary(taxa_sums(d))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#10.0     34.0     67.0   1135.8    166.5 161155.0  

d <- prune_taxa(taxa_sums(d) > 0, d) 

## Counts per Samples
summary(sample_sums(d))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#173    1036    1410    1543    2020    4276 

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
head(sort(sample_sums(d)),50)
low.counts  <- head(sort(sample_sums(d)),12)
low.names   <- attr(low.counts, "names")
low.names
#12 samples with <300 counts!!!! 

# Remove samples with low counts:
all.names   <- sample_names(d)
high.names  <- all.names[!(all.names %in% low.names)]
d.high  <- prune_samples(high.names, d)

D <- d.high 

### Rarefaction curve

## All samples 
ggrare(D, step = 1000, color = "Farm", label = "X.SampleID", se = FALSE)

## Samples with low counts
D.low <- prune_samples(sample_sums(D) < 750, D)
ggrare(D.low, step = 1000, color = "Farm", label = "X.SampleID", se = FALSE)

## Remove samples with counts below 750
D.b2k <- prune_samples(sample_sums(D) > 750, D)
D = D.b2k
#Note: It is important that a sample reaches the plateau otherwise we underestimate the diversity of the sample.

## The bad
badZOTU <- c("ZOTU307","ZOTU215","ZOTU119","ZOTU230","ZOTU301")

## All
allZOTU  <- taxa_names(D)

## The good
goodZOUT <- allZOTU[!(allZOTU %in% badZOTU)]

## Keep only the good ones
D.b2k.goodZOTU <- prune_taxa(goodZOUT, D.b2k)
D = D.b2k.goodZOTU
#Transform OTU abundance
Dt <- transform_sample_counts(D, function(OTU) OTU/sum(OTU))







#Subsetting
pr = Dt
pr = subset_samples(pr, Farm == "mel")
pru = subset_samples(Dt, Fertilization == "unfertilized")
pru = as.data.frame(sample_data(pru))

###NUTRIENT AND BIOMASS ACCUMULATION ANALYSIS###
##Database
pr = as.data.frame(sample_data(Dt))
##Subsetting
#Farm
pr = subset(pr, Farm == "mat")
#Soil
#Fertilization
pr = subset(pr, Fertilization == "unfertilized")
#Block
pr = subset(pr, Block != 1)
pr = subset(pr, Block != 2)
pr = subset(pr, Block == 9)
#Testing (ANOVA + TUKEY'S HSD)
lm = lm(pr$above.DM ~ pr$Soil)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "pr$Soil", conf.level=0.95)
#Boxplot
qplot(x=pru$nodule.number, y=pru$root.DM, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 


##INTERESTING FINDINGS
library(lme4)
library(nlme)

#pr = read.csv("p451_runrun180802_nodA_MapFile.csv", header = T)

setwd("T:/ETH/Data/Pot_experiment/Nutrients")
pr = read.csv("Final_database_ambient.csv", header = T)
pr = pr[,c("Block", "Fertilization", "Soil", "Farm", "Mg.mg.g")]
#, "root.DM", "nodule.number", "cluster.presence"
pr = na.omit(pr)
fit = lme(Mg.mg.g ~ Fertilization*Soil, random= ~1|Block/Farm, data=pr)
lmu = anova(fit)


p.adjust(lmu$`p-value`, method = "bonferroni")
library(wesanderson)
p10 <- ggplot(pr, aes(Fertilization, plant.DM, fill = Soil)) + geom_boxplot() + theme_classic()
p10 = p10 + scale_fill_manual(values=wes_palette(n=3, name="Cavalcanti1"))
p10 + geom_point(size = 3, position=position_dodge(width=0.75),aes(group=Soil, colour = cluster.presence)) + scale_color_grey()

p = qplot(x=root.DM, y=nodule.number, data = pr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
p + theme_bw() + theme_classic() + scale_color_grey()


#Mixed+Melkraal+Unfertilized has higher biomass than cultivation in Melkraal, but wild has the highest
pr = sample_data(Dt)
pr = subset(pr, Block != 1)
pr = subset(pr, Block != 2)
pr = subset(pr, Farm == "mel")
pr = subset(pr, Fertilization == "unfertilized")
lm = lm(pr$plant.DM~ pr$Soil)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "pr$Soil", conf.level=0.95)
plot(pr$plant.DM ~ pr$Fertilization, na.rm = T)                 
#Is it soil or rhizobia? !!!!!ONLY CEC AND SANDINESS EXPLAIN PLANT DM!!!!!
#Plant biomass
pr = sample_data(Dt)
pr = subset(pr, Fertilization == "unfertilized")
soil.model = lm(pr$plant.DM ~ pr$CEC.cmol.kg. + pr$WHC...10kPa. + pr$Coarse.sand... + pr$OrgC.soil + pr$pH + pr$Pot.soil.N.conc + pr$Pot.soil.P.conc + pr$Soil.Ca.mg.g. + pr$Soil.K.mg.g. + pr$Soil.Mg.mg.g. + pr$Soil.Mn.mg.g.)
anova(soil.model)
lmu = lm(pr$plant.DM ~ pr$Soil)
tukey = TukeyHSD(aov(lmu), "pr$Soil", conf.level=0.95)

#N content
pr = sample_data(Dt)
pr = subset(pr, Fertilization == "unfertilized")
soil.model = lm(pr$Leaf.N.content ~ pr$CEC.cmol.kg. + pr$WHC...10kPa. + pr$Coarse.sand... + pr$OrgC.soil + pr$pH + pr$Pot.soil.N.conc + pr$Pot.soil.P.conc + pr$Soil.Ca.mg.g. + pr$Soil.K.mg.g. + pr$Soil.Mg.mg.g. + pr$Soil.Mn.mg.g.)
anova(soil.model)
plot(pr$Leaf.N.content ~ pr$Soil*pr$Farm, na.rm = T) 
#!!!N content Melkkraal!!! NOT DUE TO SOIL N CONCENTRATION!!!
pr = sample_data(Dt)
pr = subset(pr, Farm == "mel")
pr = subset(pr, Fertilization == "fertilized")
lm = lm(pr$Leaf.N.content ~ pr$Soil)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "pr$Soil", conf.level=0.95)
plot(pr$Leaf.N.content ~ pr$Soil, na.rm = T)  
#N concentrations
pr = sample_data(Dt)
pr = subset(pr, Fertilization == "unfertilized")
soil.model = lm(pr$Leaf.Nconc.mg.g. ~ pr$CEC.cmol.kg. + pr$WHC...10kPa. + pr$Coarse.sand... + pr$OrgC.soil + pr$pH + pr$Pot.soil.N.conc + pr$Pot.soil.P.conc + pr$Soil.Ca.mg.g. + pr$Soil.K.mg.g. + pr$Soil.Mg.mg.g. + pr$Soil.Mn.mg.g.)
anova(soil.model)
plot(pr$Leaf.Nconc.mg.g. ~ pr$Soil*pr$Farm, na.rm = T) 
#15N signature NOT USEFUL
pr = sample_data(Dt)
pr = subset(pr, Fertilization == "unfertilized")
soil.model = lm(pr$Leaf.15N.standardized ~ pr$CEC.cmol.kg. + pr$WHC...10kPa. + pr$Coarse.sand... + pr$OrgC.soil + pr$pH + pr$Pot.soil.N.conc + pr$Pot.soil.P.conc + pr$Soil.Ca.mg.g. + pr$Soil.K.mg.g. + pr$Soil.Mg.mg.g. + pr$Soil.Mn.mg.g.)
anova(soil.model)
plot(pr$Leaf.15N.standardized ~ pr$Soil*pr$Farm, na.rm = T) 
lmu=lm(pr$Leaf.15N.standardized ~ pr$Farm, na.rm = T)
tukey = TukeyHSD(lmu, "pr$Soil:pr$Farm", conf.level=0.95)
#P content
pr = sample_data(Dt)
pr = subset(pr, Fertilization == "unfertilized")
soil.model = lm(pr$Leaf.Pcont ~ pr$CEC.cmol.kg. + pr$WHC...10kPa. + pr$Coarse.sand... + pr$OrgC.soil + pr$pH + pr$Pot.soil.N.conc + pr$Pot.soil.P.conc + pr$Soil.Ca.mg.g. + pr$Soil.K.mg.g. + pr$Soil.Mg.mg.g. + pr$Soil.Mn.mg.g.)
anova(soil.model)
plot(pr$Leaf.Pcont ~ pr$Soil*pr$Farm, na.rm = T) 
#!!!P content Melkkraal!!! BUT BECAUSE MORE P IN WILD!!!
pr = sample_data(Dt)
pr = subset(pr, Farm == "mel")
pr = subset(pr, Fertilization == "unfertilized")
lm = lm(pr$Leaf.Pcont ~ pr$Soil)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "pr$Soil", conf.level=0.95)
plot(pr$Leaf.Pcont ~ pr$Soil, na.rm = T)  
##Micronutrients
#Mn concentrations
pr = sample_data(Dt)
pr = subset(pr, Fertilization == "unfertilized")
soil.model = lm(pr$Leaf.Mn..ug.g. ~ pr$CEC.cmol.kg. + pr$WHC...10kPa. + pr$Coarse.sand... + pr$OrgC.soil + pr$pH + pr$Pot.soil.N.conc + pr$Pot.soil.P.conc + pr$Soil.Ca.mg.g. + pr$Soil.K.mg.g. + pr$Soil.Mg.mg.g. + pr$Soil.Mn.mg.g.)
anova(soil.model)
plot(pr$Leaf.Pcont ~ pr$Soil*pr$Farm, na.rm = T) 
#Mn Melkkraal
pr = sample_data(Dt)
pr = subset(pr, Farm == "mel")
pr = subset(pr, Fertilization == "unfertilized")
lm = lm(pr$Leaf.Mn..ug.g. ~ pr$Soil)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "pr$Soil", conf.level=0.95)
plot(pr$Leaf.Mn..ug.g. ~ pr$Soil, na.rm = T)  


###EXPLORATION OF THE RHIZOBIAL COMMUNITIES OF MELKKRAAL WILD AND MIXED
pr = Dt
pr = subset_samples(pr, Farm == "mel")
pru = subset_samples(pr, Fertilization == "unfertilized")
##Ordination !!CAP ordination constrained by plant DM and N content does not help!! => MUST COME FROM SINGLE STRAINS
bdist <- distance(pru, "bray")
wunif.dist <- distance(pru, "unifrac", weighted = T)
ordn=ordinate(pru, "MDS", "bray", weighted = T)
p6=plot_ordination(pru, ordn, color="Soil", label = rownames(sample_data(pru)))
p6 + geom_point(size=4) + ggtitle("NodA communities Melkkraal") 

cap_ord <- ordinate(physeq=pru, method = "CAP", distance =bdist, formula = ~ plant.DM + Leaf.N.content, na.action=na.omit)
cap_plot <- plot_ordination(physeq =pru, ordination = cap_ord, color = "Soil")
p = cap_plot + geom_point(size=5)
p

##PAIRWISE NMDS DIFFERENCES AND PLANT RESPONSES IN MELKKRAAL
#Ordination scores
scores.ord = ordn$vectors[,1:2]
#Plant responses
response = sample_data(pru)
#Euclidean distance matrices
bdist <- distance(pru, "unifrac", weighted = T)
bdist <- distance(pru, "bray")
dist.ordination.PCoA1 = dist(scores.ord[,1], method = "euclidean")
dist.ordination.PCoA2 = dist(scores.ord[,2], method = "euclidean")
dist.totalbiomass = dist(response[,35], method = "euclidean")
plot(dist.ordination.PCoA1 ~ dist.totalbiomass)
plot(dist.ordination.PCoA2 ~ dist.totalbiomass)
mantel(bdist, dist.totalbiomass, method = "pearson", permutations = 9999, na.rm = T)


##CAP analysis
#remove data points with missing metadata
D_not_na <- pru %>% subset_samples(!is.na(plant.DM)&!is.na(nodule.number)&!is.na(Leaf.N.content)&!is.na(Leaf.d15N))
bray_not_na <- phyloseq::distance(physeq = D_not_na, method="bray")
cap_ord <- ordinate(physeq=D_not_na, method = "CAP", distance =bray_not_na, formula = ~ plant.DM + nodule.number + Leaf.N.content + Leaf.d15N)
#CAP plot
cap_plot <- plot_ordination(physeq =D_not_na, ordination = cap_ord, color = "Soil")
p = cap_plot + geom_point(size=5)
p
#plant variables as arrows
arrowmat <- vegan::scores(cap_ord, display ="bp")
#add labels
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
#arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, yend = CAP2,  x = 0, y = 0, shape = NULL, color = NULL,  label=labels)
label_map <- aes(x = 1.3 * CAP1,  y = 1.3 * CAP2, shape = NULL, color = NULL, label= labels)
arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
p + geom_segment(mapping = arrow_map, size =1, data = arrowdf, color = "black", arrow = arrowhead) + 
  geom_text(mapping = label_map, size = 5, data = arrowdf, show.legend = FALSE)



##PERMANOVA
meta = meta(pru)
adonis(distance(pru, method="bray")~Soil, data=meta, permutations = 9999)
#Regression of function and NMDS scores
scores.ord = scores(ordn, choices =c(1,2), display = "sites")
Plant.response = as.data.frame(sample_data(pru))
Plant.response.scores = cbind(scores.ord, Plant.response)
#Plant biomass and N content NOT VERY USEFUL
reg = ggplot(Plant.response.scores, aes(x=NMDS2, y=Leaf.N.content, color = Soil)) +
  geom_point(size=5)  + geom_smooth(method = "lm")
reg

##INDICATOR SPECIES
library(indicspecies)
indicspecies.soil = multipatt(as.data.frame(t(otu_table(pru))), cluster = sample_data(pru)$Soil, control = how(nperm = 9999))
indsoil = summary(indicspecies.soil)
#Indicators cultivated
ind.cult = as.data.frame(c("ZOTU12","ZOTU19","ZOTU148","ZOTU166","ZOTU79"))
ind.wild = as.data.frame(c("ZOTU21","ZOTU162"))
ind.wild.mix =as.data.frame(c("ZOTU149"))
ind.cult.mix = as.data.frame(c("ZOTU170","ZOTU64","ZOTU105","ZOTU86"))
ind.cult.wild = as.data.frame(c("ZOTU99","ZOTU121","ZOTU6")) ##!!!ZOTUS THAT BECOME DILUTED IN THE MIX!!!
all = matrix(c("ZOTU12","ZOTU19","ZOTU148","ZOTU166","ZOTU79","ZOTU21","ZOTU162","ZOTU149",
               "ZOTU170","ZOTU64","ZOTU105","ZOTU86","ZOTU99","ZOTU121","ZOTU6"),ncol=1)
#Regressions ind spp abundance-plant response
OTU_id = as.matrix(rownames(tax_table(pru)))
ind.spp.subset = subset_taxa(pru, OTU_id%in%c("ZOTU12","ZOTU19","ZOTU148","ZOTU166","ZOTU79","ZOTU21","ZOTU162","ZOTU149",
                                              "ZOTU170","ZOTU64","ZOTU105","ZOTU86","ZOTU99","ZOTU121","ZOTU6"))
abds = as.data.frame(t(otu_table(ind.spp.subset)))
resp = as.data.frame(sample_data(ind.spp.subset))
#ZOTU21 INTERESTING
qplot(x=abds$ZOTU21, y=resp$plant.DM, color=resp$Soil, geom="point") + geom_point(size=5, na.action=na.omit) + geom_smooth(method="lm", se = T) 
#ZOTU99 INTERESTING
qplot(x=abds$ZOTU99, y=resp$plant.DM, color=resp$Soil, geom="point") + geom_point(size=5, na.action=na.omit) + geom_smooth(method="lm", se = T) 


##DOMINANT OTU ANALYSIS
####Find the most abundant taxa per plant####
#Goes through a phyloseq object, picks out the most abundant taxa and gives the abundance for each
#and identifies which taxa is most abundant for which sample
# Function to find the top most abundant ASV/OTU per sample, number of ASV/OTU per sample depends on user input. This is particularly relevant for finding strain level community structure when a few genera dominates the communities, for  example "is it a single variant of Pseudomonas dominating all the samples?"
find.top.asv <- function(x,num){
  x <- pru
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
test = find.top.asv(pru,3)
write.table(as.matrix(test), "top_abundant_taxa_Melkkraal.txt")
toptax = read.csv("top_abundant_taxa_Melkkraal.csv", header = T, sep=",", row.names = 1) 
# toptax = as.data.frame(t(toptax))

#Get top 3 per plant
tax123 = toptax[,94:96] 
sample_data(pru) = cbind(sample_data(pru), tax123)
dominant.OTU <- subset(otu_table(pru), rownames(otu_table(pru)) %in% c("ZOTU7", "ZOTU36", "ZOTU27", "ZOTU33","ZOTU18", "ZOTU86","ZOTU22","ZOTU97", "ZOTU21","ZOTU42","ZOTU99","ZOTU23","ZOTU43","ZOTU64","ZOTU105","ZOTU12"))
dom.physeq <- merge_phyloseq(dominant.OTU, tax_table(pru), sample_data(pru), phy_tree(pru))

#Recalculate evenness for only dominant strains
even <- evenness(dom.physeq, "all")
pielou = even$pielou
sample_data(dom.physeq) = as.data.frame(cbind(sample_data(dom.physeq), pielou))

##Ordination
#NMDS soils
ord <- ordinate(dom.physeq, "NMDS", "bray")
ord <- ordinate(dom.physeq, "MDS", "unifrac", weighted = T)
plot_ordination(dom.physeq, ord, color = "Soil") +
  geom_point(size = 5)
##PERMANOVA
meta = meta(dom.physeq)
adonis(distance(dom.physeq, method="bray") ~ Soil, data = meta, permutations = 9999)

#Indicator species
indicspecies.soil = multipatt(as.data.frame(t(otu_table(dom.physeq))), cluster = sample_data(dom.physeq)$Soil, control = how(nperm = 9999))
indsoil = summary(indicspecies.soil)

###Correlation of ZOTU relative abundance to plant response using top abundant ZOTUs###
##Remember that tests are not possible due to non-normality of the abundances
library(ggpubr)
OTU_id = as.matrix(rownames(tax_table(dom.physeq)))
species_indicator = OTU_id[,1]
#Prepare dbb
df1 <- data.frame(sample_data(dom.physeq))
df2 <- data.frame(otu_table(dom.physeq))
df2.t <- t(df2) #transpose df2
df_corr <- cbind(df1, df2.t) #combine by column (r001...)
#Correlation plots with indicspp
#ZOTU21
qplot(x=ZOTU21, y=plant.DM, color=Soil, data=df_corr, geom="point") + geom_point(size=5) + geom_smooth(method="lm", se = T) 
#ZOTU27 IS THE MOST DOMINANT IN THE MIX
qplot(x=ZOTU27, y=plant.DM, color=Soil, data=df_corr, geom="point") + geom_point(size=5) 
ggplot(x=ZOTU27, y=Soil, color=Soil, data=df_corr) + 
##EVENNESS
plot(df_corr$pielou ~ df_corr$Soil, na.rm = T)



##PCA to determine integrated plant response to the different treatments using only Unfertilized treatment
library(ggfortify)
library(cluster)
#Subset
Dtu = subset_samples(Dt, Fertilization == "unfertilized")
response = as.data.frame(sample_data(Dt))
#Plot
resp.DM = na.omit(response[,c(4,5,27,28,29,30,31)])
PC<-prcomp(resp.DM[,3:7])
autoplot(PC, data = resp.DM, colour = "Farm", shape = "Soil", size = 5, loadings = TRUE, loadings.label = TRUE, 
        loadings.label.colour = "black")




##RICHNES-PLANT DM UNFERTILIZED
Dtu = subset_samples(Dt, Fertilization == "unfertilized")
rich <- richness(Dt)
en = evenness(Dt)
response = as.data.frame(sample_data(Dt))
response = as.data.frame(cbind(response, rich, en))

qplot(x=`0`, y=nodule.number, color=Soil, data=response, geom="point") + geom_point(size=5) + geom_smooth(method = "lm")
lmu = lm(response$`0` ~ response$nodule.number)
anova(lmu)
cor(response$nodule.number, response$`0`)

aa = as.data.frame(sample_data(Dt))
ggplot(aa, aes(x=Fertilization, y=nodule.number)) + geom_boxplot() + theme_bw() + theme_classic()


##CLUSTER ROOT COSTS AND RESPONSES
cl.root = sample_data(Dtu)
ggplot(cl.root, aes(x=Fertilization, y=cluster.number)) + geom_boxplot(notch=TRUE)
lmu = lm(cl.root$cluster.number ~ cl.root$Soil)
anova(lmu)

qplot(x=root.DM, y=cluster.number, color=Soil, data=cl.root, geom="point") + geom_point(size=5) + geom_smooth(method = "lm")
lmu = lm(response$cluster.number ~ response$root.DM)
anova(lmu)
cor(response$cluster.number, response$root.DM)
rcorr = cor.test(response$cluster.number, response$root.DM, method = "pearson")


###ABUNDANCE OF NON-MESORHIZOBIUM OTUS GYRB IN LARGE AND SMALL PLANTS
Dt <- transform_sample_counts(D, function(OTU) OTU/sum(OTU))
#Dt = filter_taxa(Dt, function(x) sum(x) > .01, TRUE)

#Dt = subset_samples(Dt, Fertilization=="fertilized")
ldt = subset_samples(Dt, plant.DM.class=="L")
ldm = subset_samples(Dt, plant.DM.class=="M")
ldh = subset_samples(Dt, plant.DM.class=="H")
ldt = subset_taxa(ldt, Genus != "Mesorhizobium")
ldm = subset_taxa(ldm, Genus != "Mesorhizobium")
ldh = subset_taxa(ldh, Genus != "Mesorhizobium")

#Low
sum.other.L = as.data.frame(sample_sums(ldt))
meso.L = 1-sum.other.L[,1]
sum.other.L$Mesorhizobium = meso.L
sum.other.L$others = sum.other.L$`sample_sums(ldt)`
sum.other.L = sum.other.L[,2:3]
Size = c(rep("A", 38)) #ALL
sum.other.L$size = Size 
root.dry.matter = as.data.frame(sample_data(ldt)$root.DM)
sum.other.L$root.dry.matter = root.dry.matter
nodulation = as.data.frame(sample_data(ldt)$nodule.number)
sum.other.L$nodulation = nodulation
Fertilization = as.data.frame(sample_data(ldt)$Fertilization)
sum.other.L$Fertilization = Fertilization

#Medium
sum.other.M = as.data.frame(sample_sums(ldm))
meso.M = 1-sum.other.M[,1]
sum.other.M$Mesorhizobium = meso.M
sum.other.M$others = sum.other.M$`sample_sums(ldm)`
sum.other.M = sum.other.M[,2:3]
Size = c(rep("B", 154)) #ALL
sum.other.M$size = Size 
root.dry.matter = as.data.frame(sample_data(ldm)$root.DM)
sum.other.M$root.dry.matter = root.dry.matter
nodulation = as.data.frame(sample_data(ldm)$nodule.number)
sum.other.M$nodulation = nodulation
Fertilization = as.data.frame(sample_data(ldm)$Fertilization)
sum.other.M$Fertilization = Fertilization


#High
sum.other.H = as.data.frame(sample_sums(ldh))
meso.H = 1-sum.other.H[,1]
sum.other.H$Mesorhizobium = meso.H
sum.other.H$others = sum.other.H$`sample_sums(ldh)`
sum.other.H = sum.other.H[,2:3]
Size = c(rep("C", 26)) #ALL AND FERT
sum.other.H$size = Size
root.dry.matter = as.data.frame(sample_data(ldh)$root.DM)
sum.other.H$root.dry.matter = root.dry.matter
nodulation = as.data.frame(sample_data(ldh)$nodule.number)
sum.other.H$nodulation = nodulation
Fertilization = as.data.frame(sample_data(ldh)$Fertilization)
sum.other.H$Fertilization = Fertilization


#aa = data.frame(t(sum.other.L), t(sum.other.M), t(sum.other.H))
aa = data.frame(t(sum.other.L), t(sum.other.M), t(sum.other.H))
all.taxa = as.data.frame(t(aa), stringsAsFactors = F)
all.taxa$Mesorhizobium = as.numeric(as.character(all.taxa$Mesorhizobium))
all.taxa$others = as.numeric(as.character(all.taxa$others))
all.taxa$root.dry.matter = as.numeric(as.character(all.taxa$root.dry.matter))
all.taxa$nodulation = as.numeric(as.character(all.taxa$nodulation))

#Plant dry matter
ggplot(all.taxa, aes(x=size, y=Mesorhizobium, color = size)) + 
  geom_boxplot()

p = qplot(x=root.dry.matter, y=others, data = all.taxa, geom="point", colour = Fertilization) + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
p + theme_bw() + theme_classic() + scale_color_grey()
rcorr = cor.test(all.taxa$root.dry.matter, all.taxa$others, method = "pearson")

p = qplot(x=nodulation, y=others, data = all.taxa, geom="point", colour = Fertilization) + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
p + theme_bw() + theme_classic() + scale_color_grey()

rcorr = cor.test(all.taxa$nodulation, all.taxa$others, method = "pearson")


lm = lm(all.taxa$others ~ all.taxa$plant.dry.matter)
anova(lm)
#By size groups
lm = lm(all.taxa$others ~ all.taxa$size)
lmu = anova(lm)
tukey = TukeyHSD(lmu, "all.taxa$size", conf.level=0.95)

#Nodule number
ggplot(all.taxa, aes(x=size, y=nodulation, color = size)) + 
  geom_boxplot()
lm = lm(all.taxa$nodulation ~ all.taxa$size)
lmu = aov(lm)
tukey = TukeyHSD(lmu, "all.taxa$size", conf.level=0.95)
qplot(x=all.taxa$nodulation, y=all.taxa$others, geom="point") + geom_point(size=3) + geom_smooth(method="lm", se = T, colour = 1)
cor(all.taxa$nodulation, all.taxa$others)
lm = lm(all.taxa$others ~ all.taxa$nodulation)
anova(lm)

###ALPHA DIVERSITY LINKS TO NON-MESORHIZOBIUM ABUNDANCE
#write.csv(as.data.frame(all.taxa), "NonMesoAlphaDatabase.csv")
alpha.non = read.csv("NonMesoAlphaDatabase.csv", header = T)

p = qplot(x=others, y=Richness, data = alpha.non, geom="point", colour = Fertilization) + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
p + theme_bw() + theme_classic() + scale_color_grey()
rcorr = cor.test(alpha.non$others, alpha.non$Richness, method = "pearson")


##NON-MESO RICHNESS

Dt <- transform_sample_counts(D, function(OTU) OTU/sum(OTU))
Dt = subset_taxa(Dt, Genus != "Mesorhizobium")

gyrB_richness = richness(Dt)
gyrB_richness = cbind(sample_data(Dt), gyrB_richness)

gyrB_simpson = evenness(Dt, index = "all")
gyrB_simpson = cbind(sample_data(Dt), gyrB_simpson)

write.csv(gyrB_simpson, "gyrB_simpson_nonmeso.csv")
write.csv(gyrB_richness, "gyrB_richness_nonmeso.csv")

alpha.non = read.csv("NonMesoAlphaDatabase.csv", header = T)

p = qplot(x=NonMeso_richness, y=Richness, data = alpha.non, geom="point", colour = Fertilization) + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
p + theme_bw() + theme_classic() + scale_color_grey()
rcorr = cor.test(alpha.non$NonMeso_richness, alpha.non$Richness, method = "pearson")


##MESORHIZOBIUM RICHNESS
Dt <- transform_sample_counts(D, function(OTU) OTU/sum(OTU))
Dt = subset_taxa(Dt, Genus == "Mesorhizobium")

meso_richness = richness(Dt)
meso_richness = cbind(sample_data(Dt), meso_richness)

meso_simpson = evenness(Dt, index = "all")
meso_simpson = cbind(sample_data(Dt), meso_simpson)

write.csv(meso_simpson, "gyrB_simpson_meso.csv")
write.csv(meso_richness, "gyrB_richness_meso.csv")

alpha.meso = read.csv("MesoVsNonMesoAlphaDatabase.csv", header = T)

p = qplot(x=Meso_richness, y=Richness, data = alpha.meso, geom="point", colour = Fertilization) + geom_point(size=5) + geom_smooth(method="lm", se = T, colour = 1)
p + theme_bw() + theme_classic() + scale_color_grey()
rcorr = cor.test(alpha.meso$Meso_richness, alpha.meso$Richness, method = "pearson")



##TESTING THE ROOIBOS NUTRITIONAL STATUS ACCORDING TO ROOT SYSTEM BIOMASS CLASS WITHIN THE FERTILIZATION TREATMENT
##REFERENCE VARIABLES: LEAF N CONC, LEAF P CONC, LEAF NP RATIO, LEAF 15N
library(lme4)
library(nlme)

Dt = subset_samples(Dt,Fertilization == "fertilized")
pr = meta(Dt)
pr = pr[,c("Block", "Farm", "root.DM", "Root.15N")]
pr = na.omit(pr)
fit = lme(Leaf.d15N ~ root.DM, random= ~1|Block/Farm, data=pr)
lmu = anova(fit)

fit.tuk = lm(pr$Leaf.d15N ~ pr$plant.DM.class)
lmu.tuk = aov(fit.tuk)
tukey = TukeyHSD(lmu.tuk, "pr$plant.DM.class", conf.level=0.95)
pr$plant.DM.class<-factor(pr$plant.DM.class, levels=c("L", "M", "H"))
boxplot(Leaf.d15N ~ plant.DM.class, data=pr)

##RICHNESS AND EVENNESS AGAINST PLANT NUTRITIONAL STATUS
rich.gyrB = read.csv("richness_iNEXT_gyrB_medianendpoint_OK.csv", sep = ",", header = T)
simp.gyrB = read.csv("simpson_iNEXT_gyrB_medianendpoint_new.csv", sep = ",", header = T)
#alpha.div = as.data.frame(cbind(rich.gyrB, simp.gyrB$Simpson))
rich.gyrB$Richness = rich.gyrB$qD
simp.gyrB$Simpson = simp.gyrB$qD

rich.gyrB = subset(rich.gyrB, Fertilization == "fertilized")
p = qplot(x=Root.15N, y=root.DM, data = pr, geom="point") + geom_point(size=5) 
p + geom_smooth(method="lm", se = T, colour = 1)+ theme_bw() + theme_classic() + scale_color_grey()
 
rcorr = cor.test(pr$Leaf.d15N, pr$root.DM, method = "pearson")




