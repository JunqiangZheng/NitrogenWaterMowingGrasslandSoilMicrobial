rm(list=ls())
# setwd("~/Dropbox/FAST/Additional file 4/R_Input/")
# output <- "~/Dropbox/FAST/Additional file 4/R_Output/"
# setwd("~/GitHub/MicrobiomeStatPlot/252NetworkModule/40168_2017_389_MOESM4_ESM/R_Input")
setwd("E:/工作资料集/数据/多伦 草原站/郑俊强老师 YF156-M201708043/MPL 20170976 and 77 32个样 （多因子一 GCB）/发表用/network analysis/Hartman et al Microbiome_N/Figure_occu_03/R_input/")
output <- "E:/工作资料集/数据/多伦 草原站/郑俊强老师 YF156-M201708043/MPL 20170976 and 77 32个样 （多因子一 GCB）/发表用/network analysis/Hartman et al Microbiome_N/Figure_occu_03/R_output/"
###################################
##### prepare packages #####
###################################
library(BiodiversityR) #package ???tcltk??? is not available (for R version 3.2.2)
library(ggplot2)
library(vegan)
library(reshape2)
library(Hmisc)
library(plotrix)
library(phyloseq)
library(MASS)
# BiocManager::install("bioDist")
library(bioDist) #package ???bioDist??? is not available (for R version 3.2.2)
library(igraph)
library(car)
library(coin)
library(edgeR)
library(formatR)
library(gridExtra)
library(gplots)
library(indicspecies)
library(sciplot)
library(ape)
library(grid)
library(RVAideMemoire)
library(gridBase)
library(TukeyC)
library(corrplot)
library(userfriendlyscience)
library(caret)
library(multcompView)
source("vennDia.R")
source("CorrDF.R")
source("plotOTU.R")
source("cor.mtest.R")
source("variance_functions.R")
source("maPalette.R")
source("triangle_shape.R")
source("star_shape.R")
debuggingState(on=FALSE)
options(scipen=10)

###################################
##### Import and prepare data #####
###################################

#######################
##### 16S #####
#######################

##### Import Data #####
otu_16s <- read.csv("bacterial_16s_otu.csv", row.names=1,sep=",", header=T, blank.lines.skip=F, check.names=F)
otu_16s <- as.matrix(otu_16s)
rownames(otu_16s) <- paste("b", rownames(otu_16s), sep="") # otu name 闁??? add "b"

##### Import design file #####
design_16s <- read.table("map.txt", header=T, row.names=1, stringsAsFactors=F, na.strings="NA")
design_16s$Clipping <- factor(design_16s$Clipping,c("CL","uC"))
design_16s$Fertilizer <-factor(design_16s$Fertilizer,c("aN","eN"))
design_16s$Precipitation <- factor(design_16s$Precipitation,c("aP","eP"))
design_16s$Treatments <- factor(design_16s$Treatments,c("uCaNaP","uCaNeP","uCeNaP", "uCeNeP","CLaNaP","CLaNeP","CLeNaP","CLeNeP"))
design_16s$CL_N <- factor(design_16s$CL_N ,c("uCaN","uCeN", "CLaN","CLeN"))
design_16s$CL_P <- factor(design_16s$CL_P,c("uCaP","uCeP","CLaP","CLeP"))
design_16s$N_P <- factor(design_16s$N_P, c("aNaP", "aNeP", "eNaP", "eNeP"))
str(design_16s)

##### Import Taxonomy #####
tax_16s <- read.table("bacterial_tax.txt", row.names=1, sep="\t", header=F ,stringsAsFactors=F,quote="")
rownames(tax_16s) <- paste("b",rownames(tax_16s),sep="")
tax_16s[tax_16s==""] <- "unassigned"
tax_16s[tax_16s=="uncultured_bacterium"] <- "unassigned"
tax_16s[tax_16s=="unidentified"] <- "unassigned"
colnames(tax_16s) <- c("Kingdom","Phylum","Class","Order", "Family", "Genus", "Species","C")

# create separate taxonomy label specifying classes of Proteobacteria
tax_16s$labels <- tax_16s$Phylum
tax_16s[ rownames(tax_16s)[tax_16s$Class=="Alphaproteobacteria" ], ]$labels <- "Alphaproteobacteria"
tax_16s[ rownames(tax_16s)[tax_16s$Class=="Betaproteobacteria" ], ]$labels <- "Betaproteobacteria"
tax_16s[ rownames(tax_16s)[tax_16s$Class=="Gammaproteobacteria" ], ]$labels <- "Gammaproteobacteria"
tax_16s[ rownames(tax_16s)[tax_16s$Class=="Deltaproteobacteria" ], ]$labels <- "Deltaproteobacteria"
table(tax_16s$labels)


##### Store Cyanobacteria and mitrochrondia sequeces #####
unique(tax_16s$Kingdom)
table(tax_16s$Kingdom)

table(tax_16s$Kingdom)
r1 <- rownames(tax_16s[tax_16s$Kingdom=="Archaea",])

unique(tax_16s$Kingdom)
r2 <- c(rownames(tax_16s[tax_16s$Kingdom=="No_blast_hit",]))

table(tax_16s$Phylum)
r3 <- rownames(tax_16s[tax_16s$Phylum=="Cyanobacteria",])

unique(tax_16s$Family)
r4 <- c(rownames(tax_16s[tax_16s$Family=="mitochondria",]))

otus_remove_16s <- c(r1,r2, r3, r4)

## Remove these from otu table, tax table
otu_filter_16s <- otu_16s[-which(rownames(otu_16s) %in% otus_remove_16s),]
tax_filter_16s <- tax_16s[rownames(otu_filter_16s),]
design_filter_16s <- droplevels(design_16s[rownames(design_16s) %in% colnames(otu_filter_16s),])
design_filter_16s <- design_filter_16s[colnames(otu_filter_16s),]

dim(otu_filter_16s)
dim(tax_filter_16s)
dim(design_filter_16s)

# check the removed results
unique(tax_filter_16s$Kingdom)
table(tax_filter_16s$Phylum)
unique(tax_filter_16s$Family)

###### 16S sequence and OTU counts ######
sum(colSums(otu_filter_16s))
sort(colSums(otu_filter_16s))
median(colSums(otu_filter_16s))

nrow(tax_filter_16s)
table(tax_filter_16s$Kingdom)

## Export filtered OTU table for multiple rarefactions in QIIME
write.table(otu_filter_16s,paste(output,"otu_filter_16s.txt",sep=""),sep="\t",quote=F)

## Order taxonmy file by OTU
otu_order_16s <- match(rownames(otu_filter_16s), rownames(tax_filter_16s))
tax_filter_16s <- tax_filter_16s[otu_order_16s,]

##### Calculate % sequences removed #####
archa_otu <- otu_16s[r1,]; dim(archa_otu); length(r1)
sum(colSums(archa_otu))/sum(colSums(otu_16s))*100

noblast_otu <- otu_16s[r2,]; dim(noblast_otu); length(r2)
sum(sum(noblast_otu))/sum(colSums(otu_16s))*100

cyano_otu <- otu_16s[r3,]; dim(cyano_otu); length(r3)
sum(colSums(cyano_otu))/sum(colSums(otu_16s))*100

mito_otu <- otu_16s[r4,]; dim(mito_otu); length(r4)
sum(sum(mito_otu))/sum(colSums(otu_16s))*100




##### Define sample Fertilizers #####
aN_samples <- rownames(design_filter_16s)[which(design_filter_16s$Fertilizer == "aN")]
eN_samples <- rownames(design_filter_16s)[which(design_filter_16s$Fertilizer == "eN")]

##### Total number of aN / eN OTUs #####
aN_16s_otu <- otu_filter_16s[,aN_samples]
aN_16s_otu <- aN_16s_otu[rowSums(aN_16s_otu) > 0,]
nrow(aN_16s_otu)

eN_16s_otu <- otu_filter_16s[,eN_samples]
eN_16s_otu <- eN_16s_otu[rowSums(eN_16s_otu) > 0,]
nrow(eN_16s_otu)

##### Define sample types #####
aNsamples <-  rownames(design_filter_16s)[which(design_filter_16s$Fertilizer == "aN")]
eNsamples <- rownames(design_filter_16s)[which(design_filter_16s$Fertilizer == "eN")]

#######################
##### ITS #####
#######################

#### Import OTU table #####
otu_its <- read.csv("Fungi_otu.csv",row.names=1,sep=",",header=T,
                      blank.lines.skip=F,check.names=F)
otu_its <- as.matrix(otu_its)
rownames(otu_its) <- paste("f",rownames(otu_its),sep="")


##### Import design file #####

design_its <- read.table("map.txt", header=T, row.names=1, stringsAsFactors=F, na.strings="NA")
design_its$Clipping <- factor(design_its$Clipping,c("CL","uC"))
design_its$Fertilizer <-factor(design_its$Fertilizer,c("aN","eN"))
design_its$Precipitation <- factor(design_its$Precipitation,c("aP","eP"))
design_its$Treatments <- factor(design_its$Treatments,c("uCaNaP","uCaNeP","uCeNaP", "uCeNeP","CLaNaP","CLaNeP","CLeNaP","CLeNeP"))
design_its$CL_N <- factor(design_its$CL_N ,c("uCaN","uCeN", "CLaN","CLeN"))
design_its$CL_P <- factor(design_its$CL_P,c("uCaP","uCeP","CLaP","CLeP"))
design_its$N_P <- factor(design_its$N_P, c("aNaP", "aNeP", "eNaP", "eNeP"))
str(design_its)

##### Import Taxonomy #####
tax_its <- read.csv("Fungi_tax.csv",row.names=1, sep=",", header=F,stringsAsFactors=F,quote="")
rownames(tax_its) <- paste("f",rownames(tax_its),sep="")
colnames(tax_its) <- c("Kingdom","Phylum","Class","Order", "Family", "Genus", "Species","C")
tax_its[tax_its==""] <- "unassigned"

##### Store Plantae, Protists and sequences of unknown origin and remove #####
#OTUs classified in kingdom Plantae or Protista are removed, as well as those whose kingdom is unassigned

unique(tax_its$Kingdom)
table(tax_its$Kingdom)
r3 <- rownames(tax_its[tax_its$Kingdom == "No blast hit",])
r4 <- rownames(tax_its[tax_its$Kingdom == "Protista",])
r5 <- rownames(tax_its[tax_its$Kingdom == "unidentified",])
r6 <- rownames(tax_its[tax_its$Phylum=="unassigned",])
r7 <- rownames(tax_its[tax_its$Phylum=="unidentified",])
r8 <- rownames(tax_its[tax_its$Phylum=="Cercozoa",])
r9 <- rownames(tax_its[tax_its$Phylum=="Ciliophora",])
r10 <- rownames(tax_its[tax_its$Phylum=="Incertae sedis",])
#r6 <- rownames(tax_its[tax_its$Phylum == "unassigned",])

otus_remove_its <- c(r3,r4,r5,r6,r7,r8,r9,r10)

## Remove these from otu table, tax table
otu_filter_its <- otu_its[-which(rownames(otu_its) %in% otus_remove_its),]
tax_filter_its <- tax_its[rownames(otu_filter_its),]
design_filter_its <- droplevels(design_its[rownames(design_its) %in% colnames(otu_filter_its),])
design_filter_its <- design_filter_its[colnames(otu_filter_its),]

unique(tax_filter_its$Kingdom)
unique(tax_filter_its$Phylum)
table(tax_filter_its$Phylum)


##### Calculate % sequences removed #####
noblast_otu <- otu_its[r3,]; dim(noblast_otu); length(r3)
sum(colSums(noblast_otu))/sum(colSums(otu_its))*100

protist_otu <- otu_its[r4,]; dim(protist_otu); length(r4)
sum(colSums(protist_otu))/sum(colSums(otu_its))*100

unassign_otu <- otu_its[rownames(otu_its) %in% r5,]; dim(unassign_otu); length(r5)
sum(colSums(unassign_otu))/sum(colSums(otu_its))*100

unassigned_otu <- otu_its[rownames(otu_its) %in% r6,]; dim(unassigned_otu); length(r6)
sum(colSums(unassigned_otu))/sum(colSums(otu_its))*100

unidentified_otu <- otu_its[rownames(otu_its) %in% r7,]; dim(unidentified_otu); length(r7)
sum(colSums(unidentified_otu))/sum(colSums(otu_its))*100

Cercozoa_otu <- otu_its[rownames(otu_its) %in% r8,]; dim(Cercozoa_otu); length(r8)
sum(colSums(Cercozoa_otu))/sum(colSums(otu_its))*100

Ciliophora_otu <- otu_its[rownames(otu_its) %in% r9,]; dim(Ciliophora_otu); length(r9)
sum(colSums(Ciliophora_otu))/sum(colSums(otu_its))*100

Incertae_sedis_otu <- otu_its[rownames(otu_its) %in% r10,]; dim(Incertae_sedis_otu); length(r10)
sum(colSums(Incertae_sedis_otu))/sum(colSums(otu_its))*100

dim(otu_filter_its)
# Number of OTUs
dim(tax_filter_its)
# Number of samples
dim(design_filter_its)

###### ITS sequence and OTU counts ######
sum(colSums(otu_filter_its))
sort(colSums(otu_filter_its))
# sort(colSums(otu_its))
median(colSums(otu_filter_its))

nrow(tax_filter_its)
table(tax_filter_its$Kingdom)

## Export filtered OTU table for multiple rarefactions in QIIME
write.table(otu_filter_its,paste(output,"otu_filter_its.txt",sep=""),sep="\t",quote=F)

## Order taxonmy file by OTU
otu_order_its <- match(rownames(otu_filter_its), rownames(tax_filter_its))
tax_filter_its <- tax_filter_its[otu_order_its, ]

##### Define sample Fertilizers #####
aN_its_samples <- rownames(design_filter_its)[which(design_filter_its$Fertilizer == "aN")]
eN_its_samples <- rownames(design_filter_its)[which(design_filter_its$Fertilizer == "eN")]

##### Total number of aN / eN OTUs #####
aN_its_otu <- otu_filter_its[,aN_its_samples]
aN_its_otu <- aN_its_otu[rowSums(aN_its_otu) > 0,]
nrow(aN_its_otu)

eN_its_otu <- otu_filter_its[,eN_its_samples]
eN_its_otu <- eN_its_otu[rowSums(eN_its_otu) > 0,]
nrow(eN_its_otu)

##### Define sample types #####
aNsamples <- rownames(design_filter_its)[which(design_filter_its$Fertilizer == "aN")]
eNsamples <- rownames(design_filter_its)[which(design_filter_its$Fertilizer == "eN")]

#######################
##### 18S #####
#######################

#### Import OTU table #####
otu_18s <- read.csv("Protist_otu.csv",row.names=1,sep=",",header=T,
                    blank.lines.skip=F,check.names=F)
otu_18s <- as.matrix(otu_18s)
rownames(otu_18s) <- paste("p",rownames(otu_18s),sep="")


##### Import design file #####

design_18s <- read.table("map.txt", header=T, row.names=1, stringsAsFactors=F, na.strings="NA")
design_18s$Clipping <- factor(design_18s$Clipping,c("uC","CL"))
design_18s$Fertilizer <-factor(design_18s$Fertilizer,c("aN","eN"))
design_18s$Precipitation <- factor(design_18s$Precipitation,c("aP","eP"))
design_18s$Treatments <- factor(design_18s$Treatments,c("uCaNaP","uCaNeP","uCeNaP", "uCeNeP","CLaNaP","CLaNeP","CLeNaP","CLeNeP"))
design_18s$CL_N <- factor(design_18s$CL_N ,c("uCaN","uCeN", "CLaN","CLeN"))
design_18s$CL_P <- factor(design_18s$CL_P,c("uCaP","uCeP","CLaP","CLeP"))
design_18s$N_P <- factor(design_18s$N_P, c("aNaP", "aNeP", "eNaP", "eNeP"))
str(design_18s)

##### Import Taxonomy #####
tax_18s <- read.csv("Protist_tax.csv",row.names=1, sep=",", header=F,stringsAsFactors=F,quote="")
rownames(tax_18s) <- paste("p",rownames(tax_18s),sep="")
colnames(tax_18s) <- c("Kingdom","Phylum","Class","Order", "Family", "Genus", "Species","C")
tax_18s[tax_18s==""] <- "unassigned"

##### Store Plantae, Protists and sequences of unknown origin and remove #####
#OTUs classified in kingdom Plantae or Protista are removed, as well as those whose kingdom is unassigned

unique(tax_18s$Kingdom)
table(tax_18s$Kingdom)
#r3 <- rownames(tax_18s[tax_18s$Kingdom == "No blast hit",])
#r4 <- rownames(tax_18s[tax_18s$Kingdom == "fungi",])
#r5 <- rownames(tax_18s[tax_18s$Kingdom == "unidentified",])

r6 <- rownames(tax_18s[tax_18s$Phylum == "Metazoa",])

otus_remove_18s <- c(r6)

## Remove these from otu table, tax table
otu_filter_18s <- otu_18s[-which(rownames(otu_18s) %in% otus_remove_18s),]
tax_filter_18s <- tax_18s[rownames(otu_filter_18s),]
design_filter_18s <- droplevels(design_18s[rownames(design_18s) %in% colnames(otu_filter_18s),])
design_filter_18s <- design_filter_18s[colnames(otu_filter_18s),]

##### Calculate % sequences removed #####
#plant_otu <- otu_18s[r3,]; dim(plant_otu); length(r3)
#sum(colSums(plant_otu))/sum(colSums(otu_18s))*100

#protist_otu <- otu_18s[r4,]; dim(protist_otu); length(r4)
#sum(colSums(protist_otu))/sum(colSums(otu_18s))*100

#unassign_otu <- otu_18s[rownames(otu_18s) %in% r5,]; dim(unassign_otu); length(r5)
#sum(colSums(unassign_otu))/sum(colSums(otu_18s))*100

Metazoa_otu <- otu_18s[rownames(otu_18s) %in% r6,]; dim(Metazoa_otu); length(r6)
sum(colSums(Metazoa_otu))/sum(colSums(otu_18s))*100

dim(otu_filter_18s)
# Number of OTUs
dim(tax_filter_18s)
# Number of samples
dim(design_filter_18s)

###### 18s sequence and OTU counts ######
sum(colSums(otu_filter_18s))
sort(colSums(otu_filter_18s))
median(colSums(otu_filter_18s))

nrow(tax_filter_18s)
table(tax_filter_18s$Kingdom)

## Export filtered OTU table for multiple rarefactions in QIIME
write.table(otu_filter_18s,paste(output,"otu_filter_18s.txt",sep=""),sep="\t",quote=F)

## Order taxonmy file by OTU
otu_order_18s <- match(rownames(otu_filter_18s), rownames(tax_filter_18s))
tax_filter_18s <- tax_filter_18s[otu_order_18s, ]

##### Define sample Fertilizers #####
aN_18s_samples <-  rownames(design_filter_18s)[which(design_filter_18s$Fertilizer == "aN")]
eN_18s_samples <- rownames(design_filter_18s)[which(design_filter_18s$Fertilizer == "eN")]

##### Total number of aN / eN OTUs #####
aN_18s_otu <- otu_filter_18s[,aN_18s_samples]
aN_18s_otu <- aN_18s_otu[rowSums(aN_18s_otu) > 0,]
nrow(aN_18s_otu)

eN_18s_otu <- otu_filter_18s[,eN_18s_samples]
eN_18s_otu <- eN_18s_otu[rowSums(eN_18s_otu) > 0,]
nrow(eN_18s_otu)

##### Define sample types #####
aNsamples <-  rownames(design_filter_18s)[which(design_filter_18s$Fertilizer == "aN")]
eNsamples <- rownames(design_filter_18s)[which(design_filter_18s$Fertilizer == "eN")]
############==============================================================================
############==============================================================================






##### Defining  aN bacteria and fungi communities for beta diversity analysis #####

## Apply sequence count threshold to bacteria aN community and TMM normalize counts
otu_16s_aN <- otu_filter_16s[, aNsamples ]
otu_16s_aN <- otu_16s_aN[rowSums(otu_16s_aN) > 0,]
dim(otu_16s_aN)
keep_OTUs_aN_16s <- which(rowSums(otu_16s_aN >= 8) >= 4)

otu_16s_aN <- otu_16s_aN[keep_OTUs_aN_16s,]

dim(otu_16s_aN)

tax_aN_16s <- tax_filter_16s[rownames(otu_16s_aN),]
design_16s_aN <- droplevels(design_filter_16s[aNsamples,])

edgeR_16s_aN <- DGEList(counts=otu_16s_aN, 
                          group=design_16s_aN$CL_P,
                          genes=tax_aN_16s)

edgeR_16s_aN <- calcNormFactors(edgeR_16s_aN)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_aN_16s <- cpm(edgeR_16s_aN, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk aN bacteria community into phyloseq objects
## for further analysis
phy_16s_aN <- otu_table(otu_norm_aN_16s,taxa_are_rows=T)
phy_tax_aN_16s <-tax_table(as.matrix(tax_aN_16s))
phy_design_aN_16s <- sample_data(design_16s_aN)
physeq_aN_norm_16s <- phyloseq(phy_16s_aN,phy_tax_aN_16s,phy_design_aN_16s)
sample_data(physeq_aN_norm_16s)$CL_P <- factor(sample_data(physeq_aN_norm_16s)$CL_P,levels=c("CLaP", "CLeP","uCaP","uCeP"))


## Apply sequence count threshold to fungi aN community and TMM normalize counts
otu_its_aN <- otu_filter_its[,aNsamples]
otu_its_aN <- otu_its_aN[rowSums(otu_its_aN) > 0,]

keep_OTUs_aN_its <- which(rowSums(otu_its_aN >= 8) >= 4)

otu_its_aN <- otu_its_aN[keep_OTUs_aN_its,]

dim(otu_its_aN)

tax_aN_its <- tax_filter_its[rownames(otu_its_aN),]
design_its_aN <- droplevels(design_filter_its[aNsamples,])

edgeR_its_aN <- DGEList(counts=otu_its_aN, 
                          group=design_its_aN$CL_P,
                          genes=tax_aN_its)

edgeR_its_aN <- calcNormFactors(edgeR_its_aN)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_aN_its <- cpm(edgeR_its_aN,normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk aN bacteria community into phyloseq objects
## for further analysis
phy_its_aN <- otu_table(otu_norm_aN_its,taxa_are_rows=T)
phy_tax_aN_its <-tax_table(as.matrix(tax_aN_its))
phy_design_aN_its <- sample_data(design_its_aN)
physeq_aN_norm_its <- phyloseq(phy_its_aN,phy_tax_aN_its,phy_design_aN_its)
sample_data(physeq_aN_norm_its)$CL_P <- factor(sample_data(physeq_aN_norm_its)$CL_P,levels=c("CLaP", "CLeP","uCaP","uCeP"))

## Apply sequence count threshold to Protist aN community and TMM normalize counts
otu_18s_aN <- otu_filter_18s[, aNsamples ]
otu_18s_aN <- otu_18s_aN[rowSums(otu_18s_aN) > 0,]
dim(otu_18s_aN)
keep_OTUs_aN_18s <- which(rowSums(otu_18s_aN >= 8) >= 4)

otu_18s_aN <- otu_18s_aN[keep_OTUs_aN_18s,]

dim(otu_18s_aN)

tax_aN_18s <- tax_filter_18s[rownames(otu_18s_aN),]
design_18s_aN <- droplevels(design_filter_18s[aNsamples,])

edgeR_18s_aN <- DGEList(counts=otu_18s_aN, 
                          group=design_18s_aN$CL_P,
                          genes=tax_aN_18s)

edgeR_18s_aN <- calcNormFactors(edgeR_18s_aN)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_aN_18s <- cpm(edgeR_18s_aN, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk aN bacteria community into phyloseq objects
## for further analysis
phy_18s_aN <- otu_table(otu_norm_aN_18s,taxa_are_rows=T)
phy_tax_aN_18s <-tax_table(as.matrix(tax_aN_18s))
phy_design_aN_18s <- sample_data(design_18s_aN)
physeq_aN_norm_18s <- phyloseq(phy_18s_aN,phy_tax_aN_18s,phy_design_aN_18s)
sample_data(physeq_aN_norm_18s)$CL_P <- factor(sample_data(physeq_aN_norm_18s)$CL_P,levels=c("CLaP", "CLeP","uCaP","uCeP"))


##### Defining  eN bacteria and fungi communities for beta diversity analysis #####

## Apply sequence count threshold to bacteria eN community and TMM normalize counts
otu_16s_eN <- otu_filter_16s[, eNsamples ]
otu_16s_eN <- otu_16s_eN[rowSums(otu_16s_eN) > 0,]
dim(otu_16s_eN)
keep_OTUs_eN_16s <- which(rowSums(otu_16s_eN >= 8) >= 4)

otu_16s_eN <- otu_16s_eN[keep_OTUs_eN_16s,]

dim(otu_16s_eN)

tax_eN_16s <- tax_filter_16s[rownames(otu_16s_eN),]
design_16s_eN <- droplevels(design_filter_16s[eNsamples,])

edgeR_16s_eN <- DGEList(counts=otu_16s_eN, 
                          group=design_16s_eN$CL_P,
                          genes=tax_eN_16s)

edgeR_16s_eN <- calcNormFactors(edgeR_16s_eN)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_eN_16s <- cpm(edgeR_16s_eN, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk eN bacteria community into phyloseq objects
## for further analysis
phy_16s_eN <- otu_table(otu_norm_eN_16s,taxa_are_rows=T)
phy_tax_eN_16s <-tax_table(as.matrix(tax_eN_16s))
phy_design_eN_16s <- sample_data(design_16s_eN)
physeq_eN_norm_16s <- phyloseq(phy_16s_eN,phy_tax_eN_16s,phy_design_eN_16s)
sample_data(physeq_eN_norm_16s)$CL_P <- factor(sample_data(physeq_eN_norm_16s)$CL_P,levels=c("CLaP", "CLeP","uCaP","uCeP"))


## Apply sequence count threshold to fungi eN community and TMM normalize counts
otu_its_eN <- otu_filter_its[,eNsamples]
otu_its_eN <- otu_its_eN[rowSums(otu_its_eN) > 0,]

keep_OTUs_eN_its <- which(rowSums(otu_its_eN >= 8) >= 4)

otu_its_eN <- otu_its_eN[keep_OTUs_eN_its,]

dim(otu_its_eN)

tax_eN_its <- tax_filter_its[rownames(otu_its_eN),]
design_its_eN <- droplevels(design_filter_its[eNsamples,])

edgeR_its_eN <- DGEList(counts=otu_its_eN, 
                          group=design_its_eN$CL_P,
                          genes=tax_eN_its)

edgeR_its_eN <- calcNormFactors(edgeR_its_eN)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_eN_its <- cpm(edgeR_its_eN,normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk eN bacteria community into phyloseq objects
## for further analysis
phy_its_eN <- otu_table(otu_norm_eN_its,taxa_are_rows=T)
phy_tax_eN_its <-tax_table(as.matrix(tax_eN_its))
phy_design_eN_its <- sample_data(design_its_eN)
physeq_eN_norm_its <- phyloseq(phy_its_eN,phy_tax_eN_its,phy_design_eN_its)
sample_data(physeq_eN_norm_its)$CL_P <- factor(sample_data(physeq_eN_norm_its)$CL_P,levels=c("CLaP", "CLeP","uCaP","uCeP"))

## Apply sequence count threshold to Protist eN community and TMM normalize counts
otu_18s_eN <- otu_filter_18s[, eNsamples ]
otu_18s_eN <- otu_18s_eN[rowSums(otu_18s_eN) > 0,]
dim(otu_18s_eN)
keep_OTUs_eN_18s <- which(rowSums(otu_18s_eN >= 8) >= 4)

otu_18s_eN <- otu_18s_eN[keep_OTUs_eN_18s,]

dim(otu_18s_eN)

tax_eN_18s <- tax_filter_18s[rownames(otu_18s_eN),]
design_18s_eN <- droplevels(design_filter_18s[eNsamples,])

edgeR_18s_eN <- DGEList(counts=otu_18s_eN, 
                          group=design_18s_eN$CL_P,
                          genes=tax_eN_18s)

edgeR_18s_eN <- calcNormFactors(edgeR_18s_eN)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_eN_18s <- cpm(edgeR_18s_eN, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk eN bacteria community into phyloseq objects
## for further analysis
phy_18s_eN <- otu_table(otu_norm_eN_18s,taxa_are_rows=T)
phy_tax_eN_18s <-tax_table(as.matrix(tax_eN_18s))
phy_design_eN_18s <- sample_data(design_18s_eN)
physeq_eN_norm_18s <- phyloseq(phy_18s_eN,phy_tax_eN_18s,phy_design_eN_18s)
sample_data(physeq_eN_norm_18s)$CL_P <- factor(sample_data(physeq_eN_norm_18s)$CL_P,levels=c("CLaP", "CLeP","uCaP","uCeP"))



##### Identifiying cropping system responsive OTUs with indicator species analysis #####

## Identify indicator species in bulk aN bacteria communities
indic_aN_16s <- as.data.frame(t(otu_norm_aN_16s))
indic_aN_groups_16s <- design_16s_aN$CL_P
length(unique(indic_aN_groups_16s))

## Define indicator species for aN bacteria community.
## Note: These calculations can be time and processor intensive
set.seed(8046)
indicatorsp_aN_16s <- multipatt(indic_aN_16s,indic_aN_groups_16s,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_aN_16s,alpha=1,indvalcomp=T)
indic_aN_df_16s <- indicatorsp_aN_16s$sign
write.table(indic_aN_df_16s,paste0(output,"indicsp_aN_16s.txt"),sep="\t",quote=F)

## Import data frame of indicator species to save time
indic_aN_df_16s <- read.table("indicsp_aN_16s.txt", header=T, sep="\t")

uCaP_aN_16s <- as.matrix(indic_aN_df_16s[which(indic_aN_df_16s$s.uCaP == 1 & indic_aN_df_16s$p.value < 0.05),])
uCeP_aN_16s <- as.matrix(indic_aN_df_16s[which(indic_aN_df_16s$s.uCeP == 1 & indic_aN_df_16s$p.value < 0.05),])
CLaP_aN_16s <- as.matrix(indic_aN_df_16s[which(indic_aN_df_16s$s.CLaP == 1 & indic_aN_df_16s$p.value < 0.05),])
CLeP_aN_16s <- as.matrix(indic_aN_df_16s[which(indic_aN_df_16s$s.CLeP == 1 & indic_aN_df_16s$p.value < 0.05),])

aN_r_values_16s <- rbind(CLaP_aN_16s,CLeP_aN_16s,uCaP_aN_16s,uCeP_aN_16s)
colnames(aN_r_values_16s)[1:4] <-c("CLaP","CLeP","uCaP","uCeP")

## Range of correlation coefficients
range(aN_r_values_16s[,"stat"])

## Total number of indicator OTUS
length(unique(rownames(aN_r_values_16s)))

## Proportion of aN bacteria OTUs responding to cropping system
length(unique(rownames(aN_r_values_16s)))/nrow(otu_norm_aN_16s)

## Proportion of aN bacteria sequences responding to cropping system
aN_16s_ra <- t(t(otu_16s_aN)/colSums(otu_16s_aN)) * 100
sum(colSums(aN_16s_ra[unique(rownames(aN_r_values_16s)),]))/sum(colSums(aN_16s_ra))

## Identify indicator species in bulk aN bacteria communities
indic_aN_its <- as.data.frame(t(otu_norm_aN_its))
indic_aN_groups_its <- design_its_aN$CL_P
length(unique(indic_aN_groups_its))

## Define indicator species for aN bacteria community.
## Note: These calculations can be time and processor intensive
set.seed(8046)
indicatorsp_aN_its <- multipatt(indic_aN_its,indic_aN_groups_its,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_aN_its,alpha=1,indvalcomp=T)
indic_aN_df_its <- indicatorsp_aN_its$sign
write.table(indic_aN_df_its,paste0(output,"indicsp_aN_its.txt"),sep="\t",quote=F)

## Import data frame of indicator species to save time
indic_aN_df_its <- read.table("indicsp_aN_its.txt", header=T, sep="\t")

CLaP_aN_its <- as.matrix(indic_aN_df_its[which(indic_aN_df_its$s.CLaP == 1 & indic_aN_df_its$p.value < 0.05),])
CLeP_aN_its <- as.matrix(indic_aN_df_its[which(indic_aN_df_its$s.CLeP == 1 & indic_aN_df_its$p.value < 0.05),])
uCaP_aN_its <- as.matrix(indic_aN_df_its[which(indic_aN_df_its$s.uCaP == 1 & indic_aN_df_its$p.value < 0.05),])
uCeP_aN_its <- as.matrix(indic_aN_df_its[which(indic_aN_df_its$s.uCeP == 1 & indic_aN_df_its$p.value < 0.05),])

aN_r_values_its <- rbind(CLaP_aN_its,CLeP_aN_its,uCaP_aN_its,uCeP_aN_its)
colnames(aN_r_values_its)[1:4] <-c("CLaP","CLeP","uCaP","uCeP")

## Range of correlation coefficients
range(aN_r_values_its[,"stat"])

## Total number of indicator OTUS
length(unique(rownames(aN_r_values_its)))

## Proportion of aN bacteria OTUs responding to cropping system
length(unique(rownames(aN_r_values_its)))/nrow(otu_norm_aN_its)

## Proportion of aN bacteria sequences responding to cropping system
aN_its_ra <- t(t(otu_its_aN)/colSums(otu_its_aN)) * 100
sum(colSums(aN_its_ra[unique(rownames(aN_r_values_its)),]))/sum(colSums(aN_its_ra))


## Identify indicator species in bulk aN bacteria communities
indic_aN_18s <- as.data.frame(t(otu_norm_aN_18s))
indic_aN_groups_18s <- design_18s_aN$CL_P
length(unique(indic_aN_groups_18s))

## Define indicator species for aN bacteria community.
## Note: These calculations can be time and processor intensive
set.seed(8046)
indicatorsp_aN_18s <- multipatt(indic_aN_18s,indic_aN_groups_18s,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_aN_18s,alpha=1,indvalcomp=T)
indic_aN_df_18s <- indicatorsp_aN_18s$sign
write.table(indic_aN_df_18s,paste0(output,"indicsp_aN_18s.txt"),sep="\t",quote=F)

## Import data frame of indicator species to save time
indic_aN_df_18s <- read.table("indicsp_aN_18s.txt", header=T, sep="\t")

CLaP_aN_18s <- as.matrix(indic_aN_df_18s[which(indic_aN_df_18s$s.CLaP == 1 & indic_aN_df_18s$p.value < 0.05),])
CLeP_aN_18s <- as.matrix(indic_aN_df_18s[which(indic_aN_df_18s$s.CLeP == 1 & indic_aN_df_18s$p.value < 0.05),])
uCaP_aN_18s <- as.matrix(indic_aN_df_18s[which(indic_aN_df_18s$s.uCaP == 1 & indic_aN_df_18s$p.value < 0.05),])
uCeP_aN_18s <- as.matrix(indic_aN_df_18s[which(indic_aN_df_18s$s.uCeP == 1 & indic_aN_df_18s$p.value < 0.05),])

aN_r_values_18s <- rbind(CLaP_aN_18s,CLeP_aN_18s,uCaP_aN_18s,uCeP_aN_18s)
colnames(aN_r_values_18s)[1:4] <-c("CLaP","CLeP","uCaP","uCeP")

## Range of correlation coefficients
range(aN_r_values_18s[,"stat"])

## Total number of indicator OTUS
length(unique(rownames(aN_r_values_18s)))

## Proportion of aN bacteria OTUs responding to cropping system
length(unique(rownames(aN_r_values_18s)))/nrow(otu_norm_aN_18s)

## Proportion of aN bacteria sequences responding to cropping system
aN_18s_ra <- t(t(otu_18s_aN)/colSums(otu_18s_aN)) * 100
sum(colSums(aN_18s_ra[unique(rownames(aN_r_values_18s)),]))/sum(colSums(aN_18s_ra))


## Identify indicator species in bulk eN bacteria communities
indic_eN_16s <- as.data.frame(t(otu_norm_eN_16s))
indic_eN_groups_16s <- design_16s_eN$CL_P
length(unique(indic_eN_groups_16s))

## Define indicator species for eN bacteria community.
## Note: These calculations can be time and processor intensive
set.seed(8046)
indicatorsp_eN_16s <- multipatt(indic_eN_16s,indic_eN_groups_16s,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_eN_16s,alpha=1,indvalcomp=T)
indic_eN_df_16s <- indicatorsp_eN_16s$sign
write.table(indic_eN_df_16s,paste0(output,"indicsp_eN_16s.txt"),sep="\t",quote=F)

## Import data frame of indicator species to save time
indic_eN_df_16s <- read.table("indicsp_eN_16s.txt", header=T, sep="\t")

CLaP_eN_16s <- as.matrix(indic_eN_df_16s[which(indic_eN_df_16s$s.CLaP == 1 & indic_eN_df_16s$p.value < 0.05),])
CLeP_eN_16s <- as.matrix(indic_eN_df_16s[which(indic_eN_df_16s$s.CLeP == 1 & indic_eN_df_16s$p.value < 0.05),])
uCaP_eN_16s <- as.matrix(indic_eN_df_16s[which(indic_eN_df_16s$s.uCaP == 1 & indic_eN_df_16s$p.value < 0.05),])
uCeP_eN_16s <- as.matrix(indic_eN_df_16s[which(indic_eN_df_16s$s.uCeP == 1 & indic_eN_df_16s$p.value < 0.05),])

eN_r_values_16s <- rbind(CLaP_eN_16s,CLeP_eN_16s,uCaP_eN_16s,uCeP_eN_16s)
colnames(eN_r_values_16s)[1:4] <-c("CLaP","CLeP","uCaP","uCeP")

## Range of correlation coefficients
range(eN_r_values_16s[,"stat"])

## Total number of indicator OTUS
length(unique(rownames(eN_r_values_16s)))

## Proportion of eN bacteria OTUs responding to cropping system
length(unique(rownames(eN_r_values_16s)))/nrow(otu_norm_eN_16s)

## Proportion of eN bacteria sequences responding to cropping system
eN_16s_ra <- t(t(otu_16s_eN)/colSums(otu_16s_eN)) * 100
sum(colSums(eN_16s_ra[unique(rownames(eN_r_values_16s)),]))/sum(colSums(eN_16s_ra))

## Identify indicator species in bulk eN bacteria communities
indic_eN_its <- as.data.frame(t(otu_norm_eN_its))
indic_eN_groups_its <- design_its_eN$CL_P
length(unique(indic_eN_groups_its))

## Define indicator species for eN bacteria community.
## Note: These calculations can be time and processor intensive
set.seed(8046)
indicatorsp_eN_its <- multipatt(indic_eN_its,indic_eN_groups_its,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_eN_its,alpha=1,indvalcomp=T)
indic_eN_df_its <- indicatorsp_eN_its$sign
write.table(indic_eN_df_its,paste0(output,"indicsp_eN_its.txt"),sep="\t",quote=F)

## Import data frame of indicator species to save time
indic_eN_df_its <- read.table("indicsp_eN_its.txt", header=T, sep="\t")

CLaP_eN_its <- as.matrix(indic_eN_df_its[which(indic_eN_df_its$s.CLaP == 1 & indic_eN_df_its$p.value < 0.05),])
CLeP_eN_its <- as.matrix(indic_eN_df_its[which(indic_eN_df_its$s.CLeP == 1 & indic_eN_df_its$p.value < 0.05),])
uCaP_eN_its <- as.matrix(indic_eN_df_its[which(indic_eN_df_its$s.uCaP == 1 & indic_eN_df_its$p.value < 0.05),])
uCeP_eN_its <- as.matrix(indic_eN_df_its[which(indic_eN_df_its$s.uCeP == 1 & indic_eN_df_its$p.value < 0.05),])

eN_r_values_its <- rbind(CLaP_eN_its,CLeP_eN_its,uCaP_eN_its,uCeP_eN_its)
colnames(eN_r_values_its)[1:4] <-c("CLaP","CLeP","uCaP","uCeP")

## Range of correlation coefficients
range(eN_r_values_its[,"stat"])

## Total number of indicator OTUS
length(unique(rownames(eN_r_values_its)))

## Proportion of eN bacteria OTUs responding to cropping system
length(unique(rownames(eN_r_values_its)))/nrow(otu_norm_eN_its)

## Proportion of eN bacteria sequences responding to cropping system
eN_its_ra <- t(t(otu_its_eN)/colSums(otu_its_eN)) * 100
sum(colSums(eN_its_ra[unique(rownames(eN_r_values_its)),]))/sum(colSums(eN_its_ra))


## Identify indicator species in bulk eN bacteria communities
indic_eN_18s <- as.data.frame(t(otu_norm_eN_18s))
indic_eN_groups_18s <- design_18s_eN$CL_P
length(unique(indic_eN_groups_18s))

## Define indicator species for eN bacteria community.
## Note: These calculations can be time and processor intensive
set.seed(8046)
indicatorsp_eN_18s <- multipatt(indic_eN_18s,indic_eN_groups_18s,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_eN_18s,alpha=1,indvalcomp=T)
indic_eN_df_18s <- indicatorsp_eN_18s$sign
write.table(indic_eN_df_18s,paste0(output,"indicsp_eN_18s.txt"),sep="\t",quote=F)

## Import data frame of indicator species to save time
indic_eN_df_18s <- read.table("indicsp_eN_18s.txt", header=T, sep="\t")

CLaP_eN_18s <- as.matrix(indic_eN_df_18s[which(indic_eN_df_18s$s.CLaP == 1 & indic_eN_df_18s$p.value < 0.05),])
CLeP_eN_18s <- as.matrix(indic_eN_df_18s[which(indic_eN_df_18s$s.CLeP == 1 & indic_eN_df_18s$p.value < 0.05),])
uCaP_eN_18s <- as.matrix(indic_eN_df_18s[which(indic_eN_df_18s$s.uCaP == 1 & indic_eN_df_18s$p.value < 0.05),])
uCeP_eN_18s <- as.matrix(indic_eN_df_18s[which(indic_eN_df_18s$s.uCeP == 1 & indic_eN_df_18s$p.value < 0.05),])

eN_r_values_18s <- rbind(CLaP_eN_18s,CLeP_eN_18s,uCaP_eN_18s,uCeP_eN_18s)
colnames(eN_r_values_18s)[1:4] <-c("CLaP","CLeP","uCaP","uCeP")

## Range of correlation coefficients
range(eN_r_values_18s[,"stat"])

## Total number of indicator OTUS
length(unique(rownames(eN_r_values_18s)))

## Proportion of eN bacteria OTUs responding to cropping system
length(unique(rownames(eN_r_values_18s)))/nrow(otu_norm_eN_18s)

## Proportion of eN bacteria sequences responding to cropping system
eN_18s_ra <- t(t(otu_18s_eN)/colSums(otu_18s_eN)) * 100
sum(colSums(eN_18s_ra[unique(rownames(eN_r_values_18s)),]))/sum(colSums(eN_18s_ra))








##########===============================================================================
#########=================================================================================
##### Bacteria
##### TMM normalize 16S/ITS counts for whole community beta diversity analysis #####

## Apply TMM normalization to entire 16s data set and create phyloseq objects for later analysis
group_16s <- design_filter_16s$Treatments
edgeR_16s <- DGEList(counts=otu_filter_16s, 
                     group=design_filter_16s$Treatments, 
                     genes=tax_filter_16s)

edgeR_16s <- calcNormFactors(edgeR_16s)

## Extract normalized counts
otu_norm_16s <- cpm(edgeR_16s, normalized.lib.sizes=T, log=F)
## Create phyloseq objects
physeq_16s_norm <- phyloseq(otu_table(otu_norm_16s, taxa_are_rows=T),
                            tax_table(as.matrix(tax_filter_16s)),
                            sample_data(design_filter_16s))

## Create bray-curtis dissimiliartiy matrix
all_dis_16s <- vegdist(t(otu_table(physeq_16s_norm)),method="bray")

## Express 16S OTU counts as relative abunance percent
otu_16s_RA <- t(t(otu_filter_16s)/colSums(otu_filter_16s)) * 100
colSums(otu_16s_RA)
nrow(otu_16s_RA)
## Get names of bacteria phyla present (use 'labels' as this specifies class within Proteobacteria)
PHYLAnames_16s <- names(sort(table(tax_filter_16s[,"labels"]), decr=T))
length(PHYLAnames_16s)
sort(table(tax_filter_16s[,"labels"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(otu_16s_RA)
for (i in PHYLAnames_16s){
  x <- array(colSums(otu_16s_RA[rownames(tax_filter_16s)[which(tax_filter_16s$labels == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}
## Create matrix
rownames(y) <- paste(PHYLAnames_16s)
colnames(y) <- paste(colnames(otu_16s_RA))
PHYLUM_mat_16s <- y
PHYLUM_mat_16s[,1:5]
colSums(PHYLUM_mat_16s)
PHYLUM_mat_16s_mean <- sort(apply(PHYLUM_mat_16s,1,mean),decr=T)
PHYLUM_mat_16s <- PHYLUM_mat_16s[names(PHYLUM_mat_16s_mean),]

## Use hierarchical clustering to order samples
norm_clu_16s <- hclust(all_dis_16s, "average")
norm_clu_16s$height
Beta_labs_16s <- norm_clu_16s$labels
Beta_order_16s <- norm_clu_16s$order
Beta_labs_16s <- Beta_labs_16s[Beta_order_16s]
PHYLUM_mat_16s <- PHYLUM_mat_16s[ ,Beta_labs_16s]

## aN phyla abundances
PHYLUM_mat_16s_aN <- PHYLUM_mat_16s[,aNsamples]
colSums(PHYLUM_mat_16s_aN)
PHYLUM_mat_16s_aN_mean <- sort(apply(PHYLUM_mat_16s_aN,1,mean),decr=T)
PHYLUM_mat_16s_aN <- PHYLUM_mat_16s_aN[names(PHYLUM_mat_16s_aN_mean),]
PHYLUM_mat_16s_aN_se <- apply(PHYLUM_mat_16s_aN,1,se)[names(PHYLUM_mat_16s_aN_mean)]

length(PHYLUM_mat_16s_aN_mean[PHYLUM_mat_16s_aN_mean > 0])

## eN phyla abundances
PHYLUM_mat_16s_eN <- PHYLUM_mat_16s[,eNsamples]
colSums(PHYLUM_mat_16s_eN)
PHYLUM_mat_16s_eN_mean <- apply(PHYLUM_mat_16s_eN,1,mean)[names(PHYLUM_mat_16s_aN_mean)]
PHYLUM_mat_16s_eN_se <- apply(PHYLUM_mat_16s_eN,1,se)[names(PHYLUM_mat_16s_aN_mean)]

length(PHYLUM_mat_16s_eN_mean[PHYLUM_mat_16s_eN_mean > 0])


### Defining bOTU colors by phylum (using the taxonomy file)
tax_filter_16s$cols <- tax_filter_16s$labels
table(tax_filter_16s$cols)

# Phyla with MEAN abundances lower than 1% relative abundances
table(apply(PHYLUM_mat_16s, 1, mean) < 1)
low_count_phyla_16s <- rownames(PHYLUM_mat_16s)[sort(apply(PHYLUM_mat_16s, 1, mean), decr=T) < 1]
# attribute grey color
for(i in low_count_phyla_16s){
  tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$Phylum==paste(i) ], ]$cols <- "lightgrey"
}
table(tax_filter_16s$cols)

# Phyla with MEAN abundances higher than 1% relative abundances
abundant_phyla_16s <- rownames(PHYLUM_mat_16s)[sort(apply(PHYLUM_mat_16s, 1, mean), decr=T) > 1]
abundant_phyla_16s
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Alphaproteobacteria" ], ]$cols <- "palegreen1"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Betaproteobacteria" ], ]$cols <- "palegreen3"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Gammaproteobacteria" ], ]$cols <- "palegreen2"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Deltaproteobacteria" ], ]$cols <- "palegreen4"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Actinobacteria" ], ]$cols <- "indianred2"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Bacteroidetes" ], ]$cols <- "steelblue1"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Firmicutes" ], ]$cols <- "tan1"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Acidobacteria" ], ]$cols <- "lightsalmon4"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Chloroflexi" ], ]$cols <- "gold1"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Verrucomicrobia" ], ]$cols <- "orchid3"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Nitrospirae" ], ]$cols <- "palevioletred2"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Gemmatimonadetes" ], ]$cols <- "peachpuff3"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Planctomycetes" ], ]$cols <- "peachpuff4"

## collaps OTU colors to prepare Phylum level colors
label_cols_16s <- tax_filter_16s[, c("labels", "cols") ]
library(plyr)
PHYLA_label_cols_16s <- ddply(label_cols_16s, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_16s) <- PHYLA_label_cols_16s[,1]
PHYLA_label_cols_16s <- PHYLA_label_cols_16s[c(abundant_phyla_16s, low_count_phyla_16s),]
PHYLA_label_cols_16s

## Legend for Phylum colors
PHYLA_label_cols_16s_legend <- PHYLA_label_cols_16s[1:13,]
PHYLA_label_cols_16s_legend[13,1] <- "other"
rownames(PHYLA_label_cols_16s_legend)[13] <- "other"
PHYLA_label_cols_16s_legend

sort(table(tax_filter_16s[,"labels"]), decr=T)


##### Fungi
##### TMM normalize ITS counts for whole community beta diversity analysis #####

## Apply TMM normalization to entire its data set and create phyloseq objects for later analysis
group_its <- design_filter_its$Treatments
edgeR_its <- DGEList(counts=otu_filter_its, 
                     group=design_filter_its$Treatments, 
                     genes=tax_filter_its)

edgeR_its <- calcNormFactors(edgeR_its)

## Extract normalized counts
otu_norm_its <- cpm(edgeR_its, normalized.lib.sizes=T, log=F)
## Create phyloseq objects
physeq_its_norm <- phyloseq(otu_table(otu_norm_its, taxa_are_rows=T),
                            tax_table(as.matrix(tax_filter_its)),
                            sample_data(design_filter_its))

## Create bray-curtis dissimiliartiy matrix
all_dis_its <- vegdist(t(otu_table(physeq_its_norm)),method="bray")

## Express its OTU counts as relative abunance percent
otu_its_RA <- t(t(otu_filter_its)/colSums(otu_filter_its)) * 100
colSums(otu_its_RA)
nrow(otu_its_RA)
## Get names of fungi phyla present
PHYLAnames_its <- names(sort(table(tax_filter_its[,"Phylum"]), decr=T))
length(PHYLAnames_its)
sort(table(tax_filter_its[,"Phylum"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(otu_its_RA)
for (i in PHYLAnames_its){
  x <- array(colSums(otu_its_RA[rownames(tax_filter_its)[which(tax_filter_its$Phylum == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}
## Create matrix
rownames(y) <- paste(PHYLAnames_its)
colnames(y) <- paste(colnames(otu_its_RA))
PHYLUM_mat_its <- y
PHYLUM_mat_its[,1:5]
colSums(PHYLUM_mat_its)
PHYLUM_mat_its_mean <- sort(apply(PHYLUM_mat_its,1,mean),decr=T)
PHYLUM_mat_its <- PHYLUM_mat_its[names(PHYLUM_mat_its_mean),]

## Use hierarchical clustering to order samples
norm_clu_its <- hclust(all_dis_its, "average")
norm_clu_its$height
Beta_labs_its <- norm_clu_its$labels
Beta_order_its <- norm_clu_its$order
Beta_labs_its <- Beta_labs_its[Beta_order_its]
PHYLUM_mat_its <- PHYLUM_mat_its[ ,Beta_labs_its]


## Soil phyla abundances
PHYLUM_mat_its_aN <- PHYLUM_mat_its[,aNsamples]
colSums(PHYLUM_mat_its_aN)
PHYLUM_mat_its_aN_mean <- sort(apply(PHYLUM_mat_its_aN,1,mean),decr=T)
PHYLUM_mat_its_aN <- PHYLUM_mat_its_aN[names(PHYLUM_mat_its_aN_mean),]
PHYLUM_mat_its_aN_se <- apply(PHYLUM_mat_its_aN,1,se)[names(PHYLUM_mat_its_aN_mean)]

length(PHYLUM_mat_its_aN_mean[PHYLUM_mat_its_aN_mean > 0])

## eN phyla abundances
PHYLUM_mat_its_eN <- PHYLUM_mat_its[,eNsamples]
colSums(PHYLUM_mat_its_eN)
PHYLUM_mat_its_eN_mean <- apply(PHYLUM_mat_its_eN,1,mean)[names(PHYLUM_mat_its_eN_mean)]
PHYLUM_mat_its_eN_se <- apply(PHYLUM_mat_its_eN,1,se)[names(PHYLUM_mat_its_eN_mean)]

length(PHYLUM_mat_its_eN_mean[PHYLUM_mat_its_eN_mean > 0])


### Defining bOTU colors by phylum (using the taxonomy file)
tax_filter_its$cols <- tax_filter_its$Phylum
table(tax_filter_its$cols)

# Phyla with MEAN abundances lower than 1% relative abundances
table(apply(PHYLUM_mat_its, 1, mean) < 1)
low_count_phyla_its <- rownames(PHYLUM_mat_its)[sort(apply(PHYLUM_mat_its, 1, mean), decr=T) < 1]
# attribute grey color
for(i in low_count_phyla_its){
  tax_filter_its[ rownames(tax_filter_its)[tax_filter_its$Phylum==paste(i) ], ]$cols <- "lightgrey"
}
table(tax_filter_its$cols)

# Phyla with MEAN abundances higher than 1% relative abundances
abundant_phyla_its <- rownames(PHYLUM_mat_its)[sort(apply(PHYLUM_mat_its, 1, mean), decr=T) > 1]
abundant_phyla_its
tax_filter_its[ rownames(tax_filter_its)[tax_filter_its$Phylum=="Ascomycota" ], ]$cols <- "indianred2"
tax_filter_its[ rownames(tax_filter_its)[tax_filter_its$Phylum=="Basidiomycota" ], ]$cols <- "steelblue1"
tax_filter_its[ rownames(tax_filter_its)[tax_filter_its$Phylum=="Chytridiomycota" ], ]$cols <- "palegreen1"
tax_filter_its[ rownames(tax_filter_its)[tax_filter_its$Phylum=="Glomeromycota" ], ]$cols <- "blue2"
tax_filter_its[ rownames(tax_filter_its)[tax_filter_its$Phylum=="Rozellomycota" ], ]$cols <- "deeppink2"
tax_filter_its[ rownames(tax_filter_its)[tax_filter_its$Phylum=="Zygomycota" ], ]$cols <- "coral1"

## collaps OTU colors to prepare Phylum level colors
label_cols_its <- tax_filter_its[, c("Phylum", "cols") ]
library(plyr)
PHYLA_label_cols_its <- ddply(label_cols_its, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_its) <- PHYLA_label_cols_its[,1]
PHYLA_label_cols_its <- PHYLA_label_cols_its[c(abundant_phyla_its, low_count_phyla_its),]
PHYLA_label_cols_its

## Legend for Phylum colors
PHYLA_label_cols_its_legend <- PHYLA_label_cols_its[1:6,]
#PHYLA_label_cols_its_legend[13,1] <- "other"
#rownames(PHYLA_label_cols_its_legend)[13] <- "other"
PHYLA_label_cols_its_legend







##### Protist
##### TMM normalize 18S counts for whole community beta diversity analysis #####

## Apply TMM normalization to entire 18s data set and create phyloseq objects for later analysis
group_18s <- design_filter_18s$Treatments
edgeR_18s <- DGEList(counts=otu_filter_18s, 
                     group=design_filter_18s$Treatments, 
                     genes=tax_filter_18s)

edgeR_18s <- calcNormFactors(edgeR_18s)

## Extract normalized counts
otu_norm_18s <- cpm(edgeR_18s, normalized.lib.sizes=T, log=F)
## Create phyloseq objects
physeq_18s_norm <- phyloseq(otu_table(otu_norm_18s, taxa_are_rows=T),
                            tax_table(as.matrix(tax_filter_18s)),
                            sample_data(design_filter_18s))

## Create bray-curtis dissimiliartiy matrix
all_dis_18s <- vegdist(t(otu_table(physeq_18s_norm)),method="bray")

## Express 18S OTU counts as relative abunance percent
otu_18s_RA <- t(t(otu_filter_18s)/colSums(otu_filter_18s)) * 100
colSums(otu_18s_RA)
nrow(otu_18s_RA)

## Get names of protist Kingdom present
PHYLAnames_18s <- names(sort(table(tax_filter_18s[,"Kingdom"]), decr=T))
length(PHYLAnames_18s)
sort(table(tax_filter_18s[,"Kingdom"]), decr=T)

## Preparation of matrix with relative abundance by Kingdom
y <- NULL
otunames <- rownames(otu_18s_RA)
for (i in PHYLAnames_18s){
  x <- array(colSums(otu_18s_RA[rownames(tax_filter_18s)[which(tax_filter_18s$Kingdom == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}
## Create matrix
rownames(y) <- paste(PHYLAnames_18s)
colnames(y) <- paste(colnames(otu_18s_RA))
PHYLUM_mat_18s <- y
PHYLUM_mat_18s[,1:5]
colSums(PHYLUM_mat_18s)
PHYLUM_mat_18s_mean <- sort(apply(PHYLUM_mat_18s,1,mean),decr=T)
PHYLUM_mat_18s <- PHYLUM_mat_18s[names(PHYLUM_mat_18s_mean),]

## Use hierarchical clustering to order samples
norm_clu_18s <- hclust(all_dis_18s, "average")
norm_clu_18s$height
Beta_labs_18s <- norm_clu_18s$labels
Beta_order_18s <- norm_clu_18s$order
Beta_labs_18s <- Beta_labs_18s[Beta_order_18s]
PHYLUM_mat_18s <- PHYLUM_mat_18s[ ,Beta_labs_18s]


## Soil phyla abundances
PHYLUM_mat_18s_aN <- PHYLUM_mat_18s[,aNsamples]
colSums(PHYLUM_mat_18s_aN)
PHYLUM_mat_18s_aN_mean <- sort(apply(PHYLUM_mat_18s_aN,1,mean),decr=T)
PHYLUM_mat_18s_aN <- PHYLUM_mat_18s_aN[names(PHYLUM_mat_18s_aN_mean),]
PHYLUM_mat_18s_aN_se <- apply(PHYLUM_mat_18s_aN,1,se)[names(PHYLUM_mat_18s_aN_mean)]

length(PHYLUM_mat_18s_aN_mean[PHYLUM_mat_18s_aN_mean > 0])

## Root phyla abundances
PHYLUM_mat_18s_eN <- PHYLUM_mat_18s[,eNsamples]
colSums(PHYLUM_mat_18s_eN)
PHYLUM_mat_18s_eN_mean <- apply(PHYLUM_mat_18s_eN,1,mean)[names(PHYLUM_mat_18s_aN_mean)]
PHYLUM_mat_18s_eN_se <- apply(PHYLUM_mat_18s_eN,1,se)[names(PHYLUM_mat_18s_aN_mean)]

length(PHYLUM_mat_18s_eN_mean[PHYLUM_mat_18s_eN_mean > 0])


### Defining bOTU colors by phylum (using the taxonomy file)
tax_filter_18s$cols <- tax_filter_18s$Kingdom
table(tax_filter_18s$cols)

# Phyla with MEAN abundances lower than 1% relative abundances
table(apply(PHYLUM_mat_18s, 1, mean) < 1)
low_count_phyla_18s <- rownames(PHYLUM_mat_18s)[sort(apply(PHYLUM_mat_18s, 1, mean), decr=T) < 1]
# attribute grey color
for(i in low_count_phyla_18s){
  tax_filter_18s[rownames(tax_filter_18s)[tax_filter_18s$Kingdom==paste(i) ], ]$cols <- "lightgrey"
}
table(tax_filter_18s$cols)

# Phyla with MEAN abundances higher than 1% relative abundances
abundant_phyla_18s <- rownames(PHYLUM_mat_18s)[sort(apply(PHYLUM_mat_18s, 1, mean), decr=T) > 1]
abundant_phyla_18s
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Alveolata" ], ]$cols <- "indianred2"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Rhizaria" ], ]$cols <- "gold2"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Amoebozoa" ], ]$cols <- "steelblue1"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Stramenopiles" ], ]$cols <- "blue2"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Hacrobia" ], ]$cols <- "palegreen1"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Opisthokonta" ], ]$cols <- "darkorchid2"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Excavata" ], ]$cols <- "tan1"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Apusozoa" ], ]$cols <- "lightsalmon4"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Protalveolata" ], ]$cols <- "deeppink2"


## collaps OTU colors to prepare Phylum level colors
label_cols_18s <- tax_filter_18s[, c("Kingdom", "cols") ]
library(plyr)
PHYLA_label_cols_18s <- ddply(label_cols_18s, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_18s) <- PHYLA_label_cols_18s[,1]
PHYLA_label_cols_18s <- PHYLA_label_cols_18s[c(abundant_phyla_18s, low_count_phyla_18s),]
PHYLA_label_cols_18s

## Legend for Phylum colors
PHYLA_label_cols_18s_legend <- PHYLA_label_cols_18s[1:9,]
PHYLA_label_cols_18s_legend[9,1] <- "Protalveolata"
rownames(PHYLA_label_cols_18s_legend)[9] <- "Protalveolata"
PHYLA_label_cols_18s_legend












##### Figure 3: Bipartite networks of OTUs associated with cropping systems #####

### Bulk aN BACTERIA community
## Construct node table for bulk aN bacteria communities from indicator species data
aN_bipartite_16s <- data.frame(from= c(rep("CLaP",length(which(aN_r_values_16s[,"CLaP"]==1))),
                                         rep("CLeP",length(which(aN_r_values_16s[,"CLeP"]==1))),
                                         rep("uCaP",length(which(aN_r_values_16s[,"uCaP"]==1))),
                                         rep("uCeP",length(which(aN_r_values_16s[,"uCeP"]==1)))),
                                 to= c(rownames(aN_r_values_16s)[which(aN_r_values_16s[,"CLaP"]==1)],
                                       rownames(aN_r_values_16s)[which(aN_r_values_16s[,"CLeP"]==1)],
                                       rownames(aN_r_values_16s)[which(aN_r_values_16s[,"uCaP"]==1)],
                                       rownames(aN_r_values_16s)[which(aN_r_values_16s[,"uCeP"]==1)]),
                                 r= c(aN_r_values_16s[which(aN_r_values_16s[,"CLaP"]==1),"stat"],
                                      aN_r_values_16s[which(aN_r_values_16s[,"CLeP"]==1),"stat"],
                                      aN_r_values_16s[which(aN_r_values_16s[,"uCaP"]==1),"stat"],
                                      aN_r_values_16s[which(aN_r_values_16s[,"uCeP"]==1),"stat"]))

## make node attribute table for each OTU
aN_bipartite_attrib_16s <- data.frame(node=unique(rownames(aN_r_values_16s)),indicgroup=0)

for (i in as.character(aN_bipartite_attrib_16s$node))
{
  aN_bipartite_attrib_16s[aN_bipartite_attrib_16s$node==i,"indicgroup"] <- paste(colnames(aN_r_values_16s)[which(aN_r_values_16s[i,1:4]==1)],collapse = "_")
}

aN_bipartite_attrib_16s <- cbind(aN_bipartite_attrib_16s,tax_aN_16s[as.character(aN_bipartite_attrib_16s$node),])

## Create bipartite network with igraph
aN_bi_16s <- graph.data.frame(aN_bipartite_16s,directed=F)
V(aN_bi_16s)$type <- V(aN_bi_16s)$name %in% aN_bipartite_16s[,1]
aN_bi_16s <- simplify(aN_bi_16s, remove.multiple=T, remove.loops=T)

## Set labels for nodes
V(aN_bi_16s)$label <- V(aN_bi_16s)$name
V(aN_bi_16s)$label <- gsub("bOTU*",NA,V(aN_bi_16s)$label)      ### note bOTU or botu capital letter or lower letter

## Set node sizes
V(aN_bi_16s)$size <- c(rep(15,4),rep(3,438))      #  438 = 442-4

## Set node shapes
V(aN_bi_16s)$shape <- c(rep("circle",4),rep("circle",438))

## Define node colors based upon phylum/class taxonomy assignment
V(aN_bi_16s)$color <- V(aN_bi_16s)$name
V(aN_bi_16s)$color[1:4] <- "white"
V(aN_bi_16s)$color <- tax_filter_16s[V(aN_bi_16s)$name, ]$cols
V(aN_bi_16s)$frame.color <- V(aN_bi_16s)$color



### Bulk aN FUNGI community
## Construct node table for bulk aN fungi communities from indicator species data
aN_bipartite_its <- data.frame(from= c(rep("CLaP",length(which(aN_r_values_its[,"CLaP"]==1))),
                                         rep("CLeP",length(which(aN_r_values_its[,"CLeP"]==1))),
                                         rep("uCaP",length(which(aN_r_values_its[,"uCaP"]==1))),
                                         rep("uCeP",length(which(aN_r_values_its[,"uCeP"]==1)))),
                                 to= c(rownames(aN_r_values_its)[which(aN_r_values_its[,"CLaP"]==1)],
                                       rownames(aN_r_values_its)[which(aN_r_values_its[,"CLeP"]==1)],
                                       rownames(aN_r_values_its)[which(aN_r_values_its[,"uCaP"]==1)],
                                       rownames(aN_r_values_its)[which(aN_r_values_its[,"uCeP"]==1)]),
                                 r= c(aN_r_values_its[which(aN_r_values_its[,"CLaP"]==1),"stat"],
                                      aN_r_values_its[which(aN_r_values_its[,"CLeP"]==1),"stat"],
                                      aN_r_values_its[which(aN_r_values_its[,"uCaP"]==1),"stat"],
                                      aN_r_values_its[which(aN_r_values_its[,"uCeP"]==1),"stat"]))

## make node attribute table for each OTU
aN_bipartite_attrib_its <- data.frame(node=unique(rownames(aN_r_values_its)),indicgroup=0)

for (i in as.character(aN_bipartite_attrib_its$node))
{
  aN_bipartite_attrib_its[aN_bipartite_attrib_its$node==i,"indicgroup"] <- paste(colnames(aN_r_values_its)[which(aN_r_values_its[i,1:4]==1)],collapse = "_")
}

aN_bipartite_attrib_its <- cbind(aN_bipartite_attrib_its,tax_aN_its[as.character(aN_bipartite_attrib_its$node),])

## Create bipartite network with igraph
aN_bi_its <- graph.data.frame(aN_bipartite_its,directed=F)
V(aN_bi_its)$type <- V(aN_bi_its)$name %in% aN_bipartite_its[,1]
aN_bi_its <- simplify(aN_bi_its, remove.multiple=T, remove.loops=T)

## Set labels for nodes
V(aN_bi_its)$label <- V(aN_bi_its)$name
V(aN_bi_its)$label <- gsub("fotu*",NA,V(aN_bi_its)$label)

## Set sizes of nodes
V(aN_bi_its)$size <- c(rep(15,4),rep(3,84))

## Set node shape
V(aN_bi_its)$shape <- c(rep("circle",4),rep("triangle",84))

## Set node colors by taxonomy assignment as phylum level
V(aN_bi_its)$color <- V(aN_bi_its)$name
V(aN_bi_its)$color[1:4] <- "white"
V(aN_bi_its)$color <- tax_filter_its[ V(aN_bi_its)$name, ]$cols
V(aN_bi_its)$frame.color <- V(aN_bi_its)$color


### Bulk aN Protist community
## Construct node table for bulk aN protist communities from indicator species data
aN_bipartite_18s <- data.frame(from= c(rep("CLaP",length(which(aN_r_values_18s[,"CLaP"]==1))),
                                         rep("CLeP",length(which(aN_r_values_18s[,"CLeP"]==1))),
                                         rep("uCaP",length(which(aN_r_values_18s[,"uCaP"]==1))),
                                         rep("uCeP",length(which(aN_r_values_18s[,"uCeP"]==1)))),
                                 to= c(rownames(aN_r_values_18s)[which(aN_r_values_18s[,"CLaP"]==1)],
                                       rownames(aN_r_values_18s)[which(aN_r_values_18s[,"CLeP"]==1)],
                                       rownames(aN_r_values_18s)[which(aN_r_values_18s[,"uCaP"]==1)],
                                       rownames(aN_r_values_18s)[which(aN_r_values_18s[,"uCeP"]==1)]),
                                 r= c(aN_r_values_18s[which(aN_r_values_18s[,"CLaP"]==1),"stat"],
                                      aN_r_values_18s[which(aN_r_values_18s[,"CLeP"]==1),"stat"],
                                      aN_r_values_18s[which(aN_r_values_18s[,"uCaP"]==1),"stat"],
                                      aN_r_values_18s[which(aN_r_values_18s[,"uCeP"]==1),"stat"]))

## make node attribute table for each OTU
aN_bipartite_attrib_18s <- data.frame(node=unique(rownames(aN_r_values_18s)),indicgroup=0)

for (i in as.character(aN_bipartite_attrib_18s$node))
{
  aN_bipartite_attrib_18s[aN_bipartite_attrib_18s$node==i,"indicgroup"] <- paste(colnames(aN_r_values_18s)[which(aN_r_values_18s[i,1:4]==1)],collapse = "_")
}

aN_bipartite_attrib_18s <- cbind(aN_bipartite_attrib_18s,tax_aN_18s[as.character(aN_bipartite_attrib_18s$node),])

## Create bipartite network with igraph
aN_bi_18s <- graph.data.frame(aN_bipartite_18s,directed=F)
V(aN_bi_18s)$type <- V(aN_bi_18s)$name %in% aN_bipartite_18s[,1]
aN_bi_18s <- simplify(aN_bi_18s, remove.multiple=T, remove.loops=T)

## Set labels for nodes
V(aN_bi_18s)$label <- V(aN_bi_18s)$name
V(aN_bi_18s)$label <- gsub("potu*",NA,V(aN_bi_18s)$label)            ###   potu or pOTU?

## Set node sizes
V(aN_bi_18s)$size <- c(rep(15,4),rep(3,145))

## Set node shapes
V(aN_bi_18s)$shape <- c(rep("circle",4),rep("square",145))

## Define node colors based upon phylum/class taxonomy assignment
V(aN_bi_18s)$color <- V(aN_bi_18s)$name
V(aN_bi_18s)$color[1:4] <- "white"
V(aN_bi_18s)$color <- tax_filter_18s[ V(aN_bi_18s)$name, ]$cols
V(aN_bi_18s)$frame.color <- V(aN_bi_18s)$color


### eN associated BACTERIA community
## Construct node table for eN-associated bacteria communities from indicator species data
eN_bipartite_16s <- data.frame(from= c(rep("CLaP",length(which(eN_r_values_16s[,"CLaP"]==1))),
                                         rep("CLeP",length(which(eN_r_values_16s[,"CLeP"]==1))),
                                         rep("uCaP",length(which(eN_r_values_16s[,"uCaP"]==1))),
                                         rep("uCeP",length(which(eN_r_values_16s[,"uCeP"]==1)))),
                                 to= c(rownames(eN_r_values_16s)[which(eN_r_values_16s[,"CLaP"]==1)],
                                       rownames(eN_r_values_16s)[which(eN_r_values_16s[,"CLeP"]==1)],
                                       rownames(eN_r_values_16s)[which(eN_r_values_16s[,"uCaP"]==1)],
                                       rownames(eN_r_values_16s)[which(eN_r_values_16s[,"uCeP"]==1)]),
                                 r= c(eN_r_values_16s[which(eN_r_values_16s[,"CLaP"]==1),"stat"],
                                      eN_r_values_16s[which(eN_r_values_16s[,"CLeP"]==1),"stat"],
                                      eN_r_values_16s[which(eN_r_values_16s[,"uCaP"]==1),"stat"],
                                      eN_r_values_16s[which(eN_r_values_16s[,"uCeP"]==1),"stat"]))

## make node attribute table for each OTU
eN_bipartite_attrib_16s <- data.frame(node=unique(rownames(eN_r_values_16s)),indicgroup=0)

for (i in as.character(eN_bipartite_attrib_16s$node))
{
  eN_bipartite_attrib_16s[eN_bipartite_attrib_16s$node==i,"indicgroup"] <- paste(colnames(eN_r_values_16s)[which(eN_r_values_16s[i,1:4]==1)],collapse = "_")
}

eN_bipartite_attrib_16s <- cbind(eN_bipartite_attrib_16s,tax_eN_16s[as.character(eN_bipartite_attrib_16s$node),])

## Create bipartite network with igraph
eN_bi_16s <- graph.data.frame(eN_bipartite_16s,directed=F)
V(eN_bi_16s)$type <- V(eN_bi_16s)$name %in% eN_bipartite_16s[,1]
eN_bi_16s <- simplify(eN_bi_16s, remove.multiple=T, remove.loops=T)

## Set node labels
V(eN_bi_16s)$label <- V(eN_bi_16s)$name
V(eN_bi_16s)$label <- gsub("bOTU*",NA,V(eN_bi_16s)$label)

## Set node sizes
V(eN_bi_16s)$size <- c(rep(15,4),rep(3,483))

## Set node shape
V(eN_bi_16s)$shape <- c(rep("circle",4), rep("circle",483))

## Define node colors based upon phylum/class taxonomy assignment
V(eN_bi_16s)$color <- V(eN_bi_16s)$name
V(eN_bi_16s)$color[1:4] <- "white"
V(eN_bi_16s)$color <- tax_filter_16s[ V(eN_bi_16s)$name, ]$cols
V(eN_bi_16s)$frame.color <- V(eN_bi_16s)$color



### eN associated FUNGI community
## Construct node table for eN-associated fungi communities from indicator species data
eN_bipartite_its <- data.frame(from= c(rep("CLaP",length(which(eN_r_values_its[,"CLaP"]==1))),
                                         rep("CLeP",length(which(eN_r_values_its[,"CLeP"]==1))),
                                         rep("uCaP",length(which(eN_r_values_its[,"uCaP"]==1))),
                                         rep("uCeP",length(which(eN_r_values_its[,"uCeP"]==1)))),
                                 to= c(rownames(eN_r_values_its)[which(eN_r_values_its[,"CLaP"]==1)],
                                       rownames(eN_r_values_its)[which(eN_r_values_its[,"CLeP"]==1)],
                                       rownames(eN_r_values_its)[which(eN_r_values_its[,"uCaP"]==1)],
                                       rownames(eN_r_values_its)[which(eN_r_values_its[,"uCeP"]==1)]),
                                 r= c(eN_r_values_its[which(eN_r_values_its[,"CLaP"]==1),"stat"],
                                      eN_r_values_its[which(eN_r_values_its[,"CLeP"]==1),"stat"],
                                      eN_r_values_its[which(eN_r_values_its[,"uCaP"]==1),"stat"],
                                      eN_r_values_its[which(eN_r_values_its[,"uCeP"]==1),"stat"]))

## make node attribute table for each OTU
eN_bipartite_attrib_its <- data.frame(node=unique(rownames(eN_r_values_its)),indicgroup=0)

for (i in as.character(eN_bipartite_attrib_its$node))
{
  eN_bipartite_attrib_its[eN_bipartite_attrib_its$node==i,"indicgroup"] <- paste(colnames(eN_r_values_its)[which(eN_r_values_its[i,1:4]==1)],collapse = "_")
}

eN_bipartite_attrib_its <- cbind(eN_bipartite_attrib_its,tax_eN_its[as.character(eN_bipartite_attrib_its$node),])

## Create bipartite network with igraph
eN_bi_its <- graph.data.frame(eN_bipartite_its,directed=F)
V(eN_bi_its)$type <- V(eN_bi_its)$name %in% eN_bipartite_its[,1]
eN_bi_its <- simplify(eN_bi_its, remove.multiple=T, remove.loops=T)

## Set node labels
V(eN_bi_its)$label <- V(eN_bi_its)$name
V(eN_bi_its)$label <- gsub("fotu*",NA,V(eN_bi_its)$label)

## Set node size
length(V(eN_bi_its)$label)
V(eN_bi_its)$size <- c(rep(15,4), rep(3,70))      # eN_bi_its  examine the parameter variable 

## Set node shape
V(eN_bi_its)$shape <- c(rep("circle",4),rep("triangle",70))

## Set node colors by taxonomy assignment as phylum level
V(eN_bi_its)$color <- V(eN_bi_its)$name
V(eN_bi_its)$color[1:4] <- "white"
V(eN_bi_its)$color <- tax_filter_its[ V(eN_bi_its)$name, ]$cols
V(eN_bi_its)$frame.color <- V(eN_bi_its)$color

V(eN_bi_its)$width <- rep(0.2, length(V(eN_bi_its)$name))

### eN associated PROTIST community
## Construct node table for eN-associated protist communities from indicator species data
eN_bipartite_18s <- data.frame(from= c(rep("CLaP",length(which(eN_r_values_18s[,"CLaP"]==1))),
                                         rep("CLeP",length(which(eN_r_values_18s[,"CLeP"]==1))),
                                         rep("uCaP",length(which(eN_r_values_18s[,"uCaP"]==1))),
                                         rep("uCeP",length(which(eN_r_values_18s[,"uCeP"]==1)))),
                                 to= c(rownames(eN_r_values_18s)[which(eN_r_values_18s[,"CLaP"]==1)],
                                       rownames(eN_r_values_18s)[which(eN_r_values_18s[,"CLeP"]==1)],
                                       rownames(eN_r_values_18s)[which(eN_r_values_18s[,"uCaP"]==1)],
                                       rownames(eN_r_values_18s)[which(eN_r_values_18s[,"uCeP"]==1)]),
                                 r= c(eN_r_values_18s[which(eN_r_values_18s[,"CLaP"]==1),"stat"],
                                      eN_r_values_18s[which(eN_r_values_18s[,"CLeP"]==1),"stat"],
                                      eN_r_values_18s[which(eN_r_values_18s[,"uCaP"]==1),"stat"],
                                      eN_r_values_18s[which(eN_r_values_18s[,"uCeP"]==1),"stat"]))

## make node attribute table for each OTU
eN_bipartite_attrib_18s <- data.frame(node=unique(rownames(eN_r_values_18s)),indicgroup=0)

for (i in as.character(eN_bipartite_attrib_18s$node))
{
  eN_bipartite_attrib_18s[eN_bipartite_attrib_18s$node==i,"indicgroup"] <- paste(colnames(eN_r_values_18s)[which(eN_r_values_18s[i,1:4]==1)],collapse = "_")
}

eN_bipartite_attrib_18s <- cbind(eN_bipartite_attrib_18s,tax_eN_18s[as.character(eN_bipartite_attrib_18s$node),])

## Create bipartite network with igraph
eN_bi_18s <- graph.data.frame(eN_bipartite_18s,directed=F)
V(eN_bi_18s)$type <- V(eN_bi_18s)$name %in% eN_bipartite_18s[,1]
eN_bi_18s <- simplify(eN_bi_18s, remove.multiple=T, remove.loops=T)

## Set node labels
V(eN_bi_18s)$label <- V(eN_bi_18s)$name
V(eN_bi_18s)$label <- gsub("potu*",NA,V(eN_bi_18s)$label)

## Set node sizes
V(eN_bi_18s)$size <- c(rep(15,4),rep(3,126))

## Set node shape
V(eN_bi_18s)$shape <- c(rep("circle",4), rep("square",126))

## Define node colors based upon phylum/class taxonomy assignment
V(eN_bi_18s)$color <- V(eN_bi_18s)$name
V(eN_bi_18s)$color[1:4] <- "white"
V(eN_bi_18s)$color <- tax_filter_18s[ V(eN_bi_18s)$name, ]$cols
V(eN_bi_18s)$frame.color <- V(eN_bi_18s)$color


##### Plot Figure 3 : aN/eN Bipartite Networks
set.seed(8046)
aN_layout_16s <- layout_with_fr(aN_bi_16s, niter=9999)
aN_layout_its <- layout_with_fr(aN_bi_its, niter=9999)
aN_layout_18s <- layout_with_fr(aN_bi_18s, niter=9999)

eN_layout_16s <- layout_with_fr(eN_bi_16s, niter=9999)
eN_layout_its <- layout_with_fr(eN_bi_its, niter=9999)
eN_layout_18s <- layout_with_fr(eN_bi_18s, niter=9999)

pdf(paste0(output,"Figure3.pdf"), encoding="MacRoman", paper = "a4")
layout(matrix(c(1,2,3,4,5,6), nrow=3, byrow=T), c(3,3,3), c(4,4,4))

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(aN_bi_16s, vertex.label.cex=0, layout=aN_layout_16s, asp=0)
plot(eN_bi_16s, vertex.label.cex=0, layout=eN_layout_16s, asp=0)
plot(aN_bi_its, vertex.label.cex=0, layout=aN_layout_its, asp=0)
plot(eN_bi_its, vertex.label.cex=0, layout=eN_layout_its, asp=0)
plot(aN_bi_18s, vertex.label.cex=0, layout=aN_layout_18s, asp=0)
plot(eN_bi_18s, vertex.label.cex=0, layout=eN_layout_18s, asp=0)

plot(1, type="n", ann=F, axes=F)
legend("center", pch=19, bty="n", ncol=2, horiz=F, x.intersp=0.4, 
       legend=PHYLA_label_cols_16s_legend$labels, 
       col=PHYLA_label_cols_16s_legend$cols)

plot(1,type="n",ann=F,axes=F)
legend("center", pch=19, bty="n", ncol=1, horiz=F, x.intersp=0.4, 
       legend=PHYLA_label_cols_its_legend$Phylum, 
       col=PHYLA_label_cols_its_legend$cols)
	   
plot(1, type="n", ann=F, axes=F)
legend("center", pch=19, bty="n", ncol=1, horiz=F, x.intersp=0.4, 
       legend=PHYLA_label_cols_18s_legend$Kingdom, 
       col=PHYLA_label_cols_18s_legend$cols)	   

dev.off()

