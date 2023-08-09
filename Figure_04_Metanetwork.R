
rm(list=ls())
setwd("E:/工作资料集/数据/多伦 草原站/郑俊强老师 YF156-M201708043/MPL 20170976 and 77 32个样 （多因子一 GCB）/发表用/network analysis/Hartman et al Microbiome_N/Figure_4_MetaNetwork/R_input")
output <- "E:/工作资料集/数据/多伦 草原站/郑俊强老师 YF156-M201708043/MPL 20170976 and 77 32个样 （多因子一 GCB）/发表用/network analysis/Hartman et al Microbiome_N/Figure_4_MetaNetwork/R_output/"

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

################################################################################
#####   MODEL 01 Import and prepare data                                    ####
################################################################################
#######################
##### 16S #####
#######################


##### Import Data #####
otu_16s <- read.csv("bacterial_16s_otu.csv", row.names=1,sep=",", header=T, blank.lines.skip=F, check.names=F)
otu_16s <- as.matrix(otu_16s)
rownames(otu_16s) <- paste("b", rownames(otu_16s), sep="") # otu name 闂????? add "b"

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



####################################################################################
#####                 MODEL: Analysis and Figures                              ######
####################################################################################
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

##### 16s overall PERMANOVA #####
## Perform PERMANVOA testing for sample type and Nitrogen addition system effects 
paov_all_16s <- adonis(all_dis_16s ~Clipping*Fertilizer*Precipitation, data=design_filter_16s, permutations=9999)


## Apply TMM normalization to entire ITS data set and create phyloseq objects for later analysis
group_its <- design_filter_its$Treatments
edgeR_its<- DGEList(counts= otu_filter_its, 
                    group=design_filter_its$Treatments,
                    genes=tax_filter_its)

edgeR_its <- calcNormFactors(edgeR_its)

## Extract normalized counts
otu_norm_its <- cpm(edgeR_its, normalized.lib.sizes=T,log=F)

## Create phyloseq objects
phy_its_norm <- otu_table(otu_norm_its,taxa_are_rows=T)
phy_tax_norm_its <- tax_table(as.matrix(tax_filter_its))
phy_design_its <- sample_data(design_filter_its)
physeq_its_norm <- phyloseq(phy_its_norm,phy_tax_norm_its,phy_design_its)

## Create bray-curtis dissimiliartiy matrix
all_dis_its <- vegdist(t(otu_table(physeq_its_norm)),method="bray")

##### ITS overall PERMANOVA #####
## Perform PERMANVOA testing for sample type and Nitrogen addition system effects
paov_all_its <- adonis(all_dis_its ~Clipping*Fertilizer*Precipitation, data=design_filter_its, permutations=9999)

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

##### 18s overall PERMANOVA #####
## Perform PERMANVOA testing for sample type and Nitrogen addition system effects 
paov_all_18s <- adonis(all_dis_18s ~Clipping*Fertilizer*Precipitation, data=design_filter_18s, permutations=9999)
##### Supplementary Table S2: global PERMANOVA #####
paov_all_16s  
paov_all_its
paov_all_18s


###################################################################################
############       MODEL 02     Supplementary Figure S3: phyla relative abundance plots      #####
###################################################################################

##### 
##### Bacteria
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
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Planctomycetes" ], ]$cols <- "darkred"


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



##### Fungi
## Express ITS OTU counts as relative abunance percent
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

PHYLUM_mat_its <- PHYLUM_mat_its[,Beta_labs_its]

## aN phyla abundances
PHYLUM_mat_its_aN <- PHYLUM_mat_its[,aNsamples]
colSums(PHYLUM_mat_its_aN)
PHYLUM_mat_its_aN_mean <- sort(apply(PHYLUM_mat_its_aN,1,mean),decr=T)
PHYLUM_mat_its_aN <- PHYLUM_mat_its_aN[names(PHYLUM_mat_its_aN_mean),]
PHYLUM_mat_its_aN <- PHYLUM_mat_its_aN[names(PHYLUM_mat_its_aN_mean),]
PHYLUM_mat_its_aN_se <- apply(PHYLUM_mat_its_aN,1,se)[names(PHYLUM_mat_its_aN_mean)]

length(PHYLUM_mat_its_aN_mean[PHYLUM_mat_its_aN_mean > 0])

## eN phyla abundances
PHYLUM_mat_its_eN <- PHYLUM_mat_its[,eNsamples]
colSums(PHYLUM_mat_its_eN)
PHYLUM_mat_its_eN_mean <- sort(apply(PHYLUM_mat_its_eN,1,mean),decr=T)
PHYLUM_mat_its_eN <- PHYLUM_mat_its_eN[names(PHYLUM_mat_its_eN_mean),]
PHYLUM_mat_its_eN <- PHYLUM_mat_its_eN[names(PHYLUM_mat_its_aN_mean),]
PHYLUM_mat_its_eN_se <- apply(PHYLUM_mat_its_eN,1,se)[names(PHYLUM_mat_its_aN_mean)]

length(PHYLUM_mat_its_eN_mean[PHYLUM_mat_its_eN_mean > 0])



### Defining fOTU colors by phylum (using the taxonomy file)
tax_filter_its$cols <- tax_filter_its$Phylum
table(tax_filter_its$cols)

# # Phyla with MEAN abundances lower than 1% relative abundances
# table(apply(PHYLUM_mat_its, 1, mean) < 1)
# low_count_phyla_its <- rownames(PHYLUM_mat_its)[apply(PHYLUM_mat_its, 1, mean) < 1]
# # attribute grey color
# for(i in low_count_phyla_its){
#   tax_filter_its[ rownames(tax_filter_its)[tax_filter_its$Phylum==paste(i) ], ]$cols <- "lightgrey"
# }

# Phyla with MEAN abundances higher than 1% relative abundances
phyla_its <- names(sort(apply(PHYLUM_mat_its, 1, mean), decr=T) )
phyla_its

tax_filter_its[ rownames(tax_filter_its)[tax_filter_its$Phylum=="Ascomycota" ], ]$cols <- "indianred2"
tax_filter_its[ rownames(tax_filter_its)[tax_filter_its$Phylum=="Basidiomycota" ], ]$cols <- "steelblue1"
tax_filter_its[ rownames(tax_filter_its)[tax_filter_its$Phylum=="Chytridiomycota" ], ]$cols <- "palegreen1"
tax_filter_its[ rownames(tax_filter_its)[tax_filter_its$Phylum=="Glomeromycota" ], ]$cols <- "blue2"
tax_filter_its[ rownames(tax_filter_its)[tax_filter_its$Phylum=="Rozellomycota" ], ]$cols <- "deeppink2"
tax_filter_its[ rownames(tax_filter_its)[tax_filter_its$Phylum=="Zygomycota" ], ]$cols <- "coral1"

## collaps OTU colors to prepare Phylum colors
PHYLA_label_cols_its <- tax_filter_its[,c("Phylum", "cols")]
library(plyr)
PHYLA_label_cols_its <- ddply(PHYLA_label_cols_its, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_its) <- PHYLA_label_cols_its[,1]
PHYLA_label_cols_its <- PHYLA_label_cols_its[phyla_its, ]
PHYLA_label_cols_its

## Legend for Phylum colors
PHYLA_label_cols_its_legend <- PHYLA_label_cols_its[, ]
# PHYLA_label_cols_its_legend[5, 1] <- "other"
# rownames(PHYLA_label_cols_its_legend)[5] <- "other"
PHYLA_label_cols_its_legend




##### Protist
## Express ITS OTU counts as relative abunance percent
otu_18s_RA <- t(t(otu_filter_18s)/colSums(otu_filter_18s)) * 100
colSums(otu_18s_RA)
nrow(otu_18s_RA)

## Get names of Protist phyla present 
PHYLAnames_18s <- names(sort(table(tax_filter_18s[,"Kingdom"]), decr=T))
length(PHYLAnames_18s)
sort(table(tax_filter_18s[,"Kingdom"]), decr=T)

## Preparation of matrix with relative abundance by kingdom
y <- NULL
otunames <- rownames(otu_18s_RA)
for (i in PHYLAnames_18s){
  x <- array(colSums(otu_18s_RA[rownames(tax_filter_18s)[which(tax_filter_18s$Kingdom == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

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

PHYLUM_mat_18s <- PHYLUM_mat_18s[,Beta_labs_18s]

## aN phyla abundances
PHYLUM_mat_18s_aN <- PHYLUM_mat_18s[,aNsamples]
colSums(PHYLUM_mat_18s_aN)
PHYLUM_mat_18s_aN_mean <- sort(apply(PHYLUM_mat_18s_aN,1,mean),decr=T)
PHYLUM_mat_18s_aN <- PHYLUM_mat_18s_aN[names(PHYLUM_mat_18s_aN_mean),]
PHYLUM_mat_18s_aN <- PHYLUM_mat_18s_aN[names(PHYLUM_mat_18s_aN_mean),]
PHYLUM_mat_18s_aN_se <- apply(PHYLUM_mat_18s_aN,1,se)[names(PHYLUM_mat_18s_aN_mean)]

length(PHYLUM_mat_18s_aN_mean[PHYLUM_mat_18s_aN_mean > 0])

## eN phyla abundances
PHYLUM_mat_18s_eN <- PHYLUM_mat_18s[,eNsamples]
colSums(PHYLUM_mat_18s_eN)
PHYLUM_mat_18s_eN_mean <- sort(apply(PHYLUM_mat_18s_eN,1,mean),decr=T)
PHYLUM_mat_18s_eN <- PHYLUM_mat_18s_eN[names(PHYLUM_mat_18s_eN_mean),]
PHYLUM_mat_18s_eN <- PHYLUM_mat_18s_eN[names(PHYLUM_mat_18s_aN_mean),]
PHYLUM_mat_18s_eN_se <- apply(PHYLUM_mat_18s_eN,1,se)[names(PHYLUM_mat_18s_aN_mean)]

length(PHYLUM_mat_18s_eN_mean[PHYLUM_mat_18s_eN_mean > 0])



### Defining fOTU colors by kingdom (using the taxonomy file)
tax_filter_18s$cols <- tax_filter_18s$Kingdom
table(tax_filter_18s$cols)

# # Kingdom with MEAN abundances lower than 1% relative abundances
# table(apply(PHYLUM_mat_18s, 1, mean) < 1)
# low_count_phyla_18s <- rownames(PHYLUM_mat_18s)[apply(PHYLUM_mat_18s, 1, mean) < 1]
# # attribute grey color
# for(i in low_count_phyla_18s){
#   tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom==paste(i) ], ]$cols <- "lightgrey"
# }

# Kingdom with MEAN abundances higher than 1% relative abundances
phyla_18s <- names(sort(apply(PHYLUM_mat_18s, 1, mean), decr=T) )
phyla_18s

tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Alveolata" ], ]$cols <- "indianred2"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Rhizaria" ], ]$cols <- "gold2"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Amoebozoa" ], ]$cols <- "steelblue1"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Stramenopiles" ], ]$cols <- "blue2"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Hacrobia" ], ]$cols <- "palegreen1"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Opisthokonta" ], ]$cols <- "darkorchid2"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Excavata" ], ]$cols <- "tan1"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Apusozoa" ], ]$cols <- "lightsalmon4"
tax_filter_18s[ rownames(tax_filter_18s)[tax_filter_18s$Kingdom=="Protalveolata" ], ]$cols <- "deeppink2"

## collaps OTU colors to prepare Kingdom colors
PHYLA_label_cols_18s <- tax_filter_18s[,c("Kingdom", "cols")]
library(plyr)
PHYLA_label_cols_18s <- ddply(PHYLA_label_cols_18s, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_18s) <- PHYLA_label_cols_18s[,1]
PHYLA_label_cols_18s <- PHYLA_label_cols_18s[phyla_18s, ]
PHYLA_label_cols_18s

## Legend for Kingdom colors
PHYLA_label_cols_18s_legend <- PHYLA_label_cols_18s[, ]
# PHYLA_label_cols_18s_legend[5, 1] <- "other"
# rownames(PHYLA_label_cols_18s_legend)[5] <- "other"
PHYLA_label_cols_18s_legend




##### Plot Supplementary Figure S3
pdf(paste0(output,"FigureS3.pdf"), encoding="MacRoman", height=6, width=10, paper="a4r")

layout(matrix(c(1,2),1,2, byrow=F))

par(oma=c(0,0,0,0), mar=c(6,4,1,5), xpd=NA)
phylum_bar_16s <- barplot(as.matrix(PHYLUM_mat_16s), col=PHYLA_label_cols_16s[rownames(PHYLUM_mat_16s),]$cols,
                          ylim=c(0,100), xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_16s, labels=design_16s$Treatments, col.axis="black", las=2, cex.axis=0.6)
title(ylab="Relative abundance (%)")
title(main="Bacteria Community")
legend(38, 100, bty="n", cex=0.7, x.intersp=0.1, y.intersp=1,
       legend=rev(PHYLA_label_cols_16s_legend$labels), 
       fill=rev(PHYLA_label_cols_16s_legend$cols), 
       border=rev(PHYLA_label_cols_16s_legend$cols) )

par(mar=c(6,4,1,5), xpd=NA)
phylum_bar_its <- barplot(as.matrix(PHYLUM_mat_its),col=PHYLA_label_cols_its[rownames(PHYLUM_mat_its),]$cols,
                          ylim=c(0,100), xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_its, labels=design_its$Treatments, col.axis="black", las=2, cex.axis=0.6)
title(main="Fungi Community")
legend(38, 100, bty="n", cex=0.7, x.intersp=0.1, y.intersp=1, 
       legend=rev(PHYLA_label_cols_its_legend$Phylum),
       fill=rev(PHYLA_label_cols_its_legend$cols),
       border=rev(PHYLA_label_cols_its_legend$cols) )
	   
par(oma=c(0,0,0,0), mar=c(6,4,1,5), xpd=NA)
phylum_bar_18s <- barplot(as.matrix(PHYLUM_mat_18s), col=PHYLA_label_cols_18s[rownames(PHYLUM_mat_18s),]$cols,
                          ylim=c(0,100), xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_18s, labels=design_18s$Treatments, col.axis="black", las=2, cex.axis=0.6)
title(ylab="Relative abundance (%)")
title(main="Protist Community")
legend(38, 100, bty="n", cex=0.7, x.intersp=0.1, y.intersp=1,
       legend=rev(PHYLA_label_cols_18s_legend$Kingdom), 
       fill=rev(PHYLA_label_cols_18s_legend$cols), 
       border=rev(PHYLA_label_cols_18s_legend$cols))

dev.off()





#######################################################################################################
#####  MODEL    Defining  aN bacteria,fungi and protist communities for beta diversity analysis   #####
#######################################################################################################

## Apply sequence count threshold to bacteria aN community and TMM normalize counts
otu_16s_aN <- otu_filter_16s[, aNsamples ]
otu_16s_aN <- otu_16s_aN[rowSums(otu_16s_aN) > 0,]
dim(otu_16s_aN)

## Apply threshold of at least 8 sequences in at least 4 samples
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

## apply threshold of at least 8 sequences in at least 4 samples
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


## Apply sequence count threshold to protist aN community and TMM normalize counts
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

## Input TMM normalized counts, taxonomy, and design of bulk aN protist community into phyloseq objects
## for further analysis
phy_18s_aN <- otu_table(otu_norm_aN_18s,taxa_are_rows=T)
phy_tax_aN_18s <-tax_table(as.matrix(tax_aN_18s))
phy_design_aN_18s <- sample_data(design_18s_aN)
physeq_aN_norm_18s <- phyloseq(phy_18s_aN,phy_tax_aN_18s,phy_design_aN_18s)
sample_data(physeq_aN_norm_18s)$CL_P <- factor(sample_data(physeq_aN_norm_18s)$CL_P,levels=c("CLaP", "CLeP","uCaP","uCeP"))




##### Defining eN bacteria, fungi and protist communities for beta diversity analysis #####
## Apply sequence count threshold to bacteria eN community and TMM normalize counts
## This defines the "eN-associated" bacteria community
otu_16s_eN <- otu_filter_16s[,eNsamples]
otu_16s_eN <- otu_16s_eN[rowSums(otu_16s_eN) > 0,]
dim(otu_16s_eN)

## Apply threshold of at least 8 sequences in at least 4 samples
keep_OTUs_eN_16s <- which(rowSums(otu_16s_eN >= 8) >= 4)
otu_16s_eN <- otu_16s_eN[keep_OTUs_eN_16s,]

## Number of OTUs in the eN-associated bacteria community
dim(otu_16s_eN)

tax_eN_16s <- tax_filter_16s[rownames(otu_16s_eN),]
design_16s_eN <- droplevels(design_filter_16s[eNsamples,])

edgeR_16s_eN <- DGEList(counts=otu_16s_eN, group=design_16s_eN$CL_P, genes=tax_eN_16s)
edgeR_16s_eN <- calcNormFactors(edgeR_16s_eN)

## Get TMM normalized counts expressed as relative abundance counts per million
otu_norm_eN_16s <- cpm(edgeR_16s_eN, normalized.lib.sizes=T, log=F)

## Create phyloseq objects for later analysis
phy_16s_eN <- otu_table(otu_norm_eN_16s,taxa_are_rows=T)
phy_tax_eN_16s <- tax_table(as.matrix(tax_eN_16s))
phy_design_eN_16s <- sample_data(design_16s_eN)
physeq_eN_norm_16s <- phyloseq(phy_16s_eN, phy_tax_eN_16s, phy_design_eN_16s)
sample_data(physeq_eN_norm_16s)$CL_P <- factor(sample_data(physeq_eN_norm_16s)$CL_P, levels=c("CLaP", "CLeP","uCaP","uCeP"))


## Apply sequence count threshold to fungi community and TMM normalize counts
## This defines the "eN-associated" bacteria community
otu_its_eN <- otu_filter_its[,eNsamples]
otu_its_eN <- otu_its_eN[rowSums(otu_its_eN) > 0,]

## apply threshold of at least 8 sequences in at least 4 samples
keep_OTUs_eN_its <- which(rowSums(otu_its_eN >= 8) >= 4)

otu_its_eN <- otu_its_eN[keep_OTUs_eN_its,]

## Number of OTUs in the eN-associated fungi community
dim(otu_its_eN)
tax_eN_its <- tax_filter_its[rownames(otu_its_eN),]
design_its_eN <- droplevels(design_filter_its[eNsamples,])

edgeR_its_eN <- DGEList(counts=otu_its_eN, group=design_its_eN$CL_P ,genes=tax_eN_its)
edgeR_its_eN <- calcNormFactors(edgeR_its_eN)

## Get TMM normalized counts expressed as relative abundance counts per million
otu_norm_eN_its <- cpm(edgeR_its_eN,normalized.lib.sizes=T,log=F)

## Create phyloseq objects for later anaylsis
phy_its_eN <- otu_table(otu_norm_eN_its,taxa_are_rows=T)
phy_tax_eN_its <-tax_table(as.matrix(tax_eN_its))
phy_design_eN_its <- sample_data(design_its_eN)
physeq_eN_norm_its <- phyloseq(phy_its_eN,phy_tax_eN_its,phy_design_eN_its)
sample_data(physeq_eN_norm_its)$CL_P <- factor(sample_data(physeq_eN_norm_its)$CL_P,levels=c("CLaP", "CLeP","uCaP","uCeP"))


## Apply sequence count threshold to protist eN community and TMM normalize counts
## This defines the "eN-associated" protist community
otu_18s_eN <- otu_filter_18s[,eNsamples]
otu_18s_eN <- otu_18s_eN[rowSums(otu_18s_eN) > 0,]
dim(otu_18s_eN)

## Apply threshold of at least 2 sequences in at least 4 samples
keep_OTUs_eN_18s <- which(rowSums(otu_18s_eN >= 8) >= 4)

otu_18s_eN <- otu_18s_eN[keep_OTUs_eN_18s,]

## Number of OTUs in the eN-associated protist community
dim(otu_18s_eN)

tax_eN_18s <- tax_filter_18s[rownames(otu_18s_eN),]
design_18s_eN <- droplevels(design_filter_18s[eNsamples,])

edgeR_18s_eN <- DGEList(counts=otu_18s_eN, group=design_18s_eN$CL_P, genes=tax_eN_18s)
edgeR_18s_eN <- calcNormFactors(edgeR_18s_eN)

## Get TMM normalized counts expressed as relative abundance counts per million
otu_norm_eN_18s <- cpm(edgeR_18s_eN, normalized.lib.sizes=T, log=F)

## Create phyloseq objects for later analysis
phy_18s_eN <- otu_table(otu_norm_eN_18s,taxa_are_rows=T)
phy_tax_eN_18s <- tax_table(as.matrix(tax_eN_18s))
phy_design_eN_18s <- sample_data(design_18s_eN)
physeq_eN_norm_18s <- phyloseq(phy_18s_eN, phy_tax_eN_18s, phy_design_eN_18s)
sample_data(physeq_eN_norm_18s)$CL_P <- factor(sample_data(physeq_eN_norm_18s)$CL_P, levels=c("CLaP", "CLeP","uCaP","uCeP"))






##########################################################################################################
##### MODEL 4 Identifiying Nitrogen addition system responsive OTUs with indicator species analysis ######
##########################################################################################################
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

CLaP_aN_16s <- as.matrix(indic_aN_df_16s[which(indic_aN_df_16s$s.CLaP == 1 & indic_aN_df_16s$p.value < 0.05),])
CLeP_aN_16s <- as.matrix(indic_aN_df_16s[which(indic_aN_df_16s$s.CLeP == 1 & indic_aN_df_16s$p.value < 0.05),])
uCaP_aN_16s <- as.matrix(indic_aN_df_16s[which(indic_aN_df_16s$s.uCaP == 1 & indic_aN_df_16s$p.value < 0.05),])
uCeP_aN_16s <- as.matrix(indic_aN_df_16s[which(indic_aN_df_16s$s.uCeP == 1 & indic_aN_df_16s$p.value < 0.05),])

aN_r_values_16s <- rbind(CLaP_aN_16s,CLeP_aN_16s,uCaP_aN_16s,uCeP_aN_16s)
colnames(aN_r_values_16s)[1:4] <-c("CLaP","CLeP","uCaP","uCeP")

## Range of correlation coefficients
range(aN_r_values_16s[,"stat"])

## Total number of indicator OTUS
length(unique(rownames(aN_r_values_16s)))

## Proportion of aN bacteria OTUs responding to Nitrogen addition system
length(unique(rownames(aN_r_values_16s)))/nrow(otu_norm_aN_16s)

## Proportion of aN bacteria sequences responding to Nitrogen addition system
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

## Proportion of aN bacteria OTUs responding to Nitrogen addition system
length(unique(rownames(aN_r_values_its)))/nrow(otu_norm_aN_its)

## Proportion of aN bacteria sequences responding to Nitrogen addition system
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

## Proportion of aN bacteria OTUs responding to Nitrogen addition system
length(unique(rownames(aN_r_values_18s)))/nrow(otu_norm_aN_18s)

## Proportion of aN bacteria sequences responding to Nitrogen addition system
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

## Proportion of eN bacteria OTUs responding to Nitrogen addition system
length(unique(rownames(eN_r_values_16s)))/nrow(otu_norm_eN_16s)

## Proportion of eN bacteria sequences responding to Nitrogen addition system
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

## Proportion of eN bacteria OTUs responding to Nitrogen addition system
length(unique(rownames(eN_r_values_its)))/nrow(otu_norm_eN_its)

## Proportion of eN bacteria sequences responding to Nitrogen addition system
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

## Proportion of eN bacteria OTUs responding to Nitrogen addition system
length(unique(rownames(eN_r_values_18s)))/nrow(otu_norm_eN_18s)

## Proportion of eN bacteria sequences responding to Nitrogen addition system
eN_18s_ra <- t(t(otu_18s_eN)/colSums(otu_18s_eN)) * 100
sum(colSums(eN_18s_ra[unique(rownames(eN_r_values_18s)),]))/sum(colSums(eN_18s_ra))    #Error in h(simpleError(msg, call)) : 
                                                                                       #在为'colSums'函数选择方法时评估'x'参数出了错: subscript out of bounds 







################################ edgeR ##############################################
##### MODEL 5  Nitrogen addition system responsive OTUs with edgeR              #####
#####################################################################################

## Get Nitrogen addition system responsive bulk aN bacteria OTUs with likelihood ratio testing in edgeR
model_mataN_16s <- model.matrix(~CL_P, data=design_16s_aN)
edgeR_16s_aN_CL_P <- DGEList(counts=otu_16s_aN, group=design_16s_aN$CL_P, genes=tax_aN_16s)
edgeR_16s_aN_CL_P <- calcNormFactors(edgeR_16s_aN_CL_P)

dge_aNCL_P_16s <- estimateGLMRobustDisp(edgeR_16s_aN_CL_P, design=model_mataN_16s)

fit_aNCL_P_16s <- glmFit(dge_aNCL_P_16s, design = model_mataN_16s)
lrt_aNCL_P_16s <- glmLRT(fit_aNCL_P_16s, coef=2:4)
CL_P_aN_16s <- topTags(lrt_aNCL_P_16s, n=Inf, p.value=0.05)
CL_P_aN_16s <- CL_P_aN_16s$table


## Get Nitrogen addition system responsive bulk aN fungi OTUs with likelihood ratio testing in edgeR
model_mataN_its <- model.matrix(~CL_P, data=design_its_aN)
edgeR_its_aN_CL_P <- DGEList(counts=otu_its_aN, group=design_its_aN$CL_P, genes=tax_aN_its)
edgeR_its_aN_CL_P <- calcNormFactors(edgeR_its_aN_CL_P)

dge_aNCL_P_its <- estimateGLMRobustDisp(edgeR_its_aN_CL_P, design=model_mataN_its)

fit_aNCL_P_its <- glmFit(dge_aNCL_P_its, design=model_mataN_its)
lrt_aNCL_P_its <- glmLRT(fit_aNCL_P_its, coef=2:4)
CL_P_aN_its <- topTags(lrt_aNCL_P_its, n=Inf, p.value=0.05)
CL_P_aN_its <- CL_P_aN_its$table


## Get Nitrogen addition system responsive bulk aN protist OTUs with likelihood ratio testing in edgeR
model_mataN_18s <- model.matrix(~CL_P, data=design_18s_aN)
edgeR_18s_aN_CL_P <- DGEList(counts=otu_18s_aN, group=design_18s_aN$CL_P, genes=tax_aN_18s)
edgeR_18s_aN_CL_P <- calcNormFactors(edgeR_18s_aN_CL_P)

dge_aNCL_P_18s <- estimateGLMRobustDisp(edgeR_18s_aN_CL_P, design=model_mataN_18s)

fit_aNCL_P_18s <- glmFit(dge_aNCL_P_18s, design = model_mataN_18s)
lrt_aNCL_P_18s <- glmLRT(fit_aNCL_P_18s, coef=2:4)
CL_P_aN_18s <- topTags(lrt_aNCL_P_18s, n=Inf, p.value=0.05)
CL_P_aN_18s <- CL_P_aN_18s$table




## Get Nitrogen addition system responsive eN-associated bacteria OTUs with likelihood ratio testing in edgeR
model_mateN_16s <- model.matrix(~CL_P, data=design_16s_eN)
edgeR_16s_eN_CL_P <- DGEList(counts=otu_16s_eN, group=design_16s_eN$CL_P, genes=tax_eN_16s)
edgeR_16s_eN_CL_P <- calcNormFactors(edgeR_16s_eN_CL_P)

dge_eNCL_P_16s <- estimateGLMRobustDisp(edgeR_16s_eN_CL_P, design=model_mateN_16s)

fit_eNCL_P_16s <- glmFit(dge_eNCL_P_16s, design=model_mateN_16s)
lrt_eNCL_P_16s <- glmLRT(fit_eNCL_P_16s, coef=2:4)
CL_P_eN_16s <- topTags(lrt_eNCL_P_16s, n=Inf, p.value=0.05)
CL_P_eN_16s <- CL_P_eN_16s$table

## Get Nitrogen addition system responsive eN-associated fungi OTUs with likelihood ratio testing in edgeR
model_mateN_its <- model.matrix(~CL_P, data=design_its_eN)
edgeR_its_eN_CL_P <- DGEList(counts=otu_its_eN, group=design_its_eN$CL_P, genes=tax_eN_its)
edgeR_its_eN_CL_P <- calcNormFactors(edgeR_its_eN_CL_P)

dge_eNCL_P_its <- estimateGLMRobustDisp(edgeR_its_eN_CL_P, design=model_mateN_its)

fit_eNCL_P_its <- glmFit(dge_eNCL_P_its, design=model_mateN_its)
lrt_eNCL_P_its <- glmLRT(fit_eNCL_P_its, coef=2:4)
CL_P_eN_its <- topTags(lrt_eNCL_P_its, n=Inf, p.value=0.05)
CL_P_eN_its <- CL_P_eN_its$table


## Get Nitrogen addition system responsive eN-associated protist OTUs with likelihood ratio testing in edgeR
model_mateN_18s <- model.matrix(~CL_P, data=design_18s_eN)
edgeR_18s_eN_CL_P <- DGEList(counts=otu_18s_eN, group=design_18s_eN$CL_P, genes=tax_eN_18s)
edgeR_18s_eN_CL_P <- calcNormFactors(edgeR_18s_eN_CL_P)

dge_eNCL_P_18s <- estimateGLMRobustDisp(edgeR_18s_eN_CL_P, design=model_mateN_18s)

fit_eNCL_P_18s <- glmFit(dge_eNCL_P_18s, design=model_mateN_18s)
lrt_eNCL_P_18s <- glmLRT(fit_eNCL_P_18s, coef=2:4)
CL_P_eN_18s <- topTags(lrt_eNCL_P_18s, n=Inf, p.value=0.05)
CL_P_eN_18s <- CL_P_eN_18s$table

##########################################################################################################################
####################################################   sensitive OTUS ####################################################

## Define Nitrogen addition senstive bulk aN bacteria OTUs (Nitrogen addition response validated by both statisical methods)
## This is important for highlighting their position in the subsequent figures
indic_edge_16s_aN <- intersect(rownames(aN_r_values_16s), rownames(CL_P_aN_16s))    # aN sensitive OTUs, aN_r_values_16s <- rbind(CLaP_aN_16s,CLeP_aN_16s,uCaP_aN_16s,uCeP_aN_16s)
write.csv(indic_edge_16s_aN, "indic_edge_16s_aN.csv")      # 
indic_edge_16s_CLaP_aN <- intersect(rownames(CLaP_aN_16s), rownames(CL_P_aN_16s))    # CLaP sensitive OTUs under aN treatment level.
write.csv(indic_edge_16s_CLaP_aN, "indic_edge_16s_CLaP_aN.csv")
indic_edge_16s_CLeP_aN <- intersect(rownames(CLeP_aN_16s), rownames(CL_P_aN_16s))    # CLeP sensitive OTUs under aN treatment level.
write.csv(indic_edge_16s_CLeP_aN, "indic_edge_16s_CLeP_aN.csv")
indic_edge_16s_uCaP_aN <- intersect(rownames(uCaP_aN_16s), rownames(CL_P_aN_16s))    # uCaP sensitive OTUs under aN treatment level.
write.csv(indic_edge_16s_uCaP_aN, "indic_edge_16s_uCaP_aN.csv")
indic_edge_16s_uCeP_aN <- intersect(rownames(uCeP_aN_16s), rownames(CL_P_aN_16s))    # uCeP sensitive OTUs under aN treatment level.
write.csv(indic_edge_16s_uCeP_aN, "indic_edge_16s_uCeP_aN.csv")


## Define Nitrogen addition senstive bulk aN fungi OTUs (Nitrogen addition response validated by both statisical methods)
## This is important for highlighting their position in the subsequent figures
indic_edge_its_aN <- intersect(rownames(aN_r_values_its),rownames(CL_P_aN_its))
write.csv(indic_edge_its_aN, "indic_edge_its_aN.csv")      # 
indic_edge_its_CLaP_aN <- intersect(rownames(CLaP_aN_its), rownames(CL_P_aN_its))    # CLaP sensitive OTUs under aN treatment level.
write.csv(indic_edge_its_CLaP_aN, "indic_edge_its_CLaP_aN.csv")
indic_edge_its_CLeP_aN <- intersect(rownames(CLeP_aN_its), rownames(CL_P_aN_its))    # CLeP sensitive OTUs under aN treatment level.
write.csv(indic_edge_its_CLeP_aN, "indic_edge_its_CLeP_aN.csv")
indic_edge_its_uCaP_aN <- intersect(rownames(uCaP_aN_its), rownames(CL_P_aN_its))    # uCaP sensitive OTUs under aN treatment level.
write.csv(indic_edge_its_uCaP_aN, "indic_edge_its_uCaP_aN.csv")
indic_edge_its_uCeP_aN <- intersect(rownames(uCeP_aN_its), rownames(CL_P_aN_its))    # uCeP sensitive OTUs under aN treatment level.
write.csv(indic_edge_its_uCeP_aN, "indic_edge_its_uCeP_aN.csv")

## Define Nitrogen addition senstive bulk aN protist OTUs (Nitrogen addition response validated by both statisical methods)
## This is important for highlighting their position in the subsequent figures
indic_edge_18s_aN <- intersect(rownames(aN_r_values_18s), rownames(CL_P_aN_18s))
write.csv(indic_edge_18s_aN, "indic_edge_18s_aN.csv")      # 
indic_edge_18s_CLaP_aN <- intersect(rownames(CLaP_aN_18s), rownames(CL_P_aN_18s))    # CLaP sensitive OTUs under aN treatment level.
write.csv(indic_edge_18s_CLaP_aN, "indic_edge_18s_CLaP_aN.csv")
indic_edge_18s_CLeP_aN <- intersect(rownames(CLeP_aN_18s), rownames(CL_P_aN_18s))    # CLeP sensitive OTUs under aN treatment level.
write.csv(indic_edge_18s_CLeP_aN, "indic_edge_18s_CLeP_aN.csv")
indic_edge_18s_uCaP_aN <- intersect(rownames(uCaP_aN_18s), rownames(CL_P_aN_18s))    # uCaP sensitive OTUs under aN treatment level.
write.csv(indic_edge_18s_uCaP_aN, "indic_edge_18s_uCaP_aN.csv")
indic_edge_18s_uCeP_aN <- intersect(rownames(uCeP_aN_18s), rownames(CL_P_aN_18s))    # uCeP sensitive OTUs under aN treatment level.
write.csv(indic_edge_18s_uCeP_aN, "indic_edge_18s_uCeP_aN.csv")

## Define Nitrogen addition senstive bulk eN bacteria OTUs (Nitrogen addition response validated by both statisical methods)
## This is important for highlighting their position in the subsequent figures
indic_edge_16s_eN <- intersect(rownames(eN_r_values_16s), rownames(CL_P_eN_16s))
write.csv(indic_edge_16s_eN, "indic_edge_16s_eN.csv")      # 
indic_edge_16s_CLaP_eN <- intersect(rownames(CLaP_eN_16s), rownames(CL_P_eN_16s))    # CLaP sensitive OTUs under eN treatment level.
write.csv(indic_edge_16s_CLaP_eN, "indic_edge_16s_CLaP_eN.csv")
indic_edge_16s_CLeP_eN <- intersect(rownames(CLeP_eN_16s), rownames(CL_P_eN_16s))    # CLeP sensitive OTUs under eN treatment level.
write.csv(indic_edge_16s_CLeP_eN, "indic_edge_16s_CLeP_eN.csv")
indic_edge_16s_uCaP_eN <- intersect(rownames(uCaP_eN_16s), rownames(CL_P_eN_16s))    # uCaP sensitive OTUs under eN treatment level.
write.csv(indic_edge_16s_uCaP_eN, "indic_edge_16s_uCaP_eN.csv")
indic_edge_16s_uCeP_eN <- intersect(rownames(uCeP_eN_16s), rownames(CL_P_eN_16s))    # uCeP sensitive OTUs under eN treatment level.
write.csv(indic_edge_16s_uCeP_eN, "indic_edge_16s_uCeP_eN.csv")

## Define Nitrogen addition senstive bulk eN fungi OTUs (Nitrogen addition response validated by both statisical methods)
## This is important for highlighting their position in the subsequent figures
indic_edge_its_eN <- intersect(rownames(eN_r_values_its),rownames(CL_P_eN_its))
write.csv(indic_edge_its_eN, "indic_edge_its_eN.csv")      # 
indic_edge_its_CLaP_eN <- intersect(rownames(CLaP_eN_its), rownames(CL_P_eN_its))    # CLaP sensitive OTUs under eN treatment level.
write.csv(indic_edge_its_CLaP_eN, "indic_edge_its_CLaP_eN.csv")
indic_edge_its_CLeP_eN <- intersect(rownames(CLeP_eN_its), rownames(CL_P_eN_its))    # CLeP sensitive OTUs under eN treatment level.
write.csv(indic_edge_its_CLeP_eN, "indic_edge_its_CLeP_eN.csv")
indic_edge_its_uCaP_eN <- intersect(rownames(uCaP_eN_its), rownames(CL_P_eN_its))    # uCaP sensitive OTUs under eN treatment level.
write.csv(indic_edge_its_uCaP_eN, "indic_edge_its_uCaP_eN.csv")
indic_edge_its_uCeP_eN <- intersect(rownames(uCeP_eN_its), rownames(CL_P_eN_its))    # uCeP sensitive OTUs under eN treatment level.
write.csv(indic_edge_its_uCeP_eN, "indic_edge_its_uCeP_eN.csv")
## Define Nitrogen addition senstive bulk eN protist OTUs (Nitrogen addition response validated by both statisical methods)
## This is important for highlighting their position in the subsequent figures
indic_edge_18s_eN <- intersect(rownames(eN_r_values_18s), rownames(CL_P_eN_18s))
write.csv(indic_edge_18s_eN, "indic_edge_18s_eN.csv")      # 
indic_edge_18s_CLaP_eN <- intersect(rownames(CLaP_eN_18s), rownames(CL_P_eN_18s))    # CLaP sensitive OTUs under eN treatment level.
write.csv(indic_edge_18s_CLaP_eN, "indic_edge_18s_CLaP_eN.csv")
indic_edge_18s_CLeP_eN <- intersect(rownames(CLeP_eN_18s), rownames(CL_P_eN_18s))    # CLeP sensitive OTUs under eN treatment level.
write.csv(indic_edge_18s_CLeP_eN, "indic_edge_18s_CLeP_eN.csv")
indic_edge_18s_uCaP_eN <- intersect(rownames(uCaP_eN_18s), rownames(CL_P_eN_18s))    # uCaP sensitive OTUs under eN treatment level.
write.csv(indic_edge_18s_uCaP_eN, "indic_edge_18s_uCaP_eN.csv")
indic_edge_18s_uCeP_eN <- intersect(rownames(uCeP_eN_18s), rownames(CL_P_eN_18s))    # uCeP sensitive OTUs under eN treatment level.
write.csv(indic_edge_18s_uCeP_eN, "indic_edge_18s_uCeP_eN.csv")


####################################### aN meta and keystone OTUs  ##################################
##### MODEL 6  aN meta co-occurrence network creation and analysis and defining keystone OTUs  ######
#####################################################################################################
## Combine OTU counts of both kingdoms together
otu_norm_aN_combine <- rbind(otu_norm_aN_16s, otu_norm_aN_its, otu_norm_aN_18s)     # note the sample order must be same

## Perform Spearman correlation of all OTU pairs
all_aN_cor <- rcorr(t(otu_norm_aN_combine), type=c("spearman"))

## Create data frame of co-occurring OTUs
all_cor_aN_df <- CorrDF(all_aN_cor$r, all_aN_cor$P)
all_cor_aN_df$padj <- p.adjust(all_cor_aN_df$p, method="none")

all_cor_aN_df_padj <- all_cor_aN_df[which(all_cor_aN_df$cor > 0.7),]
all_cor_aN_df_padj <- all_cor_aN_df_padj[which(all_cor_aN_df_padj$padj < 0.001),]

## Make node attribute table
indic_edge_aN_combine <- c(indic_edge_16s_aN, indic_edge_its_aN, indic_edge_18s_aN)

aN_r_values_combine <- rbind(aN_r_values_16s, aN_r_values_its, aN_r_values_18s)

nodeattrib_aN_combine <- data.frame(node=union(all_cor_aN_df_padj$from,all_cor_aN_df_padj$to))
nodeattrib_aN_combine$indicgroup <- 0

for (i in as.character(nodeattrib_aN_combine$node))
{
  if (i %in% indic_edge_aN_combine == TRUE)
  {nodeattrib_aN_combine[nodeattrib_aN_combine$node==i,"indicgroup"] <- paste(colnames(aN_r_values_combine)[which(aN_r_values_combine[i,1:4]==1)],collapse = "_")}
  else
  {nodeattrib_aN_combine[nodeattrib_aN_combine$node==i,"indicgroup"]<- "NA"}
}

rownames(nodeattrib_aN_combine) <- as.character(nodeattrib_aN_combine$node)

all_aN_net <- graph_from_data_frame(all_cor_aN_df_padj,direct=F,vertices=nodeattrib_aN_combine)

## Number of nodes
length(V(all_aN_net))

## Number of bacteria and fungi nodes
length(grep("bOTU_*",names(V(all_aN_net))))
length(grep("fotu*",names(V(all_aN_net))))
length(grep("potu*",names(V(all_aN_net))))

## Connections 
bb_occur_aN <- droplevels(all_cor_aN_df_padj[with(all_cor_aN_df_padj, grepl("bOTU_*",from) & grepl("bOTU_*",to)),])
nrow(bb_occur_aN)

ff_occur_aN <- droplevels(all_cor_aN_df_padj[with(all_cor_aN_df_padj, grepl("fotu*",from) & grepl("fotu*",to)),])
nrow(ff_occur_aN)

pp_occur_aN <- droplevels(all_cor_aN_df_padj[with(all_cor_aN_df_padj, grepl("potu*",from) & grepl("potu*",to)),])
nrow(pp_occur_aN)

bf_occur_aN <- droplevels(all_cor_aN_df_padj[with(all_cor_aN_df_padj, grepl("fotu*",from) & grepl("bOTU_*",to)),])
nrow(bf_occur_aN)

bp_occur_aN <- droplevels(all_cor_aN_df_padj[with(all_cor_aN_df_padj, grepl("potu*",from) & grepl("bOTU_*",to)),])
nrow(bp_occur_aN)

fp_occur_aN <- droplevels(all_cor_aN_df_padj[with(all_cor_aN_df_padj, grepl("potu*",from) & grepl("fotu*",to)),])
nrow(fp_occur_aN)


## Node degree values
aN_all_deg <- sort(degree(all_aN_net,mode="all"),decr=T)
max(aN_all_deg)
mean(aN_all_deg)

## Define keystone species based on top 1% of node degree values
n <- 1
aN_all_keystone <- aN_all_deg[aN_all_deg > quantile(aN_all_deg,prob=1-n/100)]
length(aN_all_keystone)

## Get bacteria and fungi keystone OTUs
aN_net_keystone_bac <- names(aN_all_keystone)[grep("bOTU_*",names(aN_all_keystone))]
length(aN_net_keystone_bac)

aN_net_keystone_fun <- names(aN_all_keystone)[grep("fotu*",names(aN_all_keystone))]
length(aN_net_keystone_fun)

aN_net_keystone_pro <- names(aN_all_keystone)[grep("potu*",names(aN_all_keystone))]
length(aN_net_keystone_pro)

## get csOTUs present in network
cs_aN_net <- names(V(all_aN_net))[names(V(all_aN_net)) %in% indic_edge_aN_combine]

cs_aN_net_bac <- cs_aN_net[grep("bOTU_*",cs_aN_net)]
length(cs_aN_net_bac)

cs_aN_net_fun <- cs_aN_net[grep("fotu*",cs_aN_net)]
length(cs_aN_net_fun)

cs_aN_net_pro <- cs_aN_net[grep("potu*",cs_aN_net)]
length(cs_aN_net_pro)

## Keystone OTUs present in csOTUs
intersect(aN_net_keystone_bac,cs_aN_net_bac)        # keystone, sensitive OTU. 
intersect(aN_net_keystone_fun,cs_aN_net_fun)        
intersect(aN_net_keystone_pro,cs_aN_net_pro)

########################### keystone OTUs and sensitive OTUs ###############################################
write.csv(aN_net_keystone_bac, "aN_net_keystone_bac.csv")
write.csv(aN_net_keystone_fun, "aN_net_keystone_fun.csv")
write.csv(aN_net_keystone_pro, "aN_net_keystone_pro.csv")
write.csv(cs_aN_net_bac,"cs_aN_net_bac.csv")
write.csv(cs_aN_net_fun,"cs_aN_net_fun.csv")
write.csv(cs_aN_net_pro,"cs_aN_net_pro.csv")
############################################################################################################
## Set node shape
V(all_aN_net)$shape <- V(all_aN_net)$name
V(all_aN_net)$shape[V(all_aN_net)$shape %in% names(aN_all_keystone)] <- "star"
V(all_aN_net)$shape[V(all_aN_net)$shape %in% rownames(otu_norm_aN_16s)] <- "circle"
V(all_aN_net)$shape[V(all_aN_net)$shape %in% rownames(otu_norm_aN_its)] <- "triangle"
V(all_aN_net)$shape[V(all_aN_net)$shape %in% rownames(otu_norm_aN_18s)] <- "square"

cs <- c("CLaP","CLaP_CLeP","CLeP","CLaP_uCaP","CLeP_uCeP","uCaP","uCaP_uCeP","uCeP")

unique(V(all_aN_net)$indicgroup)
V(all_aN_net)$color <- V(all_aN_net)$indicgroup
V(all_aN_net)$color[!V(all_aN_net)$color %in% cs] <- "gray30"
V(all_aN_net)$color[V(all_aN_net)$color == "CLaP"] <- "dodgerblue4"
V(all_aN_net)$color[V(all_aN_net)$color == "CLaP_CLeP"] <- "dodgerblue3"
V(all_aN_net)$color[V(all_aN_net)$color == "CLeP"] <- "dodgerblue1"
V(all_aN_net)$color[V(all_aN_net)$color == "CLaP_uCaP"] <- "sienna4"
V(all_aN_net)$color[V(all_aN_net)$color == "CLeP_uCeP"] <- "sienna3"
V(all_aN_net)$color[V(all_aN_net)$color == "uCaP"] <- "firebrick4"
V(all_aN_net)$color[V(all_aN_net)$color == "uCaP_uCeP"] <- "firebrick3"
V(all_aN_net)$color[V(all_aN_net)$color == "uCeP"] <- "firebrick2"
V(all_aN_net)$frame.color <- V(all_aN_net)$color

aN_all_net_csnodes <- rownames(nodeattrib_aN_combine[nodeattrib_aN_combine$indicgroup %in% cs,])
aN_all_net_csnodes_bac <- aN_all_net_csnodes[grep("bOTU*",aN_all_net_csnodes)]
aN_all_net_csnodes_fun <- aN_all_net_csnodes[grep("fotu*",aN_all_net_csnodes)]
aN_all_net_csnodes_pro <- aN_all_net_csnodes[grep("potu*",aN_all_net_csnodes)]

V(all_aN_net)$size <- V(all_aN_net)$name
V(all_aN_net)$size[V(all_aN_net)$size %in% names(aN_all_keystone)] <- 4
V(all_aN_net)$size[V(all_aN_net)$size %in% aN_all_net_csnodes_bac] <- 4
V(all_aN_net)$size[V(all_aN_net)$size %in% aN_all_net_csnodes_fun] <- 6
V(all_aN_net)$size[V(all_aN_net)$size %in% aN_all_net_csnodes_pro] <- 6

V(all_aN_net)$size[V(all_aN_net)$size %in% rownames(otu_norm_aN_16s)] <- 2
V(all_aN_net)$size[V(all_aN_net)$size %in% rownames(otu_norm_aN_its)] <- 3
V(all_aN_net)$size[V(all_aN_net)$size %in% rownames(otu_norm_aN_18s)] <- 3
aNcombine_nodesizes <- as.numeric(V(all_aN_net)$size)



####################################### modules ###########################################
##### MODEL 7 Explore community structure of aN meta-network by defining modules      #####
###########################################################################################
## Make vectors of network nodes responding to different Nitrogen addition systems
CIT_nodes_aN <- rownames(nodeattrib_aN_combine[nodeattrib_aN_combine$indicgroup=="CLaP",])
CNT_nodes_aN <- rownames(nodeattrib_aN_combine[nodeattrib_aN_combine$indicgroup=="CLeP",])
CON_nodes_aN <- rownames(nodeattrib_aN_combine[nodeattrib_aN_combine$indicgroup=="CLaP_CLeP",])
IT_nodes_aN <- rownames(nodeattrib_aN_combine[nodeattrib_aN_combine$indicgroup=="CLaP_uCaP",])
OIT_nodes_aN <- rownames(nodeattrib_aN_combine[nodeattrib_aN_combine$indicgroup=="uCaP",])
ORT_nodes_aN <- rownames(nodeattrib_aN_combine[nodeattrib_aN_combine$indicgroup=="uCeP",])
ORG_nodes_aN <- rownames(nodeattrib_aN_combine[nodeattrib_aN_combine$indicgroup=="uCaP_uCeP",])
NTRT_nodes_aN <- rownames(nodeattrib_aN_combine[nodeattrib_aN_combine$indicgroup=="CLeP_uCeP",])

cs_nodes_aN_all <- c(CIT_nodes_aN,CNT_nodes_aN,CON_nodes_aN,IT_nodes_aN,OIT_nodes_aN,ORT_nodes_aN,ORG_nodes_aN,NTRT_nodes_aN)
Bcs_nodes_aN_all <- cs_nodes_aN_all[grep("bOTU*",cs_nodes_aN_all)]
Fcs_nodes_aN_all <- cs_nodes_aN_all[grep("fotu*",cs_nodes_aN_all)]
Pcs_nodes_aN_all <- cs_nodes_aN_all[grep("potu*",cs_nodes_aN_all)]

## Perform cluster analysis using greedy clustering algorithm 
cfg_aN <- cluster_fast_greedy(as.undirected(all_aN_net))

## Subset for top 20 biggest nodes
aN_modules <- sort(table(membership(cfg_aN)),decr=T)
aN_modules_20 <- aN_modules[1:20]

sum(aN_modules_20)/sum(aN_modules)
sm20_plot <- aN_modules_20
names(sm20_plot) <- as.factor(1:20)
aN_modules_cs <- table(factor(membership(cfg_aN)[cs_nodes_aN_all],levels=names(aN_modules)))
aN_modules_cs_20 <- aN_modules_cs[names(aN_modules_20)]
smcs20_plot <- aN_modules_cs_20
names(smcs20_plot) <- as.factor(1:20)

## Make vector of nodes in top 20 modules
aN_modules_points <- membership(cfg_aN)[membership(cfg_aN) %in% names(aN_modules_20)]

aN_points <- NULL
for(i in aN_modules_points){
  aNx <- which(names(aN_modules_20)==i)
  aN_points <- c(aN_points, aNx)
}
names(aN_points) <- names(aN_modules_points)

## Set node colors by Nitrogen addition system sensitivity 
aN_all_cols <- sort(aN_points)
aN_all_cols[!names(aN_all_cols) %in% cs] <- "gray30"
aN_all_cols[names(aN_all_cols) %in% CIT_nodes_aN] <- "dodgerblue4"
aN_all_cols[names(aN_all_cols) %in% CON_nodes_aN] <- "dodgerblue3"
aN_all_cols[names(aN_all_cols) %in% CNT_nodes_aN] <- "dodgerblue1"
aN_all_cols[names(aN_all_cols) %in% IT_nodes_aN] <- "sienna4"
aN_all_cols[names(aN_all_cols) %in% NTRT_nodes_aN]<- "sienna3"
aN_all_cols[names(aN_all_cols) %in% OIT_nodes_aN]<- "firebrick4"
aN_all_cols[names(aN_all_cols) %in% ORG_nodes_aN] <- "firebrick3"
aN_all_cols[names(aN_all_cols) %in% ORT_nodes_aN] <- "firebrick2"

aN_all_pch <- sort(aN_points)
aN_all_pch[names(aN_all_pch) %in% rownames(otu_norm_aN_16s)] <- 1
aN_all_pch[names(aN_all_pch) %in% rownames(otu_norm_aN_its)] <- 2
aN_all_pch[names(aN_all_pch) %in% rownames(otu_norm_aN_18s)] <- 2
aN_all_pch[names(aN_all_pch) %in% intersect(rownames(otu_norm_aN_16s),cs_nodes_aN_all)] <- 16
aN_all_pch[names(aN_all_pch) %in% intersect(rownames(otu_norm_aN_its),cs_nodes_aN_all)] <- 17
aN_all_pch[names(aN_all_pch) %in% intersect(rownames(otu_norm_aN_18s),cs_nodes_aN_all)] <- 17
aN_all_pch[names(aN_all_pch) %in% names(aN_all_keystone)] <- 8

aN_all_cex <- sort(aN_points)
aN_all_cex[!names(aN_all_cex) %in% cs_nodes_aN_all] <- 1
aN_all_cex[names(aN_all_cex) %in% cs_nodes_aN_all] <- 2

aN_mods_list_cs <- list()
for (i in names(aN_modules_cs_20)){
  x1 <- names(membership(cfg_aN)[membership(cfg_aN)==i])
  x2 <- x1[x1 %in% cs_nodes_aN_all]
  aN_mods_list_cs[[i]] <- as.numeric(V(all_aN_net)[x2])
}





#########################################################################################################
##### MODEL 8 eN meta co-occurrence network creation and analysis and defining keystone OTUs        #####
#########################################################################################################
## Combine OTU tables from three kingdoms together
otu_norm_eN_combine <- rbind(otu_norm_eN_16s, otu_norm_eN_its, otu_norm_eN_18s)

## Perform Spearman correlation of all OTU pairs
all_eN_cor <- rcorr(t(otu_norm_eN_combine), type=c("spearman"))

## Create data frame of co-occurring OTUs
all_cor_eN_df <- CorrDF(all_eN_cor$r, all_eN_cor$P)
all_cor_eN_df$padj <- p.adjust(all_cor_eN_df$p, method="none")

all_cor_eN_df_padj <- all_cor_eN_df[which(all_cor_eN_df$cor > 0.7),]
all_cor_eN_df_padj <- all_cor_eN_df_padj[which(all_cor_eN_df_padj$padj < 0.001),]

## Make node attribute table
indic_edge_eN_combine <- c(indic_edge_16s_eN, indic_edge_its_eN, indic_edge_18s_eN)

eN_r_values_combine <- rbind(eN_r_values_16s, eN_r_values_its, eN_r_values_18s)

nodeattrib_eN_combine <- data.frame(node=union(all_cor_eN_df_padj$from, all_cor_eN_df_padj$to))
nodeattrib_eN_combine$indicgroup <- 0

for (i in as.character(nodeattrib_eN_combine$node))
{
  if (i %in% indic_edge_eN_combine == TRUE)
  {nodeattrib_eN_combine[nodeattrib_eN_combine$node==i,"indicgroup"] <- paste(colnames(eN_r_values_combine)[which(eN_r_values_combine[i,1:4]==1)],collapse = "_")}
  else
  {nodeattrib_eN_combine[nodeattrib_eN_combine$node==i,"indicgroup"]<- "NA"}
}

rownames(nodeattrib_eN_combine) <- as.character(nodeattrib_eN_combine$node)

all_eN_net <- graph_from_data_frame(all_cor_eN_df_padj, direct=F, vertices=nodeattrib_eN_combine)

## Number of nodes
length(V(all_eN_net))

## Number of bacteria and fungi nodes
length(grep("bOTU_*",names(V(all_eN_net))))
length(grep("fotu*",names(V(all_eN_net))))
length(grep("potu*",names(V(all_eN_net))))

## Connections 
bb_occur_eN <- droplevels(all_cor_eN_df_padj[with(all_cor_eN_df_padj, grepl("bOTU_*",from) & grepl("bOTU_*",to)),])
nrow(bb_occur_eN)

ff_occur_eN <- droplevels(all_cor_eN_df_padj[with(all_cor_eN_df_padj, grepl("fotu*",from) & grepl("fotu*",to)),])
nrow(ff_occur_eN)

pp_occur_eN <- droplevels(all_cor_eN_df_padj[with(all_cor_eN_df_padj, grepl("potu*",from) & grepl("potu*",to)),])
nrow(pp_occur_eN)

bf_occur_eN <- droplevels(all_cor_eN_df_padj[with(all_cor_eN_df_padj, grepl("fotu*",from) & grepl("bOTU_*",to)),])
nrow(bf_occur_eN)

bp_occur_eN <- droplevels(all_cor_eN_df_padj[with(all_cor_eN_df_padj, grepl("potu*",from) & grepl("bOTU_*",to)),])
nrow(bp_occur_eN)

fp_occur_eN <- droplevels(all_cor_eN_df_padj[with(all_cor_eN_df_padj, grepl("potu*",from) & grepl("fotu*",to)),])
nrow(fp_occur_eN)

## Node degree values
eN_all_deg <- sort(degree(all_eN_net, mode="all"), decr=T)
max(eN_all_deg)
mean(eN_all_deg)

## Define keystone species based on top 1% of node degree values
n <- 1
eN_all_keystone <- eN_all_deg[eN_all_deg > quantile(eN_all_deg, prob=1-n/100)]
length(eN_all_keystone)

## Get bacteria and fungi keystone OTUs
eN_net_keystone_bac <- names(eN_all_keystone)[grep("bOTU_*", names(eN_all_keystone))]
length(eN_net_keystone_bac)

eN_net_keystone_fun <- names(eN_all_keystone)[grep("fotu*", names(eN_all_keystone))]
length(eN_net_keystone_fun)

eN_net_keystone_pro <- names(eN_all_keystone)[grep("potu*", names(eN_all_keystone))]
length(eN_net_keystone_pro)


## get csOTUs present in network
cs_eN_net <- names(V(all_eN_net))[names(V(all_eN_net)) %in% indic_edge_eN_combine]

cs_eN_net_bac <- cs_eN_net[grep("bOTU_*", cs_eN_net)]
length(cs_eN_net_bac)

cs_eN_net_fun <- cs_eN_net[grep("fotu*", cs_eN_net)]
length(cs_eN_net_fun)

cs_eN_net_pro <- cs_eN_net[grep("potu*", cs_eN_net)]
length(cs_eN_net_pro)

## Keystone OTUs present in csOTUs
intersect(eN_net_keystone_bac, cs_eN_net_bac)
intersect(eN_net_keystone_fun, cs_eN_net_fun)
intersect(eN_net_keystone_pro, cs_eN_net_pro)

write.csv(eN_net_keystone_bac, "eN_net_keystone_bac.csv")   # export keystone OTUs under eN
write.csv(eN_net_keystone_fun, "eN_net_keystone_fun.csv")   # export keystone OTUs under eN
write.csv(eN_net_keystone_pro, "eN_net_keystone_pro.csv")   # export keystone OTUs under eN
write.csv(cs_eN_net_bac, "cs_eN_net_bac.csv")               # export sensitive OTUs under eN
write.csv(cs_eN_net_fun, "cs_eN_net_fun.csv")               # export sensitive OTUs under eN
write.csv(cs_eN_net_pro, "cs_eN_net_pro.csv")               # export sensitive OTUs under eN



## Set node shape
V(all_eN_net)$shape <- V(all_eN_net)$name
V(all_eN_net)$shape[V(all_eN_net)$shape %in% names(eN_all_keystone)] <- "star"
V(all_eN_net)$shape[V(all_eN_net)$shape %in% rownames(otu_norm_eN_16s)] <- "circle"
V(all_eN_net)$shape[V(all_eN_net)$shape %in% rownames(otu_norm_eN_its)] <- "triangle"
V(all_eN_net)$shape[V(all_eN_net)$shape %in% rownames(otu_norm_eN_18s)] <- "square"

cs <- c("CLaP","CLaP_CLeP","CLeP","CLaP_uCaP","CLeP_uCeP","uCaP","uCaP_uCeP","uCeP")
unique(V(all_eN_net)$indicgroup)
V(all_eN_net)$color <- V(all_eN_net)$indicgroup
V(all_eN_net)$color[!V(all_eN_net)$color %in% cs] <- "gray30"
V(all_eN_net)$color[V(all_eN_net)$color == "CLaP"] <- "dodgerblue4"
V(all_eN_net)$color[V(all_eN_net)$color == "CLaP_CLeP"] <- "dodgerblue3"
V(all_eN_net)$color[V(all_eN_net)$color == "CLeP"] <- "dodgerblue1"
V(all_eN_net)$color[V(all_eN_net)$color == "CLaP_uCaP"] <- "sienna4"
V(all_eN_net)$color[V(all_eN_net)$color == "CLeP_uCeP"] <- "sienna3"
V(all_eN_net)$color[V(all_eN_net)$color == "uCaP"] <- "firebrick4"
V(all_eN_net)$color[V(all_eN_net)$color == "uCaP_uCeP"] <- "firebrick3"
V(all_eN_net)$color[V(all_eN_net)$color == "uCeP"] <- "firebrick2"
V(all_eN_net)$frame.color <- V(all_eN_net)$color

eN_all_net_csnodes <- rownames(nodeattrib_eN_combine[nodeattrib_eN_combine$indicgroup %in% cs,])
eN_all_net_csnodes_bac <- eN_all_net_csnodes[grep("bOTU*",eN_all_net_csnodes)]
eN_all_net_csnodes_fun <- eN_all_net_csnodes[grep("fotu*",eN_all_net_csnodes)]
eN_all_net_csnodes_pro <- eN_all_net_csnodes[grep("potu*",eN_all_net_csnodes)]

V(all_eN_net)$size <- V(all_eN_net)$name
V(all_eN_net)$size[V(all_eN_net)$size %in% names(eN_all_keystone)] <- 4
V(all_eN_net)$size[V(all_eN_net)$size %in% eN_all_net_csnodes_bac] <- 4
V(all_eN_net)$size[V(all_eN_net)$size %in% eN_all_net_csnodes_fun] <- 6
V(all_eN_net)$size[V(all_eN_net)$size %in% eN_all_net_csnodes_pro] <- 6
V(all_eN_net)$size[V(all_eN_net)$size %in% rownames(otu_norm_eN_16s)] <- 2
V(all_eN_net)$size[V(all_eN_net)$size %in% rownames(otu_norm_eN_its)] <- 3
V(all_eN_net)$size[V(all_eN_net)$size %in% rownames(otu_norm_eN_18s)] <- 3
eNcombine_nodesizes <- as.numeric(V(all_eN_net)$size)




#############################################################################################
##### MODEL 9  Explore community structure of eN meta-network by defining modules       #####
#############################################################################################
## Make vectors of network nodes responding to different Nitrogen addition systems
CIT_nodes_eN <- rownames(nodeattrib_eN_combine[nodeattrib_eN_combine$indicgroup=="CLaP",])
CNT_nodes_eN <- rownames(nodeattrib_eN_combine[nodeattrib_eN_combine$indicgroup=="CLeP",])
CON_nodes_eN <- rownames(nodeattrib_eN_combine[nodeattrib_eN_combine$indicgroup=="CLaP_CLeP",])
IT_nodes_eN <- rownames(nodeattrib_eN_combine[nodeattrib_eN_combine$indicgroup=="CLaP_uCaP",])
OIT_nodes_eN <- rownames(nodeattrib_eN_combine[nodeattrib_eN_combine$indicgroup=="uCaP",])
ORT_nodes_eN <- rownames(nodeattrib_eN_combine[nodeattrib_eN_combine$indicgroup=="uCeP",])
ORG_nodes_eN <- rownames(nodeattrib_eN_combine[nodeattrib_eN_combine$indicgroup=="uCaP_uCeP",])
NTRT_nodes_eN <- rownames(nodeattrib_eN_combine[nodeattrib_eN_combine$indicgroup=="CLeP_uCeP",])

cs_nodes_eN_all <- c(CIT_nodes_eN, CNT_nodes_eN, CON_nodes_eN, IT_nodes_eN, OIT_nodes_eN, ORT_nodes_eN, ORG_nodes_eN, NTRT_nodes_eN)
Bcs_nodes_eN_all <- cs_nodes_eN_all[grep("bOTU*", cs_nodes_eN_all)]
Fcs_nodes_eN_all <- cs_nodes_eN_all[grep("fotu*", cs_nodes_eN_all)]
Pcs_nodes_eN_all <- cs_nodes_eN_all[grep("potu*", cs_nodes_eN_all)]
## Perform cluster analysis using greedy clustering algorithm 
cfg_eN <- cluster_fast_greedy(as.undirected(all_eN_net))

## Subset for 20 biggest nodes
eN_modules <- sort(table(membership(cfg_eN)),decr=T)
eN_modules_20 <- eN_modules[1:20]
sum(eN_modules_20)/sum(eN_modules)
rm20_plot <- eN_modules_20
names(rm20_plot) <- as.factor(1:20)
eN_modules_cs <- table(factor(membership(cfg_eN)[cs_nodes_eN_all],levels=names(eN_modules)))
eN_modules_cs_20 <- eN_modules_cs[names(eN_modules_20)]
rmcs20_plot <- eN_modules_cs_20
names(rmcs20_plot) <- as.factor(1:20)

## Make vector of nodes in top 20 modules
eN_modules_points <- membership(cfg_eN)[membership(cfg_eN) %in% names(eN_modules_20)]

eN_points <- NULL
for(i in eN_modules_points){
  eNx <- which(names(eN_modules_20)==i)
  eN_points <- c(eN_points,eNx)
}
names(eN_points) <- names(eN_modules_points)

## Set node colors by Nitrogen addition system sensitivity 
eN_all_cols <- sort(eN_points)
eN_all_cols[!names(eN_all_cols) %in% cs] <- "gray30"
eN_all_cols[names(eN_all_cols) %in% CIT_nodes_eN] <- "dodgerblue4"
eN_all_cols[names(eN_all_cols) %in% CON_nodes_eN] <- "dodgerblue3"
eN_all_cols[names(eN_all_cols) %in% CNT_nodes_eN] <- "dodgerblue1"
eN_all_cols[names(eN_all_cols) %in% IT_nodes_eN] <- "sienna4"
eN_all_cols[names(eN_all_cols) %in% NTRT_nodes_eN]<- "sienna3"
eN_all_cols[names(eN_all_cols) %in% OIT_nodes_eN]<- "firebrick4"
eN_all_cols[names(eN_all_cols) %in% ORG_nodes_eN] <- "firebrick3"
eN_all_cols[names(eN_all_cols) %in% ORT_nodes_eN] <- "firebrick2"

eN_all_pch <- sort(eN_points)
eN_all_pch[names(eN_all_pch) %in% rownames(otu_norm_eN_16s)] <- 1
eN_all_pch[names(eN_all_pch) %in% rownames(otu_norm_eN_its)] <- 2
eN_all_pch[names(eN_all_pch) %in% rownames(otu_norm_eN_18s)] <- 2
eN_all_pch[names(eN_all_pch) %in% intersect(rownames(otu_norm_eN_16s),cs_nodes_eN_all)] <- 16
eN_all_pch[names(eN_all_pch) %in% intersect(rownames(otu_norm_eN_its),cs_nodes_eN_all)] <- 17
eN_all_pch[names(eN_all_pch) %in% intersect(rownames(otu_norm_eN_18s),cs_nodes_eN_all)] <- 17
eN_all_pch[names(eN_all_pch) %in% names(eN_all_keystone)] <- 8

eN_all_cex <- sort(eN_points)
eN_all_cex[!names(eN_all_cex) %in% cs_nodes_eN_all] <- 1
eN_all_cex[names(eN_all_cex) %in% cs_nodes_eN_all] <- 2

eN_mods_list_cs <- list()
for (i in names(eN_modules_cs_20)){
  x1 <- names(membership(cfg_eN)[membership(cfg_eN) == i])
  x2 <- x1[x1 %in% cs_nodes_eN_all]
  eN_mods_list_cs[[i]] <- as.numeric(V(all_eN_net)[x2])
}



############################################################################################################
#####   MODEL 10    Figure 4: aN and eN community meta co-occurrence networks                         #####
############################################################################################################
## Note: the permuations for the layout of the network can be very time consuming and processor intensive
set.seed(8051)
coords_aN_all <- layout_(all_aN_net, with_fr(niter=9999, grid="nogrid"))
write.table(coords_aN_all,paste0(output,"coords_aN_all.txt"),sep="\t",row.names=F,col.names=F,quote=F)

## Note: the permuations for the layout of the network can be very time consuming and processor intensive
set.seed(8051)
coords_eN_all <- layout_(all_eN_net, with_fr(niter=9999, grid="nogrid"))
write.table(coords_eN_all,paste0(output,"coords_eN_all.txt"),sep="\t",row.names=F,col.names=F,quote=F)

## Import pre-calculated FR layout coordinates to save time 
#coords_aN_all <- as.matrix(read.table("coords_aN_all.txt"))
#dimnames(coords_aN_all) <-  NULL

## Import pre-calculated FR layout coordinates to save time 
#coords_eN_all <- as.matrix(read.table("coords_eN_all.txt"))
#dimnames(coords_eN_all) <-  NULL


###############################################################################################
######################                 Plot Figure 4               #############################################
pdf(paste0(output,"Figure4a1.pdf"),width=36,height=18)
par(mfrow=c(1,2), mar=c(0,0,0,0))


aN_cols <- c("burlywood1","DarkOliveGreen1","DarkOrchid1", "MediumPurple1")
plot(all_aN_net,vertex.label=NA, edge.width=0.01, vertex.size=aNcombine_nodesizes, layout=coords_aN_all,
     mark.groups=list(aN_mods_list_cs$`3`,aN_mods_list_cs$`6`,aN_mods_list_cs$`5`, aN_mods_list_cs$`1`),         # Firstly, check "aN_mods_list_cs", module 
     mark.col=aN_cols, mark.border=aN_cols)
legend("bottomleft",legend=c("Module 3", "Module 6", "Module 5", "Module 1"),col=aN_cols,
       bty="n",fill=aN_cols,border=aN_cols, x.intersp=0.5, y.intersp=1,cex=2)

eN_cols <- c("OliveDrab1", "darkseagreen1", "CadetBlue1", "bisque", "aquamarine")
plot(all_eN_net,vertex.label=NA,edge.width=0.01, vertex.size=eNcombine_nodesizes,layout=coords_eN_all,
     mark.groups=list(eN_mods_list_cs$`1`,eN_mods_list_cs$`4`,eN_mods_list_cs$`3`, eN_mods_list_cs$`2`, eN_mods_list_cs$`7`),       # check eN_mods_list_cs
     mark.col=eN_cols, mark.border=eN_cols)
legend("bottomright",legend=c("Module 1", "Module 4", "Module 3", "Module 2", "Module 7"), col=eN_cols,
       bty="n", fill=eN_cols, border=eN_cols, x.intersp=0.5, y.intersp=1,cex=2)

dev.off()
###############################################################################################

## Plot separate legend for Figure 4
pdf(paste0(output,"Figure4a_legend.pdf"),height=7,width=10)

mat2 <- cbind(c(1,1,1,1,2,2,2,2,3,3,3,3), c(1,1.5,2,2.5,1,1.5,2,2.5,1,1.5,2,2.5))
cols <- c("white","dodgerblue1", "dodgerblue3", "dodgerblue4", "gray30","sienna3", "white", "sienna4", "white","firebrick2", "firebrick3", "firebrick4")
names(cols) <- c("CLeP", "shared", "CLaP", "shared", "", "shared", "uCeP", "shared", "uCaP")
pchlegend <- c(19,19,19,19,8,19,19,19,19,19,19,19)
plot(mat2, col=cols, pch=pchlegend, cex=5, xlim=c(.5,5), ylim=c(.5,3.5), frame=F, axes=F, xlab=NA, ylab=NA)
points(mat2+0.2, col=cols, pch=17, cex=5)
text(1:3, rep(3, 3), label=c("CL", "Shared", "uC"), adj=0.5, cex=2)
text(rep(3.5, 2), c(1.5,2,2.5), label=c("eP", "Shared", "aP"), adj=0,cex=2)
text(2.7, 1, label="Keystone OTU",cex=2)

dev.off()



########################################################################################
###                      MODEL   11       plotting average module response                                              #####
########################################################################################
pdf(paste0(output,"Figure4b.pdf"),width=10,height=7/5)
par(mfrow=c(1,10), mar=c(0.5,3.5,2,0.5))

CS_cols <- c("dodgerblue4","dodgerblue1","firebrick4","firebrick2")
names(CS_cols) <- c("CLaP","CLeP", "uCaP","uCeP")

## aN module 3
bargraph.CI(design_16s_aN$CL_P, colSums(otu_norm_aN_combine[cfg_aN[[3]],])/1000, 
            las=2, ylab="Cumulative relative abundance", cex.lab=.8, cex.axis=.9, cex.names=.7,
            err.width=.025, main="aN M3", col=CS_cols, border=F)
## aN module 6
bargraph.CI(design_16s_aN$CL_P, colSums(otu_norm_aN_combine[cfg_aN[[6]],])/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="aN M6", col=CS_cols, border=F)
## aN module 5
bargraph.CI(design_16s_aN$CL_P, colSums(otu_norm_aN_combine[cfg_aN[[5]],])/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="aN M5", col=CS_cols, border=F)
## aN module 1
bargraph.CI(design_16s_aN$CL_P, colSums(otu_norm_aN_combine[cfg_aN[[1]],])/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="aN M1", col=CS_cols, border=F)

## eN module 1
bargraph.CI(design_16s_eN$CL_P, colSums(otu_norm_eN_combine[cfg_eN[[1]],])/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="eN M1", col=CS_cols, border=F)
## eN module 4
bargraph.CI(design_16s_eN$CL_P, colSums(otu_norm_eN_combine[cfg_eN[[4]],])/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="eN M4", col=CS_cols, border=F)
## eN module 3
bargraph.CI(design_16s_eN$CL_P, colSums(otu_norm_eN_combine[cfg_eN[[3]],])/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="eN M3", col=CS_cols, border=F)
## eN module 2
bargraph.CI(design_16s_eN$CL_P, colSums(otu_norm_eN_combine[cfg_eN[[2]],])/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="eN M2", col=CS_cols, border=F)
## eN module 7
bargraph.CI(design_16s_eN$CL_P, colSums(otu_norm_eN_combine[cfg_eN[[7]],])/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="eN M7", col=CS_cols, border=F)
## Legend
plot.new()
par(mar=c(0.5,0,2,0))
legend("left",  bty="n", cex=1, #x.intersp=0.1, y.intersp=1,
       legend=names(CS_cols), 
       fill=CS_cols, 
       border=CS_cols , xpd=T)

dev.off()


 

##########################################################################
####       MODEL 12 taxonomies of module OTUs                                                ######
##########################################################################
## defining bacteria and fungi of aN and eN modules

# aN module 3
aN_M3 <- cfg_aN[[3]]
bacteria_in_aN_M3 <- aN_M3[grep("bOTU", as.vector(aN_M3))]
fungi_in_aN_M3 <- aN_M3[grep("fotu", as.vector(aN_M3))]
protist_in_aN_M3 <- aN_M3[grep("potu", as.vector(aN_M3))]
# aN module 6
aN_M6 <- cfg_aN[[6]]
bacteria_in_aN_M6<- aN_M6[grep("bOTU", as.vector(aN_M6))]
fungi_in_aN_M6<- aN_M6[grep("fotu", as.vector(aN_M6))]                      ## note OTU or otu, capital or lower letter
protist_in_aN_M6 <- aN_M6[grep("potu", as.vector(aN_M6))]
# aN module 5
aN_M5 <- cfg_aN[[5]]
bacteria_in_aN_M5 <- aN_M5[grep("bOTU", as.vector(aN_M5))]
fungi_in_aN_M5<- aN_M5[grep("fotu", as.vector(aN_M5))]
protist_in_aN_M5 <- aN_M5[grep("potu", as.vector(aN_M5))]
# aN module 1
aN_M1 <- cfg_aN[[1]]
bacteria_in_aN_M1 <- aN_M1[grep("bOTU", as.vector(aN_M1))]
fungi_in_aN_M1 <- aN_M1[grep("fotu", as.vector(aN_M1))]
protist_in_aN_M1 <- aN_M1[grep("potu", as.vector(aN_M1))]


# eN module 1
eN_M1 <- cfg_eN[[1]]
bacteria_in_eN_M1 <- eN_M1[grep("bOTU", as.vector(eN_M1))]
fungi_in_eN_M1 <- eN_M1[grep("fotu", as.vector(eN_M1))]
protist_in_eN_M1 <- eN_M1[grep("potu", as.vector(eN_M1))]
# eN module 4
eN_M4 <- cfg_eN[[4]]
bacteria_in_eN_M4 <- eN_M4[grep("bOTU", as.vector(eN_M4))]
fungi_in_eN_M4 <- eN_M4[grep("fotu", as.vector(eN_M4))]
protist_in_eN_M4 <- eN_M4[grep("potu", as.vector(eN_M4))]
# eN module 3
eN_M3 <- cfg_eN[[3]]
bacteria_in_eN_M3 <- eN_M3[grep("bOTU", as.vector(eN_M3))]
fungi_in_eN_M3 <- eN_M3[grep("fotu", as.vector(eN_M3))]
protist_in_eN_M3 <- eN_M3[grep("potu", as.vector(eN_M3))]
# eN module 2
eN_M2 <- cfg_eN[[2]]
bacteria_in_eN_M2 <- eN_M2[grep("bOTU", as.vector(eN_M2))]
fungi_in_eN_M2 <- eN_M2[grep("fotu", as.vector(eN_M2))]
protist_in_eN_M2<- eN_M2[grep("potu", as.vector(eN_M2))]
# eN module 7
eN_M7 <- cfg_eN[[7]]
bacteria_in_eN_M7 <- eN_M7[grep("bOTU", as.vector(eN_M7))]
fungi_in_eN_M7 <- eN_M7[grep("fotu", as.vector(eN_M7))]
protist_in_eN_M7 <- eN_M7[grep("potu", as.vector(eN_M7))]


### counts of bacteria OTUs
bacteria_aN_M3 <- as.data.frame(table(tax_filter_16s[bacteria_in_aN_M3, "labels"] ) )
colnames(bacteria_aN_M3) <- c("Class", "aN_M3")
bacteria_aN_M6 <- as.data.frame(table(tax_filter_16s[bacteria_in_aN_M6, "labels"] ) )
colnames(bacteria_aN_M6) <- c("Class", "aN_M6")
bacteria_aN_M5 <- as.data.frame(table(tax_filter_16s[bacteria_in_aN_M5, "labels"] ) )
colnames(bacteria_aN_M5) <- c("Class", "aN_M5")
bacteria_aN_M1 <- as.data.frame(table(tax_filter_16s[bacteria_in_aN_M1, "labels"] ) )
colnames(bacteria_aN_M1) <- c("Class", "aN_M1")

bacteria_aN_modules <- merge(bacteria_aN_M3, bacteria_aN_M6, all=T, by="Class") 
bacteria_aN_modules <- merge(bacteria_aN_modules, bacteria_aN_M5, all=T, by="Class") 
bacteria_aN_modules <- merge(bacteria_aN_modules, bacteria_aN_M1, all=T, by="Class") 
bacteria_aN_modules

bacteria_eN_M1 <- as.data.frame(table(tax_filter_16s[bacteria_in_eN_M1, "labels"] ) )
colnames(bacteria_eN_M1) <- c("Class", "eN_M1")
bacteria_eN_M4<- as.data.frame(table(tax_filter_16s[bacteria_in_eN_M4, "labels"] ) )
colnames(bacteria_eN_M4) <- c("Class", "eN_M4")
bacteria_eN_M3<- as.data.frame(table(tax_filter_16s[bacteria_in_eN_M3, "labels"] ) )
colnames(bacteria_eN_M3) <- c("Class", "eN_M3")
bacteria_eN_M2<- as.data.frame(table(tax_filter_16s[bacteria_in_eN_M2, "labels"] ) )
colnames(bacteria_eN_M2) <- c("Class", "eN_M2")
bacteria_eN_M7<- as.data.frame(table(tax_filter_16s[bacteria_in_eN_M7, "labels"] ) )
colnames(bacteria_eN_M7) <- c("Class", "eN_M7")

bacteria_eN_modules <- merge(bacteria_eN_M1, bacteria_eN_M4, all=T, by="Class") 
bacteria_eN_modules <- merge(bacteria_eN_modules, bacteria_eN_M3, all=T, by="Class")
bacteria_eN_modules <- merge(bacteria_eN_modules, bacteria_eN_M2, all=T, by="Class")
bacteria_eN_modules <- merge(bacteria_eN_modules, bacteria_eN_M7, all=T, by="Class")
bacteria_eN_modules

bacteria_modules <- merge(bacteria_aN_modules, bacteria_eN_modules, all=T, by="Class") 
bacteria_all_OTUs <- as.data.frame(table(tax_filter_16s[, "labels"] ) )
colnames(bacteria_all_OTUs) <- c("Class", "all bOTUs")
bacteria_modules <- merge(bacteria_modules, bacteria_all_OTUs, all=T, by="Class") 
bacteria_modules

bacteria_modules_mat <- bacteria_modules[2:11]                              ## note the number in []
rownames(bacteria_modules_mat) <- bacteria_modules$Class
bacteria_modules_mat[is.na(bacteria_modules_mat)] <- 0
colSums(bacteria_modules_mat)

bacteria_modules_prop <- t(t(bacteria_modules_mat)/colSums(bacteria_modules_mat) ) * 1
bacteria_modules_prop
colSums(bacteria_modules_prop)


### counts of fungi OTUs module 3,6,5,1
fungi_aN_M3 <- as.data.frame(table(tax_filter_its[fungi_in_aN_M3, "Phylum"] ) )
colnames(fungi_aN_M3) <- c("Class", "aN_M3")
fungi_aN_M6<- as.data.frame(table(tax_filter_its[fungi_in_aN_M6, "Phylum"] ) )
colnames(fungi_aN_M6) <- c("Class", "aN_M6")
fungi_aN_M5<- as.data.frame(table(tax_filter_its[fungi_in_aN_M5, "Phylum"] ) )
colnames(fungi_aN_M5) <- c("Class", "aN_M5")
fungi_aN_M1 <- as.data.frame(table(tax_filter_its[fungi_in_aN_M1, "Phylum"] ) )
colnames(fungi_aN_M1) <- c("Class", "aN_M1")

fungi_aN_modules <- merge(fungi_aN_M3, fungi_aN_M6, all=T, by="Class") 
fungi_aN_modules <- merge(fungi_aN_modules, fungi_aN_M5, all=T, by="Class") 
fungi_aN_modules <- merge(fungi_aN_modules, fungi_aN_M1, all=T, by="Class") 
fungi_aN_modules
 #  module 1, 4, 3, 2, 7 
fungi_eN_M1 <- as.data.frame(table(tax_filter_its[fungi_in_eN_M1, "Phylum"] ) )
colnames(fungi_eN_M1) <- c("Class", "eN_M1")
fungi_eN_M4 <- as.data.frame(table(tax_filter_its[fungi_in_eN_M4, "Phylum"] ) )
colnames(fungi_eN_M4) <- c("Class", "eN_M4")
fungi_eN_M3 <- as.data.frame(table(tax_filter_its[fungi_in_eN_M3, "Phylum"] ) )
colnames(fungi_eN_M3) <- c("Class", "eN_M3")
fungi_eN_M2 <- as.data.frame(table(tax_filter_its[fungi_in_eN_M2, "Phylum"] ) )
colnames(fungi_eN_M2) <- c("Class", "eN_M2")
fungi_eN_M7 <- as.data.frame(table(tax_filter_its[fungi_in_eN_M7, "Phylum"] ) )
colnames(fungi_eN_M7) <- c("Class", "eN_M7")

fungi_eN_modules <- merge(fungi_eN_M1, fungi_eN_M4, all=T, by="Class") 
fungi_eN_modules <- merge(fungi_eN_modules, fungi_eN_M3, all=T, by="Class") 
fungi_eN_modules <- merge(fungi_eN_modules, fungi_eN_M2, all=T, by="Class") 
fungi_eN_modules <- merge(fungi_eN_modules, fungi_eN_M7, all=T, by="Class") 
fungi_eN_modules

fungi_modules <- merge(fungi_aN_modules, fungi_eN_modules, all=T, by="Class") 
fungi_all_OTUs <- as.data.frame(table(tax_filter_its[, "Phylum"] ) )
colnames(fungi_all_OTUs) <- c("Class", "all fotus")
fungi_modules <- merge(fungi_modules, fungi_all_OTUs, all=T, by="Class") 
fungi_modules

fungi_modules_mat <- fungi_modules[2:11]                      ## note the number in []
rownames(fungi_modules_mat) <- fungi_modules$Class
fungi_modules_mat[is.na(fungi_modules_mat)] <- 0
colSums(fungi_modules_mat)

fungi_modules_prop <- t(t(fungi_modules_mat)/colSums(fungi_modules_mat) ) * 1
fungi_modules_prop
colSums(fungi_modules_prop)


### counts of protist OTUs  module 3,6,5,1
protist_aN_M3 <- as.data.frame(table(tax_filter_18s[protist_in_aN_M3, "Kingdom"] ) )
colnames(protist_aN_M3) <- c("Class", "aN_M3")
protist_aN_M6 <- as.data.frame(table(tax_filter_18s[protist_in_aN_M6, "Kingdom"] ) )
colnames(protist_aN_M6) <- c("Class", "aN_M6")
protist_aN_M5 <- as.data.frame(table(tax_filter_18s[protist_in_aN_M5, "Kingdom"] ) )
colnames(protist_aN_M5) <- c("Class", "aN_M5")
protist_aN_M1<- as.data.frame(table(tax_filter_18s[protist_in_aN_M1, "Kingdom"] ) )
colnames(protist_aN_M1) <- c("Class", "aN_M1")

protist_aN_modules <- merge(protist_aN_M3, protist_aN_M6, all=T, by="Class") 
protist_aN_modules <- merge(protist_aN_modules, protist_aN_M5, all=T, by="Class") 
protist_aN_modules <- merge(protist_aN_modules, protist_aN_M1, all=T, by="Class") 
protist_aN_modules

# module 1,4,3,2,7
protist_eN_M1 <- as.data.frame(table(tax_filter_18s[protist_in_eN_M1, "Kingdom"] ) )
colnames(protist_eN_M1) <- c("Class", "eN_M1")
protist_eN_M4 <- as.data.frame(table(tax_filter_18s[protist_in_eN_M4, "Kingdom"] ) )
colnames(protist_eN_M4) <- c("Class", "eN_M4")
protist_eN_M3 <- as.data.frame(table(tax_filter_18s[protist_in_eN_M3, "Kingdom"] ) )
colnames(protist_eN_M3) <- c("Class", "eN_M3")
protist_eN_M2 <- as.data.frame(table(tax_filter_18s[protist_in_eN_M2, "Kingdom"] ) )
colnames(protist_eN_M2) <- c("Class", "eN_M2")
protist_eN_M7 <- as.data.frame(table(tax_filter_18s[protist_in_eN_M7, "Kingdom"] ) )
colnames(protist_eN_M7) <- c("Class", "eN_M7")
protist_eN_modules <- merge(protist_eN_M1, protist_eN_M4, all=T, by="Class") 
protist_eN_modules <- merge(protist_eN_modules, protist_eN_M3, all=T, by="Class") 
protist_eN_modules <- merge(protist_eN_modules, protist_eN_M2, all=T, by="Class") 
protist_eN_modules <- merge(protist_eN_modules, protist_eN_M7, all=T, by="Class") 
protist_eN_modules

protist_modules <- merge(protist_aN_modules, protist_eN_modules, all=T, by="Class") 
protist_all_OTUs <- as.data.frame(table(tax_filter_18s[, "Kingdom"] ) )
colnames(protist_all_OTUs) <- c("Class", "all potus")
protist_modules <- merge(protist_modules, protist_all_OTUs, all=T, by="Class") 
protist_modules

protist_modules_mat <- protist_modules[2:11]
rownames(protist_modules_mat) <- protist_modules$Class
protist_modules_mat[is.na(protist_modules_mat)] <- 0
colSums(protist_modules_mat)

protist_modules_prop <- t(t(protist_modules_mat)/colSums(protist_modules_mat) ) * 1
protist_modules_prop
colSums(protist_modules_prop)




##########################################################################
####       MODEL 13 PLOT                                            ######
##########################################################################
pdf(paste0(output,"Figure4c.pdf"), width = 7, height = 7/6)
par(mfrow=c(1,6), mar=c(2.5,2,1,0))
### bacteria
# PHYLA_label_cols_16s_legend
table(rownames(bacteria_modules_prop) %in% PHYLA_label_cols_16s$labels) 
bp <- barplot(cbind(bacteria_modules_prop[,1:4], NA, bacteria_modules_prop[,5:9], NA, bacteria_modules_prop[,10]),
              las=2, border=NA, axes=F, cex.names=.5,                                ## note the number in []
              col=PHYLA_label_cols_16s[rownames(bacteria_modules_prop),]$cols )
text(bp, 1.1, labels=c(colSums(bacteria_modules_mat)[1:4], NA,                       ## note the number in []
                       colSums(bacteria_modules_mat)[5:9], NA,                       ## note the number in []
                       colSums(bacteria_modules_mat)[10]), xpd=T, cex=.4, las=2)     ## note the number in []

## legend from Fig. S3
plot.new()
legend("left", bty="n", cex=0.5, x.intersp=0.4, y.intersp=0.75,
       legend=rev(PHYLA_label_cols_16s_legend$labels), 
       fill=rev(PHYLA_label_cols_16s_legend$cols), 
       border=rev(PHYLA_label_cols_16s_legend$cols))

### Fungi
# PHYLA_label_cols_its_legend
table(rownames(fungi_modules_prop) %in% PHYLA_label_cols_its$Phylum) 

fp <- barplot(cbind(fungi_modules_prop[,1:4], NA, fungi_modules_prop[,5:9], NA, fungi_modules_prop[,10]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_its[rownames(fungi_modules_prop),]$cols )
text(fp, 1.1, labels=c(colSums(fungi_modules_mat)[1:4], NA,
                       colSums(fungi_modules_mat)[5:9], NA,
                       colSums(fungi_modules_mat)[10]), xpd=T, cex=.4, las=2)

## legend from Fig. S3
plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.4, y.intersp=1, 
       legend=rev(PHYLA_label_cols_its_legend$Phylum),
       fill=rev(PHYLA_label_cols_its_legend$cols),
       border=rev(PHYLA_label_cols_its_legend$cols))

### protist
# PHYLA_label_cols_18s_legend
table(rownames(protist_modules_prop) %in% PHYLA_label_cols_18s$Kingdom) 
bp <- barplot(cbind(protist_modules_prop[,1:4], NA, protist_modules_prop[,5:9], NA, protist_modules_prop[,10]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_18s[rownames(protist_modules_prop),]$cols )
text(bp, 1.1, labels=c(colSums(protist_modules_mat)[1:4], NA,
                       colSums(protist_modules_mat)[5:9], NA,
                       colSums(protist_modules_mat)[10]), xpd=T, cex=.4, las=2)

## legend from Fig. S3
plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.4, y.intersp=0.75,
       legend=rev(PHYLA_label_cols_18s_legend$Kingdom), 
       fill=rev(PHYLA_label_cols_18s_legend$cols), 
       border=rev(PHYLA_label_cols_18s_legend$cols))
dev.off()











