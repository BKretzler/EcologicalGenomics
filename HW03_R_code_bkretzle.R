# ---------------------------------------------------------
# HW03 R code
# 22 Nov 2021
# Bailey M Kretzler
# ---------------------------------------------------------

#---------------------------------------------------------#
#                Imports and Working Directory                      
#---------------------------------------------------------#
library(ggplot2)
library(gridExtra)
library(GenomicRanges)
library(GenomicFeatures)
library(patchwork)


setwd("C:/Users/kretz/R/EcologicalGenomics/ecologicalgenomics/PopGenomics")

#---------------------------------------------------------#
#                      Load in files                      
#---------------------------------------------------------#
snps <- read.table("Chr18.kept.sites", sep = "\t", header = T)

# Read in the local ancestry frequencies from Plink
AF <- read.table("Chr18_LAI_freq.afreq", skip=1,sep="\t",header=F) 


# Get the list of admixed individuals:
Admixed <- read.table("Admixed.Inds",header=F)

# Get the meta data:
meta <- read.table("Combined_Transect_Sampling_Data_2020.txt", sep="\t",header=T)


# Read in the Admixture coefficients for KBals that we made from the K=5 file:
KBals <- read.table("Admixed_KBals", sep="\t", header=F)
names(KBals) = c("ID","KBals")

# Bring in phenotype data:
pheno <- read.table("VT_Garden_Phenotypes_2021.txt",sep="\t",header=T)
clim <- read.table("climDat.txt",sep="\t",header=T)


#---------------------------------------------------------#
#                      Merge files together                      
#---------------------------------------------------------#

# first merge - meta + admixded:
meta_admx <- merge(meta, Admixed, by.x="ID", by.y="V1")
str(meta_admx)  

# Second merge - kbals + met_admx:
meta_admx_KBals <- merge(meta_admx,KBals,by="ID")


# Merge pheno/clim data with meta and KBals:
meta_admx_KBals_pheno <- merge(meta_admx_KBals,pheno,by="ID")
meta_admx_KBals_clim <- merge(meta_admx_KBals,clim,by="ID")



#---------------------------------------------------------#
#                 Plotting + analysis                     
#---------------------------------------------------------#


#---------------------------------------------------------#
## Do phenotypes show a relationship to overall ancestry 
#---------------------------------------------------------#

#  CDD
plotCDD <- ggplot(meta_admx_KBals_clim,aes(x=KBals,y=med_DD0, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Chilling Degree Data") + 
  theme_bw()

plotCDD

# GDD
plotGDD <- ggplot(meta_admx_KBals_clim,aes(x=KBals,y=mean_cGDDfreeze, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Growing Degree Days") + 
  theme_bw()

plotGDD

# Final
plotFinal <- ggplot(meta_admx_KBals_clim,aes(x=KBals,y=mean_finalFreeze, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Final Freezing Event") + 
  theme_bw()

plotFinal

# Bud flush
plotBudflush <- ggplot(meta_admx_KBals_pheno,aes(x=KBals,y=FLUSH, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Bud flush") + 
  theme_bw()

plotBudflush

## Arrange all plots into one
combo_ancestry_traits <-plotCDD + plotGDD + plotFinal + plotBudflush + plot_layout(guides = "collect") & theme(legend.position = "right") 

#<-grid_arrange_shared_legend(plotCDD, plotGDD, plotFinal, plotBudflush, nrow = 2)
ggsave("traits_by_ancestry.png", plot = combo_ancestry_traits)

# linear models testing trait ~ genome-wide admixture association
summary(lm(med_DD0~KBals + Transect.x, data=meta_admx_KBals_clim))

summary(lm(mean_cGDDfreeze~KBals + Transect.x, data=meta_admx_KBals_clim))

summary(lm(mean_finalFreeze~KBals + Transect.x, data=meta_admx_KBals_clim))

summary(lm(FLUSH~KBals + Transect.x, data=meta_admx_KBals_pheno))


### All vars show a significant relationship with proportion of P balsamifera ancestry
# strong positive correlation with CDD
# mild positive correlation with GDD
# strong positive correlation with final freeze date
# strong negative correlation with bud flush




#---------------------------------------------------------#
#local ancestry with genome wide ancestry as a covariate:
#---------------------------------------------------------#

######  Bring in Association results from Plink   ######

#########  CDD  #########

# load in plink file
CDD <- read.table("plink2.CDD.glm.linear",skip=1,sep="\t",header=F)
names(CDD) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
CDD2 <- CDD[which(CDD$TEST=="ADD"),]


# Define association outliers as the upper 1% of p-values
CDD2 <- cbind(snps, CDD2[,-c(1:2)])
CDD2$outlier = ifelse(CDD2$P<quantile(CDD2$P,0.01),2,1)


# Plot
p1 <- ggplot(CDD2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=CDD2$outlier, color=CDD2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Chilling Degree Days")+ 
  theme_bw()



####### GDD #########
# load in plink file
GDD <- read.table("plink2.GDD.glm.linear",skip=1,sep="\t",header=F)
names(GDD) = c("CHROM",  "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
GDD2 <- GDD[which(GDD$TEST=="ADD"),]

# outliers
GDD2 <- cbind(snps, GDD2[,-c(1,2)])
GDD2$outlier = ifelse(GDD2$P<quantile(GDD2$P,0.01),2,1)


#plot
p2 <- ggplot(GDD2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=GDD2$outlier, color=GDD2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Growing Degree Days")+ 
  theme_bw()


#########  Final  #########

#read in plink file
final <- read.table("plink2.Final.glm.linear",skip=1,sep="\t",header=F)
names(final) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
final <- final[which(final$TEST=="ADD"),]

#outliers
final2 <- cbind(snps, final[,-c(1,2)])
final2$outlier = ifelse(final2$P<quantile(final2$P,0.01),2,1)


#plot
p3 <- ggplot(final2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=final2$outlier, color=final2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Final Freezing Event") + 
  theme_bw()

#########  Bud flush  #########
budflush <- read.table("plink2.FLUSH.glm.linear",skip=1,sep="\t",header=F)
names(budflush) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
budflush <- budflush[which(budflush$TEST=="ADD"),]
budflush2 <- cbind(snps, budflush[,-c(1,2)])
budflush2$outlier = ifelse(budflush2$P<quantile(budflush2$P,0.01),2,1)

p4 <- ggplot(budflush2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=budflush2$outlier, color=budflush2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Bud flush") + 
  theme_bw()


# Plot all together - indication of where QTLs may be
comboP1to4 <- grid.arrange(p1, p2, p3,p4, nrow = 2)
ggsave("outlier_regions.png", plot = comboP1to4)


# Get outliers for a given trait association:

CDD_outliers <- CDD2[which(CDD2$outlier==2),c(2,3,9)]
GDD_outliers <- GDD2[which(GDD2$outlier==2),c(2,3,9)]
final_outliers <- final2[which(final2$outlier==2),c(2,3,9)]
flush_outliers <- budflush2[which(budflush2$outlier==2),c(2,3,9)]

#---------------------------------------------------------#
#plot proportion of LAI along the chromosome                      
#---------------------------------------------------------#

#grab the LAI file
AF <- read.table("Chr18_LAI_freq.afreq", skip=1,sep="\t",header=F)
names(AF) = c("CHROM",  "ID",   "REF",  "ALT",  "ALT_FREQS",    "OBS_CT")
str(AF)

AF2 <- cbind(snps,AF)

windows <- seq(1,max(AF2$POS),5e4)
AF_windows <- numeric()

for(i in 1:length(windows)){
  tmp=AF2[which(AF2$POS>windows[i] & AF2$POS<windows[i+1]),"ALT_FREQS"]
  ancfreq=mean(tmp)
  AF_windows[i] = ancfreq
}

AF3 <- as.data.frame(cbind(windows,AF_windows))
names(AF3) = c("window","AvgAncFreq")

### define outlier regions on the chromosomes

upper = mean(AF3$AvgAncFreq,na.rm=T) + 2*sd(AF3$AvgAncFreq,na.rm=T)
lower = mean(AF3$AvgAncFreq,na.rm=T) - 2*sd(AF3$AvgAncFreq,na.rm=T)

outliers_upper = AF3[which(AF3$AvgAncFreq>upper),]
outliers_lower = AF3[which(AF3$AvgAncFreq<lower),]

# Print the outlier regions out
outliers_upper
outliers_lower

# And finally, make plot 
p5 <- ggplot(AF3[,-3],aes(x=window,y=AvgAncFreq)) +
  geom_line(size=0.8, color="blue") + 
  xlab("Position (bp) along chromosome") +
  ylab("Frequency P. trichocarpa ancestry") +
  geom_hline(yintercept=mean(AF2$ALT_FREQS), color = "red") + 
  geom_hline(yintercept=upper, linetype="dashed", color = "red") + 
  geom_hline(yintercept=lower, linetype="dashed", color = "red") +
  ggtitle("Chr18: Local ancestry") + 
  theme_bw()

p5

comboP1to5 <- grid.arrange(p1, p2, p3,p4,p5, nrow = 2)
ggsave("outlier_regions_LAI.png", plot = comboP1to5)

# grid.arrange(p1, p2, p3, p4, nrow = 4)



#---------------------------------------------------------#
###look at betas
#---------------------------------------------------------#

# Get the betas from each trait and look at pleiotropy between traits
betas <- cbind(CDD2[,c(1:3,9)],GDD2[,9],final2[,9],budflush2[,9])
names(betas) = c("CHROM","POS","ID","beta_CDD","beta_GDD","beta_final", "beta_budflush")
str(betas)

cor(betas[,4:7],betas[4:7])

 plot(betas$beta_CDD,betas$beta_GDD) #

p6 <- ggplot(betas,aes(x=beta_final,y=beta_GDD)) +
  geom_point(color="black") + 
  xlab("Beta final freeze date") +
  ylab("Beta growing degree days") +
  ggtitle("Correlation of final freeze and growing degree days effect sizes") + 
  theme_bw()

p6

p7 <- ggplot(betas,aes(x=beta_final,y=beta_budflush)) +
  geom_point(color="black") + 
  xlab("Beta final freeze date") +
  ylab("Beta bud flush") +
  ggtitle("Correlation of final freeze and bud flush effect sizes") + 
  theme_bw()

p7 

ggsave("betas_final_gdd.png", plot = p6)

ggsave("betas_final_flush.png", plot = p7)


#---------------------------------------------------------#
#                      Genomic Ranges                      
#---------------------------------------------------------#

CHR = "Chr18"

# CDD
CDDGR <- GRanges(CHR,IRanges(CDD2$POS-2.5e4,CDD2$POS+2.5e4),POS=CDD2$POS, P=CDD2$P, outlier=CDD2$outlier)

CDDGRout <- unlist(reduce(split(CDDGR, ~outlier)))
CDDGRout$outlier <- names(CDDGRout)
CDDGRCand <- subset(CDDGRout, outlier==2)

CDDGRCand # Print the candidate regions 

# GDD
GDDGR <- GRanges(CHR,IRanges(GDD2$POS-2.5e4,GDD2$POS+2.5e4),POS=GDD2$POS, P=GDD2$P, outlier=GDD2$outlier)

GDDGRout <- unlist(reduce(split(GDDGR, ~outlier)))
GDDGRout$outlier <- names(GDDGRout)
GDDGRCand <- subset(GDDGRout, outlier==2)

GDDGRCand # Print the candidate regions 

#Final
finalGR <- GRanges(CHR,IRanges(final2$POS-2.5e4,final2$POS+2.5e4),POS=final2$POS, P=final2$P, outlier=final2$outlier)

finalGRout <- unlist(reduce(split(finalGR, ~outlier)))
finalGRout$outlier <- names(finalGRout)
finalGRCand <- subset(finalGRout, outlier==2)

finalGRCand # Print the candidate regions

# Budflush
budflushGR <- GRanges(CHR,IRanges(budflush2$POS-2.5e4,budflush2$POS+2.5e4),POS=budflush2$POS, P=budflush2$P, outlier=budflush2$outlier)

budflushGRout <- unlist(reduce(split(budflushGR, ~outlier)))
budflushGRout$outlier <- names(budflushGRout)
budflushGRCand <- subset(budflushGRout, outlier==2)

budflushGRCand # Print the candidate regions



####    Check for Overlap    ####

# bud flush x final freeze
overlap_BF_final <- subsetByOverlaps(budflushGRCand, finalGRCand)
length(overlap_BF_final)

overlap_BF_final # Print the overlapping regions - 0 overlapping regions

# bud flush x GDD
overlap_BF_GDD <- subsetByOverlaps(budflushGRCand, GDDGRCand)
length(overlap_BF_GDD)

overlap_BF_GDD # Print the overlapping regions - 2
## some overlap with growing degree days

# bud flush x CDD
overlap_BF_CDD <- subsetByOverlaps(budflushGRCand, CDDGRCand)
length(overlap_BF_CDD)

overlap_BF_CDD # Print the overlapping regions -0
## none with chilling degree days

# CDD x final freeze
overlap_CDD_final <- subsetByOverlaps(CDDGRCand, finalGRCand)
length(overlap_CDD_final)

overlap_CDD_final # Print the overlapping regions - 0

# GDD x final freeze
overlap_GDD_final <- subsetByOverlaps(GDDGRCand, finalGRCand)
length(overlap_GDD_final)

overlap_GDD_final # Print the overlapping regions - 2


# bud CDD x GDD
overlap_CDD_GDD <- subsetByOverlaps(CDDGRCand, GDDGRCand)
length(overlap_CDD_GDD)

overlap_CDD_GDD # Print the overlapping regions - 1


#---------------------------------------------------------#
#           Check the genes against annotated genome                       
#---------------------------------------------------------#


# Import the GFF annotation file and make a transcript database
txdb <- makeTxDbFromGFF("Ptrichocarpa_533_v4.1.gene.gff3.gz", format="gff3")

txdb

# How many chromosomes are present?
head(seqlevels(txdb))

# Subset the database for just your chromosome of interest
seqlevels(txdb) <- CHR # subset for just your chromosome

# Reduce the transcript database to just the non-redundant gene names, instead of multiple entries for all the variant transcript types per gene
genes <- unlist(reduce(transcriptsBy(txdb, by="gene"))) 
genes$geneID <- names(genes)

# extract overlap genes
candGenes_BF_GDD <- subsetByOverlaps(genes, overlap_BF_GDD)

write.table(candGenes_BF_GDD$geneID, paste0("CandGenes_BF_GDD",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")


# extract gdd genes
candGenes_GDD <- subsetByOverlaps(genes, GDDGRCand)

write.table(candGenes_GDD$geneID, paste0("CandGenes_GDD",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")



# extract overlap genes
candGenes_CDD <- subsetByOverlaps(genes, CDDGRCand)

write.table(candGenes_CDD$geneID, paste0("CandGenes_CDD",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")



# extract overlap genes
candGenes_final <- subsetByOverlaps(genes, finalGRCand)

write.table(candGenes_final$geneID, paste0("CandGenes_final",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")








