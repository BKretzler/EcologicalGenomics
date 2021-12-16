#---------------------------------------------------------#
#                 Set working directory                      
#---------------------------------------------------------#

#set working directory
setwd("C:/Users/kretz/R/EcologicalGenomics/")

#check to make sure it is correct
getwd()
#---------------------------------------------------------#
#                      Load packages                      
#---------------------------------------------------------#
library(rlang)
library(RSQLite)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn) 
library(patchwork)
library(eulerr)
#---------------------------------------------------------#
#                      Load in the file                       
#---------------------------------------------------------#

####    Load in counts    ####

## F3 counts mat
countsTable <- read.table(file = "DE_counts_F3.txt",header = T, row.names =1)
# Round
countsTableRound <- round(countsTable)  

##F1 counts mat
countsTableF1 <- read.table(file = "DE_counts_F1.txt",header = T, row.names =1)
#round
countsTableRoundF1 <- round(countsTableF1)  

####    METADATA    ####

##F3
cons <- read.delim(file = "RT_tonsa_F3_samples.txt", 
                   header = T, 
                   stringsAsFactors = T)    

##F1
consF1 <- read.delim(file = "RT_tonsa_F1_samples.txt", 
                   header = T, 
                   stringsAsFactors = T)    

#---------------------------------------------------------#
#                      exploring data                      
#---------------------------------------------------------#

#------------------Reads per sample -----------------------#

##total number or reads for each sample
colSums(countsTableRound)

##average number of reads per sample
mean(colSums(countsTableRound))
mean(colSums(countsTableRoundF1))


#------------------Reads per gene-----------------------#

## mean of counts per gene
mean(rowSums(countsTableRound))
mean(rowSums(countsTableRoundF1))


#check min reads per gene
min(rowSums(countsTableRound))
min(rowSums(countsTableRoundF1))

#---------------------------------------------------------#
#             define DESeq object + exp design                      
#---------------------------------------------------------#

##create the DESeq object from the data we loaded in 
  #accounting for interaction between line and environment
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = cons, 
                              design = ~line + environment + line:environment)
##check dimensions
dim(dds)

## filter out genes with too few reads
  #keep genes with average > 10 reads per sample
dds <- dds[rowSums(counts(dds))>160]
dim(dds)

####    F1    ####

##create the DESeq object from the data we loaded in 
#accounting for interaction between line and environment
ddsF1 <- DESeqDataSetFromMatrix(countData = countsTableRoundF1, 
                              colData = cons, 
                              design = ~line + environment + 
                                line:environment)
##check dimensions
dim(ddsF1)

## filter out genes with too few reads
#keep genes with average > 10 reads per sample
ddsF1 <- ddsF1[rowSums(counts(dds))>160]
dim(ddsF1)


#---------------------------------------------------------#
#                      run DESeq model                      
#---------------------------------------------------------#

##this tests for differential expression
dds <-DESeq(dds)

## list the results:
resultsNames(dds)

####    F1    ####

##this tests for differential expression
ddsF1 <-DESeq(ddsF1)

## list the results:
resultsNames(ddsF1)

#---------------------------------------------------------#
#                      PCA of data                      
#---------------------------------------------------------#
##create vsd object
vsd <- vst(dds, blind = F)

##create data object for PCA
data <- plotPCA(vsd, intgroup = c("line", "environment"), returnData = TRUE)

## create percentVar object for PCA
percentVar <- round(100*attr(data, "percentVar"))

####    F1    ####

##create vsd object
vsdF1 <- vst(ddsF1, blind = F)

##create data object for PCA
dataF1 <- plotPCA(vsdF1, intgroup = c("line", "environment"), returnData = TRUE)

## create percentVar object for PCA
percentVarF1 <- round(100*attr(data, "percentVar"))

####    Plotting side by side PCAs    ####


ggplot(data, aes(PC1, PC2, color = environment, shape = line)) + 
  geom_point(size = 4, alpha = 0.85) + 
  xlab(paste0("PC1: ", percentVar[1], "%variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "%variance")) + 
  theme_bw()+ 
  ggtitle("F3") +
  scale_colour_viridis_d() +
ggplot(dataF1, aes(PC1, PC2, color = environment, shape = line)) + 
  geom_point(size = 4, alpha = 0.85) + 
  xlab(paste0("PC1: ", percentVarF1[1], "%variance")) + 
  ylab(paste0("PC2: ", percentVarF1[2], "%variance")) + 
  theme_bw()+ 
  ggtitle("F1") +
  scale_colour_viridis_d() + 
  plot_annotation(title = "PCA for F1 and F3 generations")+
  plot_layout(ncol = 2, guides = "collect")


#---------------------------------------------------------#
#           order and summarize specific contrasts                      
#---------------------------------------------------------#

#looking at results from the interactions

##create matrix of dds object results + order by p val
resInt <- results(dds, alpha = 0.05)
resInt <- resInt[order(resInt$padj),]

## summarize this
summary(resInt)

#---------------------------------------------------------#
#                 test for effect of env                      
#---------------------------------------------------------#
## run DESeq with line + environment 
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = cons, 
                              design = ~ line + environment)

##compare reduced model using likelihood ratio test
dds <- DESeq(dds, test="LRT", reduced=~line)

## List the results you've generated
resultsNames(dds) 

## Order and list and summarize results from specific contrasts
resEnv <- results(dds, alpha = 0.05)  #signif cut off of 0.05
resEnv <- resEnv[order(resEnv$padj),] #order by adj pval
head(resEnv)

## Summarize the results
summary(resEnv)

## remove NAs - from outliers and low counts
resEnv <- resEnv[!is.na(resEnv$padj),]

## Save genes that are significant
degsEnv <- row.names(resEnv[resEnv$padj < 0.05,]) 
degsEnv[1:5]


####    F1    ####
## run DESeq with line + environment 
ddsF1 <- DESeqDataSetFromMatrix(countData = countsTableRoundF1, 
                                colData = cons, 
                                design = ~ line + environment)

##compare reduced model using likelihood ratio test
ddsF1 <- DESeq(ddsF1, test="LRT", reduced=~line)

## List the results you've generated
resultsNames(ddsF1) 

## Order and list and summarize results from specific contrasts
resEnvF1 <- results(ddsF1, alpha = 0.05)  #signif cut off of 0.05
resEnvF1 <- resEnvF1[order(resEnvF1$padj),] #order by adj pval
head(resEnvF1)

## Summarize the results
summary(resEnvF1)

## remove NAs - from outliers and low counts
resEnvF1 <- resEnvF1[!is.na(resEnvF1$padj),]

## Save genes that are significant
degsEnvF1 <- row.names(resEnvF1[resEnvF1$padj < 0.05,]) 
degsEnvF1[1:5]



#---------------------------------------------------------#
#                 test for effect of line                      
#---------------------------------------------------------#
##create DESeq model focused on line
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = cons, 
                              design = ~ environment + line)

#Likelihood ratio test for removing line
dds <- DESeq(dds, test="LRT", reduced=~environment)
resultsNames(dds)

## Order and list and summarize results from specific contrasts
resLine <- results(dds, alpha = 0.05)
resLine <- resLine[order(resLine$padj),]
head(resLine)

##remove NAs
resLine <- resLine[!is.na(resLine$padj),]

##pull out significant genes
degsline <- row.names(resLine[resLine$padj < 0.05,])
degsline[1:5]

####    F1    ####

##create DESeq model focused on line
ddsF1 <- DESeqDataSetFromMatrix(countData = countsTableRoundF1, 
                              colData = cons, 
                              design = ~ environment + line)

#Likelihood ratio test for removing line
ddsF1 <- DESeq(ddsF1, test="LRT", reduced=~environment)
resultsNames(ddsF1)

## Order and list and summarize results from specific contrasts
resLineF1 <- results(ddsF1, alpha = 0.05)
resLineF1 <- resLineF1[order(resLineF1$padj),]
head(resLineF1)

#remove NAs
resLineF1 <- resLineF1[!is.na(resLineF1$padj),]

#pull out significant genes
degslineF1 <- row.names(resLineF1[resLineF1$padj < 0.05,])
degslineF1[1:5]
#---------------------------------------------------------#
#                  test for effect of int                      
#---------------------------------------------------------#

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = cons, 
                              design = ~ environment + line + environment:line)

dds <- DESeq(dds, test="LRT", reduced=~environment + line)
resultsNames(dds)

resInt <- results(dds, alpha = 0.05)
resInt <- resInt[order(resInt$padj),]
head(resInt)


summary(resInt)

resInt <- resInt[!is.na(resInt$padj),]

degsInt <- row.names(resInt[resInt$padj < 0.05,])
degsInt[1:5]


####    F1    ####
ddsF1 <- DESeqDataSetFromMatrix(countData = countsTableRoundF1, 
                              colData = cons, 
                              design = ~ environment + line +
                                environment:line)

ddsF1 <- DESeq(ddsF1, test="LRT", reduced=~environment + line)
resultsNames(dds)

resIntF1 <- results(ddsF1, alpha = 0.05)
resIntF1 <- resIntF1[order(resIntF1$padj),]
head(resIntF1)


summary(resIntF1)

resIntF1 <- resIntF1[!is.na(resIntF1$padj),]

degsIntF1 <- row.names(resIntF1[resIntF1$padj < 0.05,])
degsIntF1[1:5]

#---------------------------------------------------------#
#                      plot venn diagram                      
#---------------------------------------------------------#

# Total
length(degsEnv)  # 828 <- total # genes sig in env
length(degsline)  # 1645 <- total # genes sig in line
length(degsInt)  # 283 <- total # genes sig in int

# Intersections
length(intersect(degsEnv,degsline))  # 141 <- total # genes sig @ intersect of line + env
length(intersect(degsEnv,degsInt))  # 14 <- total # genes sig @ intersect of int + env
length(intersect(degsInt,degsline))  # 32 <- total # genes sig @ intersect of line + int

intEL <- intersect(degsEnv,degsline) 
length(intersect(degsInt,intEL)) # 7 <- total # genes sig @ intersect of line + env + int using object assigned above

# Number unique
828-14-141-7 # 666   <- total # genes sig and unique to env
1645-141-32-7 # 1465   <- total # genes sig and unique to line
283-14-32-7 # 230  <- total # genes sig and unique to int

#create a "fit" using these seven vals
fit1 <- euler(c("Env" = 666, "Line" = 1465, "Interaction" = 230, "Env&Line" = 141, "Env&Interaction" = 14, "Line&Interaction" = 32, "Env&Line&Interaction" = 7))

#plot these values in a venn diagram
plot(fit1,  lty = 1:3, quantities = TRUE)

#same but transparent
plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))

####    F1    ####
# Total
length(degsEnvF1)  # 448 <- total # genes sig in env
length(degslineF1)  # 226 <- total # genes sig in line
length(degsIntF1)  # 3854 <- total # genes sig in int

# Intersections
length(intersect(degsEnvF1,degslineF1))  # 37 <- total # genes sig @ intersect of line + env
length(intersect(degsEnvF1,degsIntF1))  # 44 <- total # genes sig @ intersect of int + env
length(intersect(degsIntF1,degslineF1))  # 34 <- total # genes sig @ intersect of line + int

intEL <- intersect(degsEnvF1,degslineF1) <- total # genes sig @ intersect of line + env
length(intersect(degsIntF1,intEL)) # 7 <- total # genes sig @ intersect of line + env + int using object assigned above

# Number unique
448-44-37-7 # 360   <- total # genes sig and unique to env
226-37-34-7 # 148   <- total # genes sig and unique to line
3854-44-34-7 # 3769  <- total # genes sig and unique to int

#create a "fit" using these seven vals
fit1 <- euler(c("Env" = 360, "Line" = 148, "Interaction" = 3769, "Env&Line" = 37, "Env&Interaction" = 44, "Line&Interaction" = 34, "Env&Line&Interaction" = 7))

#plot these values in a venn diagram
plot(fit1,  lty = 1:3, quantities = TRUE)

#same but transparent
plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))


####    which genes overlap bx generations    ####
## interaction
length(intersect(degsIntF1,degsInt)) #105
head(resInt[intersect(degsIntF1,degsInt),])
head(resIntF1[intersect(degsIntF1,degsInt),])


length(which(resInt[intersect(degsIntF1,degsInt),]$log2FoldChange <0))
length(which(resInt[intersect(degsIntF1,degsInt),]$log2FoldChange >0))

## Environment
length(intersect(degsEnvF1,degsEnv)) #54
head(resInt[intersect(degsEnvF1,degsEnv),])
head(resIntF1[intersect(degsEnvF1,degsEnv),])

## line
length(intersect(degslineF1,degsline)) #49
head(resInt[intersect(degslineF1,degsline),])
head(resIntF1[intersect(degslineF1,degsline),])

#combination
intE <- intersect(degsEnvF1,degsEnv)
intI <- intersect(degsIntF1,degsInt)
intL <- intersect(degslineF1,degsline)
intEL <- intersect(intE,intL)

length(intersect(intE,intL)) #1
length(intersect(intI,intL))#1
length(intersect(intE,intI)) #0
length(intersect(intEL,intI)) #0

# log fold change and base mean for these is too low to consider
resInt[intersect(intE,intL),]
resIntF1[intersect(intE,intL),]

resInt[intersect(intI,intL),]
resIntF1[intersect(intI,intL),]






