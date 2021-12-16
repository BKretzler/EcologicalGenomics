#!/bin/bash

# Rename value in <> to your chromosome number! Include the zero only if your chromosome # is <10
# if single digit must be 0x, must be capitalized Chrxx

myChr=Chr18 

cd /data/project_data/PopGenomics

# Run VCFtools to subset the big vcf file for just your chromosome
#opens masterfile and parses by chromosome
#Outputs VCF file for chromosome into shared file 

vcftools --gzvcf poplar_hybrids.maf05.vcf.gz \
--chr $myChr \
--out shared/$myChr \
--recode



# Extract the centromere coordinates for your chromosome so you can exclude those regions from your sweep analysis

#grep the centromeres from your chromosome
grep $myChr poplar_centromeres.txt > shared/${myChr}_centromere.txt 
# grab the centromere location for your chromosome

#cd to shared folder
cd shared/

# make a new directory for your chromosome analyses
mkdir ${myChr}_sweeps  

# clean up the space by moving all files into your the directory you just made
# ok if it says can't move directory into self
mv *${myChr}* ${myChr}_sweeps 


#cd into new folder
cd ${myChr}_sweeps


## test for selective sweeps
RAiSD -n $myChr \
-I ${myChr}.recode.vcf \
-f -t -R -P -D -A 0.99 \
-X ${myChr}_centromere.txt


# Estimate nucleotide diversity (pi) in sliding windows of 50kb

vcftools --vcf ${myChr}.recode.vcf \
--chr $myChr \
--window-pi 50000 \
--out $myChr


# First, need to subset the metadata file for just those individuals with balsamifera ancestry
# We can do this using an interactive R session at the commandline. 
# An alternative is to put these R commands in a script, save it with the ".r" extension, 
# and at the commandline type "Rscript myscript.r"

R # Opens an interactive R session within Unix...
Qscores <- read.table("../poplar_hybrids.LDpruned.5.Q", sep=" ",header=F)
names(Qscores) = c("K1","K2","K3","K4","K5")

meta <- read.table("../../Combined_Transect_Sampling_Data_2020.txt",sep="\t",header=T)

merged <- cbind(meta,Qscores)
str(merged)

Bals_Inds <- merged[which(merged$K4>0.5),1]  
length(Bals_Inds) # Should net you 188 individuals

Tricho_Inds <- merged[which(merged$K4<=0.5),1]
length(Tricho_Inds) # Should net you 388 individuals

# Write out your Bals and Tricho lists as tab-delimited text files
write.table(Bals_Inds, "Bals_Inds.txt", quote=F, row.names=F, col.names=F)

write.table(Tricho_Inds, "Tricho_Inds.txt", quote=F, row.names=F, col.names=F)

quit()

# When prompted with: "Save workspace image? [y/n/c]"  choose: n
# R files can be ran using Rscript.filename (i think) 
# must have libraries installed


# Calculate Fst between Balsam and Tricho using sliding windows of 50kb

vcftools --vcf ${myChr}.recode.vcf \
--weir-fst-pop Bals_Inds.txt \
--weir-fst-pop Tricho_Inds.txt \
--fst-window-size 50000 \
--out Bals_Tricho_All