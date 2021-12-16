#!/bin/bash

#set path and filename

VCF=/data/project_data/PopGenomics/poplar_hybrids.maf05.vcf.gz


cd~/

#the actual filtering:
	#says keep those with
	# 50 = 50,000 BP window
	# 10 = 10 bp sliding scale
	# 0.1 = maximum LD r2 value as a cut off
plink2 --vcf $VCF \
--threads 3 \
--double-id \
--allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out poplar_hybrids_noLD

#create new directory
mkdir Admixture 

FILE=poplar_hybrids

#make bed - basically converts to necessary format
plink2 --vcf $VCF \
--threads 3 \
--allow-extra-chr \
--make-bed \
--out Admixture/$FILE 

#extract files trimmed of linkage based on LD cut off above
plink2 --bfile Admixture/$FILE \
--threads 3 \
--allow-extra-chr \
--set-missing-var-ids @:# \
--extract poplar_hybrids_noLD.prune.in \
--make-bed \
--out Admixture/$FILE.LDpruned 

# Replace column 1 Chr #'s with 0's, since ADMIXTURE doesn't like them
cd Admixture


FILE2=poplar_hybrids.LDpruned

awk '{$1=0;print $0}' $FILE2.bim > $FILE2.bim.tmp
mv $FILE2.bim.tmp $FILE2.bim

# Run Admixture 

#variable of K assigned in class
# this can also be done in a loop to do multiple values of K to see impact of K
K=8


# j3 = # of cpus you will run
admixture -j3 --cv $FILE2.bed $K >log${K}.out

