library(maps)  
library(plotrix) 
library(ggplot2)

setwd("C:/Users/kretz/R/EcologicalGenomics/ecologicalgenomics")  #set the path to where your downloaded files are on your laptop

meta <- read.table("Combined_Transect_Sampling_Data_2020.txt", sep="\t",header=T)  

Qscores <- read.table("poplar_hybrids.LDpruned.8.Q", sep=" ",header=F)
names(Qscores) <- c("K1","K2", "K3","K4","K5","K6", "K7","K8")  # Customize to your level of K!

tiff("Admix.tiff", width = 10, height = 7, units = "in",
     res = 300)
map("world",xlim=c(-160,-100),ylim=c(40,70),fill=T, col="lightgrey")

title("Admixture K=8") # Change for your value of K
map.axes(cex.axis=1.5)
for(i in 1:nrow(meta)){
  floating.pie(meta$Longitude[i], meta$Latitude[i], 
               c(Qscores$K1[i],Qscores$K2[i],Qscores$K3[i],Qscores$K4[i],
                 Qscores$K5[i],Qscores$K6[i], Qscores$K7[i],
                 Qscores$K8[i]),
               col=c("yellow","green","blue", "purple", 
                     "red", "orange", "black", "grey"),
               radius=0.5)
}

# Customize the above to your level of K wherever you see the ellipses (...)!

dev.off()


#------------------load in data-----------------------#

Qscores_new <- read.table("poplar_hybrids.LDpruned.5.Q", sep=" ",header=F) # Rename to the new file!
names(Qscores_new) <- c("K1","K2", "K3", "K4", "K5")  


#------------------shannon diversity-----------------------#

# Calculate Shannon Diversity across the K different Q-scores per individual
  # high values of shannon diversity have more slices to their pie chart
  # more representation of k groups

K=5  # Change X to the level of K we're investigating

tmp=numeric()

for(i in 1:nrow(Qscores_new)){
  for(j in 1:K){
    tmp[j] = Qscores_new[i,j]*log(Qscores_new[i,j])
  }
  Qscores_new$ShDiv[i] = -1*sum(tmp)
}

#------------------read in het results-----------------------#

## in the same order as ancestry and metadata file
  #obs homozy - how many are homozygous
  #expected - expected homozygous under hardy weinberg ()
  # number of sites
  #F = fixation index, 1 - observed hetero/expected he
het <- read.table("het_maf05.het", sep="\t",header=T)

str(het) # What's in this dataframe?  


# Combine the meta data, heterozygosity, and admixture data

het2 <- cbind(meta,het,Qscores_new) #bind the het results with the meta data

# How does F vary within each transect?
  # closer to 0  = more admixture

ggplot(het2, aes(x=Transect, y=F, color=Transect)) +
  geom_dotplot(binaxis='y', binwidth=0.01, stackdir='center')

  # lots of variation in the transects
  # dots = individual trees
  # we see some variation in degree of admixture along the transect!

# Plot admixture diversity spatially
ggplot(het2, aes(x=Longitude, y=Latitude, color=ShDiv)) +
  geom_point(size=4, shape=20)
  # we see less diversity in the more "pure balsam" regions of the transect

# And finally, are more admixed individuals more heterozygous in their genomes?

#higher hetero = lower level of F
ggplot(het2, aes(x=F, y=ShDiv, color=Transect)) +
  geom_point(size=4, shape=20)
  ## higher levels of ancestry diversity ~ higher levels of heterozygousity
cor.test(het2$F,het2$ShDiv)

