library(ggplot2)


setwd("ecologicalgenomics/groupProject/")

dSuite<- read.table("dsuite_run_BBAA.txt", header = T)
dsuite_tree <- read.table("dsuite_run_filtered_tree.txt", header = T)


significant<- dsuite_tree[which(dsuite_tree$p.value < 0.05/56),] # no z scores are significant..

significant$sharedLoci = significant$BBAA + significant$ABBA + significant$BABA

write.table(significant, file = "dsuite_tree_filtered_significant.txt", row.names =F)

write.csv(significant, file = "dsuite_tree_filtered_significant.csv", row.names = F)

significant[which.max(significant$sharedLoci),]

range(significant$sharedLoci)
range(significant$Dstatistic)

significant[which.max(significant$Dstatistic),]




cor(dSuite$Dstatistic, dSuite$f4.ratio) # highly positively correlated, as d increases so does f4


significant$domPattern = NA

for(row in 1:nrow(significant)){
  if(significant[row,"ABBA"]<significant[row,"BABA"]){
    significant[row,'domPattern'] = "BABA"
  } else {
    significant[row,'domPattern'] = "ABBA"  
  }
}


hist(significant$Dstatistic)
hist(dsuite_tree$Dstatistic)

significant[which(significant$Dstatistic > 0.1),]

ggplot(significant, aes(x = Dstatistic)) + 
  geom_histogram(bins = 10, colour = "black", fill = 'darkolivegreen4') + 
  theme_bw() + 
  ggtitle("Histogram of Significant D-statistics") + 
  xlab("D-statistic") + 
  ylab("Frequency")



