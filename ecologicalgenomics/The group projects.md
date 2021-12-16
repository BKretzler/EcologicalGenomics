### The group projects:

* we worked on 3 different data sets:
  * Microbiomes of sea stars impacted and unimpacted by sea star wasting disease
    * in this module we got familiar with the command line and used Qiime
    * linux -> local machine to visualize
  * RNA seq of copepods exposed to different environments (high temp co2, ambient)
    * experimentally evolved then reciprocally transplanted 
    * cleaned in command line then moved to local machine and analyzed + visualized in R
  * Population genomics on poplar trees along climatic gradients + across hybrid zones
    * used snps that were called for us to look at chromosomal landscape
    * paired down to those not in linkage disequilibrium 
    * identified level of K - degree of population structure
    * looked at admixture - could look at impact of spatial location on admixture
    * mapped ancestry along the chromosome and related to phenotypes - in homework we will do this with climate
    * assumed model as additive and that all transects were similar in their hybridization
      * but what if they weren't?
    * poplar data did not have phylogenetic context - we could add this 
      * other poplar genomes exist
* What other databases do people have that are of interest?
  * can we use the hudsonica data? 
    * yes but risky and would need to be done off the server
* what thoughts do we have about the projects?
  * Dan found an interesting cannabis data set
* So what will we do?
  * poplar phylogenetics and inferring ancient introgression
  * so we need to start by finding other populus data sets
    * check genbank
  * can we integrate a spatial component
  * we will need:
    * populus balsamifera, trichocarpa, deltoides, and angustifolia
      * also nigra, fermontii, euph..
    * need to gather the data for this
    * some Chinese species that are part of this group:
      * tomentosa
  * US department of energy + JGI run phytozome
  * outgroup - salix purpurea 
  * The challenge with this will be that assemblies are not comprable
    * we need to do assembly of reference sequences
    * do we just want to target coding regions? or whole genome alignment
    * how to pull out coding regions?
      * no go to pipelines
      * have to look at annotations for genes + gene names, coordinates
  * for ancient introgression: 
    * need large enough sample of genes
    * should be random and unbiased 
  * we need to define what data we want before we get all kinds of data
    * start with what we want to know and how we will do it - what stats, etc. 
      * this determines what we need to do/ what data we need
    * phylonetworks? - proportion of genome/ancestral contribution
    * d stat
* What did steve mean by "spatial context"
  * estimation of migration surfaces
    * landscape genomics data to examine hot and cool spots of gene flow
      * examine barriers to gene flow and integression
  * can also integrate environment component
    * how is geographic region impacting the hybridization (n-s,e-w)
    * Niche modeling
* what populus genomes are available on NCBI:
  * tomentosa
  * euphretica
  * deltoides
  * davidiana (new)
  * trichocarpa
  * NOT balsamifera
* we need to make sure the ones we need are available - lots of asian populus 
* could also mine RNAseq data and align sequences from this
  * would help with aligning assemblies
  * SRA - is the way to search for transcriptomes
* we want fasta files for data of interest:
  * map each species to trichocarpa genome to get alignments
* whole genomes on gen bank may or may not be in chromosomal assemblies
  * so we just need to start with genome alignment to trichocarpa
  * 