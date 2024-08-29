
#load in packages 
library(adegenet)
library(gstudio)
library(tibble)
library(here)
library(vcfR)
library(utils)
library(hierfstat)
library(pegas)


# Load genind object from the saved file
SNP_genind <- readRDS("C:/Users/Rhall/OneDrive - North Dakota University System/Masters/Thesis Genetics/Practice/outputs/SNP_genind.rds")

#---- Calculating Basic Pop Gen Stats from SNP data ------
#number of alleles per locus 
nAll(SNP_genind) 

head(SNP_genind)

summary(SNP_genind)

SNP_hierfstat <- genind2hierfstat(SNP_genind)
#genetic diversity 
div <- summary(SNP_genind)
div
names(div)
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")
bartlett.test(list(div$Hexp, div$Hobs)) # a test : H0: Hexp = Hobs
basicstat <- basic.stats(SNP_hierfstat, diploid = TRUE, digits = 2) 
names(basicstat)   
boot.ppfis(SNP_hierfstat) 
x <- indpca(SNP_hierfstat) 
plot(x, cex = 0.7)
#Testing for HWE
hw.test(SNP_genind, B = 1000)


#---- Calculating genetic differentiation and clustering from SNP data -----
#observed and expected heterogeneity Fst 
basic.stats(SNP_genind) 
wc(SNP_genind)# Weir and Cockerham's estimate
  #Hierarchical Fst tests (= AMOVA for SNP)
loci <- SNP_hierfstat[, -1] # Remove the population column
varcomp.glob(levels = data.frame(population, county), loci, diploid = TRUE) 
test.g(loci, level = population) 
test.between(loci, test.lev = population, rand.unit = county, nperm = 100) 
  #Pairwise Fst 
genet.dist(SNP_genind, method = "WC84")
  #unsupervised clustering 
# using Kmeans and DAPC in adegenet 
set.seed(991) # Setting a seed for a consistent result
grp <- find.clusters(SNP_genind, max.n.clust = 10, n.pca = 20, choose.n.clust = FALSE) 
names(grp)
grp$grp
dapc1 <- dapc(SNP_genind, grp$grp, n.pca = 20, n.da = 6) 
scatter(dapc1) # plot of the group

#---- Individual Based Genetic Distance -----
library("poppr")
library("pegas")
library("ape")
library("adegenet")
library("ade4")
# Ind. genetic distance: euclidean distance 
distgenEUCL <- dist(SNP_genind, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
hist(distgenEUCL)
# Ind. genetic distance: Number of loci for which individuals differ
distgenDIFF <- dist.gene(SNP_hierfstat, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
hist(distgenDIFF)
# Get percent missing data per population
missing_data <- info_table(SNP_genind, type = "missing")
sum(missing_data["Total", 1:100] > 0)
barplot(missing_data["Total", 1:100], xlab = "Locus", ylab = "Percent Missing")
distgenDIFF <- dist.gene(SNP_hierfstat, method="pairwise", pairwise.deletion = TRUE, variance = FALSE)
hist(distgenDIFF)
# Number of allelic differences bbetween two individuals 
distgenDISS <- diss.dist(SNP_genind, percent = FALSE, mat = FALSE)
hist(distgenDISS)
#conclusions drawn from the analysis
boxplot(distgenEUCL, distgenDIFF, distgenDISS)

#---- Detection of the signal of selection from genome scan ----
library("devtools")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install_github("whitlock/OutFLANK")
install.packages("pcadapt")
# Open packages 
library("pcadapt")
library("qvalue")
library("OutFLANK")
library("ggplot2")

