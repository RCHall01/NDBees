
#packages 
require(vcfR)
require(adegenet)
require(readxl)
require(dplyr)
require(poppr)
require(PopGenReport)

#import excel with pops 
df <- read_excel("data_raw/GR Practice.xlsx")

#8000000
Tern <- read.vcfR("data_raw/variants.filt.8000000.vcf/variants.filt.8000000.vcf")

#convert to genind 
Tern.genind <- vcfR2genind(Tern)
summary(Tern.genind)

#insert info into genind (metadata has population info)
Tern.genind@pop <- as.factor(df$sample_site)

saveRDS(Tern.genind, file = "Tern_genind.rds")

######HWE#####
Tern.genind <- readRDS("outputs/Tern_genind.rds")

#check for HWE 
HWE <-round(pegas::hw.test(Tern.genind, B = 1000), digits = 3)

# Chi-squared test: p-value
HWE.test <- data.frame(sapply(seppop(Tern.genind),
                              function(ls) pegas::hw.test(ls, B=0)[,3]))
HWE.test.chisq <- t(data.matrix(HWE.test))
{cat("Chi-squared test (p-values):", "\n")
  round(HWE.test.chisq,3)}

# Monte Carlo: p-value
HWE.test <- data.frame(sapply(seppop(Tern.genind), 
                              function(ls) pegas::hw.test(ls, B=1000)[,4]))
HWE.test.MC <- t(data.matrix(HWE.test))
{cat("MC permuation test (p-values):", "\n")
  round(HWE.test.MC,3)}

alpha=0.05 # /96
Prop.loci.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 2, mean), 
                                   MC=apply(HWE.test.MC<alpha, 2, mean))
Prop.loci.out.of.HWE             

Prop.pops.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean), 
                                   MC=apply(HWE.test.MC<alpha, 1, mean))
Prop.pops.out.of.HWE             

#false discovery
Chisq.fdr <- matrix(p.adjust(HWE.test.chisq,method="fdr"), 
                    nrow=nrow(HWE.test.chisq))
MC.fdr <- matrix(p.adjust(HWE.test.MC, method="fdr"), 
                 nrow=nrow(HWE.test.MC))

Prop.pops.out.of.HWE2 <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean), 
                                   MC=apply(HWE.test.MC<alpha, 1, mean),
                                   Chisq.fdr=apply(Chisq.fdr<alpha, 1, mean),
                                   MC.fdr=apply(MC.fdr<alpha, 1, mean))
Prop.pops.out.of.HWE             


#save
write.table(HWE, file = "HWE_result.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(round(HWE.test.chisq, 3), file = "HWE_chi_squared_pvalues.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(round(HWE.test.MC, 3), file = "HWE_MC_permutation_pvalues.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(Prop.loci.out.of.HWE, file = "Prop_loci_out_of_HWE.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(Prop.pops.out.of.HWE, file = "Prop_pops_out_of_HWE.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(Prop.pops.out.of.HWE2, file = "Prop_pops_out_of_HWE_with_FDR.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)



######LD######
# Perform index of association
ia_result <- poppr::ia(Tern.genind, sample = 199)
write.table(ia_result, file = "ia_result.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

# Calculate linkage disequilibrium
LD.pair <- poppr::pair.ia(Tern.genind)
write.table(LD.pair, file = "LD_pair.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)


####Genetic Diversity####
# Assess allelic richness
Sum <- adegenet::summary(Tern.genind)
barplot(Sum$pop.n.all, las = 3, xlab = "", ylab = "Number of alleles")
plot(Sum$n.by.pop, Sum$pop.n.all, xlab = "Sample size", ylab = "Number of alleles")
abline(lm(Sum$pop.n.all ~ Sum$n.by.pop), col = "red")
write.table(Sum, file = "allelic_richness_summary.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

Richness <- PopGenReport::allel.rich(Tern.genind, min.alleles = NULL)
write.table(Richness$alleles.sampled, file = "alleles_sampled.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
barplot(Richness$mean.richness, las = 3, ylab = "Rarefied allelic richness (Ar)")
plot(colMeans(Richness$pop.sizes), Richness$mean.richness, xlab = "Valid sample size", ylab = "Rarefied allelic richness (Ar)")
abline(lm(Richness$mean.richness ~ colMeans(Richness$pop.sizes)), col = "red")
write.table(Richness$mean.richness, file = "allelic_richness_mean.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

# Observed and expected heterozygosity
Summ <- summary(Tern.genind)
barplot(Summ$Hexp, ylim = c(0, 1), ylab = "Expected heterozygosity")
barplot(Summ$Hobs, ylim = c(0, 1), ylab = "Observed heterozygosity")
Hobs <- t(sapply(seppop(Tern.genind), function(ls) summary(ls)$Hobs))
Hexp <- t(sapply(seppop(Tern.genind), function(ls) summary(ls)$Hexp))
write.table(Hexp, file = "expected_heterozygosity.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(Hobs, file = "observed_heterozygosity.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

# Combine results into a data frame
Tern.diversity <- data.frame(
  Pop = names(Hobs.pop),
  n = Sum$n.by.pop,
  Hobs = Hobs.pop,
  Hexp = Hexp.pop,
  Ar = Richness$mean.richness
)
write.table(Tern.diversity, file = "tern_diversity_summary.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

# Convert genind to genpop and calculate allele frequencies
Tern.genpop <- adegenet::genind2genpop(Tern.genind)
Freq <- adegenet::makefreq(Tern.genpop)
write.table(round(Freq, 2), file = "allele_frequencies.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
apply(Freq, MARGIN = 1, FUN = sum) # Just checking
