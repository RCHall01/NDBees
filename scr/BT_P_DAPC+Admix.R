#DAPC and Admixture Plot 

#packages 
require(vcfR)
require(adegenet)
require(readxl)


#import excel with pops 
df_BT <- read_excel("data_raw/Bombus Ternarius Prelim Samples.xlsx")

#8000000
vcf8_BT <- read.vcfR("data_raw/BOTR_800000_filt/variants.filt.vcf")

#convert to genind 
vcf8genind_BT <- vcfR2genind(vcf8_BT)
vcf8genind_BT

#insert info into genind (metadata has population info)
vcf8genind_BT@pop <- as.factor(df_BT$sample_site)

#run a DAPC 
dapc8_result_BT <- dapc(vcf8genind_BT)
7
2

scatter(dapc8_result_BT)

scatter(dapc8_result_BT, cex = (2), legend = TRUE, 
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2)


#admixture plot 
compoplot(dapc8_result_BT, posi.da = TRUE) 


# Extract the data frame of the VCF
vcf_df_BT <- as.data.frame(getFIX(vcf8_BT))

# Count SNPs where both REF and ALT alleles are single nucleotides
snp_count_BT <- sum(nchar(vcf_df_BT$REF) == 1 & nchar(vcf_df_BT$ALT) == 1)

print(paste("Number of SNPs:", snp_count_BT))
