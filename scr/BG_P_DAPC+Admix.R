#DAPC and Admixture Plot 

#packages 
require(vcfR)
require(adegenet)
require(readxl)


#import excel with pops 
df <- read_excel("data_raw/GR Practice.xlsx")

#8000000
vcf8 <- read.vcfR("data_raw/variants.filt.8000000.vcf/variants.filt.8000000.vcf")

#convert to genind 
vcf8genind <- vcfR2genind(vcf8)
vcf8genind

#insert info into genind (metadata has population info)
vcf8genind@pop <- as.factor(df$sample_site)

#run a DAPC 
dapc8_result <- dapc(vcf8genind)
7
2

scatter(dapc8_result)

scatter(dapc8_result, cex = (2), legend = TRUE, 
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2)


#admixture plot 
compoplot(dapc8_result, posi.da = TRUE) 


#16000000
#import data 
vcf16 <- read.vcfR("data_raw/variants.filt16000000.vcf/variants.filt.vcf")

#convert to genind 
vcfgenind16 <- vcfR2genind(vcf16)
vcfgenind16

#insert info into genind (metadata has population info)
vcfgenind16@pop <- as.factor(df$sample_site)

#run a DAPC 
dapc_result16 <- dapc(vcfgenind16)
2
2

scatter(dapc_result16)

scatter(dapc_result16, cex = (2), legend = TRUE, 
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2)


#admixture plot 
compoplot(dapc_result16, posi.da = TRUE) 

