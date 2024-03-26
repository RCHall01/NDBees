

library(adegenet)
library(gstudio)
library(tibble)
library(here)
library(vcfR)
library(utils)


#download data
# Specify the full path to your VCF file on your laptop
vcf_path_local <- "C:/Users/Rhall/OneDrive - North Dakota University System/Masters/Thesis Genetics/Practice/data_raw/samtools_bif_annotated.vcf"


# Read the VCF file using read.vcfR
vcf <- read.vcfR(vcf_path_local, verbose = FALSE)

#convert to genind
SNP_genind <- vcfR2genind(vcf)

#import text file 

txt_path_local <- "C:/Users/Rhall/OneDrive - North Dakota University System/Masters/Thesis Genetics/Practice/data_raw/Jackson_ME_specimen_info.txt"
txt <- read.table(txt_path_local, header = TRUE)

#isolate columns/rows needed
lat <- txt[1:383,5]
long <- txt[1:383,6]
pop <- txt[1:383,7]

#add to genind object 
pop(SNP_genind) <- pop


summary(SNP_genind)

saveRDS(SNP_genind, file = "SNP_genind.rds")

