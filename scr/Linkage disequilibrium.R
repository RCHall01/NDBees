# Load required libraries
library(readxl)
library(vcfR)
library(adegenet)

# Define a function to run PLINK with given options
runPLINK <- function(PLINKoptions = "") {
  system(paste("plink_win64_20231211/plink.exe", PLINKoptions))
}

# Define the paths
vcf_file <- "data_raw/variants.filt.8000000.vcf/variants.filt.8000000.vcf"
plink_prefix <- "data/vcf8genind"

# Convert VCF to PLINK format
plink_command <- paste("--vcf", vcf_file, "--make-bed --out", plink_prefix)
runPLINK(plink_command)

#allow extra chromosome 
plink_command <- paste("--vcf", vcf_file, "--make-bed --out", plink_prefix, "--allow-extra-chr")
runPLINK(plink_command)

# Define the output file for LD results
ld_output <- "data/vcf8genind_ld"

# Run PLINK LD analysis
plink_ld_command <- paste("--bfile", plink_prefix, "--r 0.2 --ld-window 99999 --ld-window-kb 1000 --out", ld_output)

runPLINK(plink_ld_command)

plink_ld_command <- paste("--bfile", plink_prefix, "--r2 --out", ld_output, "--allow-extra-chr")

# Read the PLINK LD results
ld_results_file <- paste(ld_output, ".ld", sep="")
if (!file.exists(ld_results_file)) {
  stop("PLINK LD results file not found. Please check if PLINK ran successfully.")
}
ld_results <- read.table(ld_results_file, header=TRUE)

# View the LD results
head(ld_results)



