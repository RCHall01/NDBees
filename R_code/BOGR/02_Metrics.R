#8/30/25 -> Genetic Metrics 
#packages
require(adegenet)
require(hierfstat)
require(PopGenReport)
require(reshape2)
require(tidyverse)

# Load the genind object
BOGR_genind <- readRDS("data/BOGR/BOGR_13.rds")

##### Filtering for Polymorphic Loci #####
# Check polymorphism status of loci
polymorphic_summary <- summary(isPoly(BOGR_genind)) # Check summary for polymorphic sites

# Identify loci that are polymorphic
poly_loci <- names(which(isPoly(BOGR_genind) == TRUE))

# Filter genind object for polymorphic loci
BOGR_genind_poly <- BOGR_genind[loc = poly_loci]

####stats####
# Calculate basic stats for the filtered genind object
Bbasic <- basic.stats(BOGR_genind_poly, diploid = TRUE)
saveRDS(Bbasic, file = "data/BOGR/basic_stats_results.rds")
basic <- readRDS("data/BOGR/basic_stats_results.rds")

#n.ind.samp
n_ind_samp <- apply(basic$n.ind.samp, MARGIN = 2, FUN = max, na.rm = TRUE)  

# Mean observed heterozygosity per site
BHo_gen <- apply(basic$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)

# Mean expected heterozygosity per site
BHe_gen <- apply(basic$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)

# Mean Fis per site
Fis_gen <- apply(basic$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)
boot <- boot.ppfis(BOGR_genind_poly)

# Create a data frame for heterozygosity data
Het_df <- data.frame(
  Site = names(BHo_gen),
  'n' = n_ind_samp,
  Ho = BHo_gen,
  He = BHe_gen,
  Fis = Fis_gen)

Het_df

# Write Het_df to an Excel file
writexl::write_xlsx(Het_df, "outputs/BOGR/BOGR_metrics.xlsx")

######FST#####
require(adegenet)
require(hierfstat)
require(PopGenReport)
require(reshape2)
require(tidyverse)

# Load the genind object
BOGR_genind <- readRDS("data/BOGR/BOGR_genind.rds")

# Check polymorphism status of loci
polymorphic_summary <- summary(isPoly(BOGR_genind)) # Check summary for polymorphic sites

# Identify loci that are polymorphic
poly_loci <- names(which(isPoly(BOGR_genind) == TRUE))

# Filter genind object for polymorphic loci
BOGR_genind_poly <- BOGR_genind[loc = poly_loci]

#Fst 
wc_result <- wc(BOGR_genind_poly)
global_Fst <- wc_result$FST
print(global_Fst) #0.0002631161
saveRDS(wc_result, file = "data/BOGR/WC_results.rds")
wc_result <- readRDS("data/BOGR/WC_results.rds")

#### Pairwise Fst calculation####
pairwise_Fst <- genet.dist(BOGR_genind_poly, method = "WC84")
saveRDS(pairwise_Fst, file = "data/BOGR/pairwise_Fst_results.rds")
print(pairwise_Fst)

pairwise_Fst <- readRDS("data/BOGR/pairwise_Fst_results.rds")
fst_matrix <- as.matrix(pairwise_Fst)
fst_table <- as.data.frame(fst_matrix)
fst_table <- cbind(Population = rownames(fst_table), fst_table)
write.csv(fst_table, "outputs/BOGR/pairwise_Fst_table_wide.csv", row.names = FALSE)
