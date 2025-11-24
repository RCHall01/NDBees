#8/30/2025 -> Reassign populations unmerged 
#packages
require(vcfR)
require(tidyverse)
require(adegenet)

#import excel with pops 
df <- read.csv("/mmfs1/home/rhiannon.hall/bee_samples/Nesting/data_raw/BOGR/BOGR_samples.csv") #read in 
BOGR_samples <- df[!(df$sample_id %in% c("BOGR_0031", "BOGR_0085", 
                                         paste0("BOGR_", sprintf("%04d", 1:8)))), ] #remove two related siblings and preliminary  
BOGR_samples$sample_id <- gsub("_", ".", BOGR_samples$sample_id)
write.csv(BOGR_samples, "data/BOGR/final_BOGR.csv", row.names = FALSE)

#import vcf
BOGR <- read.vcfR("data_raw/BOGR/BOGR_obj13.vcf")
BOGR_genind <- vcfR2genind(BOGR)

# Step 1: Access current sample names and trim
current_names <- indNames(BOGR_genind)
trimmed_names <- sub(".*BOGR-(\\d{4}).*", "BOGR.\\1", current_names)
indNames(BOGR_genind) <- trimmed_names
indNames(BOGR_genind)

#Step 3: Assign Pops
BOGR_genind@pop <- as.factor(BOGR_samples$sample_site)

#Step 4: Check your work 
individual_names <- indNames(BOGR_genind)
population_assignments <- pop(BOGR_genind) # Extract the population assignments
# Create a data frame to show each individual and their assigned population
individual_population_df <- data.frame(Sample = individual_names, 
                                       Population = as.character(population_assignments))
print(individual_population_df) # View the data frame

#save genind
saveRDS(BOGR_genind,"data/BOGR/BOGR_obj13.rds")

