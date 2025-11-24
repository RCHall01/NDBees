#9/13/2025 -> Reassign populations  
#packages
require(vcfR)
require(tidyverse)
require(adegenet)

#import excel with pops 
df <- read.csv("/mmfs1/home/rhiannon.hall/bee_samples/Nesting/data_raw/BOTE/BOTE_samples.csv")
BOTE_samples <- df[!(df$sample_id %in% c("BOTE_0078", "BOTE_0112", 
                                         paste0("BOTE_", sprintf("%04d", 1:12)))), ] #remove two related siblings and preliminary  
BOTE_samples$sample_id <- gsub("_", ".", BOTE_samples$sample_id)
write.csv(BOTE_samples, "data/BOTE/final_BOTE.csv", row.names = FALSE)

#import vcf
BOTE <- read.vcfR("data_raw/BOTE/BOTE_obj13.vcf")
BOTE_genind <- vcfR2genind(BOTE)

# Step 1: Access current sample names and trim
current_names <- indNames(BOTE_genind)
trimmed_names <- sub(".*BOTE-(\\d{4}).*", "BOTE.\\1", current_names)
indNames(BOTE_genind) <- trimmed_names
indNames(BOTE_genind)

#Step 3: Assign Pops
BOTE_genind@pop <- as.factor(BOTE_samples$sample_site)

#Step 4: Check your work 
individual_names <- indNames(BOTE_genind)
population_assignments <- pop(BOTE_genind) # Extract the population assignments
# Create a data frame to show each individual and their assigned population
individual_population_df <- data.frame(Sample = individual_names, 
                                       Population = as.character(population_assignments))
print(individual_population_df) # View the data frame

#save genind
saveRDS(BOTE_genind,"data/BOTE/BOTE_genind13.rds")

