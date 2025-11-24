#8/30/2025 -> PCA visualization STACKS 
#packages
require(adegenet)
require(tidyverse)
require(ggplot2)
require(FactoMineR)

###Create obj 2 genind ####
#pull in data
require(vcfR)

#import excel with pops 
df <-read_csv("data/BOGR/final_BOGR.csv")

#import vcf
BOGR <- read.vcfR("data_raw/BOGR/BOGR_obj2.vcf")
BOGR_obj2 <- vcfR2genind(BOGR)

# Step 1: Access current sample names and trim
current_names <- indNames(BOGR_obj2)
trimmed_names <- sub(".*BOGR-(\\d{4}).*", "BOGR.\\1", current_names)
indNames(BOGR_obj2) <- trimmed_names
indNames(BOGR_obj2)

#Step 3: Assign Pops
BOGR_obj2@pop <- as.factor(df$sample_site)

#Step 4: Check your work 
individual_names <- indNames(BOGR_obj2)
population_assignments <- pop(BOGR_obj2) # Extract the population assignments
# Create a data frame to show each individual and their assigned population
individual_population_df <- data.frame(Sample = individual_names, 
                                       Population = as.character(population_assignments))
print(individual_population_df) # View the data frame

#save genind
saveRDS(BOGR_obj2,"data/BOGR/BOGR_obj2.rds")

###Run population structure####
BOGR_genind <- readRDS("data/BOGR/BOGR_obj2.rds")
#run pca
X <- tab(BOGR_genind, NA.method = "mean")
pca_result <- dudi.pca(X, scannf = FALSE, nf = 3)

#pulling out percentage of variance explained 
pca_result$eig
summary(pca_result)

# Create a named vector for renaming populations
rename_populations <- c(
  "Adams_Blaha_Wolf_Butte" = "Adams",
  "Agnes_Marsh_WPA" = "Agnes Marsh WPA", 
  "BIL-144-099-36" = "BIL",
  "Cass_Erie_Dam_WMA" = "Cass Erie Dam WPA",
  "Cavalier_Thorson_WPA" = "Cavalier Thorson WPA", 
  "EDD-150-065-24" = "EDD", 
  "Emmons_Schiermeister_WPA" = "Emmons Schiermeister WPA",
  "Fitzsimond_Slough_WMA" = "Fitzsimond Slough WMA",
  "FOS-147-062-26" = "FOS", 
  "Griggs_Evers_WPA" = "Griggs Evers WPA", 
  "Hettinger_Hillsview_WPA" = "Hettinger Hillsview WPA", 
  "McHenry_Taylor_Central " = "McHenry Taylor Central",
  "Morton_Graner_Bottoms_WMA" = "Morton Graner Bottoms WPA", 
  "Pierce_Long_Lake_WPA" = "Pierce Long Lake WPA",
  "Sheridan_Trego_West_1" = "Sherdian Trego West 1",
  "Wells_Pipestone_WPA" = "Wells Pipestone WPA", 
  "WEL-149-069-36" = "WEL") 

pop(BOGR_genind) <- recode(pop(BOGR_genind), !!!rename_populations)

# Confirm the populations are renamed
levels(pop(BOGR_genind))

#Graph
x_limits <- c(-100, 200)  # Change to desired limits for PC1
y_limits <- c(-100, 200)  # Change to desired limits for PC2

pca_data <- data.frame(PC1 = pca_result$li[, 1], 
                       PC2 = pca_result$li[, 2], 
                       Population = as.factor(pop(BOGR_genind)))


color_scheme <- c('#E48312','#566E3D','#8E293B','#BD582C','#865640',
                  '#9B8357','#C2BC80','#94A088','#A53010',
                  '#728653','#92AA4C','#572530','#25351C','#FAD36A', 
                  '#550307','#646742','#FEC44F')
shape_list <- c(15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 15, 16, 17, 18, 19, 20, 21, 22)

PCA1 <- 0.65   
PCA2 <-  0.64

pca_data$Population <- factor(pca_data$Population, levels = sort(levels(factor(pca_data$Population))))

pca <- ggplot(data = pca_data, aes(x = PC1, y = PC2, color = Population, shape = Population)) +
  geom_point(alpha = 0.6, size = 3) +  
  labs(x = "PC1", y = "PC2") +
  scale_color_manual(name = "Population", values = color_scheme) +  
  scale_shape_manual(values = shape_list) +
  xlab(paste("PC1 ", PCA1, "%", sep = ""))+
  ylab(paste("PC2 ", PCA2, "%", sep = ""))+
  theme_minimal(base_size = 13) +  
  coord_cartesian(xlim = x_limits, ylim = y_limits)+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = "transparent"),
        legend.position = "right",
        axis.title = element_text(size = 20, family = "serif"),
        axis.text = element_text(size = 20, family = "serif"),
        legend.title = element_text(size = 20, family = "serif"),
        legend.text = element_text(size = 18, family = "serif"))

pca

#ggsave("./fig/BOGR_PCA.png", pca, width = 10, height = 10, bg = "white", dpi = "print")
#want to save as rdata to use cowplot later 
saveRDS(pca, file = "outputs/BOGR/BOGR_PCA.rds")

#####Admixture####
#packages
require(tidyverse)
require(ggplot2)

# Create the data frame
data <- data.frame(
  K = 1:10,
  CrossValidationScore = c(0.20538,
                           0.22467,
                           0.23836,
                           0.25236,
                           0.26823,
                           0.28723,
                           0.30282,
                           0.31775,
                           0.32987,
                           0.34428))


# Create the plot
admix <- ggplot(data, aes(x = K, y = CrossValidationScore)) +
  geom_line(color = "#566E3D", size = 1.2) +  
  geom_point(color = "#25351C", size = 3) + 
  scale_x_continuous(breaks = 1:10) +
  labs( x = "K",
        y = "Cross-Validation Score") +
  theme_bw()+
  theme(text = element_text(size = 18, family = "serif")) +
  theme(plot.title = element_text(hjust = 0.5, size = 30))



admix
#want to save as rdata to use cowplot later 
saveRDS(admix, file = "outputs/BOGR/BOGR_Admix.rds")
