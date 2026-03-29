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
df <-read_csv("data/BOTE/final_BOTE.csv")

#import vcf
BOTE <- read.vcfR("data_raw/BOTE/BOTE_obj2.vcf")
BOTE_obj2 <- vcfR2genind(BOTE)

# Step 1: Access current sample names and trim
current_names <- indNames(BOTE_obj2)
trimmed_names <- sub(".*BOTE-(\\d{4}).*", "BOTE.\\1", current_names)
indNames(BOTE_obj2) <- trimmed_names
indNames(BOTE_obj2)

#Step 3: Assign Pops
BOTE_obj2@pop <- as.factor(df$sample_site)

#Step 4: Check your work 
individual_names <- indNames(BOTE_obj2)
population_assignments <- pop(BOTE_obj2) # Extract the population assignments
# Create a data frame to show each individual and their assigned population
individual_population_df <- data.frame(Sample = individual_names, 
                                       Population = as.character(population_assignments))
print(individual_population_df) # View the data frame

#save genind
saveRDS(BOTE_obj2,"data/BOTE/BOTE_obj2.rds")

###Run population structure####
BOTE_genind <- readRDS("data/BOTE/BOTE_obj2.rds")

#run pca
X <- tab(BOTE_genind, NA.method = "mean")
pca_result <- dudi.pca(X, scannf = FALSE, nf = 3)

#pulling out percentage of variance explained 
pca_result$eig
summary(pca_result)

# Create a named vector for renaming populations
rename_populations <- c(
  "EDD-150-065-24" = "F",
  "Fitzsimond_Slough_WMA" = "H",
  "FOS-147-062-26" = "I",
  "GRI-148-061-16" = "J", 
  "GRI-148-061-36" = "K", 
  "Griggs_Evers_WPA" = "L", 
  "Nelson_Engen" = "P", 
  "Pierce_Long_Lake_WPA" = "Q", 
  "ROL-162-072-36" = "R", 
  "Steele_Otto_Spies_WMA" = "T",
  "Stutsman_Northwestern_Lake_WPA" = "U",
  "Towner_Berg_West" = "V", 
  "Towner_Nikolasien_WPA" = "W")

pop(BOTE_genind) <- recode(pop(BOTE_genind), !!!rename_populations)

# Confirm the populations are renamed
levels(pop(BOTE_genind))

#Graph
x_limits <- c(-100, 200)  # Change to desired limits for PC1
y_limits <- c(-100, 200)  # Change to desired limits for PC2

pca_data <- data.frame(PC1 = pca_result$li[, 1], 
                       PC2 = pca_result$li[, 2], 
                       Population = as.factor(pop(BOTE_genind)))


color_scheme <- c('#E48312','#566E3D','#8E293B','#BD582C','#865640',
                  '#9B8357','#C2BC80','#94A088','#92AA4C','#572530',
                  '#25351C','#FAD36A','#FEC44F')

shape_list <- c(15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 15, 16)

PCA1 <-  0.9408
PCA2 <- 0.8753

pca_data$Population <- factor(pca_data$Population, levels = sort(levels(factor(pca_data$Population))))

# Create the scatter plot using ggplot2
pca <- ggplot(data = pca_data, aes(x = PC1, y = PC2, color = Population, shape = Population)) +
  geom_point(alpha = 0.6, size = 3) +  
  labs( x = "PC1", y = "PC2") +
  scale_color_manual(name = "Sample Site", values = color_scheme) +  
  scale_shape_manual(name = "Sample Site", values = shape_list) +
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

#ggsave("./fig/BOTE_PCA.png", pca, width = 10, height = 10, bg = "white", dpi = "print")
#want to save as rdata to use cowplot later 
saveRDS(pca, file = "outputs/BOTE/BOTE_PCA.rds")





####Admixture####
#packages
require(tidyverse)
require(ggplot2)


# Create the data frame
data <- data.frame(
  K = 1:10,
  CrossValidationScore = c(0.24883,
                           0.27035,
                           0.28824,
                           0.30716,
                           0.32474,
                           0.34651,
                           0.36313,
                           0.38206,
                           0.39566,
                           0.41015))

# Create the plot
admix <- ggplot(data, aes(x = K, y = CrossValidationScore)) +
  geom_line(color = "#F99923", linewidth = 1.2) +  
  geom_point(color = "#C1600F", size = 3) + 
  scale_x_continuous(breaks = 1:10) +
  labs( x = "K",
        y = "Cross-Validation Score") +
  theme_bw()+
  theme(text = element_text(size = 18, family = "serif")) +
  theme(plot.title = element_text(hjust = 0.5, size = 30))

admix
#want to save as rdata to use cowplot later 
saveRDS(admix, file = "outputs/BOTE/BOTE_Admix.rds")
