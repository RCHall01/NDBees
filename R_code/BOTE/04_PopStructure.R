##11/27/2024 -> PCA visualization STACKS 
#packages
require(adegenet)
require(tidyverse)
require(ggplot2)
require(FactoMineR)

#pull in data
BOTE_genind <- readRDS("data/BOTE/BOTE_genind13.rds")

#run pca
X <- tab(BOTE_genind, NA.method = "mean")
pca_result <- dudi.pca(X, scannf = FALSE, nf = 3)

#pulling out percentage of variance explained 
pca_result$eig
summary(pca_result)

# Create a named vector for renaming populations
rename_populations <- c(
  "EDD-150-065-24" = "EDD",
  "Fitzsimond_Slough_WMA" = "Fitzsimond Slough WMA",
  "FOS-147-062-26" = "FOS",
  "GRI-148-061-16" = "GRI-16", 
  "GRI-148-061-36" = "GRI-36", 
  "Griggs_Evers_WPA" = "Griggs Evers WPA", 
  "Nelson_Engen" = "Nelson Engen", 
  "Pierce_Long_Lake_WPA" = "Pierce Long Lake WPA", 
  "ROL-162-072-36" = "ROL", 
  "Steele_Otto_Spies_WMA" = "Steele Otto Spies WMA",
  "Stutsman_Northwestern_Lake_WPA" = "Stutsman Northwestern Lake WPA",
  "Towner_Berg_West" = "Towner Berg West", 
  "Towner_Nikolasien_WPA" = "Towner Nikolasien WPA")

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

PCA1 <-  0.9275
PCA2 <- 0.8667

pca_data$Population <- factor(pca_data$Population, levels = sort(levels(factor(pca_data$Population))))

# Create the scatter plot using ggplot2
pca <- ggplot(data = pca_data, aes(x = PC1, y = PC2, color = Population, shape = Population)) +
  geom_point(alpha = 0.6, size = 3) +  
  labs( x = "PC1", y = "PC2") +
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
