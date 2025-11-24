#mantel test -> 12-9-2024
#packages 
require(tidyverse)
require(vegan)
require(geosphere)
require(adegenet)
require(sf)
require(ggplot2)

####bring in data#### 
metadata <- read_csv("data/BOGR/BOGR_POP_META.csv")
Pairwise_Fst <- readRDS("data/BOGR/pairwise_Fst_results.rds")
BOGR_genind <- readRDS("data/BOGR/BOGR_13.rds")

#set seed 
set.seed(5555)

#####Create Geographic Distance####
#convert to UTM 
coords <- data.frame(
  longitude = metadata[,2], 
  latitude = metadata[,3]) #pull coords

points_sf <- st_as_sf(coords, coords = c("X", "Y"), crs = 4326)

points_utm <- st_transform(points_sf, crs = 32614) # Transform to UTM Zone 14N

# Extract UTM coordinates into separate columns
points_utm$easting <- st_coordinates(points_utm)[,1]
points_utm$northing <- st_coordinates(points_utm)[,2]

# Extract UTM coordinates (easting and northing)
coords_utm <- st_coordinates(points_utm)  

# Compute pairwise distances (in meters)
dist_utm_m <- dist(coords_utm)

# Convert to kilometers
Dgeo <- dist_utm_m / 1000

####Mantel w/pairwise####
mantel_result <- mantel(Dgeo, Pairwise_Fst, method = "pearson")
print(mantel_result) #Mantel statistic r: -0.07922   Significance: 0.767  
saveRDS(mantel_result, file = "outputs/BOGR/mantel_Pairwise")
mantel_result <- readRDS("outputs/BOGR/mantel_Pairwise")

#visualize
fst_vec <- as.vector(Pairwise_Fst)
geo_vec <- as.vector(Dgeo)

##create ggplot##
mantel_df <- data.frame(
  geo = geo_vec,
  fst = fst_vec)

Pair_Fst <- ggplot(mantel_df, aes(x = geo, y = fst)) +
  geom_point(color = "#728653", size = 1.5) +
  annotate("text",
           x = max(mantel_df$geo, na.rm = TRUE) * 0.8,
           y = max(mantel_df$fst, na.rm = TRUE) * 0.9,
           label = "P = 0.767",
           color = "black",
           size = 7,
           hjust = 0,
           family = "serif") +
  labs(x = "Geographic Distance (km)",
       y = "Pairwise FST") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.title = element_text(size = 20, family = "serif"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, family = "serif", color = "black"),
    axis.text.y = element_text(angle = 45, hjust = 1, size = 20, family = "serif", color = "black"),
    legend.title = element_text(size = 20, family = "serif"),
    legend.text = element_text(size = 18, family = "serif")
  )

Pair_Fst

#save
saveRDS(Pair_Fst, file = "outputs/BOGR/BOGR_Pairwise.rds")

####mantel w/Jost's
#calculate
#GD.pop.Joost <- mmod::pairwise_D(BOGR_genind, linearized = FALSE)
#saveRDS(GD.pop.Joost, file = "data/BOGR/GD.pop.Joost.rds")
GD.pop.Joost <- readRDS("data/BOGR/GD.pop.Joost.rds")

#mantel_result <- mantel(Dgeo, GD.pop.Joost, method = "pearson")
saveRDS(mantel_result, file = "outputs/BOGR/JostMantel")
print(mantel_result) #Mantel statistic  r: -0.07633 Significance: 0.739
mantel_result <- readRDS("outputs/BOGR/JostMantel")

#visualize
fst_vec <- as.vector(GD.pop.Joost)
geo_vec <- as.vector(Dgeo)

##create ggplot##
mantel_df <- data.frame(
  geo = geo_vec,
  fst = fst_vec)

Josts <- ggplot(mantel_df, aes(x = geo, y = fst)) +
  geom_point(color = "#550307", size = 1.5) +
  annotate("text",
           x = max(mantel_df$geo, na.rm = TRUE) * 0.8,
           y = max(mantel_df$fst, na.rm = TRUE) * 0.9,
           label = "P = 0.739",
           color = "black",
           size = 7,
           hjust = 0,
           family = "serif") +
  labs(x = "Geographic Distance (km)",
       y = "Jost's D") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.title = element_text(size = 20, family = "serif"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, family = "serif", color = "black"),
    axis.text.y = element_text(angle = 45, hjust = 1, size = 20, family = "serif", color = "black"),
    legend.title = element_text(size = 20, family = "serif"),
    legend.text = element_text(size = 18, family = "serif"))

Josts

#save
saveRDS(Josts, file = "outputs/BOGR/BOGR_Jost.rds")

####NeiGst####
#GD.pop.NeiGst <- mmod::pairwise_Gst_Nei(BOGR_genind, linearized = FALSE)
#saveRDS(GD.pop.NeiGst, file = "data/BOGR/GD.pop.NeiGst.rds")
GD.pop.NeiGst <- readRDS("data/BOGR/GD.pop.NeiGst.rds")

#Mantel w/Nei Gst
mantel_result <- mantel(Dgeo, GD.pop.NeiGst, method = "pearson")
print(mantel_result) #Mantel statistic r: -0.07322 Significance: 0.692 
saveRDS(mantel_result, file = "outputs/BOGR/mantel_Nei")
mantel_result <- readRDS("outputs/BOGR/mantel_Nei")

#visualize
fst_vec <- as.vector(GD.pop.NeiGst)
geo_vec <- as.vector(Dgeo)

##create ggplot##
mantel_df <- data.frame(
  geo = geo_vec,
  fst = fst_vec)

Nei <- ggplot(mantel_df, aes(x = geo, y = fst)) +
  geom_point(color = "#E48312", size = 1.5) +
  annotate("text",
           x = max(mantel_df$geo, na.rm = TRUE) * 0.8,
           y = max(mantel_df$fst, na.rm = TRUE) * 0.9,
           label = "P = 0.692",
           color = "black",
           size = 7,
           hjust = 0,
           family = "serif") +
  labs(x = "Geographic Distance (km)",
       y = "Nei Gst") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.title = element_text(size = 20, family = "serif"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, family = "serif", color = "black"),
    axis.text.y = element_text(angle = 45, hjust = 1, size = 20, family = "serif", color = "black"),
    legend.title = element_text(size = 20, family = "serif"),
    legend.text = element_text(size = 18, family = "serif"))
Nei

#save
saveRDS(Nei, file = "outputs/BOGR/BOGR_Nei.rds")



###Gst/Hendrick####
#cannot perform a mantel on but can calulate global value
GD.pop.Hedrick <- mmod::Gst_Hedrick(BOGR_genind)
saveRDS(GD.pop.Hedrick, file = "data/BOGR/GD.pop.Hedrick.rds")
GD.pop.Hedrick <- readRDS("data/BOGR/GD.pop.Hedrick.rds")

#Mantel w/Hendrik
str(GD.pop.Hedrick)
GD.pop.Hedrick.matrix <- do.call(rbind, GD.pop.Hedrick)
GD.pop.Hedrick <- as.dist(GD.pop.Hedrick.matrix)

mantel_result <- mantel(geo_distances, GD.pop.Hedrick, method = "pearson")
print(mantel_result)









require(tidyverse)
require(vegan)
require(geosphere)
require(adegenet)
require(sf)
require(mmod)
BOGR_genind <- readRDS("data/BOGR/BOGR_genind.rds")

#Jost's
BGGD.pop.Joost <- mmod::pairwise_D(BOGR_genind, linearized = FALSE)
saveRDS(BGGD.pop.Joost, file = "data/BOGR/GD.pop.Joost.rds")


BGGD.pop.NeiGst <- mmod::pairwise_Gst_Nei(BOGR_genind, linearized = FALSE)
saveRDS(BGGD.pop.NeiGst, file = "data/BOGR/GD.pop.NeiGst.rds")
