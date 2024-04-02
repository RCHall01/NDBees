#packages
require(tidyverse)
require(terra)
require(tidyterra)
require(ggrepel)
require(cowplot)
require(readxl)
require(ggspatial)
require(dplyr)

#import data
RH_Finals <- read_excel("data_raw/Final Samples.xlsx")
Finals <- read_excel("data_raw/Finals.xlsx")
counties <- vect('./data/Counties Boundaries/County_Boundaries.shp')
r <- rast(counties)
r <- raster::raster(r)


#just B. Gris and both plot 
BG <- Finals %>%
  filter(Species != "Bombus ternarius") %>%
  filter(Type != "Bombus ternarius")

BG_Map <- ggplot() +
  tidyterra::geom_spatvector(data = counties, fill = NA, colour = "black", size = 0.75) +
  geom_point(data = BG, aes(x = X, y = Y, color = Species), size = 2.5) +
  geom_label_repel(data = BG, aes(x = X, y = Y, label = n_tubed, color = Type),
                   min.segment.length = 0.05, size = 3, force = 10, max.overlaps = 100) +
  scale_colour_manual(values=c('darkgreen', 'darkred')) +
  theme_void(base_size = 16)
BG_Map 
#just B. tern and both plot
BT <- Finals %>%
  filter(Species != "Bombus griseocollis") %>%
  filter(Type != "Bombus griseocollis")


BT_Map <-ggplot() +
  tidyterra::geom_spatvector(data = counties, fill = NA, colour = "black", size = 0.75) +
  geom_point(data = BT, aes(x = X, y = Y, color = Species), size = 2.5) +
  geom_label_repel(data = BT, aes(x = X, y = Y, label = n_tubed, color = Type),
                   min.segment.length = 0.05, size = 3, force = 10, max.overlaps = 100) +
  scale_colour_manual(values=c('blue', 'darkred')) +
  theme_void(base_size = 16)
BT_Map

#plot points together 
map <- ggplot() +
  tidyterra::geom_spatvector(data = counties, fill = NA, colour = "black", size = 0.75) +
  geom_point(data = Finals, aes(x = X, y = Y, color = Species), size = 2.5) +
  geom_label_repel(data = Finals, aes(x = X, y = Y, label = n_tubed, color = Type),
                   min.segment.length = 0.05, size = 3, force = 10, max.overlaps = 100) +
  scale_colour_manual(values=c('darkgreen', 'blue', 'darkred')) +
  theme_void(base_size = 16)

#create scale/compass 
scale <- ggplot() +
  tidyterra::geom_spatvector(data = counties, fill = NA, colour = "white", size = 0.75) +
  annotation_scale(location = c("br")) +
  annotation_north_arrow(location = "br", which_north = "true", pad_x = unit(0.5, "in"),
                         pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering()) +
  theme_void(base_size = 16)

#map, scale, and compass together to then edit
mapandscale <- plot_grid(map, scale, nrow = 2, align = "v", rel_heights = c(3, 3))

mapandscale


