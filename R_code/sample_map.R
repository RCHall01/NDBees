#packages
require(tidyverse)
require(terra)
require(tidyterra)
require(ggrepel)
require(cowplot)
require(ggspatial)
require(dplyr)
require(ggtext)

#import data
Finals <- read_csv("data_raw/finals.csv")
counties <- vect('./data/ED/CB/County_Boundaries.shp')
r <- rast(counties)
r <- raster::raster(r)

#plot points together 
species_labels <- c("*Bombus griseocollis*", "*Bombus ternarius*", "Both")

map2 <- ggplot() +
  geom_spatvector(data = counties, fill = "white", colour = "black", size = 0.75) +
  geom_point(data = Finals, aes(x = X, y = Y, color = Species, shape = Species), size = 2.5) +
  geom_label_repel(data = Finals, aes(x = X, y = Y, label = used, color = Type),
                   min.segment.length = 0.05, size = 4, force = 10, max.overlaps = 100, show.legend = FALSE) +
  scale_color_manual(
    values = c('#566E3D', '#F7951D', '#8E293B'),
    labels = species_labels,
    guide = guide_legend(order = 1, override.aes = list(shape = c(16, 17, 15)))) +
  scale_shape_manual(
    values = c(16, 17, 15),
    labels = species_labels,
    guide = guide_legend(order = 1, override.aes = list(color = c('#566E3D', '#F7951D', '#8E293B')))) +
  annotation_scale() + 
  annotation_north_arrow(location = "tr", which_north = "true") +
  theme_void() +
  theme(
    legend.title = element_text(size = 14, family = "serif"),
    legend.text = element_markdown(size = 12, family = "serif"))

map2


ggsave("./fig/AllAGmap.png",
       
       map2,
       
       bg = "transparent",
       
       dpi = "print",
       
       height = 10,
       
       width = 8.27)
