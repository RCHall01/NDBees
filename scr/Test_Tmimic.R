#for bringing in points, starting with just CK's ND data, but could expand to include gbif at a later date.
#test
library(rnaturalearthdata)
library(rnaturalearth)
library(dismo)
library(maps)
library(mapproj)
library(ggthemes)
library(sf)
library(tidyverse)
library(conflicted)
library(terra)
library(tidyterra)
library(ggrepel)
#note that older version of feddata has missing tile issue. dev says to pull dev version
#devtools::install_github("ropensci/FedData")
library(FedData)
library(ggspatial)

#bring in county map
counties <- vect('./data/Counties Boundaries/County_Boundaries.shp')

#get landcover
#template raster
r <- rast(counties)
r <- raster::raster(r)

#load data
nlcd <- get_nlcd(
  r,
  label = 'ND',
  year = 2019,
  dataset = c("landcover"),
  extraction.dir = "./data/FedData",
  #raster.options = c("COMPRESS=DEFLATE", "ZLEVEL=9"),
  force.redo = F)

#can also load with the main function
nlcd_plot <- rast('data/FedData/NA_NALCMS_landcover_2020_30m.tif')
nlcd_plot <- crop(nlcd_plot, nd)

#csv of the colors to use
nlcd_legend <- read.csv('./data/FedData/nlcd_legend.csv')

#quick rename
nlcd_legend$value <- nlcd_legend$ï..labels
nlcd_legend$ï..labels <- NULL

#need to add unclassified
nlcd_legend2 <- rbind(data.frame(colors = '0000ffff', labels = 'Unclassified'), nlcd_legend)

#remove pren snow and ice
nlcd_legend2 <- nlcd_legend2[c(-3),]

#fix projection issue
counties_crs <- crs(counties)
crs(nlcd_plot) <- counties_crs
nlcd_plot2 <- terra::project(nlcd_plot, counties_crs)

#although labels are correct, the colors are not quite right from that csv I found online, maybe an older version of nlcd
#pull out from second position perrenial snow/ice
#also change herbaceous to a slightly more green color original was '#edeccd'
nlcd_legend3 <- data.frame(nlcd_legend2$labels, colors = c('0000ffff', '#486da2', '#e1cdce', '#dc9881', '#f10100', '#ad0101', '#b3afa4', '#6ba966','#1d6533', '#bdcc93', '#d1bb82','#a4cc51', '#ddd83e','#ae7229','#bbd7ed', '#71a4c1'))

#quick rename
nlcd_legend3$value <- nlcd_legend3$labels


#crop the nlcd to better match the counties file
nlcd_plot3 <- mask(nlcd_plot2, counties)

#this version doesn't have 'unclassified'
nlcd_legend3 <- nlcd_legend3[c(-1),]
nlcd_plot3$NA_NALCMS_landcover_2020_30m

nlcd_plot3$category <- factor(nlcd_plot3$category)
#plot by species
ggplot() +
  geom_spatraster(data = nlcd_plot3) + scale_fill_manual(values = nlcd_legend3$colors, na.value = NA)+ 
  tidyterra::geom_spatvector(data = counties, fill = NA, colour = "black", size = 0.75)+
  theme_void(base_size = 16)



