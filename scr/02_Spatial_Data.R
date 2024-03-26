#load packages
require(landscapemetrics)
require(dplyr)
require(sf)
require(terra)
require(tibble)
require(tmap)
require(RColorBrewer) 
require(readxl)

#import data set 
Finals <- read_excel("data_raw/Finals.xlsx")
class(Finals)

sf_finals <- st_as_sf(Finals, coords = c("X", "Y"), crs = 32611)

Sites.sf.longlat <- st_transform(sf_finals, crs = 4326)
head(Sites.sf.longlat)

#map sampling sites 
tmap_mode("view")
tm_shape(sf_finals) + tm_sf(col="red", size=1)
