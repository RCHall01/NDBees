
#load packages
require(dismo)
require(maps)
require(mapproj)
require(ggthemes)
require(sf)
require(tidyverse)
require(conflicted)
require(terra)
require(corrplot)
require(caret)
require(sp)
require(readxl)
require(landscapemetrics)
require(dplyr)
require(tibble)
require(tmap)
require(RColorBrewer)
require(ggthemes)
require(paletteer)
require(ggplot2)


#####WorldClim rename#####
#preemptive rename for correlation
new_column_names <- c(
  "BIO1",  "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", 
  "BIO16", "BIO17", "BIO18", "BIO19", "BIO2", "BIO3", "BIO4", "BIO5", 
  "BIO6", "BIO7", "BIO8", "BIO9")

#preemptive rename for final selection
new_column_names_2 <- c(
  "Mean Temperature of Warmest Quarter",                                    
  "Precipitation of Wettest Month",                                          
  "Precipitation of Driest Month",                                            
  "Precipitation Seasonality (Coefficient of Variation)",                     
  "Precipitation of Wettest Quarter",                                        
  "Mean Diurnal Range (Mean of monthly (max temp - min temp))",              
  "Temperature Annual Range (BIO5-BIO6)",                                    
  "Mean Temperature of Wettest Quarter")

####WorldClim#####
#import data set 
Finals <- read_excel("data_raw/Finals.xlsx")

sf_finals <- st_as_sf(Finals, coords = c("X", "Y"), crs = 4326)

Sites.sf.longlat <- st_transform(sf_finals, crs = 4326)
head(Sites.sf.longlat)

#map sampling sites 
tmap_mode("view")
tm_shape(sf_finals) + tm_sf(col="red", size=1)

#load raster and clip
#create list of env data
bio_files <- list.files(path = './data/wc2.1_2.5m_bio', pattern = '*.tif', all.files = TRUE, full.names = TRUE)

#load in the rasters
bio_layers <- rast(bio_files)
names(bio_layers) <- new_column_names

#load ND shapefile from TIGER for extent
nd <- vect('./data/tl_2016_38_place/tl_2016_38_place.shp')
bio_layers <- crop(bio_layers, nd)

#check correlations across the state
extracted_df <- terra::extract(bio_layers, Finals[,3:4])

#remove ID column
extracted_df <- extracted_df[,-1]

#calculate the correlation among our variables at our points
mydata.cor <- cor(extracted_df, method = 'spearman', use = 'complete.obs')

#make plot
corrplot(mydata.cor, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, diag=FALSE)

#set up the correlation value cutoff you want to use to make the list of highly correlated variables
hc <- findCorrelation(mydata.cor, cutoff = 0.65)

#sort the list of highlight correlated variables
hc = sort(hc)

#remove the correlated ones from the data frame. This will create what df for
predictors_final_list = extracted_df[,-c(hc)]

#calculate the correlation among our variables at our points
mydata.cor <- cor(predictors_final_list, method = 'spearman', use = 'complete.obs')

#make plot
corrplot(mydata.cor, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, diag=FALSE)


#cut out correlated from the raster stack, then save it
predictors_final <- subset(bio_layers, names(predictors_final_list))
names(predictors_final) <- new_column_names_2

plot(predictors_final)

#names of rasters in the stack
bio_names <- names(predictors_final)

plot(predictors_final)

#vector of ND counties 
counties <- vect('./data/Counties Boundaries/County_Boundaries.shp')
WorldClimND <- mask(predictors_final, counties)

#Plot of just annual precip with data points
par(mar=c(3,3,1,2))
plot(WorldClimND, "Precipitation of Wettest Quarter" , col= rev(paletteer_c("grDevices::ag_GrnYl", 130)), 
     main = "Precipitation of Wettest Quarter")
plot(counties, add = TRUE)

ggplot() +
  tidyterra::geom_spatvector(data = counties, fill = NA, colour = "black", size = 0.75) +
  tidyterra::geom_spatraster(data = WorldClimND, aes(fill = `Precipitation of Wettest Quarter`)) +
  scale_fill_gradientn(colors = rev(paletteer_c("grDevices::ag_GrnYl", 130)), name = "Precipitation (mm)") +
  tidyterra::geom_spatvector(data = counties, fill = NA, colour = "black", size = 0.75)+
  labs(title = "Precipitation of Wettest Quarter") + 
  theme_void(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5), 
    legend.key.size = unit(2, "lines")) +
theme_void(base_size = 10) 



#Plot of just annual precip with data points
par(mar=c(3,3,1,2))
plot(predictors_final, "Temperature Annual Range (BIO5-BIO6)" , col=rev(
  paletteer_c("ggthemes::Red-Blue Diverging", 30)), main = "Temperature Annual Range")

plot(WorldClimND, "Temperature Annual Range (BIO5-BIO6)",
     col=rev(paletteer_c("ggthemes::Red-Blue Diverging", 130)), 
     main = "Range of Annual Temperature")
plot(counties, add = TRUE)

####calculation#####
#extract raster values 
Sites.terra <- terra::extract(predictors_final, Sites.sf.longlat, bind=TRUE)
Sites.sf <- sf::st_as_sf(Sites.terra)
Sites.sf

#calc Landscape Metrics 
landscapemetrics::check_landscape(predictors_final)
landscapemetrics::list_lsm(level = "landscape", type = "diversity metric")
landscapemetrics::list_lsm(level = "landscape", type = "diversity metric")

#Calculate percentage of landscape of class
percentage_class <- lsm_c_pland(landscape = predictors_final$`Temperature Annual Range (BIO5-BIO6)`)

percentage_class




#####LandCover####

nlcd_files <- rast('./data/FedData/NA_NALCMS_landcover_2020_30m.tif')
nlcd_legend <- read.csv('data/FedData/nlcd_legend.csv')
counties <- vect('./data/Counties Boundaries/County_Boundaries.shp')
r <- rast(counties)
r <- raster::raster(r)

#crop 
nlcd_layers <- crop(nd, nlcd_files)
nlcd_layers <- mask(nlcd_files, counties)

#From Travis: edit legend 
nlcd_legend$value <- nlcd_legend$ï..labels
nlcd_legend$ï..labels <- NULL
nlcd_legend2 <- rbind(data.frame(colors = '0000ffff', labels = 'Unclassified'), nlcd_legend)
nlcd_legend2 <- nlcd_legend2[c(-3),]
nlcd_legend3 <- data.frame(nlcd_legend2$labels, colors = c('0000ffff', '#486da2', '#e1cdce',
                                                           '#dc9881', '#f10100', '#ad0101', '#b3afa4', 
                                                           '#6ba966','#1d6533', '#bdcc93', '#d1bb82',
                                                           '#a4cc51', '#ddd83e','#ae7229','#bbd7ed', '#71a4c1'))

class(nlcd_layers)
nlcd_layers
counties


plot(nlcd_files)
nlcd_files


######LandCover Change ######
#load in 
nlcc_plot <- rast('./data/FedData/NA_NALCMS_landchange_2015v3_2020_30m/NA_NALCMS_landchange_2015v3_2020_30m.tif')
nlcc_plot <- crop(nlcc_plot, nd)
plot(nlcc_plot)

