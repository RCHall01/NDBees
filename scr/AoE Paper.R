#load packages 
require(readxl) #read Excel
require(sf) #working with spatial data
require(raster) #Working with raster data 
require(tidyverse) #data manipulation and visualization
require(corrplot) #correlation plots
require(terra) #spatial data manipulation and analysis
require(caret) #feature selection
require(paletteer) #color palettes 

#preemptive rename for correlation
new_column_names <- c(
  "BIO1",  "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", 
  "BIO16", "BIO17", "BIO18", "BIO19", "BIO2", "BIO3", "BIO4", "BIO5", 
  "BIO6", "BIO7", "BIO8", "BIO9")

#preemptive rename for final selection
new_column_names_2 <- c(
 "BIO10 = Mean Temperature of Warmest Quarter",
 "BIO13 = Precipitation of Wettest Month",      
 "BIO15 = Precipitation Seasonality (Coefficient of Variation)",
 "BIO16 = Precipitation of Wettest Quarter",    
 "BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))",
 "BIO7 = Temperature Annual Range (BIO5-BIO6)", 
 "BIO8 = Mean Temperature of Wettest Quarter")

####WorldClim#####
#import data set 
Finals <- read_excel("data_raw/Finals.xlsx")

#load raster and clip
#create list of env data, pulling from a folder wc2.1_2.5_bio and all files ending in .tif were included, and full.names is set to true to pull the path to said file
bio_files <- list.files(path = './data/wc2.1_2.5m_bio', pattern = '*.tif', all.files = TRUE, full.names = TRUE)

#load in the rasters
#create raster stack from the list of files 
bio_layers <- rast(bio_files)
names(bio_layers) <- new_column_names #rename from wc2 names to BIO1-19

#load ND shapefile
#clip the raster layers to the shapefile extent 
nd <- vect('./data/tl_2016_38_place/tl_2016_38_place.shp')
bio_layers <- crop(bio_layers, nd) 

#check correlations across the state
# Extract raster values at the locations specified in Finals
extracted_df <- terra::extract(bio_layers, Finals[,3:4])

#remove ID column
extracted_df <- extracted_df[,-1]

#calculate the correlation among our variables at our points
mydata.cor <- cor(extracted_df, method = 'spearman', use = 'complete.obs')

#make a correlation plot
corrplot(mydata.cor, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, diag=FALSE)

#set up the correlation value cutoff you want to use to make the list of highly correlated variables
hc <- findCorrelation(mydata.cor, cutoff = 0.60)

#sort the list of highlight correlated variables
hc = sort(hc)

# Remove highly correlated variables from the extracted data
predictors_final_list = extracted_df[,-c(hc)]

# Calculate the correlation among our variables at our points after filtering
mydata.cor <- cor(predictors_final_list, method = 'spearman', use = 'complete.obs')

# Create a new correlation plot after filtering
corrplot(mydata.cor, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, diag=FALSE)

# Cut out correlated from the raster stack, then save it
# Subset the raster stack to include only the selected variables
predictors_final <- subset(bio_layers, names(predictors_final_list))
names(predictors_final) <- new_column_names_2

plot(predictors_final)

# Get the names of the rasters in the subsetted raster stack
bio_names <- names(predictors_final)

plot(predictors_final)

# Load shapefile for ND counties
counties <- vect('./data/Counties Boundaries/County_Boundaries.shp')

# Mask raster to counties
WorldClimND <- mask(predictors_final, counties)

#Plot of just annual precip with data points and counties
plot(WorldClimND, "BIO16 = Precipitation of Wettest Quarter" , col= rev(paletteer_c("grDevices::ag_GrnYl", 130)), 
     main = "Precipitation of Wettest Quarter\n(Precipitation (mm))")
plot(counties, add = TRUE)
points(Finals$X, Finals$Y, pch = 21, col = "black", bg = "white")

#Plot of just annual precip with data points and counties
plot(WorldClimND, "BIO10 = Mean Temperature of Warmest Quarter",
     col=rev(paletteer_c("ggthemes::Red-Blue Diverging", 130)), 
     main = "Mean Temperature of Warmest Quarter\n(Temperature (°C))")
plot(counties, add = TRUE)
points(Finals$X, Finals$Y, pch = 21, col = "black", bg = "white")

