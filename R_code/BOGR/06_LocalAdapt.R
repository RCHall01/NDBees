#packages
require(vegan)    # Used to run PCA & RDA
require(lfmm)     # Used to run LFMM
require(qvalue)   # Used to post-process LFMM output
require(adegenet)
require(raster)   # For handling raster data
require(sf)       # For handling spatial points
require(terra)    # Alternative raster handling (if needed)
require(sp)       # For spatial data handling
require(dplyr)    # For data wrangling
require(tidyr)
require(landscapemetrics)
require(vcfR)
require(tidyverse)
require(ggcorrplot)
require(caret)

####Load in data####
BOGR_genind <- readRDS("data/BOGR/BOGR_obj13.rds")

#### create genotype matrix####
get_mode <- function(x) {
  uniqx <- na.omit(unique(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

# Extract genotype data
BOGR_df <- tab(BOGR_genind, NA.method = NULL)  # Don't replace NA yet

# Apply the mode function to each column (each allele)
geno_data_filled <- apply(BOGR_df, 2, function(x) {
  x[is.na(x)] <- get_mode(x)  # Replace NA with mode
  return(x)
})

# Convert back to a dataframe
geno_data_filled <- as.data.frame(geno_data_filled)

# Ensure genotype data is numeric
geno_matrix <- as.matrix(geno_data_filled) 

#### prepare environmental data####
#load rasters 
bio_files <- stack(list.files(path = './data/ED/wc2.1_2.5m_bio', pattern = '*.tif', full.names = TRUE))
worldclim_raster<- rast(bio_files)

landcover_raster <- rast("./data/NLCD/FedData/ND_NLCD_Land_Cover_2019.tif")

landcover_reproj <- project(landcover_raster, crs(worldclim_raster), method = "near")

crs(landcover_reproj) == crs(worldclim_raster)  # Should return TRUE

#clip rasters to nd
nd <- vect('./data/ED/tl_2016_38_place/tl_2016_38_place.shp')
worldclim_raster <- crop(worldclim_raster, nd)

#reduce worldclim to only the ones we want to keep
variables <- c("wc2.1_2.5m_bio_4", "wc2.1_2.5m_bio_12", "wc2.1_2.5m_bio_19", 
               "wc2.1_2.5m_bio_13", "wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_7", 
               "wc2.1_2.5m_bio_8")

worldclim_raster <- worldclim_raster[[variables]]
#plot(worldclim_raster)

#### Convert to spatial points ####
#bring in points 
finals <-  read.csv('./data/BOGR/final_BOGR.csv')

coordinates(finals) <- ~X + Y  # Adjust column names as needed
proj4string(finals) <- CRS("+proj=longlat +datum=WGS84")

# Define the UTM Zone 14N projection
utm_crs <- CRS("+proj=utm +zone=14 +datum=WGS84 +units=m +no_defs")

# Transform to UTM
finals_utm <- spTransform(finals, utm_crs)

# Extract UTM coordinates
finals_utm$easting <- coordinates(finals_utm)[,1]
finals_utm$northing <- coordinates(finals_utm)[,2]

###create buffers ###
finals_sf <- st_as_sf(finals_utm)
buffers_50km <- st_buffer(finals_sf, dist = 5000)
buffers_50km <- st_transform(buffers_50km, crs(landcover_raster))

lc_buffer <- sample_lsm(
  landcover_raster,
  y = buffers_50km,
  shape = "polygon",
  progress = TRUE,
  what = "lsm_c_pland",  # Proportion of land cover classes
  plot_id = finals_sf$plot_id,  # Ensure plot_id exists
  return_raster = TRUE
)

land_cover <- lc_buffer[, c("plot_id", "class", "value")] %>%
  pivot_wider(names_from = class, values_from = value)

land_cover[is.na(land_cover)] <- 0

#write.csv(land_cover, "data/land_cover_50km.csv", row.names = FALSE)

#plot(landcover_raster)
#plot(buffers_50km, add = TRUE, col = "purple")

#### Extract raster values at sample locations ####
finals_utm_vect <- vect(finals)
# Extract WorldClim raster values
worldclim_data <- terra::extract(worldclim_raster, finals_utm_vect, df = TRUE)

# Extract Landcover raster values
landcover_data <- terra::extract(landcover_raster, finals_utm_vect, df = TRUE)

# Merge data (excluding the first ID column from extract output)
env_data <- cbind(worldclim_data[, -1], land_cover[, -1])

# Check the updated data
head(env_data)

env_data <- env_data %>% 
  rename(
    "Temperature Seasonality"= "wc2.1_2.5m_bio_4",
    "Annual Precipitation" = "wc2.1_2.5m_bio_12",
    "Precip of Cold Quarter" =  "wc2.1_2.5m_bio_19",
    "Precip of Wettest Quarter" = "wc2.1_2.5m_bio_13",
    "Annual Mean Temp" = "wc2.1_2.5m_bio_1",
    "Temp Annual Range" = "wc2.1_2.5m_bio_7",
    "Mean Temp of Wettest Quarter" = "wc2.1_2.5m_bio_8",
    "Open Water" = "11",
    "Developed Open Space" = "21",
    "Developed Low Intensity" = "22",
    "Developed Medium Intensity" = "23",
    "Developed High Intensity" = "24",
    "Barren Land (Rock/Sand/Clay)" = "31",
    "Deciduous Forest" = "41",
    "Evergreen Forest" = "42",
    "Mixed Forest" = "43",
    "Shrub/Scrub" = "52",
    "Grassland/Herbaceous" = "71",
    "Pasture/Hay" = "81",
    "Cultivated Crops" = "82",
    "Woody Wetlands" = "90",
    "Emergent Herbaceous Wetlands" = "95")

####check for correlations####
# Create a correlation plot
cor_matrix <- cor(env_data, use = "pairwise.complete.obs")
corr_BOGR <- ggcorrplot(cor_matrix, method = "circle", type = "lower", lab = TRUE, colors = c("#566E3D", "white", "#FAD36A"))
corr_BOGR
saveRDS(corr_BOGR, "outputs/BOGR/corr.rds")

# Find highly correlated variables
cor_threshold <- 0.8
high_cor_indices <- findCorrelation(cor_matrix, cutoff = cor_threshold)
high_cor_vars <- colnames(cor_matrix)[high_cor_indices]  # Convert indices to names

print(high_cor_vars)

env_data <- env_data[,-c(1, 2, 11, 16)]

# Create a correlation plot
cor_matrix <- cor(env_data, use = "pairwise.complete.obs")
final_corr <- ggcorrplot(cor_matrix, method = "circle", type = "lower", lab = TRUE, colors = c("#566E3D", "white", "#FAD36A"))
saveRDS(final_corr, "outputs/BOGR/final_corr.rds")

# Find highly correlated variables
high_cor_indices <- findCorrelation(cor_matrix, cutoff = cor_threshold)
high_cor_vars <- colnames(cor_matrix)[high_cor_indices]# Convert indices to names
####prepare data ####
# Ensure environmental data is numeric (excluding any non-numeric columns)
env_matrix <- env_data[, sapply(env_data, is.numeric)]  

####run RDA ####
rda_model <- rda(geno_matrix ~ ., data = env_matrix, scale = T)
summary(rda_model)$cont  
rda_model

RsquareAdj(rda_model)
summary(rda_model)$concont

screeplot(rda_model)
#anova.cca(rda_model, by="axis")
vif.cca(rda_model)

plot(rda_model, scaling=3)  ## d.fault is axes 1 and 2

load.rda <- summary(rda_model)$species[,1:3]

hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x) ## f.nd loadings +/- z SD from mean loading     
  x[x < lims[1] | x > lims[2]]           # locus names in these tails
}

cand1 <- outliers(load.rda[,1], 3) ## 1322
cand2 <- outliers(load.rda[,2], 3) ## 358
cand3 <- outliers(load.rda[,3], 3) ## 1270
ncand <- length(cand1) + length(cand2) + length(cand3)
ncand ## 2950

rda.cand <- c(names(cand1), names(cand2), names(cand3)) ## j.st the names of the candidates

rda.cand <- rda.cand[!duplicated(rda.cand)] ## 1.4 unique candidates 

# Set up the color scheme for plotting:
bgcol  <- ifelse(colnames(geno_data_filled) %in% rda.cand, 'gray32', '#00000000')
snpcol <- ifelse(colnames(geno_data_filled) %in% rda.cand, '#FEC44F', '#00000000')

## a.es 1 & 2 - zooming in to just the SNPs here...
plot(rda_model, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), main="RDA: Axes 1 & 2")
points(rda_model, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3)
points(rda_model, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3)
text(rda_model, scaling=3.5, display="bp", col="#25351C", cex=1.2)

## a.es 2 & 3
plot(rda_model, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3), main="RDA, axes 2 and 3")
points(rda_model, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3, choices=c(2,3))
points(rda_model, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3, choices=c(2,3))
text(rda_model, scaling=3, display="bp", col="#25351C", cex=1.2, choices=c(2,3))

intersetcor(rda_model)[,1:3]

#identify SNPS
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=18)  # 18 columns for 18 predictors
colnames(foo) <- c("Precip of Cold Quarter", 
                   "Precip of Wettest Quarter", 
                   "Annual Mean Temp", 
                   "Temp Annual Range", 
                   "Mean Temp of Wettest Quarter",
                   "Open Water", 
                   "Developed Open Space",
                   "Developed Low Intensity", 
                   "Developed High Intensity", 
                   "Barren Land (Rock/Sand/Clay", 
                   "Deciduous Forest", 
                   "Evergreen Forest", 
                   "Shrub/Scrub",
                   "Grassland/Herbaceous", 
                   "Pasture/Hay", 
                   "Cultivated Crops", 
                   "Woody Wetlands", 
                   "Emergent Herbaceous Wetlands")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- geno_matrix[,nam]
  foo[i,] <- apply(env_matrix,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)


rda.cand <- c(names(cand1), names(cand2), names(cand3)) ## just the names of the candidates

length(rda.cand[duplicated(rda.cand)]) 

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) # 
table(foo[foo[,1]==3,2]) # 

cand <- cand[!duplicated(cand$snp),]

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,12] <- names(which.max(abs(bar[4:11]))) # gives the variable
  cand[i,13] <- max(abs(bar[4:11]))              # gives the correlation
}

colnames(cand)[12] <- "predictor"
colnames(cand)[13] <- "correlation"

table(cand$predictor) 
write.csv(table(cand$predictor), "data/BOGR/rda_bypredictor.csv")

##write RDA SNPS to an csv with snp location 
vcf <- read.table("data_raw/BOGR/BOGR_obj13.vcf", comment.char = "#", stringsAsFactors = FALSE)
head(vcf[,1:5])  # Preview: CHROM, POS, ID, REF, ALT
# Add column names for clarity
colnames(vcf)[1:5] <- c("CHROM", "POS", "ID", "REF", "ALT")

# Make a cleaned SNP column without the suffix
cand$snp_clean <- sub("\\.[0-9]+$", "", cand$snp)

# Merge cand with vcf by cleaned SNP ID
cand_merged <- merge(
  cand, 
  vcf[, c("CHROM", "POS", "ID")], 
  by.x = "snp_clean", 
  by.y = "ID",
  all.x = TRUE
)

# Preview
head(cand_merged)

write.csv(cand_merged, "data/BOGR/rda_cand.csv", row.names = TRUE)

####Making beautfied RDA####
# Assign colors to predictors
env_colors <- c(
  "Precip of Cold Quarter" = "#1f78b4",
  "Precip of Wettest Quarter" = "#a6cee3",
  "Mean Temp of Wettest Quarter" = "#6a3d9a",
  "Annual Mean Temp" = "#e31a1c",
  "Temp Annual Range" = "#33a02c",
  "Open Water" = "#ffff33",
  "Developed Open Space" = "#fb9a99",
  "Deciduous Forest" = "#b2df8a"
)

# Apply colors to candidate SNPs based on their top predictor
sel <- cand$snp
env <- cand$predictor
env_color <- env_colors[env]

# Assign color to all SNPs based on whether they are candidates
col.pred <- colnames(geno_data_filled)
col.pred <- ifelse(col.pred %in% sel, env_colors[cand$predictor[match(col.pred, cand$snp)]], '#f1eef6')

# Outline and transparency for non-candidates
empty <- col.pred
empty[empty == "#f1eef6"] <- rgb(0, 1, 0, alpha = 0)
empty.outline <- ifelse(empty == "#00FF0000", "#00FF0000", "gray32")

# Legend colors
bg <- unname(env_colors)
names(bg) <- NULL

# axes 1 & 2
# Axes 1 & 2
plot(rda_model, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), main="RDA, axes 1 and 2")
points(rda_model, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(rda_model, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(rda_model, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=names(env_colors), bty="n", col="gray32", pch=21, cex=0.9, pt.bg=bg)

# Axes 2 & 3
plot(rda_model, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3), main="RDA, axes 2 and 3")
points(rda_model, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(2,3))
points(rda_model, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(2,3))
text(rda_model, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(2,3))
legend("bottomright", legend=names(env_colors), bty="n", col="gray32", pch=21, cex=0.9, pt.bg=bg)

####LFMM####
pred.pca <- rda(env_data, scale=T)
summary(pred.pca)$cont
screeplot(pred.pca, main = "Screeplot: Eigenvalues of bee Predictor Variables")

## correlations between the PC axis and predictors:
round(scores(pred.pca, choices=1:8, display="species", scaling=0), digits=3)
pred.PC1 <- scores(pred.pca, choices=1, display="sites", scaling=0)

#determine K
k <- 1 

#run LFMM
bee.lfmm <- lfmm_ridge(Y=geno_matrix, X=pred.PC1, K=k) ## c.ange K as you see fit

#Identify LFMM candidates 
bee.pv <- lfmm_test(Y=geno_matrix, X=pred.PC1, lfmm=bee.lfmm, calibrate="gif")

names(bee.pv) # this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values

bee.pv$gif

hist(bee.pv$pvalue[,1], main="Unadjusted p-values")        
hist(bee.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

# Let's change the GIF and readjust the p-values:
zscore <- bee.pv$score[,1]   # zscores for first predictor, we only have one in our case...
(gif <- bee.pv$gif[1])       ## d.fault GIF for this predictor

new.gif1 <- 2.0               ## c.oose your new GIF

# Manual adjustment of the p-values:
adj.pv1 <- pchisq(zscore^2/new.gif1, df=1, lower = FALSE)

hist(bee.pv$pvalue[,1], main="Unadjusted p-values")        
hist(bee.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values (GIF=2.8)")

bee.qv <- qvalue(bee.pv$calibrated.pvalue)$qvalues

length(which(bee.qv < 0.1)) ## h.w many SNPs have an FDR < 10%?

(bee.FDR.1 <- colnames(geno_matrix)[which(bee.qv < 0.1)]) ## i.entify which SNPs these are

#### PCA Outlier test ####
rm(list=ls())

require(devtools)
require(pcadapt)
require(qvalue)
require(dplyr)

vcf.path = "data_raw/BOGR/BOGR_obj13.vcf"
meta.path = "data/BOGR/final_BOGR.csv"

meta <- read.csv(meta.path)
genos <- read.pcadapt(vcf.path, type=c("vcf"))

#Use pcadapt to make a PCA
x <- pcadapt(input=genos, K=5)
plot(x, option="screeplot")

#Create a manhattan plot
plot(x, option="manhattan")

#Adjust pvalues to identify statistical outliers
qval <- qvalue(x$pvalues)$qvalues
outliers <- which(qval<0.1)
length(outliers)#5

# Read the VCF, keeping CHROM, POS, ID
snp_info <- read.table(
  vcf.path,
  comment.char = "#",
  stringsAsFactors = FALSE
)[, 1:3]

colnames(snp_info) <- c("CHROM", "POS", "ID")

# Subset outliers by row index
outlier_info <- snp_info[outliers, ]

head(outlier_info)

write.csv(outlier_info, "data/BOGR/pcaAdapt.csv", row.names = FALSE)

##Fst Outlier test####
# Fst Outlier Test with OutFLANK #
require(OutFLANK)
require(vcfR)
require(devtools)
require(qvalue)
require(dplyr)

vcf.path = "data_raw/BOGR/BOGR_obj13.vcf"
meta.path = "data/BOGR/final_BOGR.csv"

meta <- read.csv(meta.path)

data <- read.vcfR(vcf.path)
geno <- extract.gt(data)
dim(geno)

head(geno[,1:10])

#Change formatting of the geno
G <- geno  #we are doing this because we will be running a lot of different 
#things with G, and if we mess up we want to be able to go back to geno
G[geno %in% c("0/0")] <- 0
G[geno  %in% c("0/1")] <- 1
G[geno %in% c("1/1")] <- 2
G[is.na(G)] <- 9
tG <- t(G)
dim(tG)

#Calculate per locus Fst
fst <- MakeDiploidFSTMat(tG, locusNames=colnames(tG), popNames=meta$sample_site)

head(fst)

# Split LocusName into Chrom, Position, Strand
fst <- fst %>%
  tidyr::separate(LocusName, into = c("Chrom", "Position", "Strand"), sep = ":", remove = FALSE) %>%
  mutate(Position = as.numeric(Position))

# Run OutFLANK (same as before)
OF <- OutFLANK(fst, LeftTrimFraction = 0.01, RightTrimFraction = 0.01,
               Hmin = 0.05, NumberOfSamples = 2, qthreshold = 0.01)

OutFLANKResultsPlotter(OF, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.005,
                       Zoom = FALSE, RightZoomFraction = 0.05, titletext = NULL)

# Find outliers
P1 <- pOutlierFinderChiSqNoCorr(fst, Fstbar = OF$FSTNoCorrbar,
                                dfInferred = OF$dfInferred, qthreshold = 0.05, Hmin = 0.1)

outliers <- P1$OutlierFlag == TRUE
table(outliers)

# Plot with Chromosome and Position
plot(P1$Position, P1$FST,
     xlab = "Position", ylab = "FST",
     col = rgb(0, 0, 0, alpha = 0.1), pch = 16)
points(P1$Position[outliers], P1$FST[outliers], pch = 21, col = "#92AA4C", bg = "#92AA4C")

loci <- fst %>%
  filter(LocusName %in% P1$LocusName[outliers]) %>%
  mutate(FST = P1$FST[match(LocusName, P1$LocusName)],
         OutlierFlag = P1$OutlierFlag[match(LocusName, P1$LocusName)])

head(loci)

# Save outliers with Chrom/Position
write.csv(loci, "BOGR_OutFLANK_outliers.csv", row.names = FALSE)


hist(fst$FST, breaks=50)

summary(fst$FST) #highest FST is higher than the mean (which is a good sign)

#Fit the data to a Chi-Square distribution to see if the tails have more SNPs
#than expected
OF <- OutFLANK(fst,LeftTrimFraction=0.01,RightTrimFraction=0.01,
               Hmin=0.05,NumberOfSamples=2,qthreshold=0.01)
OutFLANKResultsPlotter(OF,withOutliers=T,
                       NoCorr=T,Hmin=0.1,binwidth=0.005,
                       Zoom=F,RightZoomFraction=0.05,titletext=NULL)

#Find SNPs that are statistical outliers
P1 <- pOutlierFinderChiSqNoCorr(fst,Fstbar=OF$FSTNoCorrbar,
                                dfInferred=OF$dfInferred,qthreshold=0.05,Hmin=0.1) 
#Changed the Hmin to a lower threshold because this dataset has generally very low heterozygosity
outliers <- P1$OutlierFlag==TRUE #which of the SNPs are outliers?
table(outliers)


plot(P1$Position, P1$FST,
     xlab = "Position", ylab = "FST",
     col = rgb(0, 0, 0, alpha = 0.1), pch = 16)
points(P1$Position[outliers],P1$FST[outliers], pch = 21,col="#E48312", bg = "#E48312")
loci <- P1 %>%
  filter(OutlierFlag==TRUE)
head(loci)

#to create excel file with correct chromosome and posiiton information 
variant_info <- as.data.frame(getFIX(data))[, 1:3]
colnames(variant_info) <- c("CHROM", "POS", "ID")
variant_info$POS <- as.numeric(as.character(variant_info$POS))
missing_id <- is.na(variant_info$ID) | variant_info$ID == "."
if (any(missing_id)) {
  variant_info$ID[missing_id] <- paste0(variant_info$CHROM[missing_id], ":", variant_info$POS[missing_id])
}
if (nrow(P1) != nrow(variant_info)) {
  stop("Row count mismatch: nrow(P1) = ", nrow(P1), ", nrow(variant_info) = ", nrow(variant_info),
       ". Can't safely combine by row order.")
}
drop_cols <- c("CHROM", "POS", "ID", "Chrom", "Position", "Strand")
P1 <- P1[, !(names(P1) %in% drop_cols), drop = FALSE]
P1 <- cbind(P1, variant_info)
outlier_df <- P1 %>%
  filter(OutlierFlag == TRUE) %>%
  dplyr::select(ID, CHROM, POS, FST, qvalues, OutlierFlag)
print(head(outlier_df))

# Save as CSV
write.csv(outlier_df, "data/BOGR/OutFlank_cand.csv", row.names = FALSE)
