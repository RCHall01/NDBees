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

#### start here on second run ####
BOTE_genind <- readRDS("data/BOTE/BOTE_genind13.rds")

#### create genotype matrix####
get_mode <- function(x) {
  uniqx <- na.omit(unique(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

# Extract genotype data
BOTE_df <- tab(BOTE_genind, NA.method = NULL)  # Don't replace NA yet

# Apply the mode function to each column (each allele)
geno_data_filled <- apply(BOTE_df, 2, function(x) {
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
finals <-  read.csv('./data/BOTE/final_BOTE.csv')

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
corr_BOTE <- ggcorrplot(cor_matrix, method = "circle", type = "lower", lab = TRUE, colors =c("#BD582C", "white", "#FAD36A"))
corr_BOTE
saveRDS(corr_BOTE, "outputs/BOTE/corr.rds")

# Find highly correlated variables
cor_threshold <- 0.8
high_cor_indices <- findCorrelation(cor_matrix, cutoff = cor_threshold)
high_cor_vars <- colnames(cor_matrix)[high_cor_indices]  # Convert indices to names

print(high_cor_vars)

env_data <- env_data[,-c(11, 16, 17)]

# Create a correlation plot
cor_matrix <- cor(env_data, use = "pairwise.complete.obs")
final_corr <- ggcorrplot(cor_matrix, method = "circle", type = "lower", lab = TRUE, colors =c("#BD582C", "white", "#FAD36A"))
saveRDS(final_corr, "outputs/BOTE/final_corr.rds")

# Find highly correlated variables
high_cor_indices <- findCorrelation(cor_matrix, cutoff = cor_threshold)
high_cor_vars <- colnames(cor_matrix)[high_cor_indices]  # Convert indices to names

print(high_cor_vars)
print(colnames(env_data))
####prepare data ####
# Ensure environmental data is numeric (excluding any non-numeric columns)
env_matrix <- env_data[, sapply(env_data, is.numeric)]  

####run RDA ####
rda_model <- rda(geno_matrix ~ ., data = env_matrix, scale = T)
summary(rda_model)$cont  
rda_model
saveRDS(rda_model, file = "outputs/BOTE/BOTE_RDA.rds")
rda_model <- readRDS("outputs/BOTE/BOTE_RDA.rds")

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

cand1 <- outliers(load.rda[,1], 3) ## 274
cand2 <- outliers(load.rda[,2], 3) ## 218
cand3 <- outliers(load.rda[,3], 3) ## 208
ncand <- length(cand1) + length(cand2) + length(cand3)
ncand ## 4466

rda.cand <- c(names(cand1), names(cand2), names(cand3)) ## j.st the names of the candidates

length(rda.cand[duplicated(rda.cand)]) 

rda.cand <- rda.cand[!duplicated(rda.cand)] ## 700 unique candidates 

# Set up the color scheme for plotting:
bgcol  <- ifelse(colnames(geno_data_filled) %in% rda.cand, 'gray32', '#00000000')
snpcol <- ifelse(colnames(geno_data_filled) %in% rda.cand, '#F9A602', '#00000000')

## a.es 1 & 2 - zooming in to just the SNPs here...
plot(rda_model, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), main="RDA: Axes 1 & 2", cex = 2)
points(rda_model, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3)
points(rda_model, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3)
text(rda_model, scaling=2, display="bp", col="#25351C", cex=1.2)

## a.es 2 & 3
plot(rda_model, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3), main="RDA: Axes 2 & 3", cex = 2)
points(rda_model, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3, choices=c(2,3))
points(rda_model, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3, choices=c(2,3))
text(rda_model, scaling=3, display="bp", col="#25351C", cex=1.2, choices=c(2,3))

intersetcor(rda_model)[,1:3]

##### Combine all candidate SNPs
cand1 <- cbind.data.frame(rep(1, length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2, length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3, length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis", "snp", "loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# Create correlation matrix (foo) of SNP x predictors
ncand <- nrow(cand)  # <- You were using `ncand` before using it; define it here
foo <- matrix(nrow = ncand, ncol = ncol(env_matrix))
colnames(foo) <- colnames(env_matrix)

for (i in 1:ncand) {
  snp_name <- cand$snp[i]
  snp_values <- geno_matrix[, snp_name]
  foo[i, ] <- apply(env_matrix, 2, function(x) cor(x, snp_values))
}

# Combine with main dataframe
cand <- cbind.data.frame(cand, foo)

# Remove duplicates (keep only one SNP per axis)
cand <- cand[!duplicated(cand$snp), ]

# Assign predictor and correlation based on max |correlation|
cor_cols <- (ncol(cand) - ncol(env_matrix) + 1):ncol(cand)  # column range for predictors
cand$predictor <- apply(cand[, cor_cols], 1, function(row) names(which.max(abs(row))))
cand$correlation <- apply(cand[, cor_cols], 1, function(row) max(abs(row)))

# Preview
head(cand[, c("snp", "axis", "loading", "predictor", "correlation")])

table(cand$predictor)
write.csv(table(cand$predictor), "data/BOTE/rda_bypredictor.csv")

##write RDA SNPS to an csv with snp location 
vcf <- read.table("data_raw/BOTE/BOTE_obj13.vcf", comment.char = "#", stringsAsFactors = FALSE)
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
  all.x = TRUE)

# Preview
head(cand_merged)

write.csv(cand_merged, "data/BOTE/rda_cand.csv", row.names = TRUE)
          
#####LFMM####

pred.pca <- rda(env_data, scale=T)
summary(pred.pca)$cont
screeplot(pred.pca, main = "Screeplot: Eigenvalues of bee Predictor Variables")

## c.rrelations between the PC axis and predictors:
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

(bee.FDR.1 <- colnames(geno_matrix)[which(bee.qv < 0.1)])## 92


# Combine SNP IDs with environmental variable name
# Assuming you ran LFMM with multiple predictors, like pred.PC1 and pred.PC2
bee.qv.multi <- apply(bee.pv$calibrated.pvalue, 2, function(p) qvalue(p)$qvalues)

# Find SNPs with q-value < 0.1 for any predictor
significant_snps <- apply(bee.qv.multi, 1, function(row) {
  predictors <- colnames(bee.qv.multi)[which(row < 0.1)]
  if (length(predictors) > 0) paste(predictors, collapse = ",") else NA
})

lfmm_results <- data.frame(
  SNP = colnames(geno_matrix),
  qvalue_PC1 = bee.qv.multi[,1],
  associated_with = significant_snps)

lfmm_candidates <- lfmm_results[!is.na(lfmm_results$associated_with), ]

# View results
print(lfmm_candidates)

env_pca <- prcomp(env_matrix, scale. = TRUE)
summary(env_pca)  # To see variance explained
loadings <- env_pca$rotation[,1]  # Loadings for PC1
loadings_sorted <- sort(loadings, decreasing = TRUE)
print(loadings_sorted)
write.csv(loadings_sorted, "data/BOTE/LFMM_loadings.csv", row.names = TRUE)


# Read the VCF and label columns
vcf <- read.table("data_raw/BOTE/BOTE_adapt.vcf", comment.char = "#", stringsAsFactors = FALSE)
colnames(vcf)[1:5] <- c("CHROM", "POS", "ID", "REF", "ALT")

# Create a data frame of significant SNPs (FDR < 0.1)
bee_sig_snps <- data.frame(snp = bee.FDR.1)

bee_sig_snps$ID <- sub("\\.(0|1)$", "", sub("_", ".", bee_sig_snps$snp))


# If your VCF has "." in the ID column, construct IDs manually from CHROM and POS
vcf$ID[vcf$ID == "."] <- paste(vcf$CHROM[vcf$ID == "."], vcf$POS[vcf$ID == "."], sep = "_")

# Merge the SNPs with their locations
snp_locations <- merge(bee_sig_snps, vcf["ID"], by = "ID", all.x = TRUE)

# Preview result
head(snp_locations)

# Save as CSV
write.csv(snp_locations, "data/BOTE/LFMM_cand.csv", row.names = FALSE)




#### PCA Outlier test ####
rm(list=ls())

require(devtools)
require(pcadapt)
require(qvalue)
require(dplyr)

vcf.path = "data_raw/BOTE/BOTE_obj13.vcf"
meta.path = "data/BOTE/final_BOTE.csv"

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
length(outliers)#1231

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

write.csv(outlier_info, "data/BOTE/pcaAdapt.csv", row.names = FALSE)

#### FST Outlier test ####
# Fst Outlier Test with OutFLANK 
require(devtools)
require(qvalue)
require(dplyr)
require(OutFLANK)
require(vcfR)

vcf.path = "data_raw/BOTE/BOTE_obj13.vcf"
meta.path = "data/BOTE/final_BOTE.csv"

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
points(P1$Position[outliers], P1$FST[outliers], pch = 21, col = "#E48312", bg = "#E48312")

loci <- fst %>%
  filter(LocusName %in% P1$LocusName[outliers]) %>%
  mutate(FST = P1$FST[match(LocusName, P1$LocusName)],
         OutlierFlag = P1$OutlierFlag[match(LocusName, P1$LocusName)])

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
write.csv(outlier_df, "data/BOTE/OutFlank_cand.csv", row.names = FALSE)
