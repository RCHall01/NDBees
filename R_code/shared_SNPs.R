require(tidyr)
require(ggVennDiagram)

####BOGR####
# 1. List your CSV file paths
BG_Outflank <- read.csv("data/BOGR/adaptsnps/OutFlank_cand.csv")
BG_RDA <- read.csv("data/BOGR/adaptsnps/rda_cand.csv")
BG_pcadapt <- read.csv("data/BOGR/adaptsnps/pcaAdapt.csv")

# 3. For each data frame, create a unique key from CHROM + POSIT
BG_Outflank$key <- paste(BG_Outflank$CHROM, BG_Outflank$POS, sep = "_")
BG_RDA$key <- paste(BG_RDA$CHROM, BG_RDA$POS, sep = "_")
BG_pcadapt$key <- paste(BG_pcadapt$CHROM, BG_pcadapt$POS, sep = "_")

shared_all_three <- Reduce(intersect, list(BG_Outflank$key, BG_RDA$key, BG_pcadapt$key))
ggVennDiagram(x = list(BG_Outflank$key, BG_RDA$key, BG_pcadapt$key))


####BOTE####
# 1. List your CSV file paths
BT_Outflank <- read.csv("data/BOTE/adaptsnps/OutFlank_cand.csv")
BT_RDA <- read.csv("data/BOTE/adaptsnps/rda_cand.csv")
BT_pcadapt <- read.csv("data/BOTE/adaptsnps/pcaAdapt.csv")


# 3. For each data frame, create a unique key from CHROM + POSIT
BT_Outflank$key <- paste(BT_Outflank$CHROM, BT_Outflank$POS, sep = "_")
BT_RDA$key <- paste(BT_RDA$CHROM, BT_RDA$POS, sep = "_")
BT_pcadapt$key <- paste(BT_pcadapt$CHROM, BT_pcadapt$POS, sep = "_")

ggVennDiagram(x = list(BT_Outflank$key, BT_RDA$key, BT_pcadapt$key))

#####Check between species####
RDA_shared <- intersect(BT_RDA$key, BG_RDA$key)
Outflank_shared <- intersect(BT_Outflank$key, BG_Outflank$key)
pcadapt_shared <- intersect(BT_pcadapt$key, BG_pcadapt$key)

shared_all <- Reduce(intersect, list(BT_Outflank$key, BT_RDA$key, BT_pcadapt$key, BG_Outflank$key, BG_RDA$key, BG_pcadapt$key))

ggVennDiagram(x = list(BT_Outflank$key, BT_RDA$key, BT_pcadapt$key, BG_Outflank$key, BG_RDA$key, BG_pcadapt$key))

