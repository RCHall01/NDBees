#packages 
require(GenomicRanges)
require(rtracklayer)
require(clusterProfiler)
require(readr)
require(dplyr)

###bring in data ####
#reference genome 
gff <- import("data_raw/ncbi_dataset/ncbi_dataset/data/GCF_024542735.1/genomic.gff")
genes <- gff[gff$type == "gene"]


#gene ontology/Gaf file
gaf_file <- "data_raw/GCF_024542735.1_iyBomHunt1.1_gene_ontology.gaf/GCF_024542735.1_iyBomHunt1.1_gene_ontology.gaf"
gaf <- read_tsv(gaf_file, comment = "!", col_names = FALSE)

#create term2gene
term2gene <- gaf %>%
  select(gene_id = X2, go_id = X5) %>%
  distinct()

head(term2gene)

#snp list 
RDA <- read.csv("data_raw/BOTE_adapt/BT_rda_cand.csv")
RDA <- RDA[,-c(7:24)]
RDA$ID <- RDA$snp_clean

PCA <- read.csv("data_raw/BOTE_adapt/BT_pcaAdapt.csv")

OF <- read.csv("data_raw/BOTE_adapt/BT_OutFlank_cand.csv")
OF <- OF[,-c(4:6)]

snp_data <- bind_rows(RDA, PCA, OF)

snp_gr <- GRanges(seqnames = snp_data$CHROM, ranges = IRanges(start = snp_data$POS, end = snp_data$POS))

#map SNPS to Genes
overlaps <- findOverlaps(snp_gr, genes)
snp_genes <- genes[subjectHits(overlaps)]
snp_data$gene_id <- NA

# map vector of SNPs to genes
hits <- genes[overlaps@to,]
hits$description

df <- as.data.frame(cbind(hits$description, hits$Name))
df <- df %>% group_by(V1) %>%  tally()

arrange(df, desc(n))

write.csv(df, "data/BOTE_GOresults.csv")


