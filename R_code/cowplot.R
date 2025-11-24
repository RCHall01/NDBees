#this script is just used to cowplot objects together across BOGR and BOTE
require(cowplot)
require(ggplot2)
####PCAs####
BOTE_PCA <- readRDS("outputs/BOTE/BOTE_PCA.rds")
BOTE_PCA
BOGR_PCA <- readRDS("outputs/BOGR/BOGR_PCA.rds")
BOGR_PCA

plot_pca_final <- plot_grid(BOTE_PCA, BOGR_PCA,
                            labels = c("a)", "b)"), ncol = 1)  
plot_pca_final

ggsave("./fig/BOTE_PCA.png", BOTE_PCA, width = 10, height = 6, bg = "white", dpi = "print")
ggsave("./fig/BOGR_PCA.png", BOGR_PCA, width = 10, height = 6, bg = "white", dpi = "print")

ggsave("./fig/PCA.png", plot_pca_final, width = 11, height = 11 , bg = "white", dpi = "print")

####Admix####
BOTE_Admix <- readRDS("outputs/BOTE/BOTE_Admix.rds")
BOGR_Admix <- readRDS("outputs/BOGR/BOGR_Admix.rds")


plot_admix_final <- plot_grid(BOTE_Admix, BOGR_Admix,
                              labels = c("a)", "b)"), ncol = 1)  
plot_admix_final

ggsave("./fig/BOTE_Admix.png", BOTE_Admix, width = 6, height = 6, bg = "white", dpi = "print")
ggsave("./fig/BOGR_Admix.png", BOGR_Admix, width = 6, height = 6, bg = "white", dpi = "print")

ggsave("./fig/Admix.png", plot_admix_final, width = 8, height = 11, bg = "white", dpi = "print")


###Mantel####
BOTE_pairwise <- readRDS("outputs/BOTE/BOTE_Pairwise.rds")
BOTE_Nei <- readRDS("outputs/BOTE/BOTE_Nei.rds")
BOTE_Jost <- readRDS("outputs/BOTE/BOTE_Jost.rds")

BOGR_pairwise <- readRDS("outputs/BOGR/BOGR_Pairwise.rds")
BOGR_Nei <- readRDS("outputs/BOGR/BOGR_Nei.rds")
BOGR_Jost <- readRDS("outputs/BOGR/BOGR_Jost.rds")

plot_mantel_final <- plot_grid(BOTE_pairwise, BOGR_pairwise,
                               BOTE_Nei, BOGR_Nei,
                               BOTE_Jost, BOGR_Jost,
                               labels = c("a)", "b)", "c)", "d)","e)", "f)"), 
                               ncol = 2,
                               nrow = 3)  
plot_mantel_final

ggsave("./fig/mantel.png", plot_mantel_final, width = 15 , height = 17, bg = "white", dpi = "print")

ggsave("./fig/BOTE_pairwise.png", BOTE_pairwise, width = 7, height = 7, bg = "white", dpi = "print")
ggsave("./fig/BOGR_pairwise.png", BOGR_pairwise, width = 7, height = 7, bg = "white", dpi = "print")

ggsave("./fig/BOTE_Nei.png", BOTE_Nei, width = 7, height = 7, bg = "white", dpi = "print")
ggsave("./fig/BOGR_Nei.png", BOGR_Nei, width = 7, height = 7, bg = "white", dpi = "print")

ggsave("./fig/BOTE_Jost.png", BOTE_Jost, width = 7, height = 7, bg = "white", dpi = "print")
ggsave("./fig/BOGR_Jost.png", BOGR_Jost, width = 7, height = 7, bg = "white", dpi = "print")


####Corr####
BOTE_Corr <- readRDS("outputs/BOTE/corr.rds")
BOGR_Corr <- readRDS("outputs/BOGR/corr.rds")

ggsave("./fig/BOTE_Corr.png", BOTE_Corr, width = 12, height = 12, bg = "white", dpi = "print")
ggsave("./fig/BOGR_Corr.png", BOGR_Corr, width = 12, height = 12, bg = "white", dpi = "print")


plot_Corr_final <- plot_grid(BOTE_Corr, BOGR_Corr,
                             labels = c("a)", "b)"), ncol = 1, nrow = 2)  
plot_Corr_final

ggsave("./fig/Corr.png", plot_Corr_final, width = 15, height = 20, bg = "white", dpi = "print")


BOTE_Corr <- readRDS("outputs/BOTE/final_corr.rds")
BOGR_Corr <- readRDS("outputs/BOGR/final_corr.rds")

ggsave("./fig/BOTE_Final_Corr.png", BOTE_Corr, width = 12, height = 12, bg = "white", dpi = "print")
ggsave("./fig/BOGR_Final_Corr.png", BOGR_Corr, width = 12, height = 12, bg = "white", dpi = "print")


###Ne####
BOTE_He <- readRDS("outputs/BOTE/BOTE_HE.rds")
BOTE_LD <- readRDS("outputs/BOTE/BOTE_LD.rds")

BOGR_He <- readRDS("outputs/BOGR/BOGR_HE.rds")
BOGR_LD <- readRDS("outputs/BOGR/BOGR_LD.rds")

plot_Ne <- plot_grid(BOTE_He, BOGR_He,
                               BOTE_LD, BOGR_LD,
                               labels = c("a)", "b)", "c)", "d)"), 
                               ncol = 2,
                               nrow = 2)  
plot_Ne
ggsave("./fig/Ne.png", plot_Ne, width = 15, height = 15, bg = "white", dpi = "print")

