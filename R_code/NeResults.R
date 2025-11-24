rm(list=ls())
#packages
library(dplyr)
library(ggplot2)

ne <- read.csv("data/effpopsize.csv")

ne_fixed <- ne %>%
  mutate(
    plot_ceiling = case_when(
      species == "Bombus griseocollis" & method == "Linkage Disequilibrium" ~ 233800,
      species == "Bombus griseocollis" & method == "Heterozygote Excess" ~ 2000, 
      species == "Bombus ternarius" & method == "Heterozygote Excess" ~ 70,
      species == "Bombus ternarius" & method == "Linkage Disequilibrium" ~ 22000,
      TRUE ~ max(hci[is.finite(hci)], na.rm = TRUE) * 1.1),
    hci_plot = ifelse(is.infinite(hci) | is.na(hci), plot_ceiling, hci),
    inf_flag = is.infinite(hci) | is.na(hci) | is.infinite(estimate) | is.na(estimate))


##plot 
species_colors <- c("Bombus griseocollis" = "#566E3D", 
                    "Bombus ternarius" = "#E48312")

# Plot a: LD - Bombus griseocollis
p_d <- ggplot(ne_fixed %>% filter(method=="Linkage Disequilibrium", species=="Bombus griseocollis"),
              aes(x=cutoff, y=estimate)) +
  geom_point(size=5, color="#566E3D") +
  geom_errorbar(aes(ymin=lci, ymax=hci_plot), color="#566E3D", linewidth = 1.5) +
  geom_text(aes(label = round(estimate, 1)),family = "serif", hjust = -.25, color="#566E3D", size = 5) +
  geom_text(data = subset(ne_fixed %>% filter(method=="Linkage Disequilibrium", species=="Bombus griseocollis"), inf_flag),
            aes(label="∞", y=hci_plot), vjust=-0.15, color="#566E3D", family = "serif", size = 10) +
  xlab("Minor Allele Frequency Cutoff") +
  ylab("Effective Population Size") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(size = 20, family = "serif"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, family = "serif", color = "black"),
    axis.text.y = element_text(angle = 45, hjust = 1, size = 20, family = "serif", color = "black"),
    plot.title = element_text(size=22, family="serif"))

p_d
saveRDS(p_d, file = "outputs/BOGR/BOGR_LD.rds")

p_b <- ggplot(ne_fixed %>% filter(method=="Heterozygote Excess", species=="Bombus griseocollis"),
              aes(x=cutoff, y=estimate)) +
  geom_point(size=5, color="#566E3D") +
  geom_errorbar(aes(ymin=lci, ymax=hci_plot), color="#566E3D", linewidth = 1.5) +
  #geom_text(aes(label = round(estimate, 1)),family = "serif", hjust = -.25, color="#566E3D", size = 5) +
  geom_text(data = subset(ne_fixed %>% filter(method=="Heterozygote Excess", species=="Bombus griseocollis"), inf_flag),
            aes(label="∞", y=hci_plot), vjust=-0.15, color="#566E3D", family = "serif", size = 10) +
  scale_y_log10()+
  xlab("Minor Allele Frequency Cutoff") +
  ylab("Effective Breeding Population Size") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(size = 20, family = "serif"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, family = "serif", color = "black"),
    axis.text.y = element_text(angle = 45, hjust = 1, size = 20, family = "serif", color = "black"),
    plot.title = element_text(size=22, family="serif"))

p_b
saveRDS(p_b, file = "outputs/BOGR/BOGR_HE.rds")

p_c <- ggplot(ne_fixed %>% filter(method=="Linkage Disequilibrium", species=="Bombus ternarius"),
              aes(x=cutoff, y=estimate)) +
  geom_point(size=5, color="#E48312") +
  geom_errorbar(aes(ymin=lci, ymax=hci_plot), color="#E48312", linewidth = 1.5) +
  geom_text(aes(label = round(estimate, 1)),family = "serif", hjust = -.25, color="#E48312", size = 5) +
  geom_text(data = subset(ne_fixed %>% filter(method=="Linkage Disequilibrium", species=="Bombus ternarius"), inf_flag),
            aes(label="∞", y=hci_plot), vjust=-0.15, color="#E48312", family = "serif", size = 10) +
  xlab("Minor Allele Frequency Cutoff") +
  ylab("Effective Population Size") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(size = 20, family = "serif"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, family = "serif", color = "black"),
    axis.text.y = element_text(angle = 45, hjust = 1, size = 20, family = "serif", color = "black"),
    plot.title = element_text(size=22, family="serif"))

p_c
saveRDS(p_c, file = "outputs/BOTE/BOTE_LD.rds")

p_a <- ggplot(ne_fixed %>% filter(method=="Heterozygote Excess", species=="Bombus ternarius"),
              aes(x=cutoff, y=estimate)) +
  geom_point(size=5, color="#E48312") +
  geom_errorbar(aes(ymin=lci, ymax=hci_plot), color="#E48312", linewidth = 1.5) +
  geom_text(aes(label = round(estimate, 1)),family = "serif", hjust = -.25, color="#E48312", size = 5) +
  geom_text(data = subset(ne_fixed %>% filter(method=="Heterozygote Excess", species=="Bombus ternarius"), inf_flag),
            aes(label="∞", y=hci_plot), vjust=-0.15, color="#E48312", family = "serif", size = 10) +
  xlab("Minor Allele Frequency Cutoff") +
  ylab("Effective Breeding Population Size") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(size = 20, family = "serif"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, family = "serif", color = "black"),
    axis.text.y = element_text(angle = 45, hjust = 1, size = 20, family = "serif", color = "black"),
    plot.title = element_text(size=22, family="serif"))

p_a
saveRDS(p_a, file = "outputs/BOTE/BOTE_HE.rds")
