#packages
if (!require("devtools")) install.packages("devtools")
devtools::install_github("thierrygosselin/radiator")
library(radiator)

#use radiator to convert vcf to genepop
genomic_converter(
  data = "data_raw/BOTE/BOTE_Ne.vcf", strata = "data_raw/BOTE/Strata.txt",
  output = c("genepop"))


