#packages
if (!require("devtools")) install.packages("devtools")
devtools::install_github("thierrygosselin/radiator")
library(radiator)

#use radiator to convert vcf to genepop
genomic_converter(
  data = "data_raw/BOGR/BOGR_Ne.vcf", strata = "data_raw/BOGR/Strata.txt",
  output = c("genepop"))


