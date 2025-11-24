# Read first columns (assumed row names) from both Excel files
file1 <- read.csv("data/BOTE_GOresults.csv", header = TRUE)
file2 <- read.csv("data/BOGR_GOresults.csv", header = TRUE)

# Check for shared row names
shared <- intersect(file1[[1]], file2[[2]])

# Print shared row names
print(shared)
