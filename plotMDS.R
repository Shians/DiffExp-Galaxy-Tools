# This tool takes in an RNAseq feature counts table along with associated
# information to output a MDS plot

# Collects arguments from command line
argv <- commandArgs(TRUE)

################################################################################
### Function Delcaration
################################################################################
# Function to sanitise group information
sanitiseGroups <- function(string) {
  string <- gsub(" *[,] *", ",", string)
  string <- gsub("^\\s+|\\s+$", "", string)
  return(string)
}

# Function to load libaries without messages
silentLibrary <- function(...) {
  list <- c(...)
  for (package in list){
    suppressPackageStartupMessages(library(package, character.only=TRUE))
  }
}

# Load all required libraries
silentLibrary("methods", "statmod", "splines", "edgeR", "limma")

# Grab Arguments
countPath <- as.character(argv[1])
outDir <- as.character(argv[2])
factData <- unlist(strsplit(as.character(argv[3]), ","))
labelData <- unlist(strsplit(as.character(argv[4]), ","))
labelSize <- as.numeric(argv[5])

# Sanitise factor data and labels
factData <- sanitiseGroups(factData)
labelData <- sanitiseGroups(labelData)

# Load in data
counts <- read.table(countPath, header=TRUE, sep="\t")
row.names(counts) <- counts$GeneID
counts <- counts[ , !(colnames(counts)=="GeneID")]

# Extract counts
data <- list()
data$counts <- counts

# Creating naming data
samplenames <- colnames(data$counts)
factor <- factData
factor <- as.factor(factor)
sampleanno <- data.frame("sampleID"=samplenames, "factor"=factor,
                         "group"=factor)

# Attaching the naming data, library size and norming factor
data$samples <- sampleanno
data$samples$lib.size <- colSums(data$counts)
data$samples$norm.factors <- 1
row.names(data$samples) <- colnames(data$counts)
data <- new("DGEList", data)

pdf(outDir)
labels <- labelData
plotMDS(data, labels=labels, col=as.numeric(factor), cex=labelSize)
invisible(dev.off())
