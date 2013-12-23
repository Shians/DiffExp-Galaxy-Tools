# Collects arguments from command line
argv <- commandArgs(TRUE)

# Load all required libraries
library(edgeR, quietly=T, warn.conflicts=F)
library(limma, quietly=T, warn.conflicts=F)
library(statmod, quietly=T, warn.conflicts=F)
library(splines, quietly=T, warn.conflicts=F)

# Grab Arguments
countPath <- as.character(argv[1])
annoPath <- as.character(argv[2])
outDir <- as.character(argv[3])
factData <- unlist(strsplit(as.character(argv[4]), ","))
labelData <- unlist(strsplit(as.character(argv[5]), ","))
labelSize <- as.numeric(argv[6])

# Load in data
load(countPath)
load(annoPath)

# Extract counts and annotation data
data <- list()
data$counts <- counts$counts
data$genes <- geneanno

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

# Naming the sample data
row.names(data$samples) <- colnames(data$counts)

# Giving data the "DGEList" type
data <- new("DGEList", data)

pdf(outDir)
labels <- labelData
plotMDS(data, labels=labels, col=as.numeric(factor), cex=labelSize)
invisible(dev.off())
