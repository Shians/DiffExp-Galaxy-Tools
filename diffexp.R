# This tool takes in a matrix of feature counts as well as gene annotations and
# outputs a table of top expressions as well as various plots for differential
# expression analysis
#
# ARGS: countPath       -Path to RData input containing counts
#       annoPath        -Path to RData input containing gene annotations
#       outDir          -Path to PDF output file for plots
#       topOut          -Path to csv output file for top table with weights
#       topNoWeightOut  -Path to csv output file for top table without weights
#       rdaOut          -Path to RData output for all objects used for plots
#       normOpt         -String specifying type of normalisation used
#       weightOpt       -String specifying usage of weights
#       contrastData    -String containing contrasts of interest
#       cpmReq          -Float specifying cpm requirement
#       sampleReq       -Integer specifying cpm requirement
#       factorData      -String containing factor names and levels
#
# OUT:  Voom Plot
#       BCV Plot
#       MA Plot
#       Top Table

# Collects arguments from command line
argv <- commandArgs(TRUE)

# Load all required libraries
library(edgeR, quietly=TRUE, warn.conflicts=FALSE)
library(limma, quietly=TRUE, warn.conflicts=FALSE)
library(statmod, quietly=TRUE, warn.conflicts=FALSE)
library(splines, quietly=TRUE, warn.conflicts=FALSE)

# Grab arguments
countPath <- as.character(argv[1])
annoPath <- as.character(argv[2])
outDir <- as.character(argv[3])
topOut <- as.character(argv[4])
topNoWeightOut <- as.character(argv[5])
rdaOpt <- as.character(argv[6])
rdaOut <- as.character(argv[7])
normOpt <- as.character(argv[8])
weightOpt <- as.character(argv[9])
robustOpt <- as.character(argv[10])
contrastData <- as.character(argv[11])
cpmReq <- as.numeric(argv[12])
sampleReq <- as.numeric(argv[13])
factorData <- unlist(strsplit(as.character(argv[14]), "::"))

# Process arguments
if (weightOpt=="yes"){
    wantWeight <- TRUE
} else {
    wantWeight <- FALSE
}

if (rdaOpt=="yes"){
    wantRda <- TRUE
} else {
    wantRda <- FALSE
}

if (robustOpt=="yes"){
    wantRobust <- TRUE
} else {
    wantRobust <- FALSE
}

if (annoPath=="None"){
    haveAnno <- FALSE
} else {
    haveAnno <- TRUE
}

factorLevelData <- unlist(strsplit(factorData[3],","))
factorLevels <- factor(factorLevelData)

# Read in data
counts <- read.table(countPath)
if (haveAnno){
    geneanno <- read.table(annoPath)
}

# Extract counts and annotation data
data <- list()
data$counts <- counts
if (haveAnno){
    data$genes <- geneanno
} else {
    data$genes <- data.frame(GeneID=row.names(counts))
}

# Filter out genes that do not have a required cpm in a required number of
# samples
sel <- rowSums(cpm(data$counts) > cpmReq) >= sampleReq
data$counts <- data$counts[sel, ]
data$genes <- data$genes[sel, ]

# Creating naming data
samplenames <- colnames(data$counts)
sampleanno <- data.frame("sampleID"=samplenames, "factorLevels"=factorLevels,
                         "group"=factorLevels)

# Generating the DGEList object "data"
data$samples <- sampleanno
data$samples$lib.size <- colSums(data$counts)
data$samples$norm.factors <- 1
row.names(data$samples) <- colnames(data$counts)
data <- new("DGEList", data)

# Calculating normalising factor, estimating dispersion
data <- calcNormFactors(data, method=normOpt)
data <- estimateDisp(data, robust=wantRobust)

# Opening connection to output file
pdf(outDir)

plotBCV(data, main="BCV Plot")

# Generating design information
design <- model.matrix(~0+factorLevels)
colnames(design) <- gsub("factorLevels","",colnames(design))
contrasts <- makeContrasts(contrasts=contrastData, levels=design)

# Generate voom data and Mean-variance plot
vData <- voom(data, design = design, plot=TRUE)

# Generating weights
aw <- arrayWeights(vData, design)
wts <- asMatrixWeights(aw, dim(vData)) * vData$w

if (wantWeight){
    # Generating fit data and top table with weights
    voomFit <- lmFit(vData, design, weights=wts)
    voomFit <- contrasts.fit(voomFit, contrasts)
    voomFit <- eBayes(voomFit)
    top <- topTable(voomFit, coef=1, number=Inf, sort.by="P")
    write.table(top, file=topOut, row.names=FALSE, sep="\t")
    
    # Plot MA (log ratios vs mean average) using limma package on weighted data
    status = decideTests(voomFit[,1])
    limma::plotMA(voomFit, status=status,
    main=paste("MA Plot:", contrastData), col=c("black","firebrick","green"),
    xlab="Average Expression", ylab="logFC")
    abline(h=0, col="grey", lty=2)
    print(summary(status))
} else {
    # Generating fit data and top table without weights
    voomFitNoWeights <- lmFit(vData, design)
    voomFitNoWeights <- contrasts.fit(voomFitNoWeights, contrasts)
    voomFitNoWeights <- eBayes(voomFitNoWeights)
    topNoWeights <- topTable(voomFitNoWeights, coef=1, number=Inf, sort.by="P")
    write.table(topNoWeights, file=topNoWeightOut, row.names=FALSE, sep="\t")
    
    # Plot MA (log ratios vs mean average) using limma package on unweighted 
    # data
    status = decideTests(voomFitNoWeights[,1])
    limma::plotMA(voomFitNoWeights, status=status,
    main=paste("MA Plot:", contrastData), col=c("black","firebrick","green"),
    xlab="Average Expression", ylab="logFC")
    abline(h=0, col="grey", lty=2)
    print(summary(status))
}

if (wantRda){
    if (wantWeight){
        save(data, vData, labels, factorLevels, wts, voomFit, top,
             file=rdaOut, ascii=TRUE)
    } else {
        save(data, vData, labels, factorLevels, wts, voomFitNoWeights, 
             topNoWeights, 
             file=rdaOut, ascii=TRUE)
    }
}

invisible(dev.off()) 