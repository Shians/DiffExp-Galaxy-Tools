# Collects arguments from command line
argv <- commandArgs(TRUE)

# Load all required libraries
library(edgeR, quietly=T, warn.conflicts=F)
library(limma, quietly=T, warn.conflicts=F)
library(statmod, quietly=T, warn.conflicts=F)
library(splines, quietly=T, warn.conflicts=F)

# Grab arguments
countPath <- as.character(argv[1])
annoPath <- as.character(argv[2])
outDir <- as.character(argv[3])
topOut <- as.character(argv[4])
topNoWeightOut <- as.character(argv[5])
options <- unlist(strsplit(as.character(argv[6]), ","))
normOpt <- as.character(argv[7])
weightOpt <- as.character(argv[8])
cpmReq <- as.numeric(argv[9])
sampleReq <- as.numeric(argv[10])
factor <- as.character(argv[11])

# Process arguments
if (any(options=="bcv")){
    bcvOPT=T;
    cat("BCV Plot Requested\n")
} else {
    bcvOPT=F;
    cat("BCV Plot Omitted\n")

}

if (any(options=="voom")){
    voomOPT=T;
    cat("voom Plot Requested\n")
} else {
    voomOPT=F;
    cat("voom Plot Omitted\n")
}

if (any(options=="ma")){
    maOPT=T;
    cat("MA Plot Requested\n")
} else {
    maOPT=F;
    cat("MA Plot Omitted\n")
}

if (weightOpt=="yes"){
    wantWeight=T;
    cat("Weights used\n")
} else {
    wantWeight=F;
    cat("Weights not used\n")
}

factor <- unlist(strsplit(as.character(argv[11]),"::"))
factorLevelData <- unlist(strsplit(factor[3],","))

###########################################################
### BCV (Biological Coefficient of Variation) Analysis
###########################################################

# Load in data
load(countPath)
load(annoPath)

# Extract counts and annotation data
data <- list()
data$counts <- counts$counts
data$genes <- geneanno

# Select only the genes that have more than 0.5 cpm in at least 3 samples
sel <- rowSums(cpm(data$counts) > cpmReq) >= sampleReq
data$counts <- data$counts[sel, ]
data$genes <- data$genes[sel, ]

# Creating naming data
samplenames <- colnames(data$counts)
factorLevels <- factorLevelData
factorLevels <- as.factor(factorLevels)
sampleanno <- data.frame("sampleID"=samplenames, "factorLevels"=factorLevels,
                         "group"=factorLevels)

# Attaching the naming data, library size and norming factor
data$samples <- sampleanno
data$samples$lib.size <- colSums(data$counts)
data$samples$norm.factors <- 1

# Naming the sample data
row.names(data$samples) <- colnames(data$counts)

# Giving data the "DGEList" type
data <- new("DGEList", data)

# Calculating normalising factor, estimating dispersion and plotting BCV
data <- calcNormFactors(data, method=normOpt)
data <- estimateDisp(data, robust=T)
pdf(outDir)

if (bcvOPT){
    plotBCV(data)
}

###########################################################
### voom (Mean-variance modelling at the observational level)
###########################################################

# Generating design information
design <- model.matrix(~factorLevels)
colnames(design)[2] <- "WT"

# Generate voom data and Mean-variance plot
vData <- voom(data, design = design, plot=voomOPT)

###########################################################
### Toptabs
###########################################################
# Generating weights
aw <- arrayWeights(vData, design)
wts <- asMatrixWeights(aw, dim(vData)) * vData$w

# Analysis with voom and sample weights
voomFit <- lmFit(vData, design, weights=wts)
voomFit <- eBayes(voomFit)

# Generate table of the 20 top ranked genes from linear fit model
top <- topTable(voomFit, coef=2, number=20, sort.by="P")[, -6]
row.names(top)<-NULL
top$rank <- 1:20
top <- data.frame(c(top[10], top[1:9]))

# Generate a topTable for all gene entries and select for those with p-values
# lower than 0.05
top <- topTable(voomFit, coef=2, number=Inf, sort.by="P")

# Analysis without sample weights
voomFitNoWeights <- lmFit(vData, design)
voomFitNoWeights <- eBayes(voomFitNoWeights)
topNoWeights <- topTable(voomFitNoWeights, coef=2, number=Inf, sort.by="P")

# Write out tables for results with and without weights
write.csv(top, file=topOut, row.names = F)
write.csv(topNoWeights, file=topNoWeightOut, row.names = F)

if (maOPT){
    # Plot MA (log ratios vs mean average) using limma package on weighted data
    limma::plotMA(voomFit, array=2, status=decideTests(voomFit[,2]),
                  main="WT - Smchd1", xlab="Average Expression",
                  ylab="logFC: WT - Smchd1")
    abline(h=0, col="grey", lty=2)

    # Plot MA (log ratios vs mean average) using limma package on unweighted 
    # data
    limma::plotMA(voomFitNoWeights, array=2, 
                  status=decideTests(voomFitNoWeights[,2]),
                  main="WT - Smchd1", xlab="Average Expression",
                  ylab="logFC: WT - Smchd1")
    abline(h=0, col="grey", lty=2)
}
invisible(dev.off()) 