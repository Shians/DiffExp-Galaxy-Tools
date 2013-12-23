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
rdaOut <- as.character(argv[6]);
normOpt <- as.character(argv[7])
weightOpt <- as.character(argv[8])
contrastData <- as.character(argv[9])
cpmReq <- as.numeric(argv[10])
sampleReq <- as.numeric(argv[11])
factorData <- unlist(strsplit(as.character(argv[12]), "::"))

# Process arguments
if (weightOpt=="yes"){
    wantWeight=TRUE;
} else {
    wantWeight=FALSE;
}
factorLevelData <- unlist(strsplit(factorData[3],","))
factorLevels <- factor(factorLevelData)
contrasts <- makeContrasts(contrasts=contrastData, levels=levels(factorLevels))

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
data <- estimateDisp(data, robust=TRUE)
pdf(outDir)

plotBCV(data)

###########################################################
### voom (Mean-variance modelling at the observational level)
###########################################################

# Generating design information
design <- model.matrix(~0+factorLevels)
colnames(design) <- gsub("factorLevels","",colnames(design))

# Generate voom data and Mean-variance plot
vData <- voom(data, design = design, plot=TRUE)

###########################################################
### Toptabs
###########################################################
# Generating weights
aw <- arrayWeights(vData, design)
wts <- asMatrixWeights(aw, dim(vData)) * vData$w

# Analysis with voom and sample weights
voomFit <- lmFit(vData, design, weights=wts)
voomFit <- contrasts.fit(voomFit, contrasts)
voomFit <- eBayes(voomFit)

# Generate table of the 20 top ranked genes from linear fit model
top <- topTable(voomFit, coef=1, number=20, sort.by="P")[, -6]
row.names(top)<-NULL
top$rank <- 1:20
top <- data.frame(c(top[10], top[1:9]))

# Generate a topTable for all gene entries and select for those with p-values
# lower than 0.05
top <- topTable(voomFit, coef=1, number=Inf, sort.by="P")

# Analysis without sample weights
voomFitNoWeights <- lmFit(vData, design)
voomFitNoWeights <- contrasts.fit(voomFitNoWeights, contrasts)
voomFitNoWeights <- eBayes(voomFitNoWeights)
topNoWeights <- topTable(voomFitNoWeights, coef=1, number=Inf, sort.by="P")

# Write out tables for results with and without weights
write.csv(top, file=topOut, row.names=FALSE)
write.csv(topNoWeights, file=topNoWeightOut, row.names=FALSE)

if (wantWeight){
    # Plot MA (log ratios vs mean average) using limma package on weighted data
    status = decideTests(voomFit[,1])
    limma::plotMA(voomFit, status=status,
    main="MA Plot with Weights", col=c("black","firebrick","green"),
    xlab="Average Expression", ylab="logFC")
    abline(h=0, col="grey", lty=2)
    print(summary(status))
} else {
    # Plot MA (log ratios vs mean average) using limma package on unweighted 
    # data
    status = decideTests(voomFitNoWeights[,1])
    limma::plotMA(voomFitNoWeights, status=status,
    main="MA Plot without Weights", col=c("black","firebrick","green"),
    xlab="Average Expression", ylab="logFC")
    abline(h=0, col="grey", lty=2)
    print(summary(status))
}

save(data, vData, labels, factorLevels, wts, voomFit, top, voomFitNoWeights, 
     topNoWeights, file=rdaOut, ascii=TRUE)

invisible(dev.off()) 