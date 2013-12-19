# Collects arguments from command line
argv <- commandArgs(TRUE)

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

# Process options
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

if (any(options=="mds")){
    mdsOPT=T;
    cat("MDS Plot Requested\n")
} else {
    mdsOPT=F;
    cat("MDS Plot Omitted\n")
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
sel <- rowSums(cpm(data$counts) > 0.5) >= 3
data$counts <- data$counts[sel, ]
data$genes <- data$genes[sel, ]

# Creating naming data
samplenames <- colnames(data$counts)
genotype <- c("WT", "WT", "Homozygous", "Homozygous", "Homozygous", "WT") 
genotype <- as.factor(genotype)
sampleanno <- data.frame("sampleID"=samplenames, "genotype"=genotype,
                         "group"=genotype)

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
design <- model.matrix(~genotype)
colnames(design)[2] <- "WT"

# Generate voom data and Mean-variance plot
vData <- voom(data, design = design, plot=voomOPT)

# Chcking for Smch1d1 expression
ind <- match("Smchd1", vData$genes$Symbols)

###########################################################
### plotMDS (Multidimensional Scaling Plot)
###########################################################
labels <- c("WT-MD1.13", "WT-MD1.14", "Homozygoys-MD1.1", "Homozygous-MD1.2",
            "Homozygous-MD1.3", "WT-MD1.4")
if (mdsOPT){
    plotMDS(vData, labels=labels, col=as.numeric(genotype), cex=0.5, 
            xlim=c(-3.5, 5))
}

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
row.names(top) <- NULL
sel <- top$adj.P.Val<0.05
top <- top[sel, ]
summaryData <- summary(decideTests(voomFit))
rnames <- c("Decreased Expression", "No Change", "Increased Expression")
row.names(summaryData) <- rnames

# Analysis without sample weights
voomFitNoWeights <- lmFit(vData, design)
voomFitNoWeights <- eBayes(voomFitNoWeights)
topNoWeights <- topTable(voomFitNoWeights, coef=2, number=Inf, sort.by="P")
selNoWeights <- topNoWeights$adj.P.Val<0.05
topNoWeights <- top[selNoWeights, ]
summaryDataNoWeights <- summary(decideTests(voomFitNoWeights))

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

    # Make a really weird plot that is worse than the previous two
    sel <- top$adj.P.Val<0.05
    plot(top$AveExpr, top$logFC, xlab="Average Expression", 
         ylab="logFC: WT - Smchd1", pch=20, cex=0.5)
    points(top$AveExpr[sel], top$logFC[sel], pch=18, col=2)
    abline(h=0,col="grey", lty=2)
}
invisible(dev.off()) 