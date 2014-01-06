# This tool takes in a matrix of feature counts as well as gene annotations and
# outputs a table of top expressions as well as various plots for differential
# expression analysis
#
# ARGS: countPath       -Path to RData input containing counts
#       annoPath        -Path to RData input containing gene annotations
#       htmlPath        -Path to html file linking to other outputs
#       outPath         -Path to folder to write all output to
#       normOpt         -String specifying type of normalisation used
#       weightOpt       -String specifying usage of weights
#       robustOpt       -String specifying usage of robust regression
#       contrastData    -String containing contrasts of interest
#       cpmReq          -Float specifying cpm requirement
#       sampleReq       -Integer specifying cpm requirement
#       factorData      -String containing factor names and levels
#
# OUT:  Voom Plot
#       BCV Plot
#       MA Plot
#       Top Expression Table
#       HTML file linking to the ouputs
#
# Author: Shian Su - registertonysu@gmail.com - Jan 2014


################################################################################
### Input Processing
################################################################################

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
htmlPath <- as.character(argv[3])
outPath <- as.character(argv[4])
rdaOpt <- as.character(argv[5])
normOpt <- as.character(argv[6])
weightOpt <- as.character(argv[7])
robustOpt <- as.character(argv[8])
contrastData <- as.character(argv[9])
cpmReq <- as.numeric(argv[10])
sampleReq <- as.numeric(argv[11])
factorData <- unlist(strsplit(as.character(argv[12]), "::"))

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

# Generate output paths
dir.create(outPath)
pdfOut <- paste0(outPath, "/plots.pdf")
topOut <- paste0(outPath, "/toptab.txt")
rdaOut <- paste0(outPath, "/objectDump.rda")

# Initiate data for html links
linkData <- data.frame(Label=character(), Link=character(),
                       stringsAsFactors=FALSE)

# Read in data
counts <- read.table(countPath)
if (haveAnno){
  geneanno <- read.table(annoPath)
}

################################################################################
### Data Processing
################################################################################

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

# Generating design information
design <- model.matrix(~0+factorLevels)
colnames(design) <- gsub("factorLevels","",colnames(design))
contrasts <- makeContrasts(contrasts=contrastData, levels=design)

################################################################################
### Data Output
################################################################################

# Opening connection to output file
pdf(pdfOut)

linkData[1, ] <- c("Plots (.pdf)", "plots.pdf") 

plotBCV(data, main="BCV Plot")

# Generate voom data and Mean-variance plot
vData <- voom(data, design = design, plot=TRUE)

if (wantWeight){
  # Generating weights
  aw <- arrayWeights(vData, design)
  wts <- asMatrixWeights(aw, dim(vData)) * vData$w
  
  # Generating fit data and top table with weights
  voomFit <- lmFit(vData, design, weights=wts)
  voomFit <- contrasts.fit(voomFit, contrasts)
  voomFit <- eBayes(voomFit)
  top <- topTable(voomFit, coef=1, number=Inf, sort.by="P")
  write.table(top, file=topOut, row.names=FALSE, sep="\t")
  linkData <- rbind(linkData, c("Top Expressions Table with Weights (.txt)",
                                "toptab.txt"))
                    
  # Plot MA (log ratios vs mean average) using limma package on weighted data
  status = decideTests(voomFit[,1])
  limma::plotMA(voomFit, status=status,
                main=paste("MA Plot:", contrastData), 
                col=c("black","firebrick","green"),
                xlab="Average Expression", ylab="logFC")
  abline(h=0, col="grey", lty=2)
  print(summary(status))
} else {
  # Generating fit data and top table without weights
  voomFitNoWeights <- lmFit(vData, design)
  voomFitNoWeights <- contrasts.fit(voomFitNoWeights, contrasts)
  voomFitNoWeights <- eBayes(voomFitNoWeights)
  topNoWeights <- topTable(voomFitNoWeights, coef=1, number=Inf, sort.by="P")
  write.table(topNoWeights, file=topOut, row.names=FALSE, sep="\t")
  linkData <- rbind(linkData, c("Top Expressions Table without Weights (.txt)",
                                "toptab.txt"))
                    
  # Plot MA (log ratios vs mean average) using limma package on unweighted 
  # data
  status = decideTests(voomFitNoWeights[,1])
  limma::plotMA(voomFitNoWeights, status=status,
                main=paste("MA Plot:", contrastData), 
                col=c("black","firebrick","green"),
                xlab="Average Expression", ylab="logFC")
  abline(h=0, col="grey", lty=2)
  print(summary(status))
}

if (wantRda){
  if (wantWeight){
    save(data, vData, labels, factorLevels, wts, voomFit, top,
         file=rdaOut, ascii=TRUE)
  } else {
    save(data, vData, labels, factorLevels, voomFitNoWeights, 
         topNoWeights, 
         file=rdaOut, ascii=TRUE)
  }
  linkData <- rbind(linkData, c("RData (.rda)", "objectDump.rda"))
}

invisible(dev.off())

################################################################################
### HTML Generation
################################################################################

# Clear file
cat("", file=htmlPath)

# Create cat function default path set, default seperator empty and appending
# true by default
cata <- function(..., file = htmlPath, sep = "", fill = FALSE, labels = NULL, 
                 append = TRUE){
  if (is.character(file)) 
    if (file == "") 
      file <- stdout()
  else if (substring(file, 1L, 1L) == "|") {
    file <- pipe(substring(file, 2L), "w")
    on.exit(close(file))
  }
  else {
    file <- file(file, ifelse(append, "a", "w"))
    on.exit(close(file))
  }
  .Internal(cat(list(...), file, sep, fill, labels, append))
}

HtmlHead <- function(title){
  cata("<head>\n")
  cata("<title>", title, "</title>\n")
  cata("</head>\n")
}

HtmlLink <- function(address, label=address){
  cata("<a href=\"", address, "\" target=\"_blank\">", label, "</a><br />\n")
}

cata("<html>\n")
HtmlHead("EdgeR Output")
cata("<body>\n")
cata("<h3>EdgeR Analysis Output:</h3>\n")
for (i in 1:nrow(linkData)){
  HtmlLink(linkData$Link[i], linkData$Label[i])
}
cata("<p>alt-click any of the links to download</p>\n")
cata("<p>.txt files are tab seperated files that can be viewed using Excel ")
cata("or other spreadsheet programs</p>\n")
cata("</body>\n")
cata("</html>")