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
contrastData <- as.character(argv[8])
cpmReq <- as.numeric(argv[9])
sampleReq <- as.numeric(argv[10])
factorData <- list()
for (i in 11:length(argv)){
  newFact <- unlist(strsplit(as.character(argv[i]), split="::"))
  factorData <- rbind(factorData, newFact)
}

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

if (annoPath=="None"){
  haveAnno <- FALSE
} else {
  haveAnno <- TRUE
}

# Set the row names to be the name of the factor and delete first row
row.names(factorData) <- factorData[, 1]
factorData <- factorData[, -1]
factorData <- sapply(factorData, strsplit, split=",")

# Transform factor data into data frame of R factor objects
factors <- data.frame(factorData)

# Split up contrasts seperated by comma into a vector
contrastData <- unlist(strsplit(contrastData, split=","))

# Generate output paths
dir.create(outPath)
bcvOut <- paste0(outPath, "/bcvplot.png")
voomOut <- paste0(outPath, "/voomplot.png")
maOutPdf <- character()   # Initialise character vector
maOutPng <- character()
topOut <- character()
for (i in 1:length(contrastData)){
  maOutPdf[i] <- paste0(outPath, "/maplot(", contrastData[i], ").pdf")
  maOutPng[i] <- paste0(outPath, "/maplot(", contrastData[i], ").png")
  topOut[i] <- paste0(outPath, "/toptab(", contrastData[i], ").tsv")
}
rdaOut <- paste0(outPath, "/objectDump.rda")

# Initiate data for html links
linkData <- data.frame(Label=character(), Link=character(),
                       stringsAsFactors=FALSE)
imageData <- data.frame(Label=character(), Link=character(),
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
sampleanno <- data.frame("sampleID"=samplenames, "group"=factors[[1]])

# Generating the DGEList object "data"
data$samples <- sampleanno
data$samples$lib.size <- colSums(data$counts)
data$samples$norm.factors <- 1
row.names(data$samples) <- colnames(data$counts)
data <- new("DGEList", data)

# Calculating normalising factor, estimating dispersion
data <- calcNormFactors(data, method=normOpt)
data <- estimateDisp(data, robust=TRUE)

# Generating design information
pasteListName <- function(string){
  return(paste0("factors$",string))
}
factorList <- sapply(names(factors), pasteListName)
formula <- "~0"
for (i in 1:length(factorList)){
  formula <- paste(formula, factorList[i], sep="+")
}
formula <- formula(formula)
design <- model.matrix(formula)
for (i in 1:length(factorList)){
  colnames(design) <- gsub(factorList[i], "", colnames(design), fixed=TRUE)
}

# Generate contrasts information
contrasts <- makeContrasts(contrasts=contrastData, levels=design)

################################################################################
### Data Output
################################################################################

# Opening connection to output file
png(bcvOut, width=600, height=600)
plotBCV(data, main="BCV Plot")
imageData[1, ] <- c("BCV Plot", "bcvplot.png")
invisible(dev.off())

# Generate voom data and Mean-variance plot
png(voomOut, width=600, height=600)
vData <- voom(data, design = design, plot=TRUE)
imageData <- rbind(imageData, c("Voom Plot", "voomplot.png"))
invisible(dev.off())

for (i in 1:length(contrastData)){
  png(maOutPng[i], height=600, width=600)
  if (wantWeight){
    # Generating weights
    aw <- arrayWeights(vData, design)
    wts <- asMatrixWeights(aw, dim(vData)) * vData$w
    
    # Generating fit data and top table with weights
    voomFit <- lmFit(vData, design, weights=wts)
    voomFit <- contrasts.fit(voomFit, contrasts)
    voomFit <- eBayes(voomFit)
    top <- topTable(voomFit, coef=i, number=Inf, sort.by="P")
    write.table(top, file=topOut[i], row.names=FALSE, sep="\t")
    
    linkName <- paste0("Top Expressions Table(", contrastData[i],
                       ") (.tsv)")
    linkAddr <- paste0("toptab(", contrastData[i], ").tsv")
    if (nrow(linkData)==0){
      linkData[1, ] <- c(linkName, linkAddr)
    } else {
      linkData <- rbind(linkData, c(linkName, linkAddr))
    }
    
    # Plot MA (log ratios vs mean average) using limma package on weighted data
    status = decideTests(voomFit[, i])
    
    limma::plotMA(voomFit, status=status, array=i,
                  main=paste("MA Plot:", contrastData[i]), 
                  col=c("black","firebrick","green"),
                  xlab="Average Expression", ylab="logFC")
    abline(h=0, col="grey", lty=2)
    
    imgName <- paste0("MA Plot(", contrastData[i], ")")
    imgAddr <- paste0("maplot(", contrastData[i], ").png")
    imageData <- rbind(imageData, c(imgName, imgAddr))
    
    #print(summary(status))
  } else {
    # Generating fit data and top table without weights
    voomFitNoWeights <- lmFit(vData, design)
    voomFitNoWeights <- contrasts.fit(voomFitNoWeights, contrasts)
    voomFitNoWeights <- eBayes(voomFitNoWeights)
    topNoWeights <- topTable(voomFitNoWeights, coef=i, number=Inf, sort.by="P")
    write.table(topNoWeights, file=topOut[i], row.names=FALSE, sep="\t")
    
    linkName <- paste0("Top Expressions Table(", 
                       contrastData[i], ") (.tsv)")
    linkAddr <- paste0("toptab(", contrastData[i], ").tsv")
    if (nrow(linkData)==0){
      linkData[1, ] <- c(linkName, linkAddr)
    } else {
      linkData <- rbind(linkData, c(linkName, linkAddr))
    }    
    # Plot MA (log ratios vs mean average) using limma package on unweighted 
    # data
    status = decideTests(voomFitNoWeights[, i])
    
    limma::plotMA(voomFitNoWeights, status=status, array=i,
                  main=paste("MA Plot:", contrastData[i]), 
                  col=c("black","firebrick","green"),
                  xlab="Average Expression", ylab="logFC")
    abline(h=0, col="grey", lty=2)
    
    imgName <- paste0("MA Plot(", contrastData[i], ")")
    imgAddr <- paste0("maplot(", contrastData[i], ").png")
    imageData <- rbind(imageData, c(imgName, imgAddr))
    
    #print(summary(status))
  }
  invisible(dev.off())
}

if (wantRda){
  if (wantWeight){
    save(data, vData, labels, factors, wts, voomFit, top,
         file=rdaOut, ascii=TRUE)
  } else {
    save(data, vData, labels, factors, voomFitNoWeights, 
         topNoWeights, 
         file=rdaOut, ascii=TRUE)
  }
  linkData <- rbind(linkData, c("RData (.rda)", "objectDump.rda"))
}

################################################################################
### HTML Generation
################################################################################

# Clear file
cat("", file=htmlPath)

# Create cat function default path set, default seperator empty and appending
# true by default (Ripped straight from the cat function with altered argument
# defaults)
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

HtmlImage <- function(source, label=source, height=600, width=600){
  cata("<img src=\"", source, "\" alt=\"", label, "\" height=\"", height)
  cata("\" width=\"", width, "\"/>\n")
}

ListItem <- function(item){
  cata("<li>", item, "</li>\n")
}

cata("<html>\n")
HtmlHead("EdgeR Output")
cata("<body>\n")
cata("<h3>EdgeR Analysis Output:</h3>\n")

for (i in 1:nrow(imageData)){
  HtmlImage(imageData$Link[i], imageData$Label[i])
}

cata("<h4>Top Expressions Table(s):</h4>\n")

for (i in 1:nrow(linkData)){
  HtmlLink(linkData$Link[i], linkData$Label[i])
}
cata("<p>alt-click any of the links to download the file, or click the name ")
cata("of this task in the galaxy history panel and click on the floppy ")
cata("disk icon to download all files in a zip file.</p>\n")
cata("<p>.tsv files are tab seperated files that can be viewed using Excel ")
cata("or other spreadsheet programs</p>\n")

cata("<h4>Extra Information</h4>\n")
cata("<ul>\n")
if(wantWeight){
  ListItem("Weights were applied to samples")
} else {
  ListItem("Weights were not applied to samples")
}
cata("<li>", "Experiment was analysed using ", ncol(factors), " factor(s): ")
cata(names(factors), sep=", ", "</li>\n")
cata("</ul>\n")

cata("</body>\n")
cata("</html>")