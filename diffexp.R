# This tool takes in a matrix of feature counts as well as gene annotations and
# outputs a table of top expressions as well as various plots for differential
# expression analysis
#
# ARGS: 1.countPath       -Path to RData input containing counts
#       2.annoPath        -Path to RData input containing gene annotations
#       3.htmlPath        -Path to html file linking to other outputs
#       4.outPath         -Path to folder to write all output to
#       5.rdaOpt          -String specifying if RData should be saved
#       6.normOpt         -String specifying type of normalisation used
#       7.weightOpt       -String specifying usage of weights
#       8.contrastData    -String containing contrasts of interest
#       9.cpmReq          -Float specifying cpm requirement
#       10.sampleReq      -Integer specifying cpm requirement
#       11.factorData     -String containing factor names and levels
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

# Generate output folder and paths
dir.create(outPath)
makeOut <- function(filename){
  return(paste0(outPath, "/", filename))
}

bcvOutPdf <- makeOut("bcvplot.pdf")
bcvOutPng <- makeOut("bcvplot.png")
voomOutPdf <- makeOut("voomplot.pdf")
voomOutPng <- makeOut("voomplot.png") 
maOutPdf <- character()   # Initialise character vector
maOutPng <- character()
topOut <- character()
for (i in 1:length(contrastData)){
  maOutPdf[i] <- makeOut(paste0("maplot(", contrastData[i], ").pdf"))
  maOutPng[i] <- makeOut(paste0("maplot(", contrastData[i], ").png"))
  topOut[i] <- makeOut(paste0("toptab(", contrastData[i], ").tsv"))
}                         # Save output paths for each contrast as vectors
rdaOut <- makeOut("objectDump.rda")

# Initialise data for html links and images, data frame with columns Label and 
# Link
linkData <- data.frame(Label=character(), Link=character(),
                       stringsAsFactors=FALSE)
imageData <- data.frame(Label=character(), Link=character(),
                        stringsAsFactors=FALSE)

# Read in counts and geneanno data
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
sampleanno <- data.frame("sampleID"=samplenames, factors)

# Generating the DGEList object "data"
data$samples <- sampleanno
data$samples$lib.size <- colSums(data$counts)
data$samples$norm.factors <- 1
row.names(data$samples) <- colnames(data$counts)
data <- new("DGEList", data)

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

# Calculating normalising factor, estimating dispersion
data <- calcNormFactors(data, method=normOpt)
#data <- estimateDisp(data, design=design, robust=TRUE)

# Generate contrasts information
contrasts <- makeContrasts(contrasts=contrastData, levels=design)

# Name rows of factors according to their sample
row.names(factors) <- names(data$counts)

################################################################################
### Data Output
################################################################################

# BCV Plot
#png(bcvOutPng, width=600, height=600)
#plotBCV(data, main="BCV Plot")
#imageData[1, ] <- c("BCV Plot", "bcvplot.png")
#invisible(dev.off())

#pdf(bcvOutPdf)
#plotBCV(data, main="BCV Plot")
#invisible(dev.off())


for (i in 1:length(contrastData)){
  if (wantWeight){
    # Generating weights
    vData <- voom(data, design = design)
    aw <- arrayWeights(vData, design)
    
    # Creating voom data object and plot
    png(voomOutPng, width=600, height=600)
    vData <- voom(data, design=design, weights=aw, plot=TRUE)
    imageData[1, ] <- c("Voom Plot", "voomplot.png")
    invisible(dev.off())
    
    pdf(voomOutPdf)
    vData <- voom(data, design=design, weights=aw, plot=TRUE)
    linkData[1, ] <- c("Voom Plot (.pdf)", voomOutPdf)
    invisible(dev.off())
    
    # Generating fit data and top table with weights
    wts <- asMatrixWeights(aw, dim(vData)) * vData$w
    voomFit <- lmFit(vData, design, weights=wts)
    
    
    #print(summary(status))
  } else {
    # Creating voom data object and plot
    png(voomOutPng, width=600, height=600)
    vData <- voom(data, design=design, plot=TRUE)
    imageData[1, ] <- c("Voom Plot", "voomplot.png")
    invisible(dev.off())
    
    pdf(voomOutPdf)
    vData <- voom(data, design=design, plot=TRUE)
    linkData[1, ] <- c("Voom Plot (.pdf)", voomOutPdf)
    invisible(dev.off())
    
    # Generate voom fit
    voomFit <- lmFit(vData, design)
    
  }
  
  # Fit linear model and estimate dispersion with eBayes
  voomFit <- contrasts.fit(voomFit, contrasts)
  voomFit <- eBayes(voomFit)
  
  status = decideTests(voomFit[, i])
  
  # Write top expressions table
  top <- topTable(voomFit, coef=i, number=Inf, sort.by="P")
  write.table(top, file=topOut[i], row.names=FALSE, sep="\t")
  
  linkName <- paste0("Top Expressions Table(", contrastData[i], ") (.tsv)")
  linkAddr <- paste0("toptab(", contrastData[i], ").tsv")
  linkData <- rbind(linkData, c(linkName, linkAddr))
  
  # Plot MA (log ratios vs mean average) using limma package on weighted data
  pdf(maOutPdf[i])
  limma::plotMA(voomFit, status=status, array=i,
                main=paste("MA Plot:", contrastData[i]), 
                col=c("black","firebrick","green"), value=c("0", "-1", "1"),
                xlab="Average Expression", ylab="logFC")
  
  abline(h=0, col="grey", lty=2)
  
  linkName <- paste0("MA Plot(", contrastData[i], ")", " (.pdf)")
  linkAddr <- paste0("maplot(", contrastData[i], ").pdf")
  linkData <- rbind(linkData, c(linkName, linkAddr))
  invisible(dev.off())
  
  png(maOutPng[i], height=600, width=600)
  limma::plotMA(voomFit, status=status, array=i,
                main=paste("MA Plot:", contrastData[i]), 
                col=c("black","firebrick","green"), value=c("0", "-1", "1"),
                xlab="Average Expression", ylab="logFC")
  
  abline(h=0, col="grey", lty=2)
  
  imgName <- paste0("MA Plot(", contrastData[i], ")")
  imgAddr <- paste0("maplot(", contrastData[i], ").png")
  imageData <- rbind(imageData, c(imgName, imgAddr))
  invisible(dev.off())
}

# Save relevant items as rda object
if (wantRda){
  if (wantWeight){
    save(data, vData, labels, factors, wts, voomFit, top, contrasts, design,
         file=rdaOut, ascii=TRUE)
  } else {
    save(data, vData, labels, factors, voomFit, top, contrasts, design,
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

# Function to write code for html head and title
HtmlHead <- function(title){
  cata("<head>\n")
  cata("<title>", title, "</title>\n")
  cata("</head>\n")
}

# Function to write code for html links
HtmlLink <- function(address, label=address){
  cata("<a href=\"", address, "\" target=\"_blank\">", label, "</a><br />\n")
}

# Function to write code for html images
HtmlImage <- function(source, label=source, height=600, width=600){
  cata("<img src=\"", source, "\" alt=\"", label, "\" height=\"", height)
  cata("\" width=\"", width, "\"/>\n")
}

# Function to write code for html list items
ListItem <- function(...){
  cata("<li>", ..., "</li>\n")
}

TableItem <- function(...){
  cata("<td>", ..., "</td>\n")
}

TableHeadItem <- function(...){
  cata("<th>", ..., "</th>\n")
}

cata("<html>\n")
HtmlHead("EdgeR Output")
cata("<body>\n")
cata("<h3>EdgeR Analysis Output:</h3>\n")
cata("All images displayed have PDF copy at the bottom of the page, these can ")
cata("exported in a pdf viewer to high resolution image format. <br/>\n")
for (i in 1:nrow(imageData)){
  HtmlImage(imageData$Link[i], imageData$Label[i])
}

cata("<h4>Additional Output:</h4>\n")

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
if(cpmReq!=0 && sampleReq!=0){
  filterString <- paste("Genes that do not have more than", cpmReq,
                        "CPM in at least", sampleReq, "samples are considered",
                        "unexpressed and filtered out.")
  ListItem(filterString)
}
ListItem(normOpt, " was the method used to normalise library sizes.")
if(wantWeight){
  ListItem("Weights were applied to samples.")
} else {
  ListItem("Weights were not applied to samples.")
}
#cata("<li>", "Experiment was analysed using ", ncol(factors), " factor(s): ")
#cata(names(factors), sep=", ", "</li>\n")
cata("</ul>\n")

cata("<h4>Summary of experimental data</h4>\n")

cata("<table border=\"1\">\n")
cata("<tr>\n")
TableItem()
for (i in names(factors)){
  TableHeadItem(i)
}
cata("</tr>\n")

for (i in 1:nrow(factors)){
  cata("<tr>\n")
  TableHeadItem(row.names(factors)[i])
  for (j in ncol(factors)){
    TableItem(as.character(factors[i, j]))
  }
  cata("</tr>\n")
}
cata("</table>")

cit <- character()
link <- character()
link[1] <- paste0("<a href=\"",
                  "http://www.bioconductor.org/packages/release/bioc/",
                  "vignettes/limma/inst/doc/usersguide.pdf",
                  "\">", "limma User's Guide", "</a>.")

link[2] <- paste0("<a href=\"",
                  "http://www.bioconductor.org/packages/release/bioc/",
                  "vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf",
                  "\">", "edgeR User's Guide", "</a>")

cit[1] <- paste("Please cite the paper below for the limma software itself.",
                "Please also try to cite the appropriate methodology articles",
                "that describe the statistical methods implemented in limma,",
                "depending on which limma functions you are using. The",
                "methodology articles are listed in Section 2.1 of the",
                link[1])
cit[2] <- paste("Smyth, GK (2005). Limma: linear models for microarray data.",
                "In: 'Bioinformatics and Computational Biology Solutions using",
                "R and Bioconductor'. R. Gentleman, V. Carey, S. doit,.",
                "Irizarry, W. Huber (eds), Springer, New York, pages 397-420.")
cit[3] <- paste("Please cite the first paper for the software itself and the",
                "other papers for the various original statistical methods",
                "implemented in edgeR.  See Section 1.2 in the", link[2],
                "for more detail.")
cit[4] <- paste("Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a",
                "Bioconductor package for differential expression analysis",
                "of digital gene expression data. Bioinformatics 26, 139-140")
cit[5] <- paste("Robinson MD and Smyth GK (2007). Moderated statistical tests",
                "for assessing differences in tag abundance. Bioinformatics",
                "23, 2881-2887")
cit[6] <- paste("Robinson MD and Smyth GK (2008). Small-sample estimation of",
                "negative binomial dispersion, with applications to SAGE data.",
                "Biostatistics, 9, 321-332")
cit[7] <- paste("McCarthy DJ, Chen Y and Smyth GK (2012). Differential",
                "expression analysis of multifactor RNA-Seq experiments with",
                "respect to biological variation. Nucleic Acids Research 40,",
                "4288-4297")

cata("<h3>Citations</h3>\n")

cata("<h4>limma</h4>\n")
cata(cit[1])
cata("<ul>\n")
ListItem(cit[2])
cata("</ul>\n")

cata("<h4>edgeR</h4>\n")
cata(cit[3])
cata("<ul>\n")
ListItem(cit[4])
ListItem(cit[5])
ListItem(cit[6])
ListItem(cit[7])
cata("</ul>\n")

cata("</body>\n")
cata("</html>")