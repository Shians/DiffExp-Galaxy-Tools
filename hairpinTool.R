################################################################################
### Input Processing
################################################################################

# Record starting time
timeStart <- as.character(Sys.time())

# Loading and checking required packages
library(methods, quietly=TRUE, warn.conflicts=FALSE)
library(statmod, quietly=TRUE, warn.conflicts=FALSE)
library(splines, quietly=TRUE, warn.conflicts=FALSE)
library(edgeR, quietly=TRUE, warn.conflicts=FALSE)
library(limma, quietly=TRUE, warn.conflicts=FALSE)

if(packageVersion("edgeR") < "3.5.24"){
  message("Please update 'edgeR' to version >= 3.5.24 to run this code")
}

# Grabbing arguments from command line
argv <- commandArgs(TRUE)

# Remove fastq file paths after collecting from argument vector
fastqPath <- as.character(gsub("fastq::", "", argv[grepl("fastq", argv)], 
                               fixed=TRUE))
argv <- argv[!grepl("fastq::", argv, fixed=TRUE)]
hairpinPath <- as.character(argv[1])
samplePath <- as.character(argv[2])
barStart <- as.numeric(argv[3])
barEnd <- as.numeric(argv[4])
hpStart <- as.numeric(argv[5])
hpEnd <- as.numeric(argv[6])
cpmReq <- as.numeric(argv[7])
sampleReq <- as.numeric(argv[8])
fdrThresh <- as.numeric(argv[9])
workMode <- as.character(argv[10])
htmlPath <- as.character(argv[11])
folderPath <- as.character(argv[12])
if(workMode=="classic"){
  pairData <- character()
  pairData[2] <- as.character(argv[13])
  pairData[1] <- as.character(argv[14])
} else if (workMode=="glm"){
  contrastData <- as.character(argv[13])
  roastOpt <- as.character(argv[14])
  hairpinReq <- as.numeric(argv[15])
  selectOpt <- as.character(argv[16])
  selectVals <- as.character(argv[17])
}

# Function to sanitise contrast equations so there are no whitespaces
# surrounding the arithmetic operators, leading or trailing whitespace
sanitiseEquation <- function(equation){
  equation <- gsub(" *[+] *", "+", equation)
  equation <- gsub(" *[-] *", "-", equation)
  equation <- gsub(" *[/] *", "/", equation)
  equation <- gsub(" *[*] *", "*", equation)
  equation <- gsub("^\\s+|\\s+$", "", equation)
  return(equation)
}

# Process arguments
if (workMode=="glm"){
  if (roastOpt=="yes"){
    wantRoast <- TRUE
  } else {
    wantRoast <- FALSE
  }
}

# Transform selection from string into index values
if (workMode=="glm"){
  if (selectOpt=="rank"){
    selectVals <- gsub(" ", "", selectVals, fixed=TRUE)
    selectVals <- unlist(strsplit(selectVals, ","))
    
    for (i in 1:length(selectVals)){
      if (grepl(":", selectVals[i], fixed=TRUE)){
        temp <- unlist(strsplit(selectVals[i], ":"))
        selectVals <- selectVals[-i]
        a <- as.numeric(temp[1])
        b <- as.numeric(temp[2])
        selectVals <- c(selectVals, a:b)
      }
    }
    selectVals <- as.numeric(unique(selectVals))
  } else {
    selectVals <- gsub(" ", "", selectVals, fixed=TRUE)
    selectVals <- unlist(strsplit(selectVals, " "))
  }
}

# Check if grouping variable has been specified
sampleData <- read.table(samplePath, sep="\t", header=TRUE)
if(!any(grepl("group", names(sampleData)))){
  stop("'group' column not specified in sample annotation file")
}

# Split up contrasts seperated by comma into a vector and replace spaces with
# periods
if(exists("contrastData")){
  contrastData <- unlist(strsplit(contrastData, split=","))
  contrastData <- sanitiseEquation(contrastData)
  contrastData <- gsub(" ", ".", contrastData, fixed=TRUE)
}

# Replace spaces with periods in pair data
if (exists("pairData")){
  pairData <- make.names(pairData)
}

# Generate output folder and paths
dir.create(folderPath)

# Function has string input and generates an output path string
makeOut <- function(filename){
  return(paste0(folderPath, "/", filename))
}

# Function has string input and generates both a pdf and png output strings
imgOut <- function(filename){
  assign(paste0(filename, "Png"), makeOut(paste0(filename,".png")), 
         envir = .GlobalEnv)
  assign(paste0(filename, "Pdf"), makeOut(paste0(filename,".pdf")),
         envir = .GlobalEnv)
}

imgOut("barHairpin")
imgOut("barIndex")
imgOut("mds")
imgOut("bcv")
if (workMode == "classic"){
  smearPng <- makeOut(paste0("smear(", pairData[2], "-", pairData[1],").png"))
  smearPdf <- makeOut(paste0("smear(", pairData[2], "-", pairData[1],").pdf"))
  topOut <- makeOut(paste0("toptag(", pairData[2], "-", pairData[1],").tsv"))
} else {
  roastOut <- makeOut("roast.tsv")
  smearPng <- character()
  smearPdf <- character()
  topOut <- character()
  barcodePng <- character()
  barcodePdf <- character()
  for (i in 1:length(contrastData)){
    smearPng[i] <- makeOut(paste0("smear(", contrastData[i], ").png"))
    smearPdf[i] <- makeOut(paste0("smear(", contrastData[i], ").pdf"))
    topOut[i] <- makeOut(paste0("toptag(", contrastData[i], ").tsv"))
    barcodePng[i] <- makeOut(paste0("barcode(", contrastData[i], ").png"))
    barcodePdf[i] <- makeOut(paste0("barcode(", contrastData[i], ").pdf"))
  }
}
# Initialise data for html links and images, table with the link label and
# link address
linkData <- data.frame(Label=character(), Link=character(),
                       stringsAsFactors=FALSE)
imageData <- data.frame(Label=character(), Link=character(),
                        stringsAsFactors=FALSE)

################################################################################
### Data Processing
################################################################################

# Use EdgeR hairpin process and capture outputs
hpReadout <- capture.output(
  data <- processHairpinReads(fastqPath, samplePath, hairpinPath,
                              hairpinStart=hpStart, hairpinEnd=hpEnd, 
                              verbose=TRUE)
)
closeAllConnections()
hpReadout <- hpReadout[!grepl("million", hpReadout)]
hpReadout <- hpReadout[hpReadout!=""]
hpReadout <- hpReadout[-(3:4)]
hpReadout <- gsub(" -- ", "", hpReadout, fixed=TRUE)


# Make the names of groups syntactically valid (replace spaces with periods)
data$samples$group<-make.names(data$samples$group)

# Filter hairpins with low counts
sel <- rowSums(cpm(data$counts)>cpmReq)>=sampleReq
data <- data[sel, ]

# Begin differential representation analysis
# Estimate dispersions
data <- estimateDisp(data)
#sqrt(data$common.dispersion)


# Plot number of hairpins that could be matched per sample
png(barIndexPng, width=600, height=600)
barplot(height<-colSums(data$counts), las=2, main="Counts per index", 
        cex.names=0.5, cex.axis=0.8, ylim=c(0, max(height)*1.2))
imageData[1, ] <- c("Counts per Index", "barIndex.png")
invisible(dev.off())

pdf(barIndexPdf)
barplot(height<-colSums(data$counts), las=2, main="Counts per index", 
        cex.names=0.5, cex.axis=0.8, ylim=c(0, max(height)*1.2))
linkData[1, ] <- c("Counts per Index Barplot (.pdf)", "barIndex.pdf")
invisible(dev.off())

# Plot per hairpin totals across all samples
png(barHairpinPng, width=600, height=600)
barplot(height<-rowSums(data$counts), las=2, main="Counts per hairpin",
        cex.names=0.5, cex.axis=0.8, ylim=c(0, max(height)*1.2))
imageData <- rbind(imageData, c("Counts per Hairpin", "barHairpin.png"))
invisible(dev.off())

pdf(barHairpinPdf)
barplot(height<-rowSums(data$counts), las=2, main="Counts per hairpin",
        cex.names=0.5, cex.axis=0.8, ylim=c(0, max(height)*1.2))
newEntry <- c("Counts per Hairpin Barplot (.pdf)", "barHairpin.pdf")
linkData <- rbind(linkData, newEntry)
invisible(dev.off())

# Make an MDS plot to visualise relationships between replicate samples
png(mdsPng, width=600, height=600)
plotMDS(data, labels=data$samples$group, col=as.numeric(data$samples$group),
        main="MDS Plot")
imageData <- rbind(imageData, c("MDS Plot", "mds.png"))
invisible(dev.off())

pdf(mdsPdf)
plotMDS(data, labels=data$samples$group, col=as.numeric(data$samples$group),
        main="MDS Plot")
newEntry <- c("MDS Plot (.pdf)", "mds.pdf")
linkData <- rbind(linkData, newEntry)
invisible(dev.off())

if (workMode=="classic"){
  # Assess differential representation using classic exact testing methodology 
  # in edgeR
  testData <- exactTest(data, pair=pairData)
  
  top <- topTags(testData, n=Inf)
  topIDs <- top$table[top$table$FDR < fdrThresh, 1]
  write.table(top, file=topOut, row.names=FALSE, sep="\t")
  linkName <- paste0("Top Tags Table(", pairData[2], "-", pairData[1], 
                     ") (.tsv)")
  linkAddr <- paste0("toptag(", pairData[2], "-", pairData[1], ").tsv")
  linkData <- rbind(linkData, c(linkName, linkAddr))
  
  # Select hairpins with FDR < 0.05 to highlight on plot
  png(smearPng, width=600, height=600)
  plotTitle <- gsub(".", " ", 
                    paste0("Smear Plot: ", pairData[2], "-", pairData[1]),
                    fixed = TRUE)
  plotSmear(testData, pair=c(pairData[1], pairData[2]), de.tags=topIDs, 
            pch=20, cex=0.5, main=plotTitle)
  abline(h = c(-1, 0, 1), col = c("dodgerblue", "yellow", "dodgerblue"), lty=2)
  imgName <- paste0("Smear Plot(", pairData[2], "-", pairData[1], ")")
  imgAddr <- paste0("smear(", pairData[2], "-", pairData[1],").png")
  imageData <- rbind(imageData, c(imgName, imgAddr))
  invisible(dev.off())
  
  pdf(smearPdf)
  plotTitle <- gsub(".", " ", 
                    paste0("Smear Plot: ", pairData[2], "-", pairData[1]),
                    fixed = TRUE)
  plotSmear(testData, pair=c(pairData[1], pairData[2]), de.tags=topIDs, 
            pch=20, cex=0.5, main=plotTitle)
  abline(h = c(-1, 0, 1), col = c("dodgerblue", "yellow", "dodgerblue"), lty=2)
  imgName <- paste0("Smear Plot(", pairData[2], "-", pairData[1], ") (.pdf)")
  imgAddr <- paste0("smear(", pairData[2], "-", pairData[1], ").pdf")
  linkData <- rbind(linkData, c(imgName, imgAddr))
  invisible(dev.off())
} else if (workMode=="glm"){
  # Generating design information
  factors <- factor(data$sample$group)
  design <- model.matrix(~0+factors)
  
  colnames(design) <- gsub("factors", "", colnames(design), fixed=TRUE)
  
  # Split up contrasts seperated by comma into a vector
  contrastData <- unlist(strsplit(contrastData, split=","))
  for (i in 1:length(contrastData)){
    # Generate contrasts information
    contrasts <- makeContrasts(contrasts=contrastData[i], levels=design)
    
    # Fit negative bionomial GLM
    fit = glmFit(data, design)
    # Carry out Likelihood ratio test
    testData = glmLRT(fit, contrast=contrasts)
    
    # Select hairpins with FDR < 0.05 to highlight on plot
    top <- topTags(testData, n=Inf)
    topIDs <- top$table[top$table$FDR < fdrThresh, 1]
    write.table(top, file=topOut[i], row.names=FALSE, sep="\t")
    
    linkName <- paste0("Top Tags Table(", contrastData[i], ") (.tsv)")
    linkAddr <- paste0("toptag(", contrastData[i], ").tsv")
    linkData <- rbind(linkData, c(linkName, linkAddr))
    
    # Make a plot of logFC versus logCPM
    png(smearPng[i], height=600, width=600)
    plotTitle <- paste("Smear Plot:", gsub(".", " ", contrastData[i], 
                       fixed=TRUE))
    plotSmear(testData, de.tags=topIDs, pch=20, cex=0.5, main=plotTitle)
    abline(h=c(-1, 0, 1), col=c("dodgerblue", "yellow", "dodgerblue"), lty=2)
    
    imgName <- paste0("Smear Plot(", contrastData[i], ")")
    imgAddr <- paste0("smear(", contrastData[i], ").png")
    imageData <- rbind(imageData, c(imgName, imgAddr))
    invisible(dev.off())
    
    pdf(smearPdf[i])
    plotTitle <- paste("Smear Plot:", gsub(".", " ", contrastData[i], 
                       fixed=TRUE))
    plotSmear(testData, de.tags=topIDs, pch=20, cex=0.5, main=plotTitle)
    abline(h=c(-1, 0, 1), col=c("dodgerblue", "yellow", "dodgerblue"), lty=2)
    
    linkName <- paste0("Smear Plot(", contrastData[i], ") (.pdf)")
    linkAddr <- paste0("smear(", contrastData[i], ").pdf")
    linkData <- rbind(linkData, c(linkName, linkAddr))
    invisible(dev.off())
    
    genes <- as.character(data$genes$Gene)
    unq <- unique(genes)
    unq <- unq[!is.na(unq)]
    geneList = list()
    for(gene in unq){
      if (length(which(genes==gene)) >= hairpinReq){
        geneList[[gene]] = which(genes==gene)
      }
    }
    
    if (wantRoast){
      # Input preparaton for roast
      nrot = 9999
      set.seed(602214129)
      roastData = mroast(data, index=geneList, design=design,
                         contrast=contrasts, nrot=nrot)
      
      write.table(roastData, file=roastOut, sep="\t")
      linkData <- rbind(linkData, c("Gene Level Analysis Table(.tsv)", 
                                    "roast.tsv"))
      if (selectOpt=="rank"){
        selectedGenes <- rownames(roastData)[selectVals]
      } else {
        selectedGenes <- selectVals
      }
      
      
      png(barcodePng[i], width=600, height=length(selectedGenes)*150)
      par(mfrow= c(length(selectedGenes), 1))
      for(gene in selectedGenes){
        barcodeplot(testData$table$logFC, index=geneList[[gene]],
                    main=paste("Barcode Plot for", gene, "(logFCs)", 
                               gsub(".", " ", contrastData[i])),
                    labels=c("Positive logFC", "Negative logFC"))
      }
      imgName <- paste0("Barcode Plot(", contrastData[i], ")")
      imgAddr <- paste0("barcode(", contrastData[i], ").png")
      imageData <- rbind(imageData, c(imgName, imgAddr))
      dev.off()
      
      pdf(barcodePdf[i], width=8, height=2)
      for(gene in selectedGenes){
        barcodeplot(testData$table$logFC, index=geneList[[gene]],
                    main=paste("Barcode Plot for", gene, "(logFCs)", 
                               gsub(".", " ", contrastData[i])),
                    labels=c("Positive logFC", "Negative logFC"))
      }
      linkName <- paste0("Barcode Plot(", contrastData[i], ") (.pdf)")
      linkAddr <- paste0("barcode(", contrastData[i], ").pdf")
      linkData <- rbind(linkData, c(linkName, linkAddr))
      dev.off()
    }
  }
}

# Record ending time
timeEnd <- as.character(Sys.time())

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
  if (grepl("barcode", imageData$Link[i])){
    HtmlImage(imageData$Link[i], imageData$Label[i], 
              height=length(selectedGenes)*150)
  } else {
    HtmlImage(imageData$Link[i], imageData$Label[i])
  }
}
cata("<br/>\n")

cata("<h4>Input Summary:</h4>\n")
cata("<ul>\n")
ListItem(hpReadout[1])
ListItem(hpReadout[2])
cata("</ul>\n")
cata(hpReadout[3], "<br/>\n")
cata("<ul>\n")
ListItem(hpReadout[4])
ListItem(hpReadout[7])
cata("</ul>\n")
cata(hpReadout[8:11], sep="<br/>\n")

cata("<h4>Plots:</h4>\n")
for (i in 1:nrow(linkData)){
  if (!grepl(".tsv", linkData$Link[i])){
    HtmlLink(linkData$Link[i], linkData$Label[i])
  }
}

cata("<h4>Tables:</h4>\n")
for (i in 1:nrow(linkData)){
  if (grepl(".tsv", linkData$Link[i])){
    HtmlLink(linkData$Link[i], linkData$Label[i])
  }
}

cata("<p>alt-click any of the links to download the file, or click the name ")
cata("of this task in the galaxy history panel and click on the floppy ")
cata("disk icon to download all files in a zip archive.</p>\n")
cata("<p>.tsv files are tab seperated files that can be viewed using Excel ")
cata("or other spreadsheet programs</p>\n")
cata("<table border=\"0\">\n")

cata("<tr>\n")
TableItem("Task started at:"); TableItem(timeStart)
cata("</tr>\n")
cata("<tr>\n")
TableItem("Task ended at:"); TableItem(timeEnd)
cata("</tr>\n")

cata("</body>\n")
cata("</html>")
