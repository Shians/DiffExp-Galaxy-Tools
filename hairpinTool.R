# ARGS: 1.inputType         -String specifying format of input (fastq or table)
#    IF inputType is "fastQ":
#       2*.fastqPath        -One or more strings specifying path to fastq files
#       2.annoPath        -String specifying path to hairpin annotation table
#       3.samplePath        -String specifying path to sample annotation table
#       4.barStart          -Integer specifying starting position of barcode
#       5.barEnd            -Integer specifying ending position of barcode
#       6.hpStart           -Integer specifying startins position of hairpin
#                            unique region
#       7.hpEnd             -Integer specifying ending position of hairpin
#                            unique region
#    ###   
#    IF inputType is "counts":
#       2.countPath         -String specifying path to count table
#       3.annoPath          -String specifying path to hairpin annotation table
#       4.samplePath        -String specifying path to sample annotation table
#    ###
#       8.cpmReq            -Float specifying cpm requirement
#       9.sampleReq         -Integer specifying cpm requirement
#       10.fdrThresh        -Float specifying the FDR requirement
#       11.lfcThresh        -Float specifying the log-fold-change requirement
#       12.workMode         -String specifying exact test or GLM usage
#       13.htmlPath         -String specifying path to HTML file
#       14.folderPath       -STring specifying path to folder for output
#    IF workMode is "classic" (exact test)
#       15.pairData[2]      -String specifying first group for exact test
#       16.pairData[1]      -String specifying second group for exact test
#    ###
#    IF workMode is "glm"
#       15.contrastData     -String specifying contrasts to be made
#       16.roastOpt         -String specifying usage of gene-wise tests
#       17.hairpinReq       -String specifying hairpin requirement for gene-
#                            wise test
#       18.selectOpt        -String specifying type of selection for barcode
#                            plots
#       19.selectVals       -String specifying members selected for barcode
#                            plots
#
# OUT:  Bar Plot of Counts Per Index
#       Bar Plot of Counts Per Hairpin
#       MDS Plot
#       Smear Plot
#       Barcode Plots (If Genewise testing was selected)
#       Top Expression Table
#       HTML file linking to the ouputs
#
# Author: Shian Su - registertonysu@gmail.com - Jan 2014

# Record starting time
timeStart <- as.character(Sys.time())

# Loading and checking required packages
library(methods, quietly=TRUE, warn.conflicts=FALSE)
library(statmod, quietly=TRUE, warn.conflicts=FALSE)
library(splines, quietly=TRUE, warn.conflicts=FALSE)
library(edgeR, quietly=TRUE, warn.conflicts=FALSE)
library(limma, quietly=TRUE, warn.conflicts=FALSE)

if (packageVersion("edgeR") < "3.5.23") {
  message("Please update 'edgeR' to version >= 3.5.23 to run this script")
}

################################################################################
### Function declarations
################################################################################

# Function to sanitise contrast equations so there are no whitespaces
# surrounding the arithmetic operators, leading or trailing whitespace
sanitiseEquation <- function(equation) {
  equation <- gsub(" *[+] *", "+", equation)
  equation <- gsub(" *[-] *", "-", equation)
  equation <- gsub(" *[/] *", "/", equation)
  equation <- gsub(" *[*] *", "*", equation)
  equation <- gsub("^\\s+|\\s+$", "", equation)
  return(equation)
}

# Function to sanitise group information
sanitiseGroups <- function(string) {
  string <- gsub(" *[,] *", ",", string)
  string <- gsub("^\\s+|\\s+$", "", string)
  return(string)
}

# Function to change periods to whitespace in a string
unmake.names <- function(string) {
  string <- gsub(".", " ", string, fixed=TRUE)
  return(string)
}

# Function has string input and generates an output path string
makeOut <- function(filename) {
  return(paste0(folderPath, "/", filename))
}

# Function has string input and generates both a pdf and png output strings
imgOut <- function(filename) {
  assign(paste0(filename, "Png"), makeOut(paste0(filename,".png")), 
         envir = .GlobalEnv)
  assign(paste0(filename, "Pdf"), makeOut(paste0(filename,".pdf")),
         envir = .GlobalEnv)
}

# Create cat function default path set, default seperator empty and appending
# true by default (Ripped straight from the cat function with altered argument
# defaults)
cata <- function(..., file = htmlPath, sep = "", fill = FALSE, labels = NULL, 
                 append = TRUE) {
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
HtmlHead <- function(title) {
  cata("<head>\n")
  cata("<title>", title, "</title>\n")
  cata("</head>\n")
}

# Function to write code for html links
HtmlLink <- function(address, label=address) {
  cata("<a href=\"", address, "\" target=\"_blank\">", label, "</a><br />\n")
}

# Function to write code for html images
HtmlImage <- function(source, label=source, height=600, width=600) {
  cata("<img src=\"", source, "\" alt=\"", label, "\" height=\"", height)
  cata("\" width=\"", width, "\"/>\n")
}

# Function to write code for html list items
ListItem <- function(...) {
  cata("<li>", ..., "</li>\n")
}

TableItem <- function(...) {
  cata("<td>", ..., "</td>\n")
}

TableHeadItem <- function(...) {
  cata("<th>", ..., "</th>\n")
}
################################################################################
### Input Processing
################################################################################

# Grabbing arguments from command line
argv <- commandArgs(TRUE)

# Remove fastq file paths after collecting from argument vector
inputType <- as.character(argv[1])
if (inputType=="fastq") {
  fastqPath <- as.character(gsub("fastq::", "", argv[grepl("fastq::", argv)], 
                                 fixed=TRUE))
  argv <- argv[!grepl("fastq::", argv, fixed=TRUE)]
  annoPath <- as.character(argv[2])
  samplePath <- as.character(argv[3])
  barStart <- as.numeric(argv[4])
  barEnd <- as.numeric(argv[5])
  hpStart <- as.numeric(argv[6])
  hpEnd <- as.numeric(argv[7])
} else if (inputType=="counts") {
  countPath <- as.character(argv[2])
  annoPath <- as.character(argv[3])
  samplePath <- as.character(argv[4])
}
  
cpmReq <- as.numeric(argv[8])
sampleReq <- as.numeric(argv[9])
fdrThresh <- as.numeric(argv[10])
lfcThresh <- as.numeric(argv[11])
workMode <- as.character(argv[12])
htmlPath <- as.character(argv[13])
folderPath <- as.character(argv[14])
if (workMode=="classic") {
  pairData <- character()
  pairData[2] <- as.character(argv[15])
  pairData[1] <- as.character(argv[16])
} else if (workMode=="glm") {
  contrastData <- as.character(argv[15])
  roastOpt <- as.character(argv[16])
  hairpinReq <- as.numeric(argv[17])
  selectOpt <- as.character(argv[18])
  selectVals <- as.character(argv[19])
}

# Read in inputs

samples <- read.table(samplePath, header=TRUE, sep="\t")
anno <- read.table(annoPath, header=TRUE, sep="\t")
if (inputType=="counts") {
  counts <- read.table(countPath, header=TRUE, sep="\t")
}

###################### Check inputs for correctness ############################
samples$ID <- make.names(samples$ID)

if (!any(grepl("group", names(samples)))) {
  stop("'group' column not specified in sample annotation file")
} # Check if grouping variable has been specified

if (any(table(samples$ID)>1)){
  tab <- table(samples$ID)
  offenders <- paste(names(tab[tab>1]), collapse=", ")
  offenders <- unmake.names(offenders)
  stop("'ID' column of sample annotation must have unique values, values ",
       offenders, " are repeated")
} # Check that IDs in sample annotation are unique

if (inputType=="fastq") {
  
  if (any(table(anno$ID)>1)){
    tab <- table(anno$ID)
    offenders <- paste(names(tab[tab>1]), collapse=", ")
    stop("'ID' column of hairpin annotation must have unique values, values ",
    offenders, " are repeated")
  } # Check that IDs in hairpin annotation are unique
  
} else if (inputType=="counts") {
  if (any(is.na(match(samples$ID, colnames(counts))))) {
    stop("not all samples have groups specified")
  } # Check that a group has be specifed for each sample
  
  if (any(table(counts$ID)>1)){
    tab <- table(counts$ID)
    offenders <- paste(names(tab[tab>1]), collapse=", ")
    stop("'ID' column of count table must have unique values, values ",
    offenders, " are repeated")
  } # Check that IDs in count table are unique
}
if (workMode=="glm") {
  if (roastOpt == "yes") {
    if (is.na(match("Gene", colnames(anno)))) {
      tempStr <- paste("Gene-wise tests selected but'Gene' column not",
      "specified in hairpin annotation file")
      stop(tempStr)
    }
  }
}

################################################################################

# Process arguments
if (workMode=="glm") {
  if (roastOpt=="yes") {
    wantRoast <- TRUE
  } else {
    wantRoast <- FALSE
  }
}

# Split up contrasts seperated by comma into a vector and replace spaces with
# periods
if (exists("contrastData")) {
  contrastData <- unlist(strsplit(contrastData, split=","))
  contrastData <- sanitiseEquation(contrastData)
  contrastData <- gsub(" ", ".", contrastData, fixed=TRUE)
}

# Replace spaces with periods in pair data
if (exists("pairData")) {
  pairData <- make.names(pairData)
}

# Generate output folder and paths
dir.create(folderPath, showWarnings=FALSE)

# Generate links for outputs
imgOut("barHairpin")
imgOut("barIndex")
imgOut("mds")
imgOut("bcv")
if (workMode == "classic") {
  smearPng <- makeOut(paste0("smear(", pairData[2], "-", pairData[1],").png"))
  smearPdf <- makeOut(paste0("smear(", pairData[2], "-", pairData[1],").pdf"))
  topOut <- makeOut(paste0("toptag(", pairData[2], "-", pairData[1],").tsv"))
} else if (workMode=="glm") {
  smearPng <- character()
  smearPdf <- character()
  topOut <- character()
  roastOut <- character()
  barcodePng <- character()
  barcodePdf <- character()
  for (i in 1:length(contrastData)) {
    smearPng[i] <- makeOut(paste0("smear(", contrastData[i], ").png"))
    smearPdf[i] <- makeOut(paste0("smear(", contrastData[i], ").pdf"))
    topOut[i] <- makeOut(paste0("toptag(", contrastData[i], ").tsv"))
    roastOut[i] <- makeOut(paste0("roast(", contrastData[i], ").tsv"))
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

# Transform gene selection from string into index values for mroast
if (workMode=="glm") {
  if (selectOpt=="rank") {
    selectVals <- gsub(" ", "", selectVals, fixed=TRUE)
    selectVals <- unlist(strsplit(selectVals, ","))
    
    for (i in 1:length(selectVals)) {
      if (grepl(":", selectVals[i], fixed=TRUE)) {
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
                                                  
if (inputType=="fastq") {
  # Use EdgeR hairpin process and capture outputs
  hpReadout <- capture.output(
  data <- processHairpinReads(fastqPath, samplePath, annoPath,
                              hairpinStart=hpStart, hairpinEnd=hpEnd, 
                              verbose=TRUE)
  )
  
  # Remove function output entries that show processing data or is empty
  hpReadout <- hpReadout[hpReadout!=""]
  hpReadout <- hpReadout[!grepl("Processing", hpReadout)]
  hpReadout <- hpReadout[!grepl("in file", hpReadout)]
  hpReadout <- gsub(" -- ", "", hpReadout, fixed=TRUE)
  
  # Make the names of groups syntactically valid (replace spaces with periods)
  data$samples$group <- make.names(data$samples$group)
} else if (inputType=="counts") {
  # Process counts information, set ID column to be row names
  rownames(counts) <- counts$ID
  counts <- counts[ , !(colnames(counts)=="ID")]
  countsRows <- nrow(counts)
  
  # Process group information
  factors <- samples$group[match(samples$ID, colnames(counts))]
  annoRows <- nrow(anno)
  anno <- anno[match(rownames(counts), anno$ID), ]
  annoMatched <- sum(!is.na(anno$ID))
  
  if (any(is.na(anno$ID))) {
    warningStr <- paste("count table contained more hairpins than",
                        "specified in hairpin annotation file")
    warning(warningStr)
  }
  
  # Filter out rows with zero counts
  sel <- rowSums(counts)!=0
  counts <- counts[sel, ]
  anno <- anno[sel, ]
  
  # Create DGEList
  data <- DGEList(counts=counts, lib.size=colSums(counts), 
                  norm.factors=rep(1,ncol(counts)), genes=anno, group=factors)
  
  # Make the names of groups syntactically valid (replace spaces with periods)
  data$samples$group <- make.names(data$samples$group)
}

# Filter hairpins with low counts
preFilterCount <- nrow(data)
sel <- rowSums(cpm(data$counts) > cpmReq) >= sampleReq
data <- data[sel, ]
postFilterCount <- nrow(data)
filteredCount <- preFilterCount-postFilterCount

# Estimate dispersions
data <- estimateDisp(data)
commonBCV <- sqrt(data$common.dispersion)

################################################################################
### Output Processing
################################################################################

# Plot number of hairpins that could be matched per sample
png(barIndexPng, width=600, height=600)
barplot(height<-colSums(data$counts), las=2, main="Counts per index", 
        cex.names=1.0, cex.axis=0.8, ylim=c(0, max(height)*1.2))
imageData[1, ] <- c("Counts per Index", "barIndex.png")
invisible(dev.off())

pdf(barIndexPdf)
barplot(height<-colSums(data$counts), las=2, main="Counts per index", 
        cex.names=1.0, cex.axis=0.8, ylim=c(0, max(height)*1.2))
linkData[1, ] <- c("Counts per Index Barplot (.pdf)", "barIndex.pdf")
invisible(dev.off())

# Plot per hairpin totals across all samples
png(barHairpinPng, width=600, height=600)
if (nrow(data$counts)<50) {
  barplot(height<-rowSums(data$counts), las=2, main="Counts per hairpin",
          cex.names=0.8, cex.axis=0.8, ylim=c(0, max(height)*1.2))
} else {
  barplot(height<-rowSums(data$counts), las=2, main="Counts per hairpin",
          cex.names=0.8, cex.axis=0.8, ylim=c(0, max(height)*1.2),
          names.arg=FALSE)
}
imageData <- rbind(imageData, c("Counts per Hairpin", "barHairpin.png"))
invisible(dev.off())

pdf(barHairpinPdf)
if (nrow(data$counts)<50) {
  barplot(height<-rowSums(data$counts), las=2, main="Counts per hairpin",
          cex.names=0.8, cex.axis=0.8, ylim=c(0, max(height)*1.2))
} else {
  barplot(height<-rowSums(data$counts), las=2, main="Counts per hairpin",
          cex.names=0.8, cex.axis=0.8, ylim=c(0, max(height)*1.2),
          names.arg=FALSE)
}
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

if (workMode=="classic") {
  # Assess differential representation using classic exact testing methodology 
  # in edgeR
  testData <- exactTest(data, pair=pairData)
  
  top <- topTags(testData, n=Inf)
  topIDs <- top$table[(top$table$FDR < fdrThresh) &
                      (abs(top$table$logFC) > lfcThresh), 1]
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
  plotSmear(testData, de.tags=topIDs, 
            pch=20, cex=1.0, main=plotTitle)
  abline(h = c(-1, 0, 1), col = c("dodgerblue", "yellow", "dodgerblue"), lty=2)
  imgName <- paste0("Smear Plot(", pairData[2], "-", pairData[1], ")")
  imgAddr <- paste0("smear(", pairData[2], "-", pairData[1],").png")
  imageData <- rbind(imageData, c(imgName, imgAddr))
  invisible(dev.off())
  
  pdf(smearPdf)
  plotTitle <- gsub(".", " ", 
                    paste0("Smear Plot: ", pairData[2], "-", pairData[1]),
                    fixed = TRUE)
  plotSmear(testData, de.tags=topIDs, 
            pch=20, cex=1.0, main=plotTitle)
  abline(h = c(-1, 0, 1), col = c("dodgerblue", "yellow", "dodgerblue"), lty=2)
  imgName <- paste0("Smear Plot(", pairData[2], "-", pairData[1], ") (.pdf)")
  imgAddr <- paste0("smear(", pairData[2], "-", pairData[1], ").pdf")
  linkData <- rbind(linkData, c(imgName, imgAddr))
  invisible(dev.off())
} else if (workMode=="glm") {
  # Generating design information
  factors <- factor(data$sample$group)
  design <- model.matrix(~0+factors)
  
  colnames(design) <- gsub("factors", "", colnames(design), fixed=TRUE)
  
  # Split up contrasts seperated by comma into a vector
  contrastData <- unlist(strsplit(contrastData, split=","))
  for (i in 1:length(contrastData)) {
    # Generate contrasts information
    contrasts <- makeContrasts(contrasts=contrastData[i], levels=design)
    
    # Fit negative bionomial GLM
    fit = glmFit(data, design)
    # Carry out Likelihood ratio test
    testData = glmLRT(fit, contrast=contrasts)
    
    # Select hairpins with FDR < 0.05 to highlight on plot
    top <- topTags(testData, n=Inf)
    topIDs <- top$table[(top$table$FDR < fdrThresh) &
                        (abs(top$table$logFC) > lfcThresh), 1]
    write.table(top, file=topOut[i], row.names=FALSE, sep="\t")
    
    linkName <- paste0("Top Tags Table(", contrastData[i], ") (.tsv)")
    linkAddr <- paste0("toptag(", contrastData[i], ").tsv")
    linkData <- rbind(linkData, c(linkName, linkAddr))
    
    # Make a plot of logFC versus logCPM
    png(smearPng[i], height=600, width=600)
    plotTitle <- paste("Smear Plot:", gsub(".", " ", contrastData[i], 
                       fixed=TRUE))
    plotSmear(testData, de.tags=topIDs, pch=20, cex=0.8, main=plotTitle)
    abline(h=c(-1, 0, 1), col=c("dodgerblue", "yellow", "dodgerblue"), lty=2)
    
    imgName <- paste0("Smear Plot(", contrastData[i], ")")
    imgAddr <- paste0("smear(", contrastData[i], ").png")
    imageData <- rbind(imageData, c(imgName, imgAddr))
    invisible(dev.off())
    
    pdf(smearPdf[i])
    plotTitle <- paste("Smear Plot:", gsub(".", " ", contrastData[i], 
                       fixed=TRUE))
    plotSmear(testData, de.tags=topIDs, pch=20, cex=0.8, main=plotTitle)
    abline(h=c(-1, 0, 1), col=c("dodgerblue", "yellow", "dodgerblue"), lty=2)
    
    linkName <- paste0("Smear Plot(", contrastData[i], ") (.pdf)")
    linkAddr <- paste0("smear(", contrastData[i], ").pdf")
    linkData <- rbind(linkData, c(linkName, linkAddr))
    invisible(dev.off())
    
    genes <- as.character(data$genes$Gene)
    unq <- unique(genes)
    unq <- unq[!is.na(unq)]
    geneList <- list()
    for (gene in unq) {
      if (length(which(genes==gene)) >= hairpinReq) {
        geneList[[gene]] <- which(genes==gene)
      }
    }
    
    if (wantRoast) {
      # Input preparaton for roast
      nrot = 9999
      set.seed(602214129)
      roastData <- mroast(data, index=geneList, design=design,
                         contrast=contrasts, nrot=nrot)
      roastData <- cbind(GeneID=rownames(roastData), roastData)
      write.table(roastData, file=roastOut[i], row.names=FALSE, sep="\t")
      linkName <- paste0("Gene Level Analysis Table(", contrastData[i], 
                         ") (.tsv)")
      linkAddr <- paste0("roast(", contrastData[i], ").tsv")
      linkData <- rbind(linkData, c(linkName, linkAddr))
      if (selectOpt=="rank") {
        selectedGenes <- rownames(roastData)[selectVals]
      } else {
        selectedGenes <- selectVals
      }
      
      if (packageVersion("limma")<"3.19.19") {
        png(barcodePng[i], width=600, height=length(selectedGenes)*150)
      } else {
        png(barcodePng[i], width=600, height=length(selectedGenes)*300)
      }
      par(mfrow=c(length(selectedGenes), 1))
      for (gene in selectedGenes) {
        barcodeplot(testData$table$logFC, index=geneList[[gene]],
                    main=paste("Barcode Plot for", gene, "(logFCs)", 
                               gsub(".", " ", contrastData[i])),
                    labels=c("Positive logFC", "Negative logFC"))
      }
      imgName <- paste0("Barcode Plot(", contrastData[i], ")")
      imgAddr <- paste0("barcode(", contrastData[i], ").png")
      imageData <- rbind(imageData, c(imgName, imgAddr))
      dev.off()
      if (packageVersion("limma")<"3.19.19") {
        pdf(barcodePdf[i], width=8, height=2)
      } else {
        pdf(barcodePdf[i], width=8, height=4)
      }
      for (gene in selectedGenes) {
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

# Record ending time and calculate total run time
timeEnd <- as.character(Sys.time())
timeTaken <- capture.output(round(difftime(timeEnd,timeStart), digits=3))
timeTaken <- gsub("Time difference of ", "", timeTaken, fixed=TRUE)
################################################################################
### HTML Generation
################################################################################
# Clear file
cat("", file=htmlPath)

cata("<html>\n")
HtmlHead("EdgeR Output")

cata("<body>\n")
cata("<h3>EdgeR Analysis Output:</h3>\n")
cata("<h4>Input Summary:</h4>\n")
if (inputType=="fastq") {
  cata("<ul>\n")
  ListItem(hpReadout[1])
  ListItem(hpReadout[2])
  cata("</ul>\n")
  cata(hpReadout[3], "<br />\n")
  cata("<ul>\n")
  ListItem(hpReadout[4])
  ListItem(hpReadout[7])
  cata("</ul>\n")
  cata(hpReadout[8:11], sep="<br />\n")
  cata("<br />\n")
  cata("<b>Please check that read percentages are consistent with ")
  cata("expectations.</b><br >\n")
} else if (inputType=="counts") {
  cata("<ul>\n")
  ListItem("Number of Samples: ", ncol(data$counts))
  ListItem("Number of Hairpins: ", countsRows)
  ListItem("Number of annotations provided: ", annoRows)
  ListItem("Number of annotations matched to hairpin: ", annoMatched)
  cata("</ul>\n")
}

cata("The estimated common biological coefficient of variation (BCV) is: ", 
     commonBCV, "<br />\n")

cata("<h4>Output:</h4>\n")
cata("All images displayed have PDF copy at the bottom of the page, these can ")
cata("exported in a pdf viewer to high resolution image format. <br />\n")
for (i in 1:nrow(imageData)) {
  if (grepl("barcode", imageData$Link[i])) {
    if (packageVersion("limma")<"3.19.19") {
      HtmlImage(imageData$Link[i], imageData$Label[i], 
                height=length(selectedGenes)*150)
    } else {
      HtmlImage(imageData$Link[i], imageData$Label[i], 
                height=length(selectedGenes)*300)
    }
  } else {
    HtmlImage(imageData$Link[i], imageData$Label[i])
  }
}
cata("<br />\n")

cata("<h4>Plots:</h4>\n")
for (i in 1:nrow(linkData)) {
  if (!grepl(".tsv", linkData$Link[i])) {
    HtmlLink(linkData$Link[i], linkData$Label[i])
  }
}

cata("<h4>Tables:</h4>\n")
for (i in 1:nrow(linkData)) {
  if (grepl(".tsv", linkData$Link[i])) {
    HtmlLink(linkData$Link[i], linkData$Label[i])
  }
}

cata("<p>alt-click any of the links to download the file, or click the name ")
cata("of this task in the galaxy history panel and click on the floppy ")
cata("disk icon to download all files in a zip archive.</p>\n")
cata("<p>.tsv files are tab seperated files that can be viewed using Excel ")
cata("or other spreadsheet programs</p>\n")

cata("<h4>Additional Information:</h4>\n")

if (inputType == "fastq") {
  ListItem("Data was gathered from fastq raw read file(s).")
} else if (inputType == "counts") {
  ListItem("Data was gathered from a table of counts.")
}

if (cpmReq!=0 && sampleReq!=0) {
  tempStr <- paste("Hairpins that do not have more than", cpmReq,
                   "CPM in at least", sampleReq, "samples are considered",
                   "insignificant and filtered out.")
  ListItem(tempStr)
  filterProp <- round(filteredCount/preFilterCount*100, digits=2)
  tempStr <- paste0(filteredCount, " of ", preFilterCount," (", filterProp,
                   "%) hairpins were filtered out for low count-per-million.")
  ListItem(tempStr)
}

if (workMode == "classic") {
  ListItem("An exact test was performed on each hairpin.")
} else if (workMode == "glm") {
  ListItem("A generalised linear model was fitted to each hairpin.")
}



cit <- character()
link <-character()
link[1] <- paste0("<a href=\"",
                  "http://www.bioconductor.org/packages/release/bioc/",
                  "vignettes/limma/inst/doc/usersguide.pdf",
                  "\">", "limma User's Guide", "</a>.")
link[2] <- paste0("<a href=\"",
                  "http://www.bioconductor.org/packages/release/bioc/",
                  "vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf",
                  "\">", "edgeR User's Guide", "</a>")
                  
cit[1] <- paste("Robinson MD, McCarthy DJ and Smyth GK (2010).",
                "edgeR: a Bioconductor package for differential",
                "expression analysis of digital gene expression",
                "data. Bioinformatics 26, 139-140")
cit[2] <- paste("Robinson MD and Smyth GK (2007). Moderated statistical tests",
                "for assessing differences in tag abundance. Bioinformatics",
                "23, 2881-2887")
cit[3] <- paste("Robinson MD and Smyth GK (2008). Small-sample estimation of",
                "negative binomial dispersion, with applications to SAGE data.",
                "Biostatistics, 9, 321-332")

cit[4] <- paste("McCarthy DJ, Chen Y and Smyth GK (2012). Differential",
                "expression analysis of multifactor RNA-Seq experiments with",
                "respect to biological variation. Nucleic Acids Research 40,",
                "4288-4297")

cata("<h4>Citations</h4>")
cata("<ol>\n")
ListItem(cit[1])
ListItem(cit[2])
ListItem(cit[3])
ListItem(cit[4])
cata("</ol>\n")

cata("<table border=\"0\">\n")
cata("<tr>\n")
TableItem("Task started at:"); TableItem(timeStart)
cata("</tr>\n")
cata("<tr>\n")
TableItem("Task ended at:"); TableItem(timeEnd)
cata("</tr>\n")
cata("<tr>\n")
TableItem("Task run time:"); TableItem(timeTaken)
cata("<tr>\n")
cata("</table>\n")

cata("</body>\n")
cata("</html>")
