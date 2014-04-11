################################################################################
### Input Processing
################################################################################
# Record starting time
timeStart <- as.character(Sys.time())

################################################################################
### Function Delcaration
################################################################################
# Function to load libaries without messages
silentLibrary <- function(...) {
  list <- c(...)
  for (package in list){
    suppressPackageStartupMessages(library(package, character.only=TRUE))
  }
}

# Function to add counts from multiple colums into a single column, removes
# all but the first specified column which holds the sum of all specified
# columns
addCols <- function(dataframe, indices){
  if (!is.vector(counts$counts[, indices])){
    dataframe[, indices[1]] <- rowSums(dataframe[, indices])
    dataframe <- dataframe[, -indices[2:length(indices)]]
  }
  return(dataframe)
}

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

# Generate output folder and paths
makeOut <- function(filename) {
  return(paste0(outPath, "/", filename))
}

# Generating design information
pasteListName <- function(string) {
  return(paste0("factors$",string))
}

# Create cata function: default path set, default seperator empty and appending
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

HtmlSpacing <- function(string) {
  string <- gsub(" ", "&nbsp;", string, fixed=TRUE)
  string <- gsub("\t", "&nbsp;&nbsp;&nbsp;&nbsp;", string, fixed=TRUE)
  return(string)
}

TrimAlignOutput <- function(alignOutput) {
  start <- which(grepl("=== Running", alignOutput, fixed=TRUE))
  end <- which(grepl("\\\\================", alignOutput, fixed=TRUE))+1
  ranges <- numeric()
  for (i in 1:length(start)) {
    ranges <- c(ranges, start[i]:end[end > start[i]][1])
  }
  banners <- which(grepl("==========   |_____/", alignOutput, fixed=TRUE))[-1]
  for (i in 1:length(banners)) {
    ranges <- c(ranges, (banners[i]-7):(banners[i]+1))
  }
  alignOutput <- alignOutput[-ranges]
  
  ranges <- numeric()
  summary <- which(grepl("== Summary ==", alignOutput, fixed=TRUE))
  for (i in 1:length(summary)) {
    ranges <- c(ranges, (summary[i]-3):(summary[i]))
  }
  alignOutput <- alignOutput[-ranges]
  return(alignOutput)
}

################################################################################
### Input Processing
################################################################################

# Grabbing arguments from command line
argv <- commandArgs(TRUE)

# Grab arguments
fastqPath <- as.character(gsub("fastq::", "", argv[grepl("fastq::", argv)], 
                          fixed=TRUE))
argv <- argv[!grepl("fastq::", argv)]
workMode <- as.character(argv[1])
inputFormat <- as.character(argv[2])
refIndexSource <- as.character(argv[3])
refIndex <- as.character(argv[4])
annoOpt <- as.character(argv[5])
pairMode <- as.character(argv[6])
nThreads <- as.numeric(argv[7])
countsOut <- as.character(argv[8])
geneannoOut <- as.character(argv[9])
htmlPath <- as.character(argv[10])
htmlDirPath <- as.character(argv[11])

# Collecting name, group and path values for each fastq input
fastqPath <- strsplit(fastqPath, "::")
nameVec <- character()
groupVec <- character()
pathVec1 <- character()
pathVec2 <- character()
for (i in 1:length(fastqPath)){
  nameVec <- c(nameVec, fastqPath[[i]][1])
  groupVec <- c(groupVec, fastqPath[[i]][2])
  pathVec1 <- c(pathVec1, fastqPath[[i]][3])
  if (pairMode=="paired"){
    pathVec2 <- c(pathVec2, fastqPath[[i]][4])
  }
}

# Create table of information regarding fastQ files
if (pairMode=="single"){
  fastqData <- data.frame(Name=nameVec, Group=groupVec, Path=pathVec1)
} else if (pairMode=="paired"){
  fastqData <- data.frame(Name=nameVec, Group=groupVec, Path1=pathVec1,
                          Path2=pathVec2)
}

# Create file names for intermediary bam files
filesOut <- character()
for (i in 1:nrow(fastqData)){
  filesOut[i] <- paste0(fastqData$Name[i], ".bam")
}

################################################################################
### Data Processing
################################################################################
# Load all required packags
silentLibrary("Rsubread", "org.Mm.eg.db", "org.Hs.eg.db")

# Build index from fasta file if required 
if (refIndexSource=="history"){
  buildOutput <- capture.output(
  buildindex(basename="index", reference=refIndex)
  )
  indexName <- "index"
} else if (refIndexSource=="indexed"){
  indexName <- refIndex
}

# *CURRENTLY UNUSED* Align reads using selected algorithm
alignOutput <- capture.output(
if (workMode == "subjunc"){
  for (i in 1:nrow(fastqData)){
    if (pairMode=="single"){
      subjunc(index=indexName, readfile1=fastqData$Path[i], 
              output_file=filesOut[i], output_format="BAM", nthreads=nThreads)
    } else if (pairMode=="paired"){
      subjunc(index=indexName, readfile1=fastqData$Path1[i],
              readfile2=fastqData$Path2[i], output_file=filesOut[i], 
              output_format="BAM", nthreads=nThreads)
    }
  }
} else if (workMode=="align"){
   for (i in 1:nrow(fastqData)){
    if (pairMode=="single"){
      align(index=indexName, readfile1=fastqData$Path[i], 
            output_file=filesOut[i], output_format="BAM", nthreads=nThreads)
    } else if (pairMode=="paired"){
      align(index=indexName, readfile1=fastqData$Path1[i],
            readfile2=fastqData$Path2[i], output_file=filesOut[i], 
            output_format="BAM", nthreads=nThreads)
    }
  }
}
)

featureOutput <- capture.output(
# Calculated counts from aligned read data
if (pairMode=="paired"){
  counts <- featureCounts(filesOut, annot.inbuilt=annoOpt, isPairedEnd=TRUE)
} else if (pairMode=="single"){
  counts <- featureCounts(filesOut, annot.inbuilt=annoOpt, isPairedEnd=FALSE)
}
)

# Merge columns by sample label, counts belonging to same sample are merged into
# a single column.
uniqueSamples <- unique(nameVec)
for (i in 1:length(uniqueSamples)){
  indices <- which(nameVec == uniqueSamples[i])
  counts$counts <- addCols(counts$counts, indices)
  if (length(indices)>1){
    nameVec <- nameVec[-indices[2:length(indices)]]
  }
}

write.table(counts$counts, file=countsOut, sep="\t")

# Load appropriate annotation database
if (annoOpt=="hg19"){
  geneNames <- toTable(org.Hs.egGENENAME)
  symbs <- toTable(org.Hs.egSYMBOL)
  chr <- toTable(org.Hs.egCHR)
} else if (annoOpt=="mm9" | annoOpt=="mm10"){
  geneNames <- toTable(org.Mm.egGENENAME)
  symbs <- toTable(org.Mm.egSYMBOL)
  chr <- toTable(org.Mm.egCHR)
}

# Create annotation table
egids <- counts$anno[,1]
geneanno <- cbind(counts$anno, 
                  "EntrezID"=egids,
                  "Symbols"=symbs[match(egids, symbs[,1]),2],
                  "GeneName"=geneNames[match(egids, geneNames[,1]),2],
                  "Chr"=chr[match(egids, chr[,1]),2])

write.table(geneanno, file=geneannoOut, sep="\t")

if (exists("buildOutput")) {
  buildOutput <- sapply(buildOutput, HtmlSpacing)
}
alignOutput <- TrimAlignOutput(alignOutput)
alignOutput <- sapply(alignOutput, HtmlSpacing)

featureOutput <- sapply(featureOutput, HtmlSpacing)

################################################################################
### HTML Generation
################################################################################

# Clear file
cat("", file=htmlPath)

cata("<html>\n")
cata("<head>\n")
cata("<title>", "RSubread Output", "</title>\n")
cata("</head>\n")
cata("<body>\n")
cata("<font face=\"Courier New\">\n")
if (exists("buildOutput")) {
  cata(buildOutput, sep="<br />\n")
}
cata(alignOutput, sep="<br />\n")
cata("</font>\n")
cata("</body>\n")

cata("</html>")