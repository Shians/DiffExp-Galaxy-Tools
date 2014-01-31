################################################################################
### Input Processing
################################################################################
# Record starting time
timeStart <- as.character(Sys.time())

# Load all required libraries
library("Rsubread", quietly=TRUE, warn.conflicts=FALSE)
library("org.Mm.eg.db", quietly=TRUE, warn.conflicts=FALSE)
library("org.Hs.eg.db", quietly=TRUE, warn.conflicts=FALSE)

# Grabbing arguments from command line
argv <- commandArgs(TRUE)

# Grab arguments
fastqPath <- as.character(gsub("fastq::", "", argv[grepl("fastq::", argv)], 
                          fixed=TRUE))
argv <- argv[!grepl("fastq::", argv)]
workMode <- as.character(argv[1])
inputFormat <- as.character(argv[2])
refIndex <- as.character(argv[3])
annoOpt <- as.character(argv[4])
pairMode <- as.character(argv[5])
nThreads <- as.numeric(argv[6])
countsOut <- as.character(argv[7])
geneannoOut <- as.character(argv[8])

# Process arguments
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

if (pairMode=="single"){
  fastqData <- data.frame(Name=nameVec, Group=groupVec, Path=pathVec1)
} else if (pairMode=="paired"){
  fastqData <- data.frame(Name=nameVec, Group=groupVec, Path1=pathVec1,
                          Path2=pathVec2)
}

filesOut <- character()
for (i in 1:nrow(fastqData)){
  filesOut[i] <- paste(fastqData$Name[i], ".bam")
}

################################################################################
### Data Processing
################################################################################

buildindex(basename="index", reference=refIndex)

if (workMode == "subjunc"){
  for (i in 1:nrow(fastqData)){
    if (pairMode=="single"){
      subjunc(index="index", readfile1=fastqData$Path[i], 
              output_file=filesOut[i], output_format="BAM", nthreads=nThreads)
    } else if (pairMode=="paired"){
      subjunc(index="index", readfile1=fastqData$Path1[i],
              readfile2=fastqData$Path2[i], output_file=filesOut[i], 
              output_format="BAM", nthreads=nThreads)
    }
  }
} else if (workMode=="align"){
   for (i in 1:nrow(fastqData)){
    if (pairMode=="single"){
      align(index="index", readfile1=fastqData$Path[i], 
            output_file=filesOut[i], output_format="BAM", nthreads=nThreads)
    } else if (pairMode=="paired"){
      align(index="index", readfile1=fastqData$Path1[i],
            readfile2=fastqData$Path2[i], output_file=filesOut[i], 
            output_format="BAM", nthreads=nThreads)
    }
  }
}

if (pairMode=="paired"){
  counts <- featureCounts(filesOut, annot.inbuilt=annoOpt, isPairedEnd=TRUE)
} else if (pairMode=="single"){
  counts <- featureCounts(filesOut, annot.inbuilt=annoOpt, isPairedEnd=FALSE)
}

write.table(counts$counts, file=countsOut, sep="\t")

if (annoOpt=="hg19"){
  symbs <- toTable(org.Hs.egSYMBOL)
  chr <- toTable(org.Hs.egCHR)
  geneNames <- toTable(org.Hs.egGENENAME)
} else if (annoOpt=="mm9" | annoOpt=="mm10"){
  geneNames <- toTable(org.Mm.egGENENAME)
  symbs <- toTable(org.Mm.egSYMBOL)
  chr <- toTable(org.Mm.egCHR)
}

egids <- counts$anno[,1]
geneanno <- cbind(counts$anno, 
                  "EntrezID"=egids,
                  "Symbols"=symbs[match(egids, symbs[,1]),2],
                  "GeneName"=geneNames[match(egids, geneNames[,1]),2],
                  "Chr"=chr[match(egids, chr[,1]),2])

write.table(geneanno, file=geneannoOut, sep="\t")