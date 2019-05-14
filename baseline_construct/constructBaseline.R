suppressPackageStartupMessages(library(optparse))
library(tools)
option_list <- list(
    make_option(c("--panelName"), type = "character",  help = "Panel Name"),
    make_option(c("--normalDir"), type = "character",  help = "Path ot normal/MSS samples"),
    make_option(c("--msihDir"), type = "character",  help = "Path to MSI-H samples")
    make_option(c("--negplaDir"), type = "character",  help = "Path to MSS plasma samples"),
    make_option(c("--tumorPercentageFile"), type = "character",  help = "table include mss tumor,pair and tumor percentage"),
    make_option(c("--targetDepth"), type = "numeric",  help = "Target sequencing depth"),
    make_option(c("--minSupportReads"), type = "numeric",  help = "Minimum number of reads supporting ctDNA")
)

## get command line options, if help option encountered print help and exit,
## otherwise if options not found on command line then set defaults,
parser <- OptionParser(option_list=option_list)
argsL <- parse_args(parser)

## print some progress messages to stderr if "quietly" wasn't requested
if(! (file.exists(argsL$normalDir) & (file.exists(argsL$tumorDir))) ){
    print_help(parser)
    stop("Parameters are not enough!")
}
normalDir <- argsL$normalDir
tumorDir <- argsL$msihDir
negplaDir <- args$negplaDir
targetDepth <- round(as.numeric(args$targetDepth))
minSupportReads <- round(as.numeric(args$minSupportReads))
tumorPercent <- args$tumorPercentageFile
########################################################
##                    Basic functions                 ##
########################################################
#' Read coverage file "_dis" of MSIsensor
#'
#' The format of msisensor consist of 100 row of repeat length 1 to 100, each with the reads of tumor and normal
#' T: loci coverage of tumor; N: loci coverage to normal
#'
#' @param coverageFile read coverage file "_dis" by msisensor
#' @return list a list of loci coverage, each is a data frame of 100 X 3 with column names "Length/CoverageT/CoverageN"
loadData <- function(coverageFile){
    coverageLine <- readLines(coverageFile)
    if( length(coverageLine)%%3 !=0 ){
        stop(paste0(coverageFile, ": error format for unpaired mode"))
    }
    msiData <- split(coverageLine, rep(1:(length(coverageLine)/3),each=3))
    siteCoverageList <- list()
    for(m in msiData){
        T <- m[[3]]
        N <- m[[2]]
        dat <- data.frame(Length = 1:100,
                          CoverageT = as.numeric(unlist(strsplit(T," "))[-1]),
                          CoverageN = as.numeric(unlist(strsplit(N," "))[-1]))
        siteCoverageList[[m[[1]]]] <- dat
    }
    return(siteCoverageList)
}


########################################################
##           Step I Read coverage of normals           ##
########################################################
#' @param normalDir Read coverage of normals
#' @return datN a list of read coverage data frames of normal/MSS samples
#' 
lstfiles <- list.files(normalDir, pattern = "_allSite_msi_dis",recursive = TRUE,
           full.names = TRUE)
lstDat <- lapply(lstfiles, loadData)
names(lstDat) <- sapply(strsplit(basename(lstfiles),"_"),function(x) x[1])
normalData <- lstDat[[1]]
siteList <- names(lstDat[[1]])
siteLen <- as.numeric(sapply(strsplit(sapply(strsplit(siteList," "),function(x) x[4]),"[[]"),function(x) x[1]))
names(siteLen) <- siteList
siteList <- siteList[order(siteLen,decreasing = T)]
datN <- list()
for(s in siteList){
    datN[[s]] <- sapply(lstDat,function(x) x[[s]][,3])
}


########################################################
##           Step II Read covrage of 100% MSI-H       ##
########################################################
#' @param tumorDir  Read coverage of tumor and paired samples 
#' @param tumorPercent a table of tumor pair tumor percentage
#' @return lstSPeakDat
#' @return lstSPeakRegion
#' @return lstSReads
#' @return lstStumorPercent
#' 
lstfiles <- list.files(tumorDir, pattern = "_allSite_msi_dis",recursive = TRUE,
           full.names = TRUE)
lstDat <- lapply(lstfiles, loadData)
names(lstDat) <- sapply(strsplit(basename(lstfiles),"_"),function(x) x[1])
tpInfo <- read.delim(tumorPercent, stringsAsFactors = F, header = TRUE)
tp <- as.numeric(tpInfo$tumorPercent)
if(max(tp,na.rm = T)>1){
    tp <- tp/100
}
names(tp) <- tpInfo$tumor
lstDat <- lstDat[names(lstDat) %in% tpInfo$tumor] # Select coverage data with known tumor percentage
sampleList <- names(lstDat)
lstSPeakDat <- list()
lstSReads <- list()
lstSPeakRegion <- list()
for(sp in sampleList){
    lstPeakDat <- list()
    lstReads <- list()
    lstPeakRegion <- list()
    for(s in siteList){
        dT <- lstDat[[sp]][[s]][,2]
        dN <- lstDat[[sp]][[s]][,3]
        if(max(dN) < 15){
            next
        }
        dis <- 0:100
        i <- which.max(dN)
        index <- dN[i]
        peakStart <- i
        peakEnd <- i
        cutoff <- sum(dN)*0.75
        while(index < cutoff){
        if(dN[peakStart-1]<dN[peakEnd+1]){
            peakEnd <- peakEnd + 1
        }else{
            peakStart <- peakStart - 1
        }
        index <- sum(dN[peakStart:peakEnd])
        if(peakEnd - peakStart + 1 == 5){
            break
        }
        }
        lstPeakDat[[s]] <- cbind(T=dT[peakStart : peakEnd] , N=dN[peakStart : peakEnd])
        lstPeakRegion[[s]] <- peakStart:peakEnd
        lstReads[[s]] <- c(N=sum(dN),T=sum(dT))
    }
    lstSPeakDat[[sp]] <- lstPeakDat
    lstSPeakRegion[[sp]] <- lstPeakRegion
    lstSReads[[sp]] <- lstReads
}
datT <- list()
for( s in siteList){
    dat <- do.call(cbind,lapply(sampleList,function(sp){
        dT <- lstDat[[sp]][[s]][,2]
        dN <- lstDat[[sp]][[s]][,3]
        if(max(dN) < 15){
            return(rep(0,100))
        }
        return(dT - dN * lstSReads[[sp]][[s]]["T"] / lstSReads[[sp]][[s]]["N"] * tp[sp])
    }))
    dat[dat<0] <- 0
    datT[[s]] <- dat
}

for(s in siteList){
    nCount <- sum(datN[[s]])
    tCount <- sum(datT[[s]])
    datT[[s]] <- round( datT[[s]] * nCount / tCount )
}

########################################################
##           Step III Search marker loci              ##
########################################################
## Determine the probability of reads from MSS samples show MSI-H pattern 
n <- targetDepth
for( p in seq(1e-4,1e-2,1e-4)){
    x <- (1-p)^n 
    for(i in 1:floor(minSupportReads/2)){  ## MSS with reads located in MSI-H pattern should be less than minSupportReads/2
        x <- x + choose(n,i)*(1-p)^(n-i)*p^i
    }
    if(x < 0.99){
        p <- p - 1e-4
        break
    }
    message(paste(p,x))
}
## Determine the lowest probability of reads from MSI-H samples show MSI-H pattern 
n2 <- minSupportReads
for( p2 in seq(0,1,0.01)){
    x <- 0
    for(i in ceiling(minSupportReads/2):minSupportReads){
        x <- x + choose(n2,i)*(1-p2)^(n2-i)*p2^i
    }
    if( x > 0.9 ){
        break
    }
    message(paste(p2,x))
}
## Select marker loci 
library(R.utils)
lstStat <- list()
for(s in siteList){
    cID <- which(apply(datN[[s]],1,sum) / sum(datN[[s]]) < p)
    if(length(cID) <= 0 ){
        next
    }
    x1 <- seqToIntervals(cID)[1,2]
    x2 <-sum(datT[[s]][1:x1,])/sum(datT[[s]])
    peakRegion <- names(which(table(unlist(sapply(lstSPeakRegion,function(x) x[[s]]))) / length(lstSPeakRegion) > 0.25))
    if(length(peakRegion)==0){
        next
    }
    x3 <- sum(datT[[s]][1:x1,]) / sum(datT[[s]][c(1:x1,as.numeric(peakRegion)),])
    x4 <- sum(datN[[s]][1:x1,]) / sum(datN[[s]][c(1:x1,as.numeric(peakRegion)),])
    y1 <- sum(datT[[s]][1:x1,])
    y2 <-  sum(datN[[s]][1:x1,])
    y3 <- sum(datN[[s]])
    y4 <- sum(datN[[s]][as.numeric(peakRegion),])
    lstStat[[s]] <- c(x1,x2,x3,x4,y1,y2,y3,y4)
}
markerTab <- do.call(rbind,lstStat)
markerTab <- markerTab[which(markerTab[,3] > p2),,drop = FALSE]
markerTab <- data.frame(m[order(markerTab[,2],decreasing = T),])
rownames(markerTab) <- markerTab[,1]
if(nrow(markerTab)<4){
    stop("Not enough markers have been found! Please modify the parameters and do it again. The failure may be caused by  the low target depth of sequencing, or too much reads required. Please use deep sequencing data or use reduce the minSupportReads parameter. Note that too small minSupportReads may also introduce noise to the model")
}

## Visualize markers 
library(pheatmap)
library(grid)
i <- 0
pdf(file.path(outdir, "markerDistribution.png"))
for(s in rownames(markerTab)){
    message(s)
    i <- i + 1
    mt <- datT[[s]]
    mt <- mt[,hclust(dist(t(mt)))$order]
    mn <- datN[[s]]
    mn <- mn[,hclust(dist(t(mn)))$order]
    m <- data.frame(mt,mn)
    colnames(m) <- c(rep("MSI-H-tumor-cells",ncol(mt)),
                     rep("MSS-normal-cells", ncol(mn)))
    colnames(m)[1:ncol(mt)][(1:ncol(mt))%%4 !=1] <- " "
    colnames(m)[(ncol(mt)+1):(ncol(mt)+ncol(mn))][(1:ncol(mn))%%4 !=1] <- " " 
    rownames(m) <- 1:100
    m <- log2(m)
    m[is.infinite(as.matrix(m))] <- 0
    vec <- unlist(m)[unlist(m)!=0]
    m[m>quantile(vec,0.9)] <- quantile(vec,0.9)
    ylim <- ceiling(as.numeric(strsplit(markerTab[s, "mssPattern"],",")[[1]][2])/10)*10+10
    pheatmap(m[1:ylim,],cluster_row = F, cluster_col = F,fontsize=10, main = s )
}
dev.off()

##################################################
## Run baseline of new markers for negative plasma cases
## to include mean and sd of H values to baseline
## 

lstfiles <- list.files(negPlaDir, pattern = "_allSite_msi_dis",recursive = TRUE,
           full.names = TRUE)
siteCoverageData <- lapply(lstfiles, loadData)
m <- markerTab
lstT <- list()
for(s in m$Site){
    dT <- siteCoverageData[[s]][,2]
    peakRegion <- eval(parse(text = paste0("seq(",m[s,"mssPattern"],")")))
    varRegion <- as.character(1:m[s,"msiPattern"])
    N <- m[s,"N"]
    K <- m[s,"K"]
    q <- sum(dT[as.numeric(varRegion)])
        q <- ifelse(q>K, K, q)
    tt <- sum(dT[as.numeric(peakRegion)])
    k <- q + tt
    p <- phyper(q,K,N-K,k,lower.tail = FALSE)
    if(-log(phyper(0,K,N-K,k,lower.tail = FALSE)) > 1){
        p <- NA
    }
    lstT[[s]] <- c(p,q,k)
}
    
parTable <- do.call(rbind, lstT)
colnames(parTable) <- c("pvalue","varReads","totalReads")
pVal <- sapply(lstT, function(y) y[1])
hVal <- -log(pVal)  #mean(x[,negID])+3*sd(x[,negID])
hVal[hVal> 100] <- 100
hVal[is.infinite(hVal)] <- 100

markerTab$meanH <- apply(hVal,1,function(x) mean(x[x > 0.1],na.rm = TRUE))
markerTab$sdH <- apply(hVal,1,function(x) sd(x[x>0.1], na.rm = TRUE))
baseline <- list(markerTable = markerTab, markerSite = markerTab$Site, rawCoverageData= list(normalData=normalData[markerTab$Site]))
save(baseline, file = file.path(outdir, paste0(panelName,".bmsisea.baseline.rdata"))
