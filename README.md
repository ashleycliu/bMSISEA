# bMSISEA
bMSISEA：Detection of micro-satellite instability by signiture enrichment analysis from circulating tumor DNA by targeted deep sequencing

The R package bMSISEA detect the status of MSI from ctDNA samples. Studies show a high sensitivity of 93.4% with suffiecient ctDNA content (tumor allele frequency (maxAF) > 0.2%) and larger than 99% specificity. Copyright belongs to Guangzhou Burning Rock Biotech Co., China. This project is free for use by Academic users. Please see enclosed license agreement for details and conditions of use.

##Required Input files:
- bamFile : plasma sample of interest aligned against reference genome, provided in bam format. 
- msi_markers : MSI marker site file (see example under “inst/example/markerLoci.msi”) - specifies the selected marker loci. User generates this file with baselineConstruction function.
- msi_baseline : MSI baseline file (see example under “inst/example/baseline.rdata”) - describes statistics of each marker locus, as well as mean and standard deviation of H value of each locus, as calculated from an MSI negative population (white blood cell samples or MSI negative plasma). User generates this file with baselineConstruction function (see below). NOTE: Baseline statistics vary markedly from assay-to-assay and lab-to-lab. It is CRITICAL that you prepare a baseline file that is specific for your analytic process, and for which data have been generated using the same protocols. Sensitivies is limited using low coverage data. Previous studies has shown a minimum of 5000X is required for targeted sequencing.

##Other required parameters
- coverageFilePath: the file of coverage data, the ouput file of function coverageCaller and input file for function siteCoverageData - msi_threshold : the cutoff of MS_score to determine MSI-H/MSS. sample name: the name of the plasma sample

## Output
a list of three element:
- msi_status character, MSI-H / MSS
- MS_score numeric, the ms score
- lociInfo data frame, the statistics of each loci
- Detection Code

## Runing Code
The MSI detection is fulfilled by three steps.
>> library(bmsisea)
>> msi_markers <- system.file("example", "markerLoci.msi", package = "bmsisea")
>> msi_baseline <- system.file("example", "baseline.rdata", package = "bmsisea")
>> coverageFilePath <- system.file("example", "test_msi", package = "bmsisea")
>> msi_threshold <- 15
## Not run: the bamFile is not available
## coverageCaller(bamFile, msi_markers, coverageFilePath)
>> siteCoverageData <- loadData(paste0(coverageFilePath, "_dis"))
>> msiOut <- msiDetect(siteCoverageData, msi_baseline, "testSample", msi_threshold)
>> print(msiOut)
