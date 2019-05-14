# bMSISEA
bMSISEA：Detection of micro-satellite instability by signiture enrichment analysis from circulating tumor DNA by targeted deep sequencing

The R package bMSISEA detect the status of MSI from ctDNA samples. Studies show a high sensitivity of 93.4% with suffiecient ctDNA content (tumor allele frequency (maxAF) > 0.2%) and larger than 99% specificity. Copyright belongs to Guangzhou Burning Rock Biotech Co., China. This project is free for use by Academic users. Please see enclosed license agreement for details and conditions of use.

## Required Input files:
- bamFile : plasma sample of interest aligned against reference genome, provided in bam format. 
- msi_markers : MSI marker site file (see example under “inst/example/markerLoci.msi”) - specifies the selected marker loci. User generates this file with baselineConstruction function.
- msi_baseline : MSI baseline file (see example under “inst/example/baseline.rdata”) - describes statistics of each marker locus, as well as mean and standard deviation of H value of each locus, as calculated from an MSI negative population (white blood cell samples or MSI negative plasma). User generates this file with baselineConstruction function (see below).
- coverageFilePath: the file of coverage data, the ouput file of function coverageCaller and input file for function siteCoverageData
- msi_threshold : the cutoff of MS_score to determine MSI-H/MSS. 
- sample name: the name of the plasma sample
 NOTE: Baseline statistics vary markedly from assay-to-assay and lab-to-lab. It is critial that you prepare a baseline file that is specific for your analytic process, and for which data have been generated using the same protocols. Sensitivies is limited using low coverage data. Previous studies has shown a minimum of 5000X is required for targeted sequencing. 

## Output
a list of three element:
- msi_status character, MSI-H / MSS
- MS_score numeric, the ms score
- lociInfo data frame, the statistics of each loci
## Running Code
The MSI detection is fulfilled by three steps.
```R
library(bmsisea)
msi_markers <- system.file("example", "markerLoci.msi", package = "bmsisea")
msi_baseline <- system.file("example", "baseline.rdata", package = "bmsisea")
coverageFilePath <- system.file("example", "test_msi", package = "bmsisea")
msi_threshold <- 15
#Not run: the bamFile is not available
#coverageCaller(bamFile, msi_markers, coverageFilePath)
siteCoverageData <- loadData(paste0(coverageFilePath, "_dis"))
msiOut <- msiDetect(siteCoverageData, msi_baseline, "testSample", msi_threshold)
print(msiOut)
```

## Baseline construction 
scripts to construct baseline is located in inst/baseline_construct. 
1) Catalog all homo-polymers in your host genome. Several algorithms available to do this, but one we have incorporated in our package is the MSIsensor(PMID:24371154), which is very easy to use. Microsatellites with long repeat length are recommendated.
2) Limit the list of microsatellites to those presented in your capture design. Location of BED format is required.
3) Select microsatellite marker loci and construct baseline files. No fewer than 20 neat MSI-H & MSS cell lines or tumor samples with known tumor cell percentage and normal samples are required.
4) Perform bMSISAE analysis of many normal plasma samples to construct the baseline of H values and cutoff of MS scores.
5) Incorporate the results of 4) to baseline data.

You can run the script  msi_baseline.sh to construct the baseline.
```
Usage:  msi_baseline.sh 
-r <path to reference genome .fasta/.fa>
-s <path to baseline construction scripts>
-o <output path>
-b <custom_assay_bed>
-t <output path of msi-h tumor, paired normal coverage file of all covered sites>
-n <output path of normal/mss coverage file of all covered sites >
-p <output path of negative plasma coverage file of all covered sites>
-f <msih_tumor_normal_tumorPercent_file; Header: tumor\tnormal\ttumor_bam\tnormal_bam\tumorPercent>
-g <mss_normal_list; Header: ID\tbam>
-i <neg_pla_list; Header: ID\tbam>]
```
