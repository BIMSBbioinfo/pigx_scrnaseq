#!/usr/bin/env Rscript
# Author: EU
# Date: June, 2018
# This script takes as input an RDS file containing a SingleCellExperiment object 
# and runs the singleCellTK shiny app

# 1. Collect Arguments ----------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)

## Default setting when no arguments passed
if(length(args) == 0) {
    args = c("--help")
}
help_command = "
run_singleCellTK.R: Runs the singleCellTK shiny app using a SingleCellExperiment object

Arguments:
--sceRdsFile Path to the RDS format file containing the SingleCellExperiment object

Example:
Rscript renderReport.R --sceRdsFile=<path to sce.RDS> \n"


## Help section
if("--help" %in% args) {
    cat(help_command, "\n")
    q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs = function(x) {
    myArgs = unlist(strsplit(x, "--"))
    myArgs = myArgs[myArgs != '']
    #when no values are provided for the argument
    myArgs = gsub(pattern = "=$", replacement = "= ", x = myArgs)
    myArgs = as.data.frame(do.call(rbind, strsplit(myArgs, "=")))
    myArgs$V2 = gsub(' ', '', myArgs$V2)
    return(myArgs)
}

argsDF = parseArgs(args)
argsL = as.list(as.character(argsDF$V2))
names(argsL) = argsDF$V1

if(!("sceRdsFile" %in% argsDF$V1)) {
    cat(help_command, "\n")
    stop("Missing argument: sceRdsFile. Provide the path to .Rds file 
         containing the SingleCellExperiment object.")
}

sceRdsFile  = argsL$sceRdsFile

# 2. Install missing packages ---------------------------------------------
req_pckgs <- c("singleCellTK", "BiocGenerics", 
               "SummarizedExperiment", "HDF5Array", "rhdf5") ## MORE?

missing_pckgs <- c()
for (pckg_name in req_pckgs)
    if (!pckg_name %in% installed.packages())
        missing_pckgs <- c(missing_pckgs, pckg_name)

if (!is.null(missing_pckgs)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite(missing_pckgs)
}

# 3. Update the object ----------------------------------------------------
sce <- readRDS(sceRdsFile)
sce <- updateObject(sce, verbose=TRUE)
rownames(sce) <- rowData(sce)$geneName

library(singleCellTK)
singleCellTK(sce)
