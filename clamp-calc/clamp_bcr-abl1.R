
source("clamp.R")

#####################################################################################################
##
## Author: Adam Ameur
##
## Function: clamp_bcr-abl1
##
## Description: Performs analysis of resistance mutations in a list of BCR-ABL1 amplicon samples.
##              This method is used for clinical routine screening for resistance mutations in Chronic
##              Myeloid Leukemia.
##
#####################################################################################################


## This function takes a file with sample info and primer sequences as arguments, and runs the screenMutation script on all samples in file.
screenMutationsCML <- function(sampleInfoFile="bcr-abl1_sample_list_with_primers.txt", assay="bcr-abl1", targetInfoFile="target_info.txt", analyzeIsoforms=TRUE, rerunAnalysis=FALSE, minReadlen=1000, maxReadlen=3000){

    ## Read all sample info
    data <- read.table(sampleInfoFile, sep="\t", header=TRUE)

    ## Read target sequence information
    targetData <- read.table(targetInfoFile, sep="\t", header=TRUE)

    ## Check that primer sequences are available for all samples!
    if(!all(as.character(data[,"primerId"]) %in% as.character(targetData[,"primerId"]))){
        stop(paste("Unknown primer sequences in file:",sampleInfoFile))
    }

    ## Loop through and launch CLAMP on each sample in file
    for(i in 1:nrow(data)){
        runId <- as.character(data[i,"runId"])
        primerId <- as.character(data[i,"primerId"])

        if(!(primerId %in% c("BCR-ABL1_Major","BCR-ABL1_Minor","BCR-ABL1_Micro","ABL1_E1","ABL1_F2"))){
            stop(paste("Wrong primer id for sample '",runId,"'. Must be either 'BCR-ABL1_Minor' or 'BCR-ABL1_Major'!",sep=""))
        }

        targetInfo <- getTargetInfo(primerId, assay=assay, targetInfoFile=targetInfoFile)

        screenMutations(runId, assay, targetInfo$mutTables, targetInfo$refFiles, targetInfo$primerFwd, targetInfo$primerRev, analyzeIsoforms=analyzeIsoforms, rerunAnalysis=rerunAnalysis, minReadlen=minReadlen, maxReadlen=maxReadlen, maxmmRound1=0, maxmmRound2=5)

    }

}
