
## TODO: More features/improvements to be added to the script:
##  1) Make references/primer sequences optional
##  2) Add automatic genotyping functionality
##  3) Create QC report
##  4) Allow for different mutation files for different references
##  5) De novo detection of indels

## Function for reading in all mutations in a text file and executing 'screen_mutations' in parallel
screenMutations <- function(sample, assay, mutationTableFiles, refFiles, fwdPrimer, revPrimer, minFreq=0.005, coverageCutoff=0.7, maxmmRound1=0, maxmmRound2=5, minq=0, minqn=0, minCoverage=100, minCoverageFraction=0.5, minFractionMuts=0.9, nrReadsForSNPcallingAndQC=500, minFreqDenovoSNPs=0.05, nrReadsForRefAssign=100, minReadlen=1, maxReadlen=1000000, windowSize=20, doGenotyping=FALSE, analyzeIsoforms=FALSE, rerunAnalysis=TRUE){

    cat(paste("Analyzing sample '",sample,"'",sep=""))

    ## Check that mutation tables and reference parameters are ok
    for(ref in refFiles){

        if(!file.exists(ref)){
            stop(paste("Reference file '",ref,"' missing!",sep=""))
        }

        refName <- sub(".fasta","",ref)
        refNameSplit <- strsplit(refName,"/")[[1]]
        refName <- refNameSplit[length(refNameSplit)]

        if(!(refName %in% names(mutationTableFiles))){
            stop(paste("Mutation table missing for '",refName,"'",sep=""))
        }

        if(!file.exists(mutationTableFiles[refName])){
            stop(paste("Mutation file '",mutationTableFiles[refName],"' missing!",sep=""))
        }

    }

    datadir <- paste("pacdata/",assay,"/",sample,"/",sep="")

    ## Input data: A fastq file containg reads
    ccsFASTQfile <- paste(datadir,sample,"_reads_of_insert.fastq",sep="")
    ccsFASTQfileGZ <- paste(datadir,sample,"_reads_of_insert.fastq.gz",sep="")
    ccsFilteredFASTQfile <- paste(datadir,sample,"_reads_of_insert_filtered.fastq",sep="")
    identicalReadsFile <-  paste(datadir,sample,"_identical_reads.fasta",sep="")
    candidateIsoformFile <-  paste(datadir,sample,"_candidate_isoforms.txt",sep="")
    logFile <-  paste(datadir,sample,"_analysis_log.txt",sep="")

    ##############################################################
    ##
    ## Output: All of the files below are created during analysis
    ##
    ##############################################################

    referenceFileList <- list()
    ccsFilteredFASTQfiles <- list()
    SNPcallingFiles <- list()
    SNPcallingFilesRefAssign <- list()
    deNovoVariantFiles <- list()
    lowCoverageFiles <- list()
    QCfiles <- list()
    qcFailedFiles <- list()
    mutationFilesRaw <- list()
    mutationFilesFinal <- list()
    mutationClonalityFiles <- list()
    mutationClonalityPlotFiles <- list()

    ## create one result directory for each reference sequence
    for(refFile in refFiles){
        if(length(grep(".fasta",refFile)) == 0){
            stop("Error. all reference file names must have suffix '.fasta'")
        }
        if(length(grep("_",refFile)) > 0){
            stop("Error. Reference files are not allowed to contain '_' characters")
        }
        if(!file.exists(refFile)){
            stop("Reference file '",refFile,"' does not exist!")
        }
        refFileSplit <- strsplit(refFile,"/")[[1]]
        refFileTmpName <- refFileSplit[length(refFileSplit)]
        refName <- sub(".fasta","",refFileTmpName)
        refDir <- paste(datadir,"/",refName,sep="")
        referenceFileList[[refName]] <- refFile
        ccsFilteredFASTQfiles[[refName]] <- paste(datadir,sample,"_",refName,"_reads_of_insert_filtered.fastq",sep="")
        SNPcallingFiles[[refName]] <- paste(datadir,sample,"_",refName,"_snpcalls.txt",sep="")
        SNPcallingFilesRefAssign[[refName]] <- paste(datadir,sample,"_",refName,"_snpcalls_refassign.txt",sep="")
        deNovoVariantFiles[[refName]] <- paste(datadir,sample,"_",refName,"_denovo_snps.txt",sep="")
        lowCoverageFiles[[refName]] <- paste(datadir,sample,"_",refName,"_low_coverage.txt",sep="")
        QCfiles[[refName]] <- paste(datadir,sample,"_",refName,"_QC.pdf",sep="")
        qcFailedFiles[[refName]] <- paste(datadir,"QC_failed_",sample,"_",refName,"_mutations_final.txt",sep="")
        mutationFilesRaw[[refName]] <- paste(datadir,sample,"_",refName,"_mutations_raw.txt",sep="")
        mutationFilesFinal[[refName]] <- paste(datadir,sample,"_",refName,"_mutations_final.txt",sep="")
        mutationClonalityFiles[[refName]] <- paste(datadir,sample,"_",refName,"_clonal_distribution.txt",sep="")
        mutationClonalityPlotFiles[[refName]] <- paste(datadir,sample,"_",refName,"_clonal_distribution.pdf",sep="")
    }

    refNames <- names(referenceFileList)

    ## #####################################################################
    ##
    ## Option 1: Keep old results if available, only generate missing files
    ##
    ## #####################################################################

    ## Check if all output files listed in the logfile exists. Otherwise rerun entire analysis pipeline.
    if(rerunAnalysis == FALSE){

        analysisComplete <- FALSE
        allFilesOk <- TRUE

        ## Check all files in logfile
        if(file.exists(logFile)){
            lines <- readLines(logFile)
            for(line in lines){
                if(length(grep("File created: ",line)) == 1){
                    fname <- as.character(strsplit(line, "File created: ")[[1]][2])
                    if(!file.exists(fname)){
                        allFilesOk <- FALSE
                    }
                }
                if(length(grep("CLAMP analysis complete.",line)) == 1){
                    analysisComplete <- TRUE
                }
            }
        }
        else{
            allFilesOk <- FALSE
        }

        if(allFilesOk && analysisComplete){
            cat(" - skipping [result files available] \n")
            return()
        }
    }


    ## #####################################################################
    ##
    ## Option 2: Delete all old results and rerun whole analysis pipeline
    ##
    ## #####################################################################


    ## ##########################################################################################
    ##
    ## 1. Clean up old results, parse mutation table and unpack read file (if a gzip file exists)
    ##
    ## ##########################################################################################

    cat(date()," - CLAMP analysis of sample '",sample,"' started. \n",sep="", file=logFile, append=TRUE)
    cat("\n")

    ## Clean up files in result directory...
    if(file.exists(ccsFilteredFASTQfile)){
        file.remove(ccsFilteredFASTQfile)
    }
    if(file.exists(identicalReadsFile)){
        file.remove(identicalReadsFile)
    }
    if(file.exists(candidateIsoformFile)){
        file.remove(candidateIsoformFile)
    }
    if(file.exists(logFile)){
        file.remove(logFile)
    }
    for(refName in refNames){
        if(file.exists(ccsFilteredFASTQfiles[[refName]])){
            file.remove(ccsFilteredFASTQfiles[[refName]])
        }
        if(file.exists(SNPcallingFiles[[refName]])){
            file.remove(SNPcallingFiles[[refName]])
        }
        if(file.exists(SNPcallingFilesRefAssign[[refName]])){
            file.remove(SNPcallingFilesRefAssign[[refName]])
        }
        if(file.exists(deNovoVariantFiles[[refName]])){
            file.remove(deNovoVariantFiles[[refName]])
        }
        if(file.exists(lowCoverageFiles[[refName]])){
            file.remove(lowCoverageFiles[[refName]])
        }
        if(file.exists(QCfiles[[refName]])){
            file.remove(QCfiles[[refName]])
        }
        if(file.exists(mutationFilesRaw[[refName]])){
            file.remove(mutationFilesRaw[[refName]])
        }
        if(file.exists(mutationFilesFinal[[refName]])){
            file.remove(mutationFilesFinal[[refName]])
        }
        if(file.exists(mutationClonalityFiles[[refName]])){
            file.remove(mutationClonalityFiles[[refName]])
        }
        if(file.exists(mutationClonalityPlotFiles[[refName]])){
            file.remove(mutationClonalityPlotFiles[[refName]])
        }
    }

    for(mutationTableFile in mutationTableFiles){

        mutationData <- read.table(mutationTableFile, header=TRUE, as.is=TRUE)

        for(i in 1:nrow(mutationData)){

            mut <- as.character(mutationData[i,"mutation"])

            for(refName in refNames){
                tmpfile <- sub("_final.txt",paste("_",mut,".txt",sep=""),mutationFilesFinal[[refName]])
                if(file.exists(tmpfile)){
                    file.remove(tmpfile)
                }
            }
        }
    }

    if(file.exists(ccsFASTQfileGZ)){
        if(!file.exists(ccsFASTQfile)){
            cat(" - Uncompressing FASTQ file...")
            cmd <- paste("gunzip ",ccsFASTQfileGZ,sep="")
            system(cmd)
            cat(" done\n")
        }
    }

    if(!file.exists(ccsFASTQfile)){
        stop(paste("FASTQ file with raw reads missing:",ccsFASTQfile))
    }

    ## #######################################################################
    ##
	## Step 2. Filter fastq file using primer sequence information
    ##
    ## #######################################################################

    cat(date()," - Filtering reads based on primer sequences [fwd:",fwdPrimer,", rev:",revPrimer,", minl:",minReadlen,", maxl:",maxReadlen,"]\n",sep="", file=logFile, append=TRUE)
    cat(" - Creating filtered FASTQ file...", sep="")

    primerTable <- generatePrimerPairs(fwdPrimer, revPrimer)

    if(file.exists(ccsFilteredFASTQfile)){
        file.remove(ccsFilteredFASTQfile)
    }

    for(i in 1:nrow(primerTable)){
        cmd <- paste("./filter_primers.pl -f ",ccsFASTQfile," -t fastq -fwd ",primerTable[i,1]," -rev ",primerTable[i,2]," -minl ", minReadlen, " -maxl ", maxReadlen, " >> ",ccsFilteredFASTQfile,sep="")
        system(cmd)
    }

    cat(date()," - File created: ",ccsFilteredFASTQfile,"\n", sep="", file=logFile, append=TRUE)
    cat(" done\n")


    ## ########################################
    ##
    ## Step 3 - Optional. Run isoform analysis
    ##
    ## ########################################

    if(analyzeIsoforms){
        cat(" - Creating file with identical reads...")
        cmd <- paste("./count_indentical_reads.pl -f ",ccsFilteredFASTQfile," > ",identicalReadsFile,sep="")
        system(cmd)
        cat(" done\n")
        if(file.info(identicalReadsFile)$size > 0){
            cat(" - Extracting candidate isoforms...")
            extractIsoforms(identicalReadsFile, candidateIsoformFile)
        }
        cat(" done\n")
    }

    ## #######################################################################################################
    ##
    ## Step 4a. (if doGenotyping == TRUE) Group fastq reads by reference genotypes
    ##
    ## TODO: Implement a BLAST (or similar) to group reads into correct reference files!!!
    ##
    ## #######################################################################################################

    if(doGenotyping){
        cat(date()," - Analysis type: genotyping\n", sep="", file=logFile, append=TRUE) ## ADDED!!
        stop("TODO! Genotyping option not yet available!!!\n")
    }

    ## #################################################################################
    ##
    ## Step 4b. (if doGenotyping == FALSE) Assign an optimal reference for the sample
    ##
    ## #################################################################################

    ## If genotyping is FALSE, figure out which reference is best and should be used for all further analysis
    if(!doGenotyping){

        cat(date()," - Analysis type: no genotyping\n", sep="", file=logFile, append=TRUE) ## ADDED!!

        ## Assign read file to each individual reference
        for(refName in refNames){
            ccsFilteredFASTQfiles[[refName]] <- ccsFilteredFASTQfile
        }

        cat(date()," - Assigning optimal reference based on ",nrReadsForRefAssign," reads\n", sep="", file=logFile, append=TRUE)
        cat(" - Assigning optimal reference...")

        optimalRef <- refNames[1]

        if(length(refNames) > 1){

            refCovAvg <- array(NA,length(refNames))
            names(refCovAvg) <- refNames
            refCovZero <- array(NA,length(refNames))
            names(refCovZero) <- refNames

            ## Run de novo SNP calling against each reference on very small number of reads to see which reference is best
            for(refName in refNames){
                cmd <- paste("./cava_snpcaller.pl -f ",ccsFilteredFASTQfiles[[refName]]," -r ",referenceFileList[[refName]]," -minq ",minq," -minqn ",minqn," -max_mm ",maxmmRound1," -n ",nrReadsForRefAssign," -w ",windowSize," --silent > ",SNPcallingFilesRefAssign[[refName]],sep="")

                system(cmd)

                tmpData <- read.table(SNPcallingFilesRefAssign[[refName]], sep="\t", header=TRUE)
                avgCov <- mean(as.numeric(tmpData[,"coverage"]))
                zeroCov <- length(which(as.numeric(tmpData[,"coverage"]) == 0))
                
                refCovAvg[refName] <- avgCov
                refCovZero[refName] <- zeroCov
            }

            ## The optimal reference is the one with least number of non-covered (zero) bases.
            ## If many references ave same number of zero coverage bases, the pick the one with highest average coverage
            tmpRefs <- names(refCovZero[which(min(refCovZero) == refCovZero)])
            optimalRef <- names(sort(refCovAvg[tmpRefs], decreasing=TRUE))[1]
            
        }

        refNames <- optimalRef

        cat(date()," - Reference to be used for further analysis: '",refNames,"'\n", sep="", file=logFile, append=TRUE)
        cat("[",refNames,"] done\n",sep="")
    }


    ## ##################################################
    ##
    ## Step 5. Run de novo SNP calling against reference
    ##
    ## ##################################################

    for(refName in refNames){

        cat(date()," - Running de novo SNP calling based on ",nrReadsForSNPcallingAndQC," reads\n", sep="", file=logFile, append=TRUE)

        cmd <- paste("./cava_snpcaller.pl -f ",ccsFilteredFASTQfiles[[refName]]," -r ",referenceFileList[[refName]]," -minq ",minq," -minqn ",minqn," -max_mm ",maxmmRound2," -n ",nrReadsForSNPcallingAndQC," -w ",windowSize," > ",SNPcallingFiles[[refName]],sep="")

        system(cmd)

        cat(date()," - File created: ",ccsFilteredFASTQfiles[[refName]],"\n", sep="", file=logFile, append=TRUE)
    }


    ## ############################################################
    ##
    ## Step 6. Analyze SNP calling results (and perform inital QC)
    ##
    ## ############################################################

    for(refName in refNames){

        cat(date()," - Analyzing SNP calling results and creating coverage statistics\n", sep="", file=logFile, append=TRUE)

        SNPcallData <- read.table(SNPcallingFiles[[refName]], sep="\t", header=TRUE)

        coverageTotal <- SNPcallData[,"cov_fwd"]+SNPcallData[,"cov_rev"]
        names(coverageTotal) <- SNPcallData[,"pos"]

        minCov <- min(coverageTotal)

        pdf(QCfiles[[refName]], height=5, width=10)

        mutationTableFile <- mutationTableFiles[refName]
        mutationData <- read.table(mutationTableFile, header=TRUE, as.is=TRUE)

        plot(coverageTotal, type="h", bty="n", ylim=c(0,nrReadsForSNPcallingAndQC), xlab=paste("Position in ",refName,sep=""), ylab="Coverage", main=paste("Coverage of ",refName,", based on ",nrReadsForSNPcallingAndQC," reads", sep=""))

        abline(h=minCov, col="red", lwd=2)

        dev.off()

        ## Create file with novel mutations (>X% frequency both on fwd and rev)

        ## 1. Only consider positions with relatively high coverage
        tmpData <- SNPcallData[which(coverageTotal >= (median(coverageTotal)*minCoverageFraction)),]
        tmpData <- tmpData[which(tmpData[,"cov_fwd"] > 0),]
        tmpData <- tmpData[which(tmpData[,"cov_rev"] > 0),]

        ## 2. Make sure mutations are detected on both strands, analyze all possible alternative alleles
        res1 <- NULL
        alt1FreqsFwd <- tmpData[,"alt1_fwd"]/tmpData[,"cov_fwd"]
        alt1FreqsRev <- tmpData[,"alt1_rev"]/tmpData[,"cov_rev"]
        positiveAlt1 <- tmpData[intersect(which(alt1FreqsFwd>=minFreqDenovoSNPs),which(alt1FreqsRev>=minFreqDenovoSNPs)),,drop=FALSE]
        baseChangeAlt1 <- paste(positiveAlt1[,"ref"],"/",positiveAlt1[,"alt1"],sep="")

        if(nrow(positiveAlt1)>0){
            res1 <- cbind(positiveAlt1[,c("pos","context")],baseChangeAlt1,positiveAlt1[,c("ref","alt1","cov_fwd","ref_fwd","alt1_fwd","cov_rev","ref_rev","alt1_rev")])
            names(res1) <- c("pos","context","basechange","ref","alt","cov_fwd","ref_fwd","alt_fwd","cov_rev","ref_rev","alt_rev")
        }

        res2 <- NULL
        alt2FreqsFwd <- tmpData[,"alt2_fwd"]/tmpData[,"cov_fwd"]
        alt2FreqsRev <- tmpData[,"alt2_rev"]/tmpData[,"cov_rev"]
        positiveAlt2 <- tmpData[intersect(which(alt2FreqsFwd>=minFreqDenovoSNPs),which(alt2FreqsRev>=minFreqDenovoSNPs)),,drop=FALSE]
        baseChangeAlt2 <- paste(positiveAlt2[,"ref"],"/",positiveAlt2[,"alt2"],sep="")

        if(nrow(positiveAlt2)>0){
            res2 <- cbind(positiveAlt2[,c("pos","context")],baseChangeAlt2,positiveAlt2[,c("ref","alt2","cov_fwd","ref_fwd","alt2_fwd","cov_rev","ref_rev","alt2_rev")])
            names(res2) <- c("pos","context","basechange","ref","alt","cov_fwd","ref_fwd","alt_fwd","cov_rev","ref_rev","alt_rev")
        }

        res3 <- NULL
        alt3FreqsFwd <- tmpData[,"alt3_fwd"]/tmpData[,"cov_fwd"]
        alt3FreqsRev <- tmpData[,"alt3_rev"]/tmpData[,"cov_rev"]
        positiveAlt3 <- tmpData[intersect(which(alt3FreqsFwd>=minFreqDenovoSNPs),which(alt3FreqsRev>=minFreqDenovoSNPs)),,drop=FALSE]
        baseChangeAlt3 <- paste(positiveAlt3[,"ref"],"/",positiveAlt3[,"alt3"],sep="")

        if(nrow(positiveAlt3)>0){
            res3 <- cbind(positiveAlt3[,c("pos","context")],baseChangeAlt3,positiveAlt3[,c("ref","alt3","cov_fwd","ref_fwd","alt3_fwd","cov_rev","ref_rev","alt3_rev")])
            names(res3) <- c("pos","context","basechange","ref","alt","cov_fwd","ref_fwd","alt_fwd","cov_rev","ref_rev","alt_rev")
        }


        ## 4. Combine results to a nice looking table and annotate SNPs with pre-defined mutations
        resCombined <- rbind(res1, res2, res3)

        if(!is.null(resCombined)){

            sequence <- sapply(1:nrow(resCombined), function(i){ sub("\\[N\\]",paste("\\[",resCombined[i,"basechange"],"\\]",sep=""),resCombined[i,"context"]) })

            resultTable <- cbind(sequence,resCombined[,c("pos","ref","alt","cov_fwd","ref_fwd","alt_fwd","cov_rev","ref_rev","alt_rev")])

            tmpTableNames <- colnames(resultTable)

            freqsFwd <- resultTable[,"alt_fwd"]/resultTable[,"cov_fwd"]
            freqsRev <- resultTable[,"alt_rev"]/resultTable[,"cov_rev"]
            freqs <- (resultTable[,"alt_fwd"]+resultTable[,"alt_rev"])/(resultTable[,"cov_fwd"]+resultTable[,"cov_rev"])

            annotatedMutations <- mutationData[match(resultTable[,"sequence"],mutationData[,"sequence"]),"mutation"]

            resultTable <- cbind(resultTable,freqsFwd,freqsRev,freqs, annotatedMutations)
            colnames(resultTable) <- c(tmpTableNames,"freq_fwd","freq_rev","freq","annotation")

            resultTable <- resultTable[order(resultTable[,"freq"],decreasing=TRUE),]

            write.table(resultTable,file=deNovoVariantFiles[[refName]],sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

            cat(date()," - File created: ",deNovoVariantFiles[[refName]],"\n", sep="", file=logFile, append=TRUE)
        }

        ## 5. Create file with low covered regions, could be used later for indel detection

        lowCoverageData <- SNPcallData[which(coverageTotal < (median(coverageTotal)*minCoverageFraction)),]

        write.table(lowCoverageData,file=lowCoverageFiles[[refName]],sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        cat(date()," - File created: ",lowCoverageFiles[[refName]],"\n", sep="", file=logFile, append=TRUE)

        ## TODO: Perform indel analysis based on low coverage data results!?

    }


    ####################################################
    ##
	## Step 7. Check if any novel mutations were found
	##
    ####################################################

    for(refName in refNames){

        if(file.exists(deNovoVariantFiles[[refName]])){

            deNovoVariants <- read.table(deNovoVariantFiles[[refName]], sep="\t", header=TRUE)

            deNovoIds <- which(is.na(deNovoVariants[,"annotation"]))

            cat(date()," - ",length(deNovoIds)," putative novel mutations at >",minFreqDenovoSNPs*100,"% found\n", sep="", file=logFile, append=TRUE)

            if(length(deNovoIds) > 0){
                cat(paste(" - NOTE: ",length(deNovoIds)," putative novel mutations at >=",minFreqDenovoSNPs*100,"% reported in file: ",deNovoVariantFiles[[refName]],"!\n",sep=""))
                cat(date()," - File created: ",deNovoVariantFiles[[refName]],"\n", sep="", file=logFile, append=TRUE)
            }
        }
    }

    ###################################
    ##
	## Step 8. Run mutation screening
	##
    ###################################

    for(refName in refNames){

        mutationTableFile <- mutationTableFiles[refName]

        mutationData <- read.table(mutationTableFile, header=TRUE, as.is=TRUE)

        cat(date()," - Screening all variable positions in file '",mutationTableFile,"'\n", sep="", file=logFile, append=TRUE)

        ## A first round of mutation screening, not allowing for mismatches
        cmd <- paste("./cava_tablescreener.pl -f ",ccsFilteredFASTQfiles[[refName]]," -m ",mutationTableFile," -minq ",minq," -minqn ",minqn," -max_mm ",maxmmRound1," > ",mutationFilesRaw[[refName]],sep="") ## TODO: Change later!!!

        system(cmd)

        ## Rerun mutation screening allowing for mismatches, only for low coverage mutations
        results <- read.table(mutationFilesRaw[[refName]], sep="\t", as.is=TRUE, header=TRUE)

        row.names(results) <- as.character(results[,"sequence"])

        coverageTotal <- rowSums(results[,c("wt_reads_fwd","wt_reads_rev","mut_reads_fwd","mut_reads_rev","other_reads_fwd","other_reads_rev")])

        rerunSequences <- results[which(coverageTotal < median(coverageTotal)*coverageCutoff),"sequence"]

        ## Only rerun mutations where routine screening is set to 'yes'
        sequencesToScreen <- mutationData[which(mutationData[,"routine"] == "yes"),"sequence"]
        rerunSequences <- intersect(rerunSequences, sequencesToScreen)

        if(length(rerunSequences)>0){

            currentRerunNr <- 0

            for(rerunSequence in rerunSequences){
                currentRerunNr <- currentRerunNr+1
                mutName <- as.character(results[rerunSequence,"mutation"])
                cat("\r - Reanalysis of low covered position ",mutName," (",currentRerunNr," of ",length(rerunSequences),")", sep="")
                tmpResults <- doMutationScreening(ccsFilteredFASTQfiles[[refName]],rerunSequence, maxmmRound2, minq, minqn)
                results[rerunSequence,] <- c(mutName,strsplit(tmpResults,"\t")[[1]])
            }
            cat(" done\n")
        }

        write.table(results, file=mutationFilesRaw[[refName]], row.names=FALSE, quote=FALSE, sep="\t")

        cat(date()," - File created: ",mutationFilesRaw[[refName]],"\n", sep="", file=logFile, append=TRUE)

        ## Format mutation results into nice table

        data <- read.table(mutationFilesRaw[[refName]],sep="\t", comment.char="#",as.is=TRUE, header=TRUE)

        detection <- rep("unresolved",nrow(data))

        dataOut <- cbind(data,detection)
        colnames(dataOut) <- c(colnames(data),"detection")
        dataOut[,"detection"] <- as.character(dataOut[,"detection"])

        for(i in 1:nrow(dataOut)){

            totalReads <- as.numeric(dataOut[i,"wt_reads"]+dataOut[i,"mut_reads"]+dataOut[i,"other_reads"])

            if(totalReads>=minCoverage){
                dataOut[i,"detection"] <- "negative"
                if( (dataOut[i,"freq_fwd"]>minFreq) && (dataOut[i,"freq_rev"] > minFreq) ){
                    dataOut[i,"detection"] <- "positive"
                }
            }
            else{
                dataOut[i,"freq"] <- NA
                dataOut[i,"freq_fwd"] <- NA
                dataOut[i,"freq_rev"] <- NA
            }
        }

        dataOut[,"freq_fwd"] <-  round(dataOut[,"freq_fwd"],3)
        dataOut[,"freq_rev"] <-  round(dataOut[,"freq_rev"],3)
        dataOut[,"freq"] <-  round(dataOut[,"freq"],3)

        dataOut <- dataOut[order(dataOut[,"freq"], decreasing=TRUE),]

        routineMutations <- as.character(mutationData[mutationData[,"routine"]=="yes","sequence"])

        positiveMutations <- as.character(dataOut[which(dataOut[,"detection"] == "positive"),"sequence"])
        unresolvedMutations <- as.character(dataOut[which(dataOut[,"detection"] == "unresolved"),"sequence"])
        unresolvedMutations <- intersect(unresolvedMutations, routineMutations)

        negativeRoutineMutations <- setdiff(routineMutations, c(positiveMutations, unresolvedMutations))

        mutationsInOrder <- c(positiveMutations, unresolvedMutations, negativeRoutineMutations)

        dataOut <- dataOut[match(mutationsInOrder,dataOut[,"sequence"]),]

        dataColnames <- c(colnames(dataOut),"routine")
        dataOut <- cbind(dataOut,mutationData[match(dataOut[,"sequence"], mutationData[,"sequence"]),"routine"])
        colnames(dataOut) <- dataColnames

        ## Output mutation results
        if( (length(unresolvedMutations)/nrow(dataOut)) < (1-minFractionMuts)){ ## Sufficent mutations with high enough coverage..
            write.table(dataOut,file=mutationFilesFinal[[refName]], sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
            cat(date()," - File created: ",mutationFilesFinal[[refName]],"\n", sep="", file=logFile, append=TRUE)
        }
        else{ ## Otherwise QC has failed..
            write.table(dataOut,file=qcFailedFiles[[refName]], sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
            cat(date()," - QC failed.\n", sep="", file=logFile, append=TRUE)
            cat(date()," - File created: ",qcFailedFiles[[refName]],"\n", sep="", file=logFile, append=TRUE)
            cat(date()," - CLAMP analysis complete.\n", sep="", file=logFile, append=TRUE)
            return(0)
        }
    }

    ## ##################################################
    ##
    ## Step 9. Analyze clonal distribution of mutations
    ##
    ## ##################################################

    for(refName in refNames){
        if(file.exists(mutationFilesFinal[[refName]])){
            mutationsFinal <- read.table(mutationFilesFinal[[refName]], sep="\t", header=TRUE)

            mutationData <- read.table(mutationTableFile, header=TRUE, as.is=TRUE)

            mutations <- mutationsFinal[mutationsFinal[,"detection"] == "positive",,drop=FALSE]

            ## Only look at clonality for mutations where routine screening is set to 'yes'
            sequencesToScreen <- mutationData[which(mutationData[,"routine"] == "yes"),"sequence"]
            mutations <- mutations[mutations[,"sequence"] %in% sequencesToScreen,]
            nrMutations <- nrow(mutations)

            if(!is.na(nrMutations)){
                if(nrMutations>=2){
                    cat(date()," - Analyzing clonal distribution of mutations\n", sep="", file=logFile, append=TRUE)
                    cat(" - Analyzing clonal distribution of mutations...")

                    clonalResults <- analyzeMutationClonality(ccsFilteredFASTQfiles[[refName]],  mutationFilesFinal[[refName]], mutationTableFile, maxmmRound2, minq, minqn, restrictToTheseIds=sequencesToScreen)
                    write.table(clonalResults, file=mutationClonalityFiles[[refName]], sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
                    cat(date()," - File created: ",mutationClonalityFiles[[refName]],"\n", sep="", file=logFile, append=TRUE)
                    plotClonalDistribution(mutationClonalityFiles[[refName]], mutationClonalityPlotFiles[[refName]])
                    cat(date()," - File created: ",mutationClonalityPlotFiles[[refName]],"\n", sep="", file=logFile, append=TRUE)
                    cat(" done\n")
                }
            }
        }
    }

    ## ############################
    ##
	## Step 10. Compress fastq file
	##
    ## ############################

    if(!file.exists(ccsFASTQfileGZ)){
        if(file.exists(ccsFASTQfile)){
            cat(" - Compressing FASTQ file...")
            cmd <- paste("gzip ",ccsFASTQfile,sep="")
            system(cmd)
            cat(" done\n")
        }
    }

    cat(" - Analysis complete.\n")
    cat(date()," - CLAMP analysis complete.\n", sep="", file=logFile, append=TRUE)

    return(1)
}


#########################################
##
## Help function for mutation screening
##
#########################################

doMutationScreening <- function(ccsFile, mutationSequence, maxMismatch, minq=0, minqn=0, detailed=FALSE, useQualFilter=TRUE, outfile=NA){

    filterQualFlag <- ""

    if(useQualFilter == FALSE){
        filterQualFlag <- " -nofilter"
    }

    detailedFlag <- ""

    if(detailed == TRUE){
        detailedFlag <- " -d"
    }

    outfileStr <- ""

    if(!is.na(outfile)){
        outfileStr <- paste(" > ",outfile,sep="")
    }


    cmd <- paste("./cava.pl -f ",ccsFile," -t fastq -q ", mutationSequence, " -max_mm ",maxMismatch, " -minq ",minq," -minqn ",minqn,filterQualFlag,detailedFlag," --silent ",outfileStr,sep="")

    ##print(cmd)

    result <- system(cmd, intern=TRUE)

}


##########################################
##
## Help functions for clonality analysis
##
##########################################


## Extracts detailed information about all positive mutations and analyzes clonality in reads..
analyzeMutationClonality <- function(fastafile, mutationResultfile, mutationTableFile, maxMismatch, minq=0, minqn=0, strandSpecificMatching=TRUE, restrictToTheseIds=c(), specialAnalysisForWt=TRUE){

    mutationData <- read.table(mutationResultfile, as.is=TRUE, header=TRUE)

    mutationTable <- read.table(mutationTableFile, as.is=TRUE, header=TRUE)

    mutations <- mutationData[mutationData[,"detection"] == "positive",,drop=FALSE]

    ## Order mutation results by 'pos' if present in the table. Otherwise use same order as in mutation table.
    orderedMuts <- mutationTable[,"sequence"]
    if("pos" %in% colnames(mutationTable)){
       orderedMuts <- mutationTable[order(mutationTable[,"pos"], decreasing=TRUE),"sequence"]
    }

    mutations <- mutations[match(orderedMuts[orderedMuts %in% mutations[,"sequence"]], mutations[,"sequence"]),] #reorder table

    if(length(restrictToTheseIds)>0){
        mutations <- mutations[mutations[,"sequence"] %in% restrictToTheseIds,]
    }

    nrMutations <- nrow(mutations)

    if(nrMutations<2){
        stop("At least 2 mutations required for analysis of clonality!")
    }

    readIds <- list()

    for(i in 1:nrMutations){
        mutID <- as.character(mutations[i,1])
        seq <- as.character(mutations[i,2])

        outfile <- sub("final.txt",paste(mutID,".txt",sep=""),mutationResultfile)

        doMutationScreening(fastafile, seq, maxMismatch, minq, minqn, detailed=TRUE, outfile=outfile)

        data <- read.table(outfile, sep="\t", as.is=TRUE, comment.char="", quote="")

        readIds[[mutID]] <- list()
        readIds[[mutID]]$wt_fwd <- as.character(data[which(data[,3] %in% c("ref_fwd","other_fwd")),1]) ## Is it good that 'other' is included here!?!?
        readIds[[mutID]]$wt_rev <- as.character(data[which(data[,3] %in% c("ref_rev","other_rev")),1])
		readIds[[mutID]]$mut_fwd <- as.character(data[which(data[,3] == "alt_fwd"),1])
        readIds[[mutID]]$mut_rev <- as.character(data[which(data[,3] == "alt_rev"),1])
    }

    mutationNames <- names(readIds)

    vec <- c("wt","mut")
    lst <- lapply(numeric(nrMutations), function(x) vec)

    allCombinations <- as.matrix(expand.grid(lst))
    allCombinations <- sapply(1:ncol(allCombinations), function(i){ allCombinations[,i]} )

    colnames(allCombinations) <- mutationNames

    results <- cbind(allCombinations,rep(0,nrow(allCombinations)))
    colnames(results) <- c(colnames(allCombinations),"read.count")
    rownames(results) <- sapply(1:nrow(allCombinations),function(i){paste(allCombinations[i,],collapse="|")})

    for(i in 1:nrow(results)){
        nreads <- -1
        matchingReads <- NULL
        matchingReadsFwd <- NULL
        matchingReadsRev <- NULL
        for(mut in mutationNames){
            status <- results[i,mut]
            if(is.null(matchingReads)){
                if(status == "wt"){
                    matchingReadsFwd <- as.character(readIds[[mut]]$wt_fwd)
                    matchingReadsRev <- as.character(readIds[[mut]]$wt_rev)
                    matchingReads <- c(matchingReadsFwd,matchingReadsRev)
                }
                if(status == "mut"){
                    matchingReadsFwd <- as.character(readIds[[mut]]$mut_fwd)
                    matchingReadsRev <- as.character(readIds[[mut]]$mut_rev)
                    matchingReads <- c(matchingReadsFwd,matchingReadsRev)
                }
            }
            else{
                if(status == "wt"){
                    tmpReadsFwd <- readIds[[mut]]$wt_fwd
                    tmpReadsRev <- readIds[[mut]]$wt_rev
                    tmpReads <- c(matchingReadsFwd,matchingReadsRev)
                    matchingReadsFwd <- intersect(matchingReadsFwd,tmpReadsFwd)
                    matchingReadsRev <- intersect(matchingReadsRev,tmpReadsRev)
                    matchingReads <- intersect(matchingReads,tmpReads)
                }
                if(status == "mut"){
                    tmpReadsFwd <- readIds[[mut]]$mut_fwd
                    tmpReadsRev <- readIds[[mut]]$mut_rev
                    tmpReads <- c(matchingReadsFwd,matchingReadsRev)
                    matchingReadsFwd <- intersect(matchingReadsFwd,tmpReadsFwd)
                    matchingReadsRev <- intersect(matchingReadsRev,tmpReadsRev)
                    matchingReads <- intersect(matchingReads,tmpReads)
                }
            }
        }

        if(strandSpecificMatching){
            nreads <- length(matchingReadsFwd)+length(matchingReadsRev)
        }
        else{
            nreads <- length(matchingReads)
        }

        results[i,"read.count"] <- nreads
    }


    results <- results[order(as.numeric(results[,"read.count"]), decreasing=TRUE),]
    results <- results[which(as.numeric(results[,"read.count"])>0),,drop=FALSE]

    ## Possible chimeric molecules are identified below, 'passed.filter' is set to 'no' for all such cases
    flaggedResults <- results

    if(nrow(flaggedResults)>1){
        flaggedResults <- flagChimericReads(results)
    }

    ## The steps below are performed to identify special cases of wt molecules that
    ## may have been wrongly labeled as chimeric in the previous filtering step.
    ## Note that this only works for somatic mutations and should not be performed for SNPs!
    if(specialAnalysisForWt){

        wtMolecules <- NULL

        if(nrow(flaggedResults)>=2){

            for(i in 1:nrow(flaggedResults)){
                for(j in 1:nrow(flaggedResults)){

                    if( (flaggedResults[i,"passed.filter"] == "yes") && (flaggedResults[j,"passed.filter"] == "yes") ){

                        dataMuts1 <- flaggedResults[i,-which(colnames(flaggedResults) %in% c("read.count","passed.filter"))]
                        dataMuts2 <- flaggedResults[j,-which(colnames(flaggedResults) %in% c("read.count","passed.filter"))]

                        mismatchPos <- which(!(dataMuts1 == dataMuts2))

                        if(length(mismatchPos) == 2){
                            wtMolecule <- as.character(dataMuts1)
                            wtMolecule[mismatchPos] <- rep("wt",2)
                            wtMolecules <- rbind(wtMolecules,c(paste(wtMolecule,collapse="|"),max(i,j)))
                        }
                    }
                }
            }
        }

        ## Go through all molecules previously classed as 'no', and change to 'yes' if the molecule can be produced
        ## by a combination of two molecules as identified above. Both molecules must be occuring ABOVE in the original
        ## clonal distribution table.
        if(!is.null(wtMolecules)){

            wtMolecules <- unique(wtMolecules)

            colnames(wtMolecules) <- c("mut.string","max.row")

            for(i in 1:nrow(flaggedResults)){
                if(flaggedResults[i,"passed.filter"] == "no"){
                    moleculeStr <- paste(flaggedResults[i,-which(colnames(flaggedResults) %in% c("read.count","passed.filter"))],collapse="|")
                    wtMoleculesSubset <- wtMolecules[which(as.numeric(wtMolecules[,"max.row"])<i),,drop=FALSE]
                    if(moleculeStr %in% wtMoleculesSubset[,"mut.string"]){
                        flaggedResults[i,"passed.filter"] <- "YES"
                    }
                }
            }
        }
    }

    return(flaggedResults)

}


## Flags all chimeric reads in a clonal distribution table.
## Mutations in the input data table must be sorted by genomic coordinates.
flagChimericReads <- function(clonalityTable){

    createChimericRead <- function(data, prefixRow, suffixRow, breakPosition){

        ## take first n positions from data at row i
        prefix <- data[prefixRow,1:breakPosition]

        ## take last n positions from data at row j
        suffix <- data[suffixRow,(breakPosition+1):ncol(data)]

        newString <- paste(c(prefix,suffix), collapse="|")

        return(newString)

    }

    createAllChimericReads <- function(data){

        results <- NULL

        ## loop through all combiations of rows
        for(i in 1:nrow(data)){
            for(j in 1:nrow(data)){
                ## loop through all possible break positions
                for(pos in 1:(ncol(data)-1)){
                    newString <- createChimericRead(data, i, j, pos)

                    results <- c(newString, results)

                }
            }
        }

        return(results)
    }

    data <- clonalityTable[,-which(colnames(clonalityTable) %in% c("read.count")),drop=FALSE]

    passed.filter <- rep("yes",nrow(data))

    ## Go through each row, check if it could be generated by some chimeric read
    ## created from all rows above.
    for(i in 1:nrow(data)){

        ## vector with values at row i
        currentRow <- paste(data[i,],collapse="|")

        ## table with all rows above row i
        if(i>1){
            dataRowsAbove <- data[1:(i-1),,drop=FALSE]
            allChimericStrings <- unique(createAllChimericReads(dataRowsAbove))

            if(grepl("mut",currentRow)==TRUE){
                if(currentRow %in% allChimericStrings){
                    passed.filter[i] <- "no"
                }
            }
        }

    }

    resultData <- cbind(clonalityTable,passed.filter)

    return(resultData)
}


## Creates a PDF plot of a 'flagged' clonal distribution table
plotClonalDistribution <- function(clonalDistFile, clonalDistributionPdf){
    dataOrig <- read.table(clonalDistFile,sep="\t",header=TRUE)

    dataYes <- dataOrig[which(toupper(dataOrig[,"passed.filter"])=="YES"),,drop=FALSE]

    nrRowsInPlot <- max(nrow(dataYes),10)

    pdf(clonalDistributionPdf, height=nrRowsInPlot, width=10)

    ## x=nr(sites), y=nr(reads/rows)
    plot <- plot(NULL,xlim=c(0,ncol(dataYes)+2),ylim=c(0,nrRowsInPlot),bty="n",xlab="",ylab="",axes=F)
    text(x=ncol(dataYes), y=nrRowsInPlot, label="Frequency", font=3)
    text(x=ncol(dataYes)+1.5, y=nrRowsInPlot, label="Reads", font=3)

    ## Loop through rows, plot mutations and percent reads
    for(i in 1:nrow(dataYes)){

        xvec <- dataYes[i,1:(ncol(dataYes)-2)]

        ycoord <- nrRowsInPlot-i

        segments(x0=0, y0=ycoord, x1=length(xvec)+1, y1=ycoord, lwd=7, col="grey")

        points(x=which(xvec == "mut"), y=rep(ycoord,length(which(xvec =="mut"))),pch=15,col="red")

        tmpNames <- colnames(dataYes[which(xvec == "mut")])

        if(length(tmpNames)>0){
            text(x=which(xvec == "mut"), y=ycoord, label=tmpNames, pos=3)
        }
        percentReads <- (dataYes[i,(ncol(dataYes)-1)]/sum(dataYes[,(ncol(dataYes)-1)]))*100
        approxPercent <- paste(substr(percentReads,1,4),"%",sep=" ")

        if(percentReads < 1){

            approxPercent <- "< 1 %"
        }

        text(x=ncol(dataYes), y=ycoord, label=approxPercent)

        text(x=ncol(dataYes)+1.5, y=ycoord, label=dataYes[i,ncol(dataYes)-1])
    }

    dev.off()
}



#######################################
##
## Help function for isoform analysis
##
#######################################

## Function for extracting all candidate isoforms from an 'identical read' file
extractIsoforms <- function(identicalReadsFile, outfile, maxIndelLen=5, maxSNPs=10){

    compareIsoforms <- function(iso1, iso2, maxIndelLen=5, maxSNPs=10){

        vec1 <- strsplit(iso1,"")[[1]]
        vec2 <- strsplit(iso2,"")[[1]]

        i <- 1

        lenOut <- max(length(vec1),length(vec2))

        minLen <- min(length(vec1),length(vec2))

        vec1out <- rep("-",lenOut)
        vec2out <- rep("-",lenOut)

        while((vec1[i]==vec2[i]) && (i<=minLen)){
            vec1[i] <- "-"
            vec2[i] <- "-"
            i <- i+1
        }

        i <- 1

        while((rev(vec1)[i]==rev(vec2)[i]) && (i<=minLen)){
            index1 <- length(vec1)-i+1
            index2 <- length(vec2)-i+1
            vec1[index1] <- "-"
            vec2[index2] <- "-"
            i <- i+1
        }

        mismatchPos1 <- which(vec1 != "-")
        mismatchPos2 <- which(vec2 != "-")

        if((length(mismatchPos1) == 0) && (length(mismatchPos2) == 0)){
            return("identical")
        }

        if(length(mismatchPos1) == length(mismatchPos2)){ ## SNP differences between isoforms!
            if(all(mismatchPos1 == mismatchPos2)){
                nrSNPs <- length(which(vec1 != vec2))
                if(nrSNPs <= maxSNPs){
                    return("snp")
                }
            }
        }
        else{
            if(abs(length(mismatchPos1)-length(mismatchPos2))<maxIndelLen){ ## Indel differences between isoforms!
                return("indel")
            }
            else{
                return("splice") ## Potential splicing differences between isoforms
            }
        }

        return(list(vec1,vec2))

    }

    require(seqinr)

    identicalReads <- read.fasta(identicalReadsFile, as.string=TRUE)

    isoforms <- NULL
    isoformLengths <- c()

   if(file.exists(outfile)){
       file.remove(outfile)
   }

    for(i in 1:length(identicalReads)){
        #cat(i,"/",length(identicalReads),"\r",sep="")
        ccsSeq <- toupper(identicalReads[[i]])
        ##ccsSeqCount <- ccsReadsMult[i,2]
        ##ccsSeqVec <- strsplit(ccsSeq, "")[[1]]
        ccsSeqLen <- nchar(ccsSeq)

        if(is.null(isoforms)){
            isoforms <- c(isoforms, ccsSeq)
        }

        else{
            newIsoform <- TRUE

            for(j in 1:length(isoforms)){
                ##print(ccsSeq)
                ##print(isoforms[j])
               res <- compareIsoforms(ccsSeq, isoforms[j], maxIndelLen=maxIndelLen, maxSNPs=maxSNPs)

               if(length(res) != 1){
                   newIsoform <- FALSE
               }
               else{
                   if(res != "splice"){
                       newIsoform <- FALSE
                   }
               }
           }
            if(newIsoform){
                isoforms <- c(isoforms, ccsSeq)
            }
        }
    }


    write.table(isoforms, file=outfile, sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)

    invisible(isoforms)
}



###############################
##
## General help functions...
##
###############################


## Generates all possible primer pairs from degenerate fwd, rev sequences
generatePrimerPairs <- function(fwdPrimer, revPrimer){

    IUPAC <- list("A"="A","C"="C","G"="G","T"="T","M"=c("A","C"),"R"=c("A","G"),"Y"=c("C","T"),"S"=c("G", "C"),"W"=c("A","T"),"K"=c("G","T"),"B"=c("C","G","T"),"D"=c("A","G","T"),"H"=c("A","C","T"),"V"=c("A","C","G"),"N"=c("A","G","C","T"))


    fwdVec <- unlist(strsplit(fwdPrimer, split=""))

    fwdComboTable <- expand.grid(IUPAC[fwdVec], stringsAsFactors=FALSE)

    fwdPrimerVector <- sapply(1:nrow(fwdComboTable), function(x) { paste(fwdComboTable[x,], collapse="") })


    revVec <- unlist(strsplit(revPrimer, split=""))

    revComboTable <- expand.grid(IUPAC[revVec], stringsAsFactors=FALSE)

    revPrimerVector <- sapply(1:nrow(revComboTable), function(x) { paste(revComboTable[x,], collapse="") })


    combinedPrimers <- expand.grid(fwdPrimerVector, revPrimerVector)
    return(combinedPrimers)

}

## Returns the reverse complement of a IUPAC or consensus sequence
revComp <- function(seq){

    baseComp <- array(NA,0)
    baseComp[c("A","C","G","T","R","Y","K","M","S","W","B","V","D","H","N","[","]","(",")")] <- c("T","G","C","A","Y","R","M","K","S","W","V","B","H","D","N","]","[",")","(")

    return(paste(rev(sapply(strsplit(toupper(seq),"")[[1]],function(x){
        if(!is.na(baseComp[x]))
        {
            return(baseComp[x])
        }else{
            return(x)
        }
    },USE.NAMES = FALSE)),collapse = ""))
}


###################################################
##
## Function for processing samples in a batch mode
##
###################################################

## This function takes a file with sample info and primer sequences as arguments, and runs the screenMutation script on all samples in file.
screenMutationsBatch <- function(sampleInfoFile, assay, mutationTableFile, primerFile="primer_sequences.txt", filterPrimerSeq=TRUE, analyzeIsoforms=FALSE, rerunAnalysis=FALSE){

    ## Read all sample info
    data <- read.table(sampleInfoFile, sep="\t", header=TRUE)

	primerData <- NULL

	if(filterPrimerSeq){
	   ## Read primer sequence information
	   primerData <- read.table(primerFile, sep="\t", header=TRUE)

	   ## Check that primer sequences are available for all samples!
	   if(!all(as.character(data[,"primerId"]) %in% as.character(primerData[,"primerId"]))){
           stop(paste("Unspecified primer sequences in file:",sampleFile))
	   }
	}

    for(i in 1:nrow(data)){

        runId <- as.character(data[i,"runId"])

		primerFwd <- NA
		primerRev <- NA

		if(filterPrimerSeq){
		    primerId <- as.character(data[i,"primerId"])
			primerFwd <- as.character(primerData[which(primerData[,"primerId"] == primerId),"fwd"])
			primerRev <- as.character(primerData[which(primerData[,"primerId"] == primerId),"rev"])
		}

        stop("CLAMP without reference/primer sequences not yet supported!!")

        screenMutations(runId, assay=assay, mutationTableFile=mutationTableFile, fwdPrimer=primerFwd, revPrimer=primerRev, analyzeIsoforms=analyzeIsoforms, rerunAnalysis=rerunAnalysis)
   }

}


## ##############################################################
##
## Function for reading the information in a 'target info' file
##
## ###############################################################

getTargetInfo <- function(primerId, assay, targetInfoFile="target_info.txt"){

    ## Read primer sequence information
    primerData <- read.table(targetInfoFile, sep="\t", header=TRUE)

    ## Check that primer sequences are available for all samples!
    if(!(primerId %in% as.character(primerData[,"primerId"]))){
        stop(paste("No target info available for primerId '",primerId,"' in file: ",targetInfoFile,"\n",sep=""))
    }

    primerFwd <- as.character(primerData[which(primerData[,"primerId"] == primerId),"fwd"])
    primerRev <- as.character(primerData[which(primerData[,"primerId"] == primerId),"rev"])
    referenceStr <- as.character(primerData[which(primerData[,"primerId"] == primerId),"references"])
    refNames <- unlist(strsplit(referenceStr,","))
    refFiles <- paste("references/",assay,"/",refNames,sep="")

    for(refFile in refFiles){
        if(!file.exists(refFile)){
            stop(paste("Reference file '",refFile,"' missing for sample ",runId,"!",sep=""))
        }
        if(length(grep(".fasta",refFile)) == 0){
            stop("Error. all reference file names must have suffix '.fasta'")
        }
    }

    mutTableStr <- as.character(primerData[which(primerData[,"primerId"] == primerId),"mutationTables"])
    mutTables <- unlist(strsplit(mutTableStr,","))

    ## Assign the same mutation table to all references if only one is defined
    if( (length(mutTables) == 1) && (length(refFiles) > 1) ){
        mutTables <- rep(mutTables[1], length(refFiles))
    }

    names(mutTables) <- sub(".fasta","",refNames)

    return(list("primerFwd"=primerFwd,"primerRev"=primerRev,"refFiles"=refFiles,"mutTables"=mutTables))

}

