#!/usr/bin/env Rscript

######## doATACseqQC ########
# This script will automate 
# the process of generating
# ATACseqQC reports
# user guide: https://bioconductor.org/packages/3.14/bioc/vignettes/ATACseqQC/inst/doc/ATACseqQC.html

# INPUT:

# bamfile: path to alignment bamfile to QC

# OUTPUT: 

# png figure of fragment sizes QC

# PARAMS:

# species: the species that is being processed 
#         current useable species: hg38, mm10

# Get the arguments as a character vector.
myargs <- commandArgs(trailingOnly=TRUE)
bamfile <- myargs[1]
species <- myargs[2]

print("loading packages (ATACseqQC, ggplot, etc)...")
suppressPackageStartupMessages(library(ggplot2, quietly=TRUE))
suppressPackageStartupMessages(library(Rsamtools, quietly=TRUE))
suppressPackageStartupMessages(library(ATACseqQC, quietly=TRUE))
suppressPackageStartupMessages(library(ChIPpeakAnno, quietly=TRUE))
suppressPackageStartupMessages(library("GenomicAlignments", quietly=TRUE))

if (species == "mm") {
  suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene, quietly=TRUE))
  suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10, quietly=TRUE))
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  bsgenome <- BSgenome.Mmusculus.UCSC.mm10
  genome <- Mmusculus
  print("species is 'mm' using mm10 for analysis")
  ### Note: Everything below is deprecated until I can figure out a way to port a 
  ### static/local package with snakemake
  # Note: phastCons60way was manually curated from GenomicAlignments, built, and installed as an R package
  # score was obtained according to: https://support.bioconductor.org/p/96226/
  # package was built and installed according to: https://www.bioconductor.org/packages/devel/bioc/vignettes/GenomicScores/inst/doc/GenomicScores.html
  # (section 5.1: Building an annotation package from a GScores object)
  #suppressWarnings(suppressPackageStartupMessages(library(GenomicScores, lib.loc="/users/dia6sx/snakeATAC/scripts/", quietly=TRUE)))
  #suppressWarnings(suppressPackageStartupMessages(library(phastCons60way.UCSC.mm10, lib.loc="/users/dia6sx/snakeATAC/scripts/", quietly=TRUE)))
} else if (species == "hs") {
  suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene, quietly=TRUE))
  suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38, quietly=TRUE))
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  bsgenome <- BSgenome.Hsapiens.UCSC.hg38
  genome <- Hsapiens
  print("species is 'hs' using hg38 for analysis")
} else {
  print(paste("params ERROR: ATACseqQC is not configured to use species =", species))
  print("exiting...")
  quit(status=1)
}

doATACseqQC <- function(bamfile, txdb, bsgenome, genome) {
    # Fragment size distribution
    print(paste("generating output for ",strsplit(basename(bamfile),split='\\.')[[1]][1],"...",sep=""))
    print("calculating Fragment size distribution...")
    bamfile.labels <- gsub(".bam", "", basename(bamfile))
    loc_to_save_figures <- paste(dirname(dirname(bamfile)),"/qc/ATACseqQC",sep="")
    if (file.exists(loc_to_save_figures)) {
        print("Warning: old figures will be overwritten")
    } else {
        dir.create(loc_to_save_figures)
    }
    png_file <- paste(loc_to_save_figures,"/",bamfile.labels,"_fragment_size_distribution.png",sep="")
    png(png_file)
    fragSizeDist(bamfile, bamfile.labels)
    dev.off()

    # Adjust the read start sites
    print("adjusting read start sites...")
    ## bamfile tags to be read in
    possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                    "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                    "TC", "UQ"), 
                    "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                                "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                                "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                                "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                                "U2"))
    bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 100),
                     param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
    tags <- names(bamTop100)[lengths(bamTop100)>0]
    ## files will be output into outPath
    ## shift the coordinates of 5'ends of alignments in the bam file
    outPath <- paste(dirname(dirname(bamfile)),"/alignments_shifted", sep="")
    seqinformation <- seqinfo(txdb)
    gal <- readBamFile(bamfile, tag=tags, asMates=TRUE, bigFile=TRUE)
    shiftedBamfile <- file.path(outPath, paste(bamfile.labels,"_shifted.bam",sep=""))
    # check if shifted Bam file exists from previous run
    if (file.exists(shiftedBamfile)) {
        print("Shifted Bamfile found.")
        print("Loading in...")
        gal <- readBamFile(shiftedBamfile, tag=tags, asMates=TRUE, bigFile=TRUE)
        ## This step is mostly for formating so splitBam can
        ## take in bamfile. Implementing shift of 0 bp on positive strand
        ## and 0 bp on negative strand because shifted Bamfile should
        ## already have these shifts
        gal1 <- shiftGAlignmentsList(gal, positive = 0L, negative = 0L)
    } else {
        # shifted bam file does not exist check if
        # old shifted alignments directory exists
        # if so remove and create new one
        if (file.exists(outPath)){
            unlink(outPath,recursive=TRUE)
        }
        dir.create(outPath)
        print("*** creating shifted bam file ***")
        gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)
    }

    # Promoter/Transcript body (PT) score
    print("calculating Promoter/Transcript body (PT) score...")
    txs <- transcripts(txdb)
    pt <- PTscore(gal1, txs)
    png_file <- paste(loc_to_save_figures,"/",bamfile.labels,"_ptscore.png",sep="")
    png(png_file)
    plot(pt$log2meanCoverage, pt$PT_score, 
        xlab="log2 mean coverage",
        ylab="Promoter vs Transcript",
        main=paste(bamfile.labels,"PT score"))
    dev.off()

    # Nucleosome Free Regions (NFR) score
    print("calculating Nucleosome Free Regions (NFR) score")
    nfr <- NFRscore(gal1, txs)
    png_file <- paste(loc_to_save_figures,"/",bamfile.labels,"_nfrscore.png",sep="")
    png(png_file)
    plot(nfr$log2meanCoverage, nfr$NFR_score, 
        xlab="log2 mean coverage",
        ylab="Nucleosome Free Regions score",
        main=paste(bamfile.labels,"\n","NFRscore for 200bp flanking TSSs",sep=""),
        xlim=c(-10, 0), ylim=c(-5, 5))
    dev.off()

    # Transcription Start Site (TSS) Enrichment Score
    print("calculating Transcription Start Site (TSS) Enrichment score")
    tsse <- TSSEscore(gal1, txs)
    png_file <- paste(loc_to_save_figures,"/",bamfile.labels,"_tss_enrichment_score.png",sep="")
    png(png_file)
    plot(100*(-9:10-.5), tsse$values, type="b", 
        xlab="distance to TSS",
        ylab="aggregate TSS score",
        main=paste(bamfile.labels,"\n","TSS Enrichment score",sep=""))
    dev.off()

    # Split reads, Heatmap and coverage curve for nucleosome positions
    print("splitting reads by fragment length...")
    genome <- genome
    outPath <- paste(dirname(dirname(bamfile)),"/alignments_split", sep="")
    TSS <- promoters(txs, upstream=0, downstream=1)
    TSS <- unique(TSS)
    ## estimate the library size for normalization
    librarySize <- estLibSize(bamfile)
    ## calculate the signals around TSSs.
    NTILE <- 101
    dws <- ups <- 1010
    splitBamfiles <- paste(outPath,"/",c("NucleosomeFree", 
                                             "mononucleosome",
                                             "dinucleosome",
                                             "trinucleosome"),".bam",sep="")
    # check if split Bam files exists from previous run
    if (all(file.exists(splitBamfiles))) {
        print("*** split bam files found! ***")
        print("Loading in...")
        sigs <- enrichedFragments(bamfiles=splitBamfiles,
                                    index=splitBamfiles, 
                                    TSS=TSS,
                                    librarySize=librarySize,
                                    TSS.filter=0.5,
                                    n.tile = NTILE,
                                    upstream = ups,
                                    downstream = dws)
    } else {
        # split bam files do not exist check if
        # old split alignments directory exists
        # if so remove and create new one
        if (file.exists(outPath)){
            unlink(outPath,recursive=TRUE)
        }
        print("*** creating split bam files ***")
        dir.create(outPath)
        ## split the reads into NucleosomeFree, mononucleosome, 
        ## dinucleosome and trinucleosome.
        ## and save the binned alignments into bam files.
        objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = outPath)
        #objs <- splitBam(bamfile, tags=tags, outPath=outPath,
            #        txs=txs, genome=genome,
            #       conservation=phastCons60way.UCSC.mm10,
            #      seqlev=paste0("chr", c(1:19, "X", "Y")))
        sigs <- enrichedFragments(gal=objs[c("NucleosomeFree", 
                                        "mononucleosome",
                                        "dinucleosome",
                                        "trinucleosome")], 
                                    TSS=TSS,
                                    librarySize=librarySize,
                                    TSS.filter=0.5,
                                    n.tile = NTILE,
                                    upstream = ups,
                                    downstream = dws)
    }
    ## log2 transformed signals
    sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
    ## plot heatmap
    png_file <- paste(loc_to_save_figures,"/",bamfile.labels,"_nucleosome_pos_heatmap.png",sep="")
    png(png_file)
    featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                        zeroAt=.5, n.tile=NTILE)
    dev.off()
    ## get signals normalized for nucleosome-free and nucleosome-bound regions.
    out <- featureAlignedDistribution(sigs, 
                                    reCenterPeaks(TSS, width=ups+dws),
                                    zeroAt=.5, n.tile=NTILE, type="l", 
                                    ylab="Averaged coverage")
    ## rescale the nucleosome-free and nucleosome signals to 0~1
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    out <- apply(out, 2, range01)
    png_file <- paste(loc_to_save_figures,"/",bamfile.labels,"_TSS_signal_distribution.png",sep="")
    png(png_file)
    matplot(out, type="l", xaxt="n", 
            xlab="Position (bp)", 
            ylab="Fraction of signal",
            main=paste(bamfile.labels,"\n","TSS signal distribution",sep=""))
    axis(1, at=seq(0, 100, by=10)+1, 
        labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
    abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
    dev.off()

    print("QC Finished.")
    print("Generated QC figures can be found in qc folder under ATACseQC")
    print(paste("*** removing temp files in",outPath,"***"))
    unlink(outPath,recursive=TRUE)
    outPath <- paste(dirname(dirname(bamfile)),"/alignments_shifted", sep="")
    print(paste("*** removing temp files in",outPath,"***"))
    unlink(outPath,recursive=TRUE)
}

doATACseqQC(bamfile, txdb, bsgenome, genome)
