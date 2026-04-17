#!/usr/bin/env Rscript 

library(ATACseqQC)
library(Rsamtools)
library(BiocIO)

args <- commandArgs(trailingOnly = TRUE)

## input is bamFile
bamfile <- args[1]
outbam  <- args[2]

## bamfile tags

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

## shift the bam file by the 5'ends
gal <- readBamFile(bamfile, tag=tags, asMates=FALSE,bigFile=TRUE)
gal1 <- shiftGAlignments(gal)
export(gal1, outbam)

