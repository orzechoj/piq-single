#!/usr/bin/Rscript
#Rscript script commondir bamout
#example:
#Rscript bam2rdata.r /cluster/thashim/basepiq/common.r /cluster/thashim/basepiq/tmp/bams.RData /cluster/cwo/dnase_seq/bams/D0_50-100_130801.bwa.mapq20.mm10.bam

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

commondir = args[1]
bamout = args[2]
bamname = args[3]

source(commondir)

###
# read bam

plusflags = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=F,isMinusStrand=F,isDuplicate=F,isNotPassingQualityControls=F))
plusstrand = readBamGappedAlignments(bamname,param=plusflags)

minusflags = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=F,isMinusStrand=T,isDuplicate=F,isNotPassingQualityControls=F))
minusstrand = readBamGappedAlignments(bamname,param=minusflags)

save(plusstrand,minusstrand,file=bamout)

