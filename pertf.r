#!/usr/bin/Rscript
#Rscript script commondir bamout
#example:
#Rscript pertf.r /cluster/thashim/basepiq/common.r /cluster/thashim/basepiq/tmp/pwmout.RData /cluster/thashim/basepiq/tmp/ /cluster/thashim/basepiq/output/

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

commonfile = args[1]
pwmdir = args[2]
tmpdir = args[3]
outdir = args[4]
bamdir = args[5]

source('cluster.r')
source('bindcall.r')
