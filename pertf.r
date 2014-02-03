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
bamfile = args[5]
pwmid = args[6]

source('utils/copybam.r')
load(paste0(pwmdir,pwmid,'.pwmout.RData'))
if(sum(clengths[1])>0){
source('loadbam.r')	
source('cluster.r')
source('bindcall.r')
}

