#!/usr/bin/Rscript
#Rscript script commondir pwmdir tmpdir outdir bamfile pwmid 
#example:
#Rscript pertf.r /cluster/thashim/basepiq/common.r /cluster/thashim/tmppiq/ /scratch/tmp/ /cluster/thashim/130130.mm10.d0/ /cluster/thashim/tmppiq/d0.RData 139

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

#location of common.r containing runtime parameters
commonfile = args[1]

#directory where pwm matches are stored
pwmdir = args[2]

#directory to use as fast temporary storage
tmpdir = args[3]

#location of output calls
outdir = args[4]

#location of the bam RData file made by bam2rdata.r
bamfile = args[5]

#which pwm file to use in pwmdir
pwmid = args[6]

source('utils/copybam.r')
load(paste0(pwmdir,pwmid,'.pwmout.RData'))
if(sum(clengths[1])>0){
source('loadbam.r')	
source('cluster.r')
source('bindcall.r')
}

