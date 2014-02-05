require('snow')

#####
# Below is an example of how to do DNase-seq calls using a multi-cluster machine using passwordless ssh + the snow package.
#####

###
# Generic helper to define a new cluster. rscript and snowlib should point to the R binary and linbrary locations on each machine.
getCluster<-function(n){
  setDefaultClusterOptions(port=10186)
  tstOptions <- c(replicate(2*n,list(host="pax6",rscript="/usr/bin/Rscript",snowlib="/afs/csail.mit.edu/u/t/thashim/R/x86_64-pc-linux-gnu-library/2.14"),simplify=FALSE),
                  replicate(2*n,list(host="cdx2",rscript="/usr/bin/Rscript",snowlib="/afs/csail.mit.edu/u/t/thashim/R/x86_64-pc-linux-gnu-library/2.14"),simplify=FALSE),
                  replicate(1*n,list(host="med1",rscript="/usr/bin/Rscript",snowlib="/afs/csail.mit.edu/u/t/thashim/R/x86_64-pc-linux-gnu-library/2.14"),simplify=FALSE),
                  replicate(2*n,list(host="hb9",rscript="/usr/bin/Rscript",snowlib="/afs/csail.mit.edu/u/t/thashim/R/x86_64-pc-linux-gnu-library/2.14"),simplify=FALSE),
                  replicate(1*n,list(host="rap1",rscript="/usr/bin/Rscript",snowlib="/afs/csail.mit.edu/u/t/thashim/R/x86_64-pc-linux-gnu-library/2.14"),simplify=FALSE),
                  replicate(2*n,list(host="foxp1",rscript="/usr/bin/Rscript",snowlib="/afs/csail.mit.edu/u/t/thashim/R/x86_64-pc-linux-gnu-library/2.14"),simplify=FALSE),
                  replicate(1*n,list(host="fkh2",rscript="/usr/bin/Rscript",snowlib="/afs/csail.mit.edu/u/t/thashim/R/x86_64-pc-linux-gnu-library/2.14"),simplify=FALSE),
                  replicate(1*n,list(host="phd1",rscript="/usr/bin/Rscript",snowlib="/afs/csail.mit.edu/u/t/thashim/R/x86_64-pc-linux-gnu-library/2.14"),simplify=FALSE)
                  )
  makeCluster(tstOptions,type="SOCK")
}

cl = getCluster(8)

#Location of the common.r file that does package load and defines parameters.
commonfile='/cluster/thashim/basepiq/common.mm.r'

#Directory in which the PWM files were outputted using pwmmatch
pwmdir='/cluster/thashim/tmppiq/'

#Temporary read/write dir. Should be fast and local to avoid bottlenecking
tmpdir='/scratch/tmp/'

#Output directory to which all call files are written
outdir='/cluster/thashim/130130.mm10.d0/'

#Processed bamfile created by bam2rdata script.
bamfile='/cluster/thashim/tmppiq/d0.RData'

#Number of motifs to process 
pwmnum = 1316

clusterExport(cl,c("commonfile","pwmdir","tmpdir","outdir","bamfile"))

cla=clusterApplyLB(cl,(1:pwmnum),function(i){
        print(i)
	assign("pwmid",i,envir=.GlobalEnv)
	pwmid = i
	if(file.exists(paste0(pwmdir,pwmid,'.pwmout.RData'))){
	load(paste0(pwmdir,pwmid,'.pwmout.RData'))
	if(sum(clengths[1])>0){
	source('/cluster/thashim/basepiq/loadbam.r')	
	source('/cluster/thashim/basepiq/cluster.r')
	source('/cluster/thashim/basepiq/bindcall.r')
	}}
})

