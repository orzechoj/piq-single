require('snow')
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

commonfile='/cluster/thashim/basepiq/common.mm.r'
pwmdir='/cluster/thashim/tmppiq/'
tmpdir='/scratch/tmp/'
outdir='/cluster/thashim/130130.mm10.d0/'
bamfile='/cluster/thashim/tmppiq/d0.RData'

clusterExport(cl,c("commonfile","pwmdir","tmpdir","outdir","bamfile"))

#exists = as.double(sapply(strsplit(list.files(outdir,'.pdf'),'-'),function(i){i[1]}))

cla=clusterApplyLB(cl,(1:1316),function(i){
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

