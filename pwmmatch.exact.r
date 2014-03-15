#!/usr/bin/Rscript
#Rscript script common-path jaspar-dir id output
#example:
#Rscript pwmmatch.r /cluster/thashim/basepiq/common.r /cluster/thashim/basepiq/pwms/jaspar.txt 141 /cluster/thashim/basepiq/tmp/pwmout.RData

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

commondir = args[1]
jaspardir = args[2]
id = as.double(args[3])
outdir = args[4]

if(file.exists(paste0(outdir,id,'.pwmout.RData'))){
  stop("pwm file already exists")
}

source(commondir)

####
# load PWMs
####

#pwmin = 'pwms/'


importJaspar <- function(file=myloc) {
  vec <- readLines(file)
  vec <- gsub("\\[|\\]", "", vec)
  start <- grep(">", vec); end <- grep(">", vec) - 1
  pos <- data.frame(start=start, end=c(end[-1], length(vec)))
  pwm <- sapply(seq(along=pos[,1]), function(x) vec[pos[x,1]:pos[x,2]])
  pwm <- sapply(seq(along=pwm), function(x) strsplit(pwm[[x]], " {1,}"))
  pwm <- lapply(seq(along=start), function(x) matrix(as.numeric(t(as.data.frame(pwm[(pos[x,1]+1):pos[x,2]]))[,-1]), nrow=4, dimnames=list(c("A", "C", "G", "T"), NULL)))
  names(pwm) <- gsub(">", "", vec[start])
  return(pwm)
}
pwmtable <- importJaspar(jaspardir)

pwmnum = id
pwmin = pwmtable[[pwmnum]]
pwmname = names(pwmtable)[pwmnum]

####
# end input script
# assert: existence of pwmin and pwmname
####


####
# motif match

pwmnorm=t(t(pwmin)/colSums(pwmin))
#informbase=colSums((log(pwmnorm+0.01)-log(1/4))*pwmnorm) #
#pwmnorm = pwmnorm[,(informbase > basecut)]
ipr=log(pwmnorm)-log(1/4)

if(!match.rc){
    pwuse = ipr
}else{
    pwuse = reverseComplement(ipr)
}

#chr names
chrstr = seqnames(genome)

coords.list = lapply(chrstr,function(i){
    print(i)
    mpwm=matchPWM(pwuse,genome[[i]],min.score=motifcut)
    pscore=PWMscoreStartingAt(pwuse,unmasked(genome[[i]]),start(mpwm))
    list(mpwm,pscore)
})

allpwm=do.call(c,lapply(coords.list,function(i){i[[2]]}))
pwmcut2=sort(allpwm,decreasing=T)[min(length(allpwm),maxcand)]
rm(allpwm)
print(pwmcut2)

coords=lapply(1:length(coords.list),function(i){
    as(coords.list[[i]][[1]],'IRanges')[coords.list[[i]][[2]] > pwmcut2]
})

coords.pwm=lapply(coords.list,function(i){i[[2]][i[[2]]>pwmcut2]})

#coords=lapply(coords.list,unlist)

clengths=sapply(coords,length)
coords.short=coords[clengths>0]
names(coords.short)=chrstr[clengths>0]
ncoords=chrstr[clengths>0]#names(coords)
coords2=sapply(coords.short,flank,width=wsize,both=T)

save(coords,coords.pwm,ipr,pwmin,pwmname,chrstr,clengths,coords.short,ncoords,coords2,file=paste0(outdir,id,'.pwmout.RData'))

}else{
clengths=0
save(clengths,file=paste0(outdir,id,'.pwmout.RData'))
}

#
#####
