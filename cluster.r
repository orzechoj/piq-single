#assert existence of
#commonfile
source(commonfile)
#pwmfile
#load(file.path(pwmdir,paste0(pwmid,'.pwmout.RData')))
#tmpdir
#outdir


#####
# make bg

load(file.path(tmpdir,paste0('background.tf',pwmid,'-chr1.RData')))

#require at least 30k sites to do a background estimate
if((ncol(pos.mat) > 30000) & (!fast.mode)){

source('covfun.r')

require(Matrix)

bgwsize=2
alldat=do.call(rbind,lapply(seq(1,ncol(pos.mat),by=max(bgwsize,ncol(pos.mat)/500)),function(offs){
	datin = t(rBind(pos.mat[1:bgwsize+offs,],neg.mat[1:bgwsize+offs,]))^2
	datin[datin>2]=2
	as.matrix(datin)
}))

datin=alldat

stabval=4
mustab=-3
sz=ncol(datin)
bind=makeblocks(1:bgwsize,2)
tabs = makeTable(datin,mustab,stabval)
matin=matrix(tabs[[1]][as.double(datin)+1],nrow(datin))
covbase = cov(matin)
covin=topiter(covbase,bind)
sci=solve(covin)
muin= mustab
mupri=rep(muin,sz)%*%sci
scid = solve(sci+diag(rep(tabs[[2]][1],sz)))


for(m in 1:5){
    print(m)
    musum=rep(0,ncol(covin))
    devsum=0
    muevsum=0
    covsumT=matrix(0,ncol(covin),ncol(covin))
    sqmuT=matrix(0,ncol(covin),ncol(covin))
    muT=rep(0,ncol(covin))
    tct=0
    numex=min(nrow(datin),50000)
    for(i in sample(1:nrow(datin),numex)){
        #print(i)
        ff=fastFit(datin[i,],mupri,tabs[[1]],tabs[[2]],scid)
        muall= t(ff[[1]])
        covsum = ff[[2]]
        devs = 1/(1/diag(ff[[2]])+tabs[[2]][datin[i,]+1])
        muevs= (ff[[1]]/diag(ff[[2]])-tabs[[1]][datin[i,]+1]*tabs[[2]][datin[i,]+1])*devs
        #
        sqmuT=sqmuT+t(muall)%*%muall
        muT=muT+colSums(muall)
        covsumT=covsumT+covsum
        devsum=devsum+sum(devs)
        muevsum=muevsum+sum(muevs)
        tct=tct+nrow(muall)
    }
    cch=sqmuT/tct-(t(t(muT))%*%t(muT))/(tct^2)
    cpost=(covsumT/tct+cch)
    cpost=topiter(cpost,bind)
    muhat=muT/tct
    stabval=devsum/(tct*ncol(covin))
    muin=muevsum/(tct*ncol(covin))
    print(c(stabval,muin))
    tabs = makeTable(datin,muin,stabval)
    print(tabs)
    covin=cpost
    sci=solve(covin)
    mupri=muhat%*%sci
    scid = solve(sci+diag(rep(tabs[[2]][1],sz)))
}

spp = toepvals(cpost,bind)
w.both = 100
mb=makeblocks(1:w.both,1)[1:bgwsize]
covmat.pos=solve(toeptomat(spp[1:bgwsize],mb,w.both))
covmat.neg=solve(toeptomat(spp[1:200 + (3*bgwsize)],mb,w.both))

pv=toepvals(covmat.pos,mb)
nv=toepvals(covmat.neg,mb)

rot.pos=toeptomat(pv/pv[1],mb,2*wsize)
rot.neg=toeptomat(nv/nv[1],mb,2*wsize)

}else{
	rot.pos = diag(wsize*2)	
	rot.neg = diag(wsize*2)	
}


#####
# make svlite




makesvlite <- function(filename,label,rot.pos,rot.neg,minread=5){
    print(filename)
    load(filename)
    if(!fast.mode){
        pm = t(rot.pos)%*%pos.mat
        nm = t(rot.neg)%*%neg.mat
    }else{
        pm = (pos.mat)
        nm = (neg.mat)
    }
    posind=apply(pm,2,function(j){
        rbind(which(j!=0),j[j!=0])
    })
    negind=apply(nm,2,function(j){
        rbind(which(j!=0),j[j!=0])
    })
    rct = (colSums(pos.mat>0) + colSums(neg.mat>0))
    sel = which(rct>minread)
    sapply(sel,function(i){
        tct=rct[i]+1
	if(tct > 1){
	pid = posind[[i]][1,]
	nid = negind[[i]][1,]
	pval = posind[[i]][2,]
	nval = negind[[i]][2,]
        if(length(pid) > 0){
	cpos = paste0(pid+1,':',pval)
	}else{
	cpos = character(0)
	}
	if(length(nid) > 0){
	cneg = paste0(nid+1+2*wsize,':',nval)
	}else{
	cneg = character(0)
	}}else{
	cpos = character(0)
	cneg = character(0)
	}
        paste0(c(label,paste0(1,":",tct),cpos,cneg),collapse=' ')
    })
}

validpos = list.files(tmpdir,paste0('positive.tf',pwmid,'-'),full.names=T)
vptrain=validpos
allpos=do.call(c,lapply(vptrain,makesvlite,label=1,minread=10,rot.pos=rot.pos,rot.neg=rot.neg))

validneg = list.files(tmpdir,paste0('background.tf',pwmid,'-'),full.names=T)
vntrain=validneg
allneg=do.call(c,lapply(vntrain,makesvlite,label= -1,minread= -1, rot.pos=rot.pos,rot.neg=rot.neg))

writeLines(c(allpos,allneg),file.path(tmpdir,paste0(pwmid,'-svlite.txt')))

#
#####

#####
# calc fig

sv.fit=sofia(file.path(tmpdir,paste0(pwmid,'-svlite.txt')),verbose=T,dimensionality=4*wsize+2,random_seed=1,lambda=2,iterations=5e+07,learner_type='logreg-pegasos')
sv.rotate = sv.fit$weights 

save(sv.fit,sv.rotate,file=file.path(tmpdir,paste0(pwmid,'.svout.RData')))

#
#####
