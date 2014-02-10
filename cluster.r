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

#use precomputed tables if wsize is 1000
if( (!fast.mode) & (wsize = 1000)){
    load('precalc/cov.premade.RData')
    tmat.pos = cov.premade[1:(2*wsize),]
    tmat.neg = cov.premade[((1:(2*wsize)) + (2*wsize)),]
}

#require at least 30k sites to do a background estimate
if((ncol(pos.mat) > 30000) & (!fast.mode) & (wsize != 1000)){

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

rot.pos = matrix(0,wsize*2,wsize*2)
rot.pos[lower.tri(rot.pos,diag=T)]=pv[1]/max(pv)
rot.neg = matrix(0,wsize*2,wsize*2)
rot.neg[lower.tri(rot.neg,diag=T)]=nv[1]/max(nv)

tmat.pos = cbind(rot.pos, matrix(0,2*wsize,2*wsize))
tmat.neg = cbind(matrix(0,2*wsize,2*wsize), rot.neg)

}else{
	rot.pos = diag(wsize*2)	
	rot.neg = diag(wsize*2)	
}



#####
# make svlite


makesvlite <- function(filename,label,rot.pos,rot.neg,minread=5){
    print(filename)
    load(filename)
    print('loaded')
    if(!fast.mode){
        pm = suppressMessages(t(rot.pos)%*%pos.mat)
        nm = suppressMessages(t(rot.neg)%*%neg.mat)
    }else{
        pm = pos.mat
        nm = neg.mat
    }
    print('findind')
    if(ncol(pm)>1){
    posind=lapply(1:ncol(pm),function(i){
        j=pm[,i]
        rbind(which(j!=0),j[j!=0])
    })
    }else{
    posind=list(rbind(which(pm[,1]!=0),pm[pm[,1]!=0,1]))
    }
    if(ncol(nm)>1){
    negind=lapply(1:ncol(nm),function(i){
        j=nm[,i]
        rbind(which(j!=0),j[j!=0])
    })
    }else{
    negind=list(rbind(which(nm[,1]!=0),nm[nm[,1]!=0,1]))	
    }
    print('runret')
    rct = (colSums(pos.mat>0) + colSums(neg.mat>0))
    sel = which(rct>minread)
    ret=sapply(sel,function(i){
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
    if(length(sel)==0){ret = character(0)}
    ret
}


if(fast.mode){
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

sv.fit=sofia(file.path(tmpdir,paste0(pwmid,'-svlite.txt')),verbose=T,dimensionality=4*wsize+2,random_seed=1,lambda=10000/length(allpos),iterations=5e+07,learner_type='logreg-pegasos',eta_type='basic')

vpos=suppressMessages(as.vector(rot.pos%*%(sv.fit$weights[2+(1:(2*wsize))])))
vneg=suppressMessages(as.vector(rot.neg%*%(sv.fit$weights[2+(1:(2*wsize))+(2*wsize)])))
sv.rotate = c(sv.fit$weights[1:2],vpos,vneg)

save(sv.fit,sv.rotate,file=file.path(outdir,paste0(pwmid,'.svout.RData')))


}else{

validpos = list.files(tmpdir,paste0('positive.tf',pwmid,'-'),full.names=T)
validneg = list.files(tmpdir,paste0('background.tf',pwmid,'-'),full.names=T)

allnom = c('label','rct',paste0('x',1:ncol(tmat.pos)))
alldat=rbind(
    do.call(rbind,lapply(validpos,function(i){
        print(i)
    load(i)
    rct = (colSums(pos.mat>0) + colSums(neg.mat>0))
    y=data.frame(label=1,rct[rct>10],x=t(as.matrix(t(tmat.pos)%*%pos.mat[,rct>10] + t(tmat.neg)%*%neg.mat[,rct>10])))
    colnames(y)=allnom
    y
})),
    do.call(rbind,lapply(validneg,function(i){
    print(i)
    load(i)
    rct = (colSums(pos.mat>0) + colSums(neg.mat>0))
    y=data.frame(label=-1,rct,x=t(as.matrix(t(tmat.pos)%*%pos.mat + t(tmat.neg)%*%neg.mat)))
    colnames(y)=allnom
    y
})))

write.svmlight(alldat[,1],as.matrix(alldat[,-1]),file.path(tmpdir,paste0(pwmid,'-svlite.txt')))

sv.fit=sofia(file.path(tmpdir,paste0(pwmid,'-svlite.txt')),verbose=T,dimensionality=ncol(alldat),random_seed=1,lambda=10000/nrow(alldat),iterations=5e+07,learner_type='logreg-pegasos',eta_type='basic')

sv.rotate=c(sv.fit$weights[1:2],tmat.pos%*%(sv.fit$weight[-(1:2)]),tmat.neg%*%(sv.fit$weight[-(1:2)]))

save(sv.fit,sv.rotate,file=file.path(outdir,paste0(pwmid,'.svout.RData')))

}


#
#####
