#assert existence of
#commonfile
source(commonfile)
#pwmfile
load(pwmdir)
#tmpdir
#outdir


#####
# make bg

if(!fast.mode){
source('covfun.r')

load(paste0(tmpdir,'background.chr1.RData'))
require(Matrix)

bgwsize=10

datin = t(rBind(pos.mat[1:bgwsize,],neg.mat[1:bgwsize,]))
datin[datin>2]=2

mean(neg.mat)

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
                                        #
    musum=rep(0,ncol(covin))
    devsum=0
    muevsum=0
    covsumT=matrix(0,ncol(covin),ncol(covin))
    sqmuT=matrix(0,ncol(covin),ncol(covin))
    muT=rep(0,ncol(covin))
    tct=0
    numex=min(nrow(datin),10000)
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
w.both = wsize*2
mb=makeblocks(1:w.both,1)[1:bgwsize]
covmat.pos=toeptomat(spp[1:bgwsize],mb,w.both)
covmat.neg=toeptomat(spp[1:200 + (3*bgwsize)],mb,w.both)

scmat.pos=solve(covmat.pos)
scmat.neg=solve(covmat.neg)
}else{
    covmat.pos = diag(wsize*2)
    covmat.neg = diag(wsize*2)
    scmat.neg = diag(wsize*2)
    scmat.pos = diag(wsize*2)
}


#####
# make svlite


tfun.sparse <- function(x,offsets){
    tmp=do.call(rbind,offsets)
    x[tmp[,1:2]]=tfun(tmp[,3])
}

makesvlite <- function(filename,label,rot.pos,rot.neg,minread=5){
    print(filename)
    load(paste0(tmpdir,filename))
    if(!fast.mode){
        pm = t(rot.pos)%*%(pos.mat)
        nm = t(rot.neg)%*%(neg.mat)
    }else{
        pm = (pos.mat)
        nm = (neg.mat)
    }
    posind=apply(pm,2,function(j){
        i=j
        rbind(which(i!=0),i[i!=0])
    })
    negind=apply(nm,2,function(j){
        i=j
        rbind(which(i!=0),i[i!=0])
    })
    rct = (colSums(pos.mat>0) + colSums(neg.mat>0))
    sel = which(rct>minread)
    sapply(sel,function(i){
        tct=rct[i]+1
        if(tct > 0){
            paste0(c(label,
                     paste0(1,":",tct),
                     paste0(posind[[i]][1,]+1,':',posind[[i]][2,]),
                     paste0(negind[[i]][1,]+1+2*wsize,':',negind[[i]][2,])),collapse=' ')
        }else{
            paste0(c(label,"1:1"),collapse=' ')
        }
    })
}

validpos = list.files(tmpdir,'positive')
allpos=do.call(c,lapply(validpos,makesvlite,label=1,minread=10,rot.pos=covmat.pos,rot.neg=covmat.neg))

validneg = list.files(tmpdir,'background')
allneg=do.call(c,lapply(validneg,makesvlite,label= -1,minread= -1, rot.pos=covmat.pos,rot.neg=covmat.neg))

writeLines(c(allpos,allneg),paste0(tmpdir,'svlite.txt'))

#
#####

#####
# calc fig

sv.fit=sofia(paste0(tmpdir,'svlite.txt'),verbose=T,dimensionality=4*wsize+2,random_seed=1,lambda=1,iterations=1e+08)
sv.rotate = sv.fit$weights 

svv=sv.fit$weights[((1:(4*wsize))+2)]+sv.fit$weights[2]

pdf(paste0(outdir,'prof.pdf'))
plot(svv,type='l')
dev.off()

save(sv.fit,sv.rotate,svv,file=paste0(tmpdir,'svout.RData'))

#
#####
