#assert existence of
#commonfile
source(commonfile)
#pwmfile
#load(file.path(pwmdir,paste0(pwmid,'.pwmout.RData')))
#tmpdir
#outdir


#####
# make bg

load(file.path(tmpdir,paste0('background.tf',pwmid,'-',seqnames(genome)[1],'.RData')))

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

converter.src <- '
Rcpp::NumericMatrix triplePos(rtriplePos);
Rcpp::NumericMatrix tripleNeg(rtripleNeg);
Rcpp::NumericVector sums(rsum);
int label = as<int>(rlabel);
int negoff = as<int>(rnegoff);
Rcpp::CharacterVector cv(sums.size());
int tripleindPos = 0;
int tripleindNeg = 0;
for(int i=0; i < sums.size(); i++){
  char str[65536];
  double sumct = sums[i];
  int cp = snprintf(str, 65536, "%d %d:%f", label, 1, sumct);
  while(tripleindPos < triplePos.nrow() && triplePos(tripleindPos,0)-1 == i){
    int crd = triplePos(tripleindPos,1)+1;
    double wt = triplePos(tripleindPos,2);
    cp += snprintf(str + cp, 65536-cp, " %d:%f", crd, wt);
    tripleindPos++;
  }
  while(tripleindNeg < tripleNeg.nrow() && tripleNeg(tripleindNeg,0)-1 == i){
    int crd = tripleNeg(tripleindNeg,1)+1+negoff;
    double wt = tripleNeg(tripleindNeg,2);
    cp += snprintf(str + cp, 65536-cp, " %d:%f", crd, wt);
    tripleindNeg++;
  }
  cv[i] = str;
}
return cv;
'
converter <- cxxfunction(signature(rtriplePos="numeric",rtripleNeg="numeric",rsum="numeric",rnegoff="integer",rlabel="integer"),converter.src,plugin="Rcpp",includes="#include <stdio.h>")

makesvlite <- function(filename,label,rot.pos,rot.neg,minread=5){
    print(filename)
    load(filename)
    rct=floor(colSums(pos.mat)+colSums(neg.mat))+1
    cv=converter(pos.triple,neg.triple,rct,2*wsize,label)
    cv[rct>minread]
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

itmax=min(5e+07,500*length(allpos))

sv.fit=sofia(file.path(tmpdir,paste0(pwmid,'-svlite.txt')),verbose=T,dimensionality=4*wsize+2,random_seed=1,lambda=10000/length(allpos),iterations=itmax,learner_type='logreg-pegasos',eta_type='basic',loop_type='balanced-stochastic')

vpos=suppressMessages(as.vector(rot.pos%*%(sv.fit$weights[2+(1:(2*wsize))])))
vneg=suppressMessages(as.vector(rot.neg%*%(sv.fit$weights[2+(1:(2*wsize))+(2*wsize)])))
sv.rotate = c(sv.fit$weights[1:2],vpos,vneg)

save(sv.fit,sv.rotate,file=file.path(tmpdir,paste0(pwmid,'.svout.RData')))


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

save(sv.fit,sv.rotate,file=file.path(tmpdir,paste0(pwmid,'.svout.RData')))

}


#
#####
