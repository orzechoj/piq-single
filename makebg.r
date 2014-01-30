
#####
# make bg

if(!fast.mode){
source('covfun.r')

load('tmp/background.chr1.RData')
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


#
#####
