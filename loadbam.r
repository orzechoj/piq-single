#assert existence of
#commonfile
source(commonfile)
#bamfile
load(bamfile)
#pwmfile
load(pwmdir)
#tmpdir 

tfun <- function(x){
    y = x
    y[x>0] = 2*sqrt(x[x>0] + 3/8)
    y
}

makeTFmatrix <- function(coords,prefix='',offset=0){
    cwidth = width(coords[[1]][1])
    obschrnames=levels(c(seqnames(plusstrand),seqnames(minusstrand)))
    validchr = obschrnames[which(obschrnames%in%ncoords)]
    for(chr in validchr){
        chrcoord=shift(coords[[chr]],offset)
        print(chr)
        pluscoord=start(plusstrand[seqnames(plusstrand)==chr])
        minuscoord=start(minusstrand[seqnames(minusstrand)==chr])
        irp=IRanges(start=pluscoord,width=1)
        fos=findOverlaps(chrcoord,irp)
        pos.unique.hits = unique(queryHits(fos))
        pos.offset=pluscoord[subjectHits(fos)]-start(chrcoord)[queryHits(fos)]+1
        ubd=findInterval(c(0,pos.unique.hits),queryHits(fos))
        pos.offsets=lapply(1:(length(pos.unique.hits)),function(id){
            reads=pos.offset[(ubd[id]+1):ubd[id+1]]
            reads.rle=rle(reads)
            cbind(pos.unique.hits[id],reads.rle$values,tfun(reads.rle$lengths))
        })
        pos.matrix=do.call(rbind,pos.offsets)
    #
        irn=IRanges(start=minuscoord,width=1)
        fos=findOverlaps(chrcoord,irn)
        neg.unique.hits = unique(queryHits(fos))
        neg.offset=minuscoord[subjectHits(fos)]-start(chrcoord)[queryHits(fos)]+1
        ubd=findInterval(c(0,neg.unique.hits),queryHits(fos))
        neg.offsets=lapply(1:(length(neg.unique.hits)),function(id){
            reads=neg.offset[(ubd[id]+1):ubd[id+1]]
            reads.rle=rle(reads)
            cbind(neg.unique.hits[id],reads.rle$values,tfun(reads.rle$lengths))
        })
        neg.matrix=do.call(rbind,neg.offsets)
#
        pos.mat = sparseMatrix(i=pos.matrix[,2],j=pos.matrix[,1],x=pos.matrix[,3],dims=c(cwidth,length(chrcoord)))
        neg.mat = sparseMatrix(i=neg.matrix[,2],j=neg.matrix[,1],x=neg.matrix[,3],dims=c(cwidth,length(chrcoord)))
        save(pos.mat,neg.mat,pos.offsets,neg.offsets,file=paste0(tmpdir,prefix,chr,'.RData'))
    }
}

makeTFmatrix(coords2,'positive.')
makeTFmatrix(coords2,'background.',10000)

#
#####
