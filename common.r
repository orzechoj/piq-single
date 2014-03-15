#####
# utils

old.repos <- getOption("repos")
on.exit(options(repos = old.repos))
new.repos <- old.repos
new.repos["CRAN"] <- "http://cran.stat.ucla.edu"
options(repos = new.repos)

#source("http://www.bioconductor.org/biocLite.R")
if(!suppressMessages(require("BiocInstaller",quietly=T))){
    install.packages("BiocInstaller", repos="http://www.bioconductor.org/packages/2.13/bioc")
}

ris <- function(x){if(!require(x,character.only=T,quietly=T,warn.conflicts=F)){
    install.packages(x)
    require(x,character.only=T,quietly=T,warn.conflicts=F)
}}
bis <- function(x){if(!require(x,character.only=T,quietly=T,warn.conflicts=F)){
    biocLite(x)
    require(x,character.only=T,quietly=T,warn.conflicts=F)
}}

wipetemp <- function(){
    x <- readline("this will wipe existing temp files (y/n)")
    if(x=="y"){
        print("wiping temp")
        unlink("tmp/*")
    }else{
        print("keeping temp")
    }
}

set.seed(1)

####
# dependencies

suppressMessages(require('BiocInstaller',quietly=T))
suppressMessages(bis("BSgenome.Hsapiens.UCSC.hg19"))
suppressMessages(bis("Rsamtools"))
suppressMessages(bis('Biostrings'))

suppressMessages(ris('RSofia'))
suppressMessages(ris('statmod'))
suppressMessages(ris('Rcpp'))
suppressMessages(ris('inline'))
suppressMessages(require(Matrix,quietly=T))




####
# params

#if a file indicating completion is found; exit.
overwrite=F

#number of kmer samples to draw
nkmer = 5000000
#motif score cutoff (log)
motifcut = 5
#motif informativeness cutoff (nats)
basecut = 0
#max candidate sites (default 500k)
maxcand = 500000
#match reverse complement (setting to T matches ONLY the R.C.)
match.rc = T


#window size
wsize = 1000

#fast mode
fast.mode = T

#purity of calls
purity.cut = 0.7

# load genome file

genome = Hsapiens

