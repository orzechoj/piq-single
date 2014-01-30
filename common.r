#####
# utils

old.repos <- getOption("repos")
on.exit(options(repos = old.repos))
new.repos <- old.repos
new.repos["CRAN"] <- "http://cran.stat.ucla.edu"
options(repos = new.repos)

#source("http://www.bioconductor.org/biocLite.R")
if(!require("BiocInstaller")){
    install.packages("BiocInstaller", repos="http://www.bioconductor.org/packages/2.13/bioc")
}

ris <- function(x){if(!require(x,character.only=T)){
    install.packages(x)
    require(x,character.only=T)
}}
bis <- function(x){if(!require(x,character.only=T)){
    biocLite(x)
    require(x,character.only=T)
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
require('BiocInstaller')
bis("BSgenome.Mmusculus.UCSC.mm10")
bis("Rsamtools")
bis('Biostrings')

ris('plotrix')
ris('RSofia')
ris('statmod')
require(Matrix)

####
# params

#number of kmer samples to draw
nkmer = 5000000
#motif score cutoff (log)
motifcut = 5
#motif informativeness cutoff (nats)
basecut = 0
#max candidate sites (default 500k)
maxcand = 500000


#window size
wsize = 1000

#fast mode
fast.mode = T

# load genome file

genome = Mmusculus


