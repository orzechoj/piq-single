#####
# utils

old.repos <- getOption("repos")
on.exit(options(repos = old.repos))
new.repos <- old.repos
new.repos["CRAN"] <- "http://cran.stat.ucla.edu"
options(repos = new.repos)


.libPaths(libpath)
.Library = libpath
.Library.site=''
#source("http://www.bioconductor.org/biocLite.R")
install.packages("BiocInstaller", repos="http://www.bioconductor.org/packages/2.13/bioc",lib=libpath)

ris <- function(x){if(!require(x,character.only=T,lib.loc=libpath)){
    install.packages(x,lib=libpath)
    require(x,lib.loc=libpath,character.only=T)
}}
bis <- function(x){if(!require(x,character.only=T,lib.loc=libpath)){
    biocLite(x,lib=libpath,lib.loc=libpath)
    require(x,lib.loc=libpath,character.only=T)
}}

