

#####
# Check before running

#match reverse complement (setting to T matches ONLY the R.C.)
match.rc = F
# load genome file
genome = Hsapiens
suppressMessages(bis("BSgenome.Hsapiens.UCSC.hg19"))

#####
# Run mode options

#if a file indicating completion is found; exit.
overwrite=F
#purity of calls
purity.cut = 0.7

#####
# Motif options

#motif score cutoff (log)
motifcut = 5
#max candidate sites (default 500k)
maxcand = 500000

#####
# Approx motif match options

#number of kmer samples to draw
nkmer = 5000000


