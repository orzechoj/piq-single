source('common.global.r')

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


######
# Optional: location of a whitelist directory. PWM matches will only be used if their window completely fits in the whitelist.
# whitelist should be in .bed format.

#whitelist = '/cluster/thashim/PIQ/capture.mm10.bed'
whitelist = NULL
