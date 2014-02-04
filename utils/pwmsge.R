require('snow')

####
# This script will match PWMs in parallel on one machine.
# 
# The core matching script is pwmmatch and is a randomized, approximate k-mer based matching, parameters for the run will be set in COMMON_DOT_R. 
####

#Number of cores to use at once on this machine
NUM_CORES = 16

#Location of common.r file
COMMON_DOT_R = '/cluster/thashim/basepiq/common.mm.r'

#Location of JASPAR format motif file
MOTIF_FILE = '/cluster/thashim/basepiq/pwms/jasparfix.txt'

#Location of TMP dir, this is where the PWM match outputs will go.
TMP_DIR = '/cluster/thashim/tmppiq/'

#Number of motifs in JASPAR file. If less than total number of motifs, only 1 to MAX_MOTIF_COUNT will be matched.
MAX_MOTIF_COUNT = 1316

cl <- makeCluster(NUM_CORES)
cls=clusterApplyLB(cl,(1:MAX_MOTIF_COUNT),function(i){
    print(i)
    system(paste0('./pwmmatch.r ',COMMON_DOT_R,' ',MOTIF_FILE,' ',i,' ',TMP_DIR))
})
