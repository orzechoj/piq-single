#####
# Below is an example of how to do DNase-seq calls using a sge based cluster
#####

#Location of the common.r file that does package load and defines parameters.
commonfile="/cluster/thashim/basepiq/common.r"

#Jaspar file directory path
jaspardir="/cluster/thashim/basepiq/pwms/jasparfix.txt"

#Directory in which the PWM files should be dumped
outdir="/cluster/thashim/tmppiqhg.rc/"

#Number of motifs to process 
for pwmid in {1..1316}
do
	echo "/cluster/thashim/basepiq/pwmmatch.r "$commonfile $jaspardir $pwmid $outdir | /usr/bin/qsub -wd /cluster/thashim/basepiq/ -e $outdir/errfile.txt -o /dev/null
done