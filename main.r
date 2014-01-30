###
# main runner script


source('common.r')
wipetemp()
if(!file.exists('tmp/pwmout.RData')){
    source('pwmmatch.r')
}
if(!file.exists('tmp/positive.chr1.RData')){
    source('loadbam.r')
}
source('makebg.r')
source('cluster.r')
source('bindcall.r')

