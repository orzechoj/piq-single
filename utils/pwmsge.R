ei=readPWMfile('pwms/motifsused.txt')

ei=sapply(rp,function(i){
    exp(i)-min(exp(i))
})
exportAsTransfacFile(ei,'pwms/motifused.cts.txt')

require('snow')

exists = as.double(sapply(strsplit(list.files('/cluster/thashim/tmppiq/',pattern='pwm'),'[.]'),function(i){i[1]}))

cl <- makeCluster(16)
cls=clusterApplyLB(cl,(1:1316)[-exists],function(i){
    print(i)
    system(paste0('./pwmmatch.r /cluster/thashim/basepiq/common.mm.r /cluster/thashim/basepiq/pwms/jasparfix.txt ',i,' /cluster/thashim/tmppiq/'))
})
#system('echo \"./pwmmatch.r /cluster/thashim/basepiq/common.r /cluster/thashim/basepiq/pwms/jasparfix.txt 1 /cluster/thashim/tmppiq/ \" | qsub')
