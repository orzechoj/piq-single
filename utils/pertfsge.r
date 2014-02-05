require('snow')

cl <- makeCluster(16)
cls=clusterApplyLB(cl,1:1316,function(i){
    print(i)
    system(paste0('./pertf.r /cluster/thashim/basepiq/common.mm.r /cluster/thashim/tmppiq/ /cluster/thashim/tmppiq/ /cluster/thashim/130130.mm10.d0/ /cluster/thashim/tmppiq/d0.RData ',i))
})
