# setwd('C:\\Users\\Abdelmonem\\Dropbox\\TF\\BMF\\BMF')

source('topFiberM_exp.R')

Xl=list()
ds_names=c('DBLP','Mushroom','Firewall1','Chess','Paleo')
ds_names=sort(ds_names)
for(ds in ds_names){
    print(load(sprintf("data/%s.RData",ds)))
    Xl[[ds]]=get(ds)
}



library(Matrix)
ds='DBLP'
    r=5
    tP=1
        X=Xl[[ds]]
        A=as(X,'TsparseMatrix')
        Xb=X==1
    
 res1=SBMF_topFibers_exp(Xb,r=r,tP=tP,SR=100,verbose=3)
    X_=res1$A %*% res1$B
    X_=as(X_,'TsparseMatrix')
    li=A@i[A@x==1]+1
    lj=A@j[A@x==1]+1
    tp=sum(X_[cbind(li,lj)]>0)
    fn=sum(X)-tp#sum(!X_[cbind(li,lj)])
    fp=sum(X_@x>0)-tp
    cv=1-(fp+fn)/(tp+fn)
    