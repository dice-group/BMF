setwd('C:\\Users\\Abdelmonem\\Dropbox\\TF\\BMF\\BMF')

source("GreCond.R")
Xl=list()
ds_names=c('DBLP','Mushroom','Firewall1','Chess','Paleo')
ds_names=sort(ds_names)
for(ds in ds_names){
    print(load(sprintf("data/%s.RData",ds)))
    Xl[[ds]]=get(ds)
}

library(Matrix)
all_Res=list()
all_r=c(1,2,5,10)
    ds='Mushroom'
    Res=NULL
    print(ds)
        X=Xl[[ds]]
        A=as(X,'TsparseMatrix')
        Xb=X==1
    for(r in c(1,2,5,10)){
        t0=proc.time()
        res1=GreConD(X,r)
        t2=proc.time()
            X_=res1$A %*% res1$B
            X_=as(X_,'TsparseMatrix')
            li=A@i[A@x==1]+1
            lj=A@j[A@x==1]+1
            tp=sum(X_[cbind(li,lj)]>0)
            fn=sum(X)-tp#sum(!X_[cbind(li,lj)])
            fp=sum(X_@x>0)-tp
            Res=rbind(Res,cbind(r,tp,fn,fp,ctime=(t2-t0)[3]))

}