# expand fiber

# 7/1/2019

#  Copyright (C) 2018 Abdelmoneim Amer Desouki, 
#   Data Science Group, Paderborn University, Germany.
#  All right reserved.
#  Email: desouki@mail.upb.de
#
#  This file is part of topFiberM.
#
#  topFiberM is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  topFiberM is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with topFiberM.  If not, see <http://www.gnu.org/licenses/>.

SBMF_topFibers_exp <- function(X,r=2,tP=0.5,verbose=2,SR=NULL){
	#top fibers but with excluding what is already taken in previous steps.
    # expand fibers instantly and eliminate Factors when better ones are found use a search limit SR.
    
    expand_col<-function(cix){
            Ac=X1[,cix]# use X1 in case more sparse factors is required and avoid nnz>1
			
            ctp=colSums(X1[Ac,,drop=FALSE])
		    cfp=colSums(!X[Ac,,drop=FALSE])
		    Bc = ifelse(ctp==0,FALSE,((ctp/(ctp+cfp)) >=tP))
            rtp=rowSums(X1[,Bc,drop=FALSE])
		    rfp=rowSums(!X[,Bc,drop=FALSE])
		    Ac = ifelse(rtp==0,FALSE,((rtp/(rtp+rfp)) >=tP))
            
			return(list(Ai=Ac,Bi=Bc))#,TP1=TP1,FP1=FP1
    }
	n=nrow(X)
	m=ncol(X)
    if(is.null(SR)) SR=r;
	As=Matrix::spMatrix(x=FALSE,i=1,j=1,nrow=n,ncol=min(SR,m,n))
	Bs=Matrix::spMatrix(x=FALSE,i=1,j=1,nrow=min(SR,m,n),ncol=m)
	X1=X
    Ai=rep(FALSE,n)
    Bi=rep(FALSE,m)
    
	i=1
	tf=NULL
    cv=NULL
	TP=0;
	FP=0;
    excluded_cols=rep(FALSE,m)
    excluded_rows=rep(FALSE,n)
    
	while(i<=min(SR,m,n) ){
		cs=colSums(X1)
		rs=rowSums(X1)
		if(sum(rs)==0 ||(sum(rs)==1 && i>r)) break;
        if(sum(excluded_rows)==n && sum(excluded_cols)==m) break;
		mxrv=max(rs[!excluded_rows])
		mxcv=max(cs[!excluded_cols])
        mxr=which(rs==mxrv & (!excluded_rows))[1]
        mxc=which(cs==mxcv & (!excluded_cols))[1]
		if(verbose>2) print(paste('mxr:',mxr ,'mxc:',mxc))
    	if(!is.na(mxc) && cs[mxc]>=rs[mxr]){
	         tmp=expand_col(mxc)
            
            Ai=tmp$Ai
            Bi=tmp$Bi
            ix=as.matrix(expand.grid(which(Ai),which(Bi)),ncol=2)
			TP1=sum(X1[ix])#gain depends on uncovered
            FP1=sum(!X[ix])
            tf1=cbind(i=i,f=2,ix=mxc,nnz=nrow(ix),gain=TP1-FP1,TP=TP1,FP=FP1)
		}else{# if a row is better
		    Bi=X1[mxr,]# use X1 in case more sparse factors is required and avoid nnz>1
			
            rtp=rowSums(X1[,Bi,drop=FALSE])
		    rfp=rowSums(!X[,Bi,drop=FALSE])#better to be calculated on X
		    Ai = ifelse(rtp==0,FALSE,((rtp/(rtp+rfp)) >=tP))
            ctp=colSums(X1[Ai,,drop=FALSE])#revise B
		    cfp=colSums(!X[Ai,,drop=FALSE])
		    Bi = ifelse(ctp==0,FALSE,((ctp/(ctp+cfp)) >=tP))

			ix=as.matrix(expand.grid(which(Ai),which(Bi)),ncol=2)
            TP1=sum(X1[ix])
            FP1=sum(!X[ix])

            tf1=cbind(i=i,f=1,ix=mxc,nnz=nrow(ix),gain=TP1-FP1,TP=TP1,FP=FP1)   
			}
            
            As[,i]=Ai
            Bs[i,]=Bi
            if(i>r){#replace min gain
            # excluded rows/excluded columns
                mgI=which.min(tf[,'gain'])
                if(tf[mgI,'gain']>=tf1[,'gain']){
                    if(tf1[,'f']==1){
                        excluded_cols[mxc]=TRUE
                    }else{
                        excluded_rows[mxr]=TRUE
                    }
                }else{#eliminate fiber
                    tfmg=tf[mgI,,drop=FALSE]
                    tf[mgI,]=tf1
                    print(sprintf('replacing one factor..old gain=%d, new gain=%d',tfmg[,'gain'],tf1[,'gain']))
                    ## reevaluate tf1 11/1/2019
                    ix2=as.matrix(expand.grid(which(As[,tfmg[,'i']]),which(Bs[tfmg[,'i'],])),ncol=2)
                    print(paste(sum(X1),'ix2:',nrow(ix2)))
                    X1[ix2]=X[ix2]# restore fiber with min gain
                    print(sum(X1))
                    if(tf1[,'f']==2){
                        tmp1=expand_col(tf1[,'ix'])
                        Ai=tmp1$Ai
                        Bi=tmp1$Bi
                        ix=as.matrix(expand.grid(which(Ai),which(Bi)),ncol=2)

                        TP1=sum(X1[ix])#gain depends on uncovered
                        FP1=sum(!X[ix])
                        tf1=cbind(i=i,f=2,ix=tf1[,'ix'],nnz=nrow(ix),gain=TP1-FP1,TP=TP1,FP=FP1)
                        print(sprintf('replacing one factor..old gain=%d, new ext gain=%d',tf[mgI,'gain'],tf1[,'gain']))
                        tf[mgI,]=tf1
                        X1[ix]=FALSE
                    }
                }
            }else{
                X1[ix]=FALSE
                tf=rbind(tf,tf1)
            }
            if(verbose>2) print(tf1)
            i=i+1;
	}
    A=As[,tf[,'i'],drop=FALSE]
    B=Bs[tf[,'i'],,drop=FALSE]
	return(list(A=A,B=B,X1=X1,tf=tf))
}

########-----------------------------------------------------------------
