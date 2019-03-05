#3/1/2019
## this code is converted from MATLAB code provided by the authers

 GreConD<-function( I, no_of_factors =NULL){#[ A, B ] 
# % GRECOND implements GreConD algorithm for Boolean matrix factorization 

# % usage: [A, B] = GreConD(I);
# % returns A \circ B = I (if the no. of factors is not limited)

# % if you are using this implementation please cite the following work
# % Belohlavek R., Vychodil V.: 
# % Discovery of optimal factors in binary data via a novel method of matrix decomposition. 
# % Journal of Computer and System Sciences 76(1)(2010), 3-20. 

 M = I==1#logical(I);
# [m, n] = size(M);
m=nrow(M)
n=ncol(M)
U = M;
k = 0;

# A = logical([]);
# B = logical([]);
A=Matrix::spMatrix(x=FALSE,i=1,j=1,nrow=m,ncol=1)
B=Matrix::spMatrix(x=FALSE,i=1,j=1,nrow=1,ncol=n)

while (sum(U)>0){#any(any(U))
    v = 0;
    d =rep(FALSE,n)#,false(1,n);
    d_old = d;
    d_mid = d#false(1,n);
    e = matrix(rep(TRUE,m),ncol=1)#true(m,1); #% extent for speed closure
    
    atr = which(colSums(U)>0); #% only not covered attributes
    
    while( 1==1){
        # for (j in which(!d[atr])){
            if(!d[j]){#~d(j)
                # % computes the value of the cover function for the candidate factor
                # % inline function for speed
                # % arrow down (speed version)
                a = as.vector(e & M[,j]);
                # % arrow up
                # b = all(M(a,:),1);
                
                # if(sum(a)==0){
                    # b=rep(TRUE,n)#as MATLAB all, return all ones
                # }else{
                  b = colSums(M[a,,drop=FALSE])==sum(a);
                 # }
                # % coverage
                # cost = sum(sum(U(a,b)));
                cost = sum(U[a,b]);
                # %
                
                if (cost > v){
                    v = cost;
                    d_mid = b;
                    cc = a;
                }
            #}
        }
        
        d = d_mid;
        e = cc;
        
        if (all(d==d_old)){
            break;
        }else{
            d_old = d;
        }
    }
    
    k = k + 1;
    print(k);
    # browser()
    # A = [A, cc];
    A[,k]=cc
    # B = [B; d];
    B[k,]=d;
    
    
    # % end if the no. of factors is reached
    if (!is.null(no_of_factors) && k==no_of_factors){
        break;
    }else{
        A=cbind(A,FALSE)
        B=rbind(B,FALSE)
    }
    
    # % delete already covered part
    U[cc, d] = FALSE;
}
return(list(A=A,B=B))
}