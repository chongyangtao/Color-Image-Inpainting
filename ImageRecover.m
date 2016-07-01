function [out]=ImageRecover(y,D,CoefMatrix)


% method:  Sliding
% CoefMatrix  256*146405
% D 64*256
    
[N M]=size(y); 
n=sqrt(size(D,1)); 
out=zeros(N,M); 
weight=zeros(N,M);
    
i=1; j=1;
for k=1:(N-n+1)*(M-n+1)
    patch=reshape(D*CoefMatrix(:,k),[n,n]); 
    out(i:i+n-1,j:j+n-1)=out(i:i+n-1,j:j+n-1)+patch; 
    weight(i:i+n-1,j:j+n-1)=weight(i:i+n-1,j:j+n-1)+1; 
    if i<N-n+1 
        i=i+1; 
    else
        i=1; j=j+1; 
    end;
end;
out=out./weight; 
return;
