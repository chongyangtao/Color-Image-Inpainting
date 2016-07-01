function [A]=OMP_Inpainting(D,X,Mask,sigma,rc_min, max_coeff)
% Author: Chongyang
% Date: 2016/06/05
% Perform sparse coding using a orthogonal matching pursuit(inpainting)
%
% Input
% D: (n x K) unit norm atoms
% X: (n x P) observations
% Mask: (n x P) mask denoting which observations are unknown
% sigma: residual error stopping criterion, normalized by signal norm
% rc_min: minimal residual correlation before stopping
% max_coeff: maximal number of non-zero coefficients for signal
%
% Output
% A: OMP coding



[n,P]=size(X); % 64 * 552108
[n,K]=size(D);
A = sparse(size(D,2),size(X,2));
for k=1:P
    if(mod(k,10000)==1)
        tic
    end
  
    Mpos=find(Mask(:,k));
    Dict=D(Mpos,:);
    W=1./sqrt(diag(Dict'*Dict)); 
    Dict=Dict*diag(W);
    x=X(Mpos,k); % patch exclude mask
    
    residual=x;
    indx = [];
    a = [];
    DD = [];
    currResNorm2 =  norm(residual);
    j = 0;
    rc_max = inf;
    threshold = norm(x)*sigma;
    while (norm(residual)> threshold  && rc_max > rc_min && sum(a~=0) < max_coeff )
        j = j+1;
        proj=Dict'*residual;
        [rc_max pos] = max(abs(proj));
        DD = [DD Dict(:,pos)];
        indx = [indx, pos(1)];
        a=pinv(DD)*x;
        residual=x-DD*a;
    end;
    if (~isempty(indx))
        A(indx,k)=a;
        A(:,k)=W.*A(:,k);
    end
    if(mod(k,10000)==0)
        fprintf('    iter:  %d/%d | (%2.2f%%) |%d \n',k,P,100*(k/P),toc);
    end
end;
return;