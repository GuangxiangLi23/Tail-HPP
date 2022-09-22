function [x]=TailHPPDi(y,A,maxiter)
tolerance = 1e-5;
[mm,nn]=size(A);

x=A\y;
xold=zeros(nn,1);
iter=0;
while (iter <=maxiter) && (norm(xold-x) > tolerance)
    iter=iter+1;
    xold=x;
    [~,sorted_idx]=sort(abs(xold),'descend');
    Index=min(iter,mm);
    Snew=sort(sorted_idx(1:Index));
    weights = ones(nn,1);
    weights(Snew)=0;
    omega=diag(weights);
    x= HPP(y,A,tolerance,maxiter,omega);
end