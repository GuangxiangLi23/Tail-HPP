function [x]=TailHPPRe(y,A,lambda,v,tolerance,maxiter)

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
    u=(A'*A.*(v*v')+lambda/2*omega)\(v.*(A'*y));
    v=(A'*A.*(u*u')+lambda/2*omega)\(u.*(A'*y));
    x=u.*v;
end