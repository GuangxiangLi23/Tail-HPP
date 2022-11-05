function [x]=HPP(y,A,tolerance,maxiter,omega,lambda)
[~,nn]=size(A);
v = ones(nn,1);

x=A\y;
xold=zeros(nn,1);
iter=0;
while (iter <= maxiter) && (norm(xold-x) > tolerance)
    xold=x;
    u=(A'*A.*(v*v')+lambda/2*omega)\(v.*(A'*y));
    v=(A'*A.*(u*u')+lambda/2*omega)\(u.*(A'*y));
    x=u.*v;
    iter=iter+1;
end
