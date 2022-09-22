function [x] = HPPOr(y,A,lambda,v,tolerance,maxiter)
%-----------------------------------------------------------------------------------------
% This function solve stand L1 by Hadamard method
% Input:
% y: observed signal;
% A: measurement matrix;
% lambda: regularized parameter; e.g., 1e-2, 1e-4, 1e-6;
% v: inital value for v;
% tolerance: tolerated number for convergence 
% maxiter: maximum number of iteration for convergence; 
%-----------------------------------------------------------------------------------------
[~,nn]=size(A);
omega=diag(ones(nn,1));
x=A\y;
xold=zeros(nn,1);
iter=0;
while (iter <=maxiter) && (norm(xold-x) > tolerance)
    xold=x;
    u=(A'*A.*(v*v')+lambda/2*omega)\(v.*(A'*y));
    v=(A'*A.*(u*u')+lambda/2*omega)\(u.*(A'*y));
    x=u.*v;
    iter=iter+1;
end