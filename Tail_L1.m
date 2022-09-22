function [xfinal] = Tail_L1(y, A, maxiter, tolerance)
%-----------------------------------------------------------------------------------------
% Input:
% y: observed signal;
% A: measurement matrix;
% maxiter: maximum number of iteration for convergence; 
% tolerance: tolerated number for convergence 
% S: number of nonzeros;
%-----------------------------------------------------------------------------------------
N=size(A,2);
M=size(A,1);
TC = zeros(N,1); %TC is the support of x..
iter = 0; %s is NOT the sparsity here, just an index variable so that we do not spend
%too long attempting tail minimization.
% Initial
cvx_begin quiet
   variable xs(N)
   minimize( norm(xs,1) )
   subject to
    A*xs==y;
cvx_end
x=xs;
xold=zeros(N,1);
% norm(x_old-x_new)<1e-6
while (iter <= maxiter) && (norm(xold-x) > tolerance)
    iter  =iter +1; 
    xold=x;
    [~,sorted_idx]=sort(abs(xold),'descend');
    Index=min(iter,M);
    Snew=sort(sorted_idx(1:Index));
    TC = ones(N,1);
    TC(Snew)=0;
     cvx_begin quiet
        variable x(N) 
        minimize(sum(TC.*abs(x)))
        subject to
        A*x == y;
    cvx_end  
end
xfinal = x;
