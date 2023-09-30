function [X,iter,normsiters]=constrainedDCAV2(C,R,X,proj_c,tau,q)
%   C       : (m,n) m data points in R^n.
%   R       : (k,1) Radius of the data points
%   X       : (k,n) initial centers in R^n.
%   proj_c  : function handle for projection operator 
%             that returns a (k,n) array P
%   tau     : (>0) projection penalty parameter.
%   q       : Number of constraints per center
[k,n]=size(X); m=size(C,1);
U=zeros(k,n);
flag=1;
iter=0;
normsiters = zeros(1000,1); 
while flag
    X_old = X;
    W     = h1(X,C,R);
    %W = EX_W_new(X,C,R);
    temp  = proj_c(X); 
    for j=1:k
        U(j,:) = sum(temp(2*j-1:2*j,:),1);
    end
    X   = (1/(m+tau*q))*(m*X+tau*U-W); 
    nrm = norm(X_old-X,'fro'); 
    normsiters(iter+1) = nrm; 
    iter = iter+1;
        
    flag  = nrm >= 1e-6;
end
%remove unneeded entries in the normsiters array
normsiters  = normsiters(1:iter);
end