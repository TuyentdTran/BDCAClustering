function [X,iter,normsiters]=constrainedDCAV2(A,X,projfun,tau,q)
%   A        : (m,n) m data points in R^n.
%   X        : (k,n) initial centers in R^n.
%   projfun  : function handle for projection operator 
%             that returns a (k,n) array P
%   tau      : (>0) projection penalty parameter.
%   q        : Number of constraints per center
[k,n] = size(X); 
m     = size(A,1);
S     = ones(k,m)*A;
U     = zeros(k,n);
flag=1;
iter=0;
normsiters = zeros(1000,1); 
while flag
    X_old  = X;
    W      = h2v2(X , A);
    temp   = projfun(X); 
    U(1,:) = temp(1,:)+temp(2,:); % hard coded assumption that k = 2 here
    U(2,:) = temp(3,:)+temp(4,:);
    Y      = tau*U+W;
    X      = (Y+S)/(m+tau*q);
    
    nrm   = norm(X_old-X,'fro'); 
    normsiters(iter+1) = nrm; 
    iter  = iter+1; 
    
    flag  = nrm >= 1e-6;
end

%remove unneeded entries in the normsiters array
normsiters  = normsiters(1:iter);
end
