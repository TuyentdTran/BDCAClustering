function [X,iter,normsiters,costiters]=constrainedDCAV2(A,X,projfun,tau,q,costfun)
%   A        : (m,n) m data points in R^n.
%   X        : (k,n) initial centers in R^n.
%   projfun  : function handle for projection operator 
%             that returns a (k,n) array P
%   tau      : (>0) projection penalty parameter.
%   q        : Number of constraints per center
%   costfun  : Function to return the cost given a solution X.
%              Evaluated at every iteration. Enter [] if you don't 
%              want the added cost of computing this. 
[k,~]=size(X); m=size(A,1);
S=ones(k,m)*A;
flag=1;
iter=0;
normsiters = zeros(1000,1); 
costiters  = zeros(1000,1);
while flag
    X_old = X;
    W     = h2v2(X , A);
    U     = projfun(X);
    Y     = tau*U+W;
    X     = (Y+S)/(m+tau*q);
    
    if( isa(costfun, 'function_handle'))
        costiters(iter+1) = costfun(X,tau);
    end
    
    nrm   = norm(X_old-X,'fro'); 
    normsiters(iter+1) = nrm; 
    iter  = iter+1; 
    
    flag  = nrm >= 1e-6;
end

%remove unneeded entries in the normsiters array
normsiters  = normsiters(1:iter);
if( isa(costfun, 'function_handle'))
    costiters   = costiters(1:iter);
else
    costiters   = 0; 
end
end