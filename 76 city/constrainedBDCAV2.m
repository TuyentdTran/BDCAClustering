function [X,iter,searchiters,lambdaiters,normsiters]=constrainedBDCAV2(A,X,projfun,tau,q,useselfadaptive)
%   A               : (m,n) m data points in R^n.
%   X               : (k,n) initial centers in R^n.
%   projfun         : function handle for projection operator 
%                     that returns a (k,n) array P
%   tau             : (>0) projection penalty parameter.
%   q               : Number of constraints per center
%   useselfadaptive : Flag to turn on self-adaptivity for lambda in BDCA 

%   alpha,beta and lambda are parameters that will be used in BDCA and set
%   inside of the function.

[k,n]       = size(X); 
m           = size(A,1);
alp         = 0.05;
beta        = 0.1;
lambdastart = 1;
gamma       = 2;
lambdahistorylength = 2;

S=ones(k,m)*A;
flag1 = 1;
iter  = 0;
searchiters  = zeros(1000,1);
lambdaiters  = zeros(1000,1);
normsiters   = zeros(1000,1);
lambdahist   = zeros(lambdahistorylength,1);
lambdabar    = zeros(lambdahistorylength,1);
lambda       = lambdastart;
    
U = zeros(k,n);

while flag1
    X_old  = X;
    W      = h2v2(X , A);
    temp   = projfun(X); 
    U(1,:) = temp(1,:)+temp(2,:); % hard coded assumption that k = 2 here
    U(2,:) = temp(3,:)+temp(4,:);
    Y      = tau*U+W;
    Xk     = (Y+S)/(m+tau*q); %DCA
    d_k    = Xk - X_old;
    
    if norm(d_k,'fro')<1e-5
        X = Xk;
        break;
    end
    
    if(useselfadaptive)
        %renewhistory
        lambdahist(1:end-1) = lambdahist(2:end);
        lambdahist(end)     = lambda;
        
        if(all(abs(lambdabar - lambdahist)<1e-6))
            %no search was done for the last lambdahistorylength iterations
            %increase lambda from previous        
            lambda = gamma*lambdahist(end); 
        else % use previous lambda
            lambda = lambdahist(end);
        end
        %renew history 
        lambdabar(1:end-1)  = lambdabar(2:end);
        lambdabar(end)      = lambda;
    else
        lambda = lambdastart;
    end
    
    X       = Xk + (d_k * lambda);
    phi_k   = cost_pen_opt(A,Xk,tau,q,projfun);
    costpen = cost_pen_opt(A,X,tau,q,projfun); 
    rhs     = alp * norm(d_k, 'fro')^2;
    d       = phi_k - lambda^2 * rhs;
    searchiter = 0; 
    while costpen > d
            lambda     = lambda * beta;
            X          = Xk + (d_k * lambda);
            costpen    = cost_pen_opt(A,X,tau,q,projfun);
            d          = phi_k - lambda^2 * rhs;
            searchiter = searchiter + 1;
    end
    searchiters(iter+1) = searchiter;
    lambdaiters(iter+1) = lambda;
    
    nrm  = norm(X_old-X,'fro');
    normsiters(iter+1) = nrm;
    iter = iter+1;
    
    flag1 = nrm>=1e-6;
end
%remove unneeded entries in the searchiters array
searchiters = searchiters(1:iter);
lambdaiters = lambdaiters(1:iter);
normsiters  = normsiters(1:iter);
end
