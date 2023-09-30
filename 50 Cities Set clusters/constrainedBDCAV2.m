function [X,iter,searchiters,lambdaiters,normsiters]=constrainedBDCAV2(C,R,X,proj_c,tau,q,useselfadaptive)
%   C               : (m,n) m data points in R^n.
%   R               : (k,1) Radius of the data points
%   X               : (k,n) initial centers in R^n.
%   proj_c          : function handle for projection operator that returns a (k,n) array P
%   tau             : (>0) projection penalty parameter.
%   q               : Number of constraints per center
%   useselfadaptive : Flag to turn on self-adaptivity for lambda in BDCA 

%   alpha,beta and lambda are parameters that will be used in BDCA and set
%   inside of the function.

[k,n] = size(X); 
m     = size(C,1);
U     = zeros(k,n);

alp         = 0.05;
beta        = 0.1;
lambdastart = 2;
lambdahistorylength = 2;
gamma       = 2;

flag = 1;
iter = 0;
searchiters = zeros(1000,1);
lambdaiters = zeros(1000,1);
normsiters  = zeros(1000,1);
lambdahist  = zeros(lambdahistorylength,1);
lambdabar   = zeros(lambdahistorylength,1);
lambda      = lambdastart;

while flag
    
    X_old = X;
    W     = h1(X,C,R);

    temp = proj_c(X); 
    for j = 1:k
        U(j,:) = sum(temp(2*j-1:2*j,:),1);
    end
    Xk  = (1/(m+tau*q))*(m*X+tau*U-W); 
    d_k = Xk - X_old;
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
    
    X = Xk + (d_k * lambda);
    phi_k = cost_pen_opt(C,R,Xk,tau,proj_c);
    costpen = cost_pen_opt(C,R,X,tau,proj_c); 
    rhs = alp * norm(d_k, 'fro')^2;
    d = phi_k - lambda^2 * rhs;
    searchiter = 0; 
    while costpen > d
            lambda = lambda * beta;
            X = Xk + (d_k * lambda);
            costpen = cost_pen_opt(C,R,X,tau,proj_c);
            d = phi_k - lambda^2 * rhs;
            searchiter = searchiter + 1;
    end
    searchiters(iter+1) = searchiter;
    lambdaiters(iter+1) = lambda;

    nrm  = norm(X_old-X,'fro');
    normsiters(iter+1) = nrm;
    iter = iter+1;

end
%remove unneeded entries in the searchiters array
searchiters = searchiters(1:iter);
lambdaiters = lambdaiters(1:iter);
normsiters  = normsiters(1:iter);
end