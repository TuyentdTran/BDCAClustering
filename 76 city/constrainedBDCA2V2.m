function [X,iter,iterlogs]=constrainedBDCA2V2(A,X,projfun,tau,sig,tauf,q,useselfadaptive)
%   A               : (m,n) m data points in R^n.
%   X               : (k,n) initial centers in R^n.
%   projfun         : function handle for projection operator 
%                     that returns a (k,n) array P
%   tau             : Initial penalty parameter
%   sig             : (>1) penalty multiplier.
%   tauf            : (>tau) final penalty.
%   q               : Number of constraints per center
%   useselfadaptive : Flag to turn on self-adaptivity for lambda in BDCA 
flag1    = 1; 
iter     = 0;
iterlogs = struct;

while tau<tauf && flag1
    [Xp,bdcaiter,searchiter,lambdaiter,dcanorms] = constrainedBDCAV2(A,X,projfun,tau,q,useselfadaptive);
    nrm   = norm(Xp-X,'fro')/(norm(X,'fro')+1);
    flag1 = nrm>=10e-6;
    X     = Xp;
    tau   = sig*tau;
    iter  = iter+1;
    iterlogs(iter).dcaiter    = bdcaiter; 
    iterlogs(iter).searchiter = searchiter;
    iterlogs(iter).lambdaiter = lambdaiter;
    iterlogs(iter).dcanorms   = dcanorms;
    iterlogs(iter).tau        = tau;
    iterlogs(iter).outernorm  = nrm;
end
end
