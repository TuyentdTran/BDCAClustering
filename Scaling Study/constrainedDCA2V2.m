function [X,iter,iterlogs]=constrainedDCA2V2(A,X,projfun,tau,sig,tauf,q,costfun)
%   A        : (m,n) m data points in R^n.
%   X        : (k,n) initial centers in R^n.
%   projfun  : function handle for projection operator 
%             that returns a (k,n) array P
%   tau      : (>0) projection penalty parameter.
%   q        : Number of constraints per center
%   costfun  : Function to return the cost given a solution X.
%              Evaluated at every iteration. Enter [] if you don't 
%              want the added cost of computing this. 
flag1    = 1;
iterlogs = struct;
iter     = 0;
while tau<tauf && flag1
    [Xp,dcaiter,dcanorms,cost] = constrainedDCAV2(A,X,projfun,tau,q,costfun);
    nrm   = norm(Xp-X,'fro')/(norm(X,'fro')+1);
    flag1 = nrm>=10e-6;
    X     = Xp;
    tau   = sig*tau;
    iter  = iter+1;
    iterlogs(iter).dcaiter    = dcaiter;
    iterlogs(iter).tau        = tau;
    iterlogs(iter).outernorms = nrm;
    iterlogs(iter).dcanorms   = dcanorms;
    iterlogs(iter).costfuns   = cost;
end