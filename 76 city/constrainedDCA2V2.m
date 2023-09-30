function [X,iter,iterlogs]=constrainedDCA2V2(A,X,projfun,tau,sig,tauf,q)
%   A        : (m,n) m data points in R^n.
%   X        : (k,n) initial centers in R^n.
%   projfun  : function handle for projection operator 
%             that returns a (k,n) array P
%   tau      : (>0) projection penalty parameter.
%   q        : Number of constraints per center
flag1    = 1;
iterlogs = struct;
iter     = 0;
while tau<tauf && flag1
    [Xp,dcaiter,dcanorms] = constrainedDCAV2(A,X,projfun,tau,q);
    nrm   = norm(Xp-X,'fro')/(norm(X,'fro')+1);
    flag1 = nrm>=10e-6;
    X     = Xp;
    tau   = sig*tau;
    iter  = iter+1;
    iterlogs(iter).dcaiter    = dcaiter;
    iterlogs(iter).tau        = tau;
    iterlogs(iter).outernorms = nrm;
    iterlogs(iter).dcanorms   = dcanorms;
end