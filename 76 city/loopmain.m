function [LogsDCA,LogsBDCA,LogsBDCA_adapt] = loopmain(N)
% Run a test N times

% Logfiles describing the runs


%% load data set
A = dlmread('eil76.dat');
%% projection function
%constraint balls
C     = [20 60;35 20;45 22];
q = 2;
r = [7,7,7];
% center lying in polygon
Z = {[40,40;40,60;20,60;20,40]};
projfun = @(mat) projhardcoded(mat,C,r,Z);

%% initializaions
tau  = 1; % projection penalty
sig  = 10;% penalty growth parameter
tauf = 1e8;
k    = 2; %number of constraints per center, assumed 2 elsewhere

for runiter = 1:N          
    X0 = [unifrnd(20,40,1,1) unifrnd(40,60,1,1);C(2,:)+r(2)*rand(1,2)/sqrt(2)];
    %% DCA
    tic
    [X,N1,iterlogsdca]        = constrainedDCA2V2(A,X0,projfun,tau,sig,tauf,q);
    LogsDCA(N-runiter+1).time = toc;
    LogsDCA(N-runiter+1).logs = iterlogsdca;
    c = cost(A,X);
    dcaiter = 0;
    for i = 1:N1
        dcaiter = dcaiter + iterlogsdca(i).dcaiter;
    end
    LogsDCA(N-runiter+1).totaliter = dcaiter;
    LogsDCA(N-runiter+1).outeriter = N1;
    LogsDCA(N-runiter+1).cost      = c;
    %% BDCA -- Adaptive
    tic
    [X2,outeriter,iterlogsbdca_adap]=constrainedBDCA2V2(A,X0,projfun,tau,sig,tauf,q,true);
    LogsBDCA_adapt(N-runiter+1).time = toc;
    LogsBDCA_adapt(N-runiter+1).logs      = iterlogsbdca_adap;
    c2 = cost(A,X2);
    bdcaiter_adap = 0;
    for i = 1:outeriter
        bdcaiter_adap = bdcaiter_adap + iterlogsbdca_adap(i).dcaiter;
    end
    LogsBDCA_adapt(N-runiter+1).totaliter = bdcaiter_adap; 
    LogsBDCA_adapt(N-runiter+1).outeriter = outeriter; 
    LogsBDCA_adapt(N-runiter+1).cost      = c2;
    
    %% BDCA -- Non Adaptive
    tic
    [X3,outeriter,iterlogsbdca]=constrainedBDCA2V2(A,X0,projfun,tau,sig,tauf,q,false);
    LogsBDCA(N-runiter+1).time = toc;
    LogsBDCA(N-runiter+1).logs      = iterlogsbdca;
    c3 = cost(A,X3);
    bdcaiter = 0;
    for i = 1:outeriter
        bdcaiter = bdcaiter + iterlogsbdca(i).dcaiter;
    end
    LogsBDCA(N-runiter+1).totaliter = bdcaiter; 
    LogsBDCA(N-runiter+1).outeriter = outeriter; 
    LogsBDCA(N-runiter+1).cost      = c3;
end

end

%% Projection function
function P=projhardcoded(X,C,r,Z)
    [~,n]=size(X);
    Pb=zeros(3,n);
    Pb(1,:)=projball2(X(1,:),C(1,:),r(1));
    Pb(2,:)=projball2(X(2,:),C(2,:),r(2));
    Pb(3,:)=projball2(X(2,:),C(3,:),r(3));
    Pp=projpolygon(X(1,:),Z);
    P=[Pp(1,:);Pb(1,:);Pb(2,:);Pb(3,:)];
end

