%%%%
% new_76 - Test DCA, BDCA and BDCA adaptive for a constrained
% clustering problem on the EIL76 data set. This is used to generate 
% a single run for each, and then plot the resulting solution.
% 
% Constraints are given by 2 intersectings sets with ball centers given by 
% variable CB and radius given by variable RB, polygonal contraints
% are given by variable Z and the affine constraints are given by variables 
% U and BU.
% 
% See example 7.1 in the paper. 
%
%%%%

clear;
clc;
format long
%% load data set
A=dlmread('eil76.dat');
[m,n]=size(A);

%% projection function
C=[20 60;35 20;45 22];
[b,~]=size(C);
q=2;
r=[7,7,7];

% center lying in polygon
Z={[40,40;40,60;20,60;20,40]};
projfun=@(mat)proj(mat,C,r,Z);

%% initializaions
tau=1; % projection penalty
sig=10;% penalty growth parameter
tauf=1e8;
ep=1e-6;
k=2;%number of constraints for each center

%initial center
X0=[unifrnd(20,40,1,1) unifrnd(40,60,1,1);C(2,:)+r(2)*rand(1,2)/sqrt(2)]; 


%% DCA
tic
[X,N1,iterlogsdca]=constrainedDCA2V2(A,X0,projfun,tau,sig,tauf,q); %ep,q);
toc
[c,I]=cost(A,X);
dcaiter = 0;
for i = 1:N1
    dcaiter = dcaiter + iterlogsdca(i).dcaiter;
end
fprintf('Total Outer iterations using DCA= %d\nTotal DCA iterations %d \nFinal cost = %f\n\n',N1,dcaiter,c)

%% BDCA
tic
[X2,outeriter,iterlogsbdca]=constrainedBDCA2V2(A,X0,projfun,tau,sig,tauf,q,false);
toc
[c2,I2]=cost(A,X2);
bdcaiter = 0;
for i = 1:outeriter
    bdcaiter = bdcaiter + iterlogsbdca(i).dcaiter;
end
fprintf('Total Outer iterations using BDCA= %d \nTotal BDCA iterations %d \n Final cost = %f\n\n',outeriter,bdcaiter,c2);


%% BDCA Adaptive
% profile on
tic
[X3,outeriter_adap,iterlogsbdca_adap]=constrainedBDCA2V2(A,X0,projfun,tau,sig,tauf,q,true);
toc
[c3,I3]=cost(A,X3);
bdcaiter_adap = 0;
for i = 1:outeriter_adap
    bdcaiter_adap = bdcaiter_adap + iterlogsbdca_adap(i).dcaiter;
end
fprintf('Total Outer iterations using BDCA adaptive = %d \nTotal BDCA iterations %d \n Final cost = %f\n\n',outeriter_adap,bdcaiter_adap,c3);


%% plots (we can use X, X1 or X2 here)
fig=figure();clf
hold on
wid=2;
axis equal

for j=1:b
    draw_circle(C(j,1),C(j,2),r(j),'k',[]);
end
plotpoly(Z,'k')

plotclusters(A,X2,I2)
hold off
grid on

%% Projection function
function P=proj(X,C,r,Z)
[~,n]=size(X);
Pb=zeros(3,n);
Pb(1,:)=projball2(X(1,:),C(1,:),r(1));
Pb(2,:)=projball2(X(2,:),C(2,:),r(2));
Pb(3,:)=projball2(X(2,:),C(3,:),r(3));
Pp=projpolygon(X(1,:),Z);
P=[Pp(1,:);Pb(1,:);Pb(2,:);Pb(3,:)];
end