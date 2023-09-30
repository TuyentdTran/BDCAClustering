%This is the code to visualize Example 7.2 in 2d.
%The number of points can be modified.
clear;
clc;
format long
close all

numpoints = 10000; %Number of points

%% load data set
A     = unifrnd(0,10,numpoints,2);
[m,n] = size(A);

%% projection function
C     = [1 5;6 4;8 8]; % constraint centers
[b,~] = size(C);
q     = 1;             % number of overalapping constraints 
r     = [1, 1, 1];     % constraint radius 
projfun = @(mat)proj(mat,C,r);

%% initializaions
tau= 1; % projection penalty
sig= 10; % penalty growth parameter
tauf=1e8;
ep=1e-6;
k= 3;%Number of centers

X=zeros(k,n);
    for i =1:k
    X(i,:)=C(i,:)+r(i)*rand(1,2)/sqrt(2);
    end
X0=X;

%% DCA
tic
[X,N1,iterlogsdca]=constrainedDCA2V2(A,X0,projfun,tau,sig,tauf,q,[]); %ep,q);
toc
[c,I]=cost(A,X);
dcaiter = 0;
for i = 1:N1
    dcaiter = dcaiter + iterlogsdca(i).dcaiter;
end
fprintf('Total Outer iterations using DCA= %d\nTotal DCA iterations %d \nFinal cost = %f\n\n',N1,dcaiter,c)

%% BDCA

tic
[X2,outeriter,iterlogsbdca]=constrainedBDCA2V2(A,X0,projfun,tau,sig,tauf,q,false,[]);
toc

[c2,I2]=cost(A,X2);
bdcaiter = 0;
for i = 1:outeriter
    bdcaiter = bdcaiter + iterlogsbdca(i).dcaiter;
end
fprintf('Total Outer iterations using BDCA= %d \nTotal BDCA iterations %d \n Final cost = %f\n\n',outeriter,bdcaiter,c2);

%% BDCA Adaptive

tic
[X3,outeriter,iterlogsbdca_adap]=constrainedBDCA2V2(A,X0,projfun,tau,sig,tauf,q,true,[]);
toc

[c3,I3]=cost(A,X3);
bdcaiter_adap = 0;
for i = 1:outeriter
    bdcaiter_adap = bdcaiter_adap + iterlogsbdca_adap(i).dcaiter;
end
fprintf('Total Outer iterations using Adaptive BDCA= %d \nTotal Adaptive BDCA iterations %d \n Final cost = %f\n\n',outeriter,bdcaiter,c2);

%% plots
fig=figure();clf
wid=2;
hold on
axis equal

plotclusters(A,X2,I2)
for i=1:k
    draw_circle(C(i,1),C(i,2),r(i),'black',[]);
end

% Plot centers
plot(X(:,1),X(:,2),'rh','MarkerFaceColor','r','MarkerSize',20)
% %initial guess
% for i = 1:k
% plot(X0(i,1),X0(i,2),'kp','MarkerFaceColor','k','MarkerSize',15)
% plot([X0(i,1),X2(i,1)],[X0(i,2),X2(i,2)],'k-.','Linewidth',1)
% end
hold off

%% Projection function
function P=proj(X,C,r)
[k,n]=size(X);
P = zeros(k,n);
for i = 1:k
    P(i,:) = projball2(X(i,:),C(i,:),r(i)); 
end
end