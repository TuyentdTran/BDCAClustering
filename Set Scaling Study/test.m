%This file is to plot Example 7.5 in the paper in 2D
clear;clc;
rng('default') % For reproducibility
%% load data set
N= 5000;
C = unifrnd(0,10,N,2); %N pts
R= 0.1*ones(N,1);

CB=[1 5 ; 2 6 ; 5 3 ; 4 3; 8 5 ; 8 4; 9 8; 8 8]; %centers of constraint balls 
RB = [1 1 1 1 1 1 1 1]; % Radius of constraint balls

q=2; %#constraints per center

projfun=@(mat)projball2(mat,C,R);
proj_c =@(mat)proj(mat,CB,RB,q);
% Parameters
tau=1;
tauf=1e8;
sig=10;
k=4; %# of centers
n=size(C,2);

% Starting point

X0= zeros(k,n);
j = 1;
for i = 1:2:2*k
    X0(j,:) = CB(i,:)+RB(i)*rand(1,n)/sqrt(2);
    j       = j+1; 
end

%% DCA 
tic
[X,N1,iterlogsdca]=constrainedDCA2V2(C,R,X0,proj_c,tau,sig,tauf,q);
toc
[c,I]=cost(C,X,projfun);
dcaiter = 0;
for i = 1:N1
    dcaiter = dcaiter + iterlogsdca(i).dcaiter;
end
fprintf('Total Outer iterations using DCA= %d\nTotal DCA iterations %d \nFinal cost = %f\n\n',N1,dcaiter,c)

%%  BDCA
tic
[X2,outeriter,iterlogsbdca]=constrainedBDCA2V2(C,R,X0,proj_c,tau,sig,tauf,q,false);
toc
[c2,I2]=cost(C,X2,projfun);
bdcaiter = 0;
for i = 1:outeriter
    bdcaiter = bdcaiter + iterlogsbdca(i).dcaiter;
end
fprintf('Total Outer iterations using BDCA= %d \nTotal BDCA iterations %d \n Final cost = %f\n\n',outeriter,bdcaiter,c2);

%% Adaptive BDCA
tic
[X3,outeriter,iterlogsbdca]=constrainedBDCA2V2(C,R,X0,proj_c,tau,sig,tauf,q,true);
toc
[c2,I2]=cost(C,X3,projfun);
bdcaiter = 0;
for i = 1:outeriter
    bdcaiter = bdcaiter + iterlogsbdca(i).dcaiter;
end
fprintf('Total Outer iterations using adaptive BDCA= %d \nTotal adaptive BDCA iterations %d \n Final cost = %f\n\n',outeriter,bdcaiter,c2);


%% plots
figure(),clf,hold on

axis equal
plot_figurev2(C,R,CB,RB,X,I)
hold on


%% Projection function
function P=proj(X,CB,RB,q)
[k,n]=size(X);
Pb = zeros(k*q,n);
for j=1:k
    Pb(2*j-1:2*j,:) = projball2(X(j,:),CB(2*j-1:2*j,:),RB(2*j-1:2*j));
end
P=Pb;
end


