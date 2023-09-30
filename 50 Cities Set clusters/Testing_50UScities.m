%%%%
% Testing50cities - Test DCA, BDCA and BDCA adaptive for a constrained 
% set clustering problem for the Top 50 populous US cities, where the sets 
% are balls with radius porportional to population. This is used to
% generate a single run for each, and then plot the resulting solution.
% 
% Constraints are given by 2 intersectings balls with centers given by 
% variable C and radius given by variable R. This is used to generate data
% for figure of Example 7.3 in the paper. 
%
%%%%

clear;
clc;

%% load data set
US = load('top50citySets.mat');
C  = [US.lon, US.lat];
R  = (0.1/sqrt(3.14))*sqrt(US.SqMi);


CB = [-80 34; -80 38; -92 37; -90 40; -115 45; -110 40] ;            
RB = [2 3 4 3 4 4];

q=2; %#constraints per center

projfun = @(mat)projball2(mat,C,R);
proj_c  = @(mat)proj(mat,CB,RB,q);
% Parameters
tau=1;
tauf=1e8;
sig=10;
k=3; %# of centers
n=2;

% Starting point

X0= zeros(k,n);
for i = 1:k
    X0(i,:)=CB(i,:)+RB(i)*rand(1,2)/sqrt(2);
end

%% DCA

tic
[X,N1,iterlogsdca]=constrainedDCA2V2(C,R,X0,proj_c,tau,sig,tauf,q); %ep,q);
toc
[c,I]=cost(C,X,projfun);
dcaiter = 0;
for i = 1:N1
    dcaiter = dcaiter + iterlogsdca(i).dcaiter;
end
fprintf('Total Outer iterations using DCA = %d\nTotal DCA iterations %d \nFinal cost = %f\n\n',N1,dcaiter,c)

%% BDCA
tic
[X2,outeriter,iterlogsbdca]=constrainedBDCA2V2(C,R,X0,proj_c,tau,sig,tauf,q,false);
toc
[c2,I2]=cost(C,X2,projfun);
bdcaiter = 0;
for i = 1:outeriter
    bdcaiter = bdcaiter + iterlogsbdca(i).dcaiter;
end
fprintf('Total Outer iterations using BDCA = %d \nTotal BDCA iterations %d \n Final cost = %f\n\n',outeriter,bdcaiter,c2);

%% BDCA Adaptive
tic
[X3,outeriter,iterlogsbdca_adapt]=constrainedBDCA2V2(C,R,X0,proj_c,tau,sig,tauf,q,true);
toc
[c3,I3]=cost(C,X3,projfun);
bdcaiter_adapt = 0;
for i = 1:outeriter
    bdcaiter_adapt = bdcaiter_adapt + iterlogsbdca_adapt(i).dcaiter;
end
fprintf('Total Outer iterations using BDCA Adaptive = %d \nTotal BDCA iterations %d \n Final cost = %f\n\n',outeriter,bdcaiter_adapt,c3);
%plotting 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load plotStatesData
clf; 
shg
start=1;
hold on
for i=1:48     % for each state in US
    numState=nnz(stateData(:,3)==i);  % count number of points in that state
    x=stateData(start:start+numState-1,2);  % collect x values in that state
    y=stateData(start:start+numState-1,1);  % collect y values in that state
    plot(x, y, 'color', [0.7 0.7 0.7]); % [0.7 0.7 0.7]--> gray
    start=start+numState;               % shift index to next state in list
end

axis equal
axis([-135 -65 15 55]);
plot_figurev2(C,R,CB,RB,X,I)%either X,X2 or X3 can be here
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


