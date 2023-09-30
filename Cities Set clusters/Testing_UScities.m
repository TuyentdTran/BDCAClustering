%%%%
% Testing_UsCities - Test DCA, BDCA and BDCA adaptive for a constrained 
% set clustering problem for or an adjustable number of top populous US 
% cities, where the sets are balls with radius porportional to population. 
% This is used to generate a single run for each, and then plot the 
% resulting solution.
% 
% Constraints are given by 2 intersectings sets with ball centers given by 
% variable CB and radius given by variable RB, polygonal contraints
% are given by variable Z and the affine constraints are given by variables 
% U and BU.
% 
% See example 7.4 in the paper. 
%
% The parameters for the script:
%
% NumCities - Number of cities to use as set data. 
%
%%%%
clear;
clc;

NumCities = 1500; 

%% load data set
US = readtable('uscities.csv'); 
% remove Hawaii Alaska and Puerto Rico, to give only mainland data
todelete1= (US.state_name == "Puerto Rico");
US(todelete1,:) = [];
todelete2= (US.state_name ==  "Hawaii");
US(todelete2,:) = [];
todelete3= (US.state_name ==  "Alaska");
US(todelete3,:) = [];

US = sortrows(US,'population','descend');

C=[US.lng(1:NumCities), US.lat(1:NumCities)];
R=(1e-3/sqrt(3.14))*sqrt(US.population(1:NumCities));


% centers within r(i) degrees lat/long of cities, i=1,..,4
CB=[	-116.6873596	43.6629384;%Caldwell,Idaho
    -104.7902 41.135;% Cheyenne, WY
    -87.6866  41.8375;% Chicago,Illinois
    -90.2451  38.6359; %St. Louis, MO
     -77.0363658	38.8951118]; %Washington,District of Columbia
RB=[4 2.5 3 4 4];

% centers lying in polygons
Z={[-115,42;-115,49;-125,49;-125,42],...%first rectangle constraint
    [-102.05,37;-102.05,41;-109.05,41;-109.05,37]};% state of colorado

% center lying west of -75 lon
U=[1,0];Bu=[-75,0];

q=2; %#constraints per center

projfun=@(mat)projball2(mat,C,R);
proj_c =@(mat)proj(mat,CB,RB,Z,U,Bu);
% Parameters
tau=1;
tauf=1e8;
sig=10;
k=4; %# of centers
n=2;

% Starting point
X0= zeros(k,n);
for i = 1:k-1
    X0(i,:)=CB(i,:)+RB(i)*rand(1,2)/sqrt(2);
end
X0(4,:)=CB(5,:)+RB(5)*rand(1,2)/sqrt(2);

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

%% BDCA 
tic
[X2,outeriter,iterlogsbdca]=constrainedBDCA2V2(C,R,X0,proj_c,tau,sig,tauf,q,false);
toc
[c2,I2]=cost(C,X2,projfun);
bdcaiter = 0;
for i = 1:outeriter
    bdcaiter = bdcaiter + iterlogsbdca(i).dcaiter;
end
fprintf('Total Outer iterations using BDCA= %d \nTotal BDCA iterations %d \n Final cost = %f\n\n',outeriter,bdcaiter,c2);

%% BDCA adaptive
tic
[X3,outeriter_adap,iterlogsbdca_adap]=constrainedBDCA2V2(C,R,X0,proj_c,tau,sig,tauf,q,true);
toc
[c3,I3]=cost(C,X3,projfun);
bdcaiter_adap = 0;
for i = 1:outeriter_adap
    bdcaiter_adap = bdcaiter_adap + iterlogsbdca_adap(i).dcaiter;
end
fprintf('Total Outer iterations using BDCA adaptive= %d \nTotal BDCA adaptive iterations %d \n Final cost = %f\n\n',outeriter_adap,bdcaiter_adap,c3);


%% plots
fig=figure();clf
wid=2;

ax=[-126,-66,23.7,51];
axis(ax);

hold on

plot_figurev2(C,R,CB,RB,X,I)%any X,X1 or X2 can be used here
hold on
plotpart(U,Bu,'k')
hold on
plotpoly(Z,'k')
hold on
borders('continental us','Color',[0.7 0.7 0.7],'LineWidth',1,'nomap')
hold off
axis(ax);
pbaspect([ax(2)-ax(1) ax(4)-ax(3) 1])
grid on

%% Projection function

function P=proj(X,CB,RB,Z,U,Bu)
Pb=zeros(5,2);
for i =1:3
    Pb(i,:) = projball2(X(i,:),CB(i,:),RB(i));
end
Pb(4,:)=projball2(X(3,:),CB(4,:),RB(4));
Pb(5,:)=projball2(X(4,:),CB(5,:),RB(5));
Q=X(1:2,:);
Pp=projpolygon(Q,Z);
Pt=projAffineR(X(4,:),U,Bu); 
P=[Pp(1,:);Pb(1,:);Pp(2,:);Pb(2,:);Pb(3,:);Pb(4,:);Pb(5,:);Pt];
end


