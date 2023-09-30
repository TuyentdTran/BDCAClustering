%%%%
% Testing_USCities_N - Test DCA, BDCA and BDCA adaptive for a constrained 
% set clustering problem for an adjustable number of top populous US cities, 
% where the sets are balls with radius porportional to population. This generates 3 
% variables 'LogsBDCA', 'LogsBDCA_adapt', and 'LogsDCA' which contain
% detailed information on all of the solves. These Log files are then saved
% to a file, the profile directory, which can then be be plotted be
% 'plotscript'
% 
% Constraints are given by 2 intersectings sets with ball centers given by 
% variable CB and radius given by variable RB, polygonal contraints
% are given by variable Z and the affine constraints are given by variables 
% U and BU. This is used to generate data for Example 7.4 in the paper. 
%
% The parameters for the script:
%
% N         - Number of times to rerun for a different random initial point.
% NumCities - Number of cities to use as set data. 
%
% Output: 
%
% LogsBDCA, LogsBDCA_adapt,LogsDCA   - 
% 
%       Log files containing detailed information about all of the
%       runs. A plotting script 'plotscript' can use these log files to 
%       generate figures of the same kind used in the paper. 
%%%%

clear;
clc;

N = 100; %Number of Runs
NumCities = 1500; 

% Parameters
tau  = 1;
tauf = 1e8;
sig  = 10;
k    = 4; %# of centers
n    = 2;
q    = 2; %#constraints per center

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

projfun=@(mat)projball2(mat,C,R);
proj_c =@(mat)proj(mat,CB,RB,Z,U,Bu);
X0 = zeros(k,n); 

for ii = 1:N
    % Starting point
    for i = 1:k-1
        X0(i,:) = CB(i,:)+RB(i)*rand(1,2)/sqrt(2);
    end
    X0(4,:) = CB(5,:)+RB(5)*rand(1,2)/sqrt(2);

    %% DCA 
    tic
    [X,N1,iterlogsdca]=constrainedDCA2V2(C,R,X0,proj_c,tau,sig,tauf,q);
    LogsDCA(N-ii+1).time = toc;
    LogsDCA(N-ii+1).logs = iterlogsdca;
    [c,I]=cost(C,X,projfun);
    dcaiter = 0;
    for i = 1:N1
        dcaiter = dcaiter + iterlogsdca(i).dcaiter;
    end
    LogsDCA(N-ii+1).totaliter = dcaiter;
    LogsDCA(N-ii+1).outeriter = N1;
    LogsDCA(N-ii+1).cost      = c;

    %% BDCA 
    tic
    [X2,outeriter,iterlogsbdca]=constrainedBDCA2V2(C,R,X0,proj_c,tau,sig,tauf,q,false);
    LogsBDCA(N-ii+1).time = toc;
    LogsBDCA(N-ii+1).logs = iterlogsbdca;
    [c2,I2]=cost(C,X2,projfun);
    bdcaiter = 0;
    for i = 1:outeriter
        bdcaiter = bdcaiter + iterlogsbdca(i).dcaiter;
    end
    LogsBDCA(N-ii+1).totaliter = bdcaiter; 
    LogsBDCA(N-ii+1).outeriter = outeriter; 
    LogsBDCA(N-ii+1).cost      = c2;

    %% BDCA adaptive
    tic
    [X3,outeriter_adap,iterlogsbdca_adap]=constrainedBDCA2V2(C,R,X0,proj_c,tau,sig,tauf,q,true);
    LogsBDCA_adapt(N-ii+1).time = toc;
    LogsBDCA_adapt(N-ii+1).logs = iterlogsbdca_adap;
    [c3,I3]=cost(C,X3,projfun);
    bdcaiter_adap = 0;
    for i = 1:outeriter_adap
        bdcaiter_adap = bdcaiter_adap + iterlogsbdca_adap(i).dcaiter;
    end
    LogsBDCA_adapt(N-ii+1).totaliter = bdcaiter_adap; 
    LogsBDCA_adapt(N-ii+1).outeriter = outeriter_adap; 
    LogsBDCA_adapt(N-ii+1).cost      = c3;
end
filename = CreateUniqueFileName('profile/Results');
save(filename,'LogsBDCA_adapt','LogsDCA','LogsBDCA');
%%Plots
plotscript;

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

function[FileName]=CreateUniqueFileName(FileName)
[fPath, fName, fExt] = fileparts(FileName);
if isempty(fExt)  % No '.mat' in FileName
  fExt     = '.mat';
  FileName = fullfile(fPath, [fName, fExt]);
end
if exist(FileName,'file')
    [fPath, fName, fExt] = fileparts(FileName);
    fDir = dir(fullfile(fPath, [fName,' (*)', fExt]));
    if isempty(fDir)
        FileName = fullfile(fPath, [fName,' (1)', fExt]);
    else
        pattern = "(" + digitsPattern + ")" + fExt;
        hasThePattern = endsWith(extractfield(fDir,'name'),pattern);
        Extracted = extract(extractfield(fDir(hasThePattern),'name'),pattern);
        num = max(cell2mat(cellfun(@(C) textscan(C,'(%d)') , Extracted,'UniformOutput',true)));
        num = num+1;
        FileName = fullfile(fPath, [fName,' (',num2str(num),')', fExt]);
    end
end
end