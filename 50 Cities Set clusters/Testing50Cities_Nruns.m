%%%%
% Testing50Cities_Nruns - Test DCA, BDCA and BDCA adaptive for a constrained 
% set clustering problem for the Top 50 populous US cities, where the sets 
% are balls with radius porportional to population. This generates 3 
% variables 'LogsBDCA', 'LogsBDCA_adapt', and 'LogsDCA' which contain
% detailed information on all of the solves. These Log files are then saved
% to a file, the profile directory, which can then be be plotted by
% 'plotscript'
% 
% Constraints are given by 2 intersectings balls with centers given by 
% variable C and radius given by variable R. This is used to generate data
% for Example 7.3 in the paper. 
%
% The parameters for the script
%
% N   - Number of times to rerun for a different random initial point.
%
% Output: 
%
% LogsBDCA, LogsBDCA_adapt, LogsDCA   - 
% 
%       Log files containing detailed information about all of the
%       runs. A plotting script 'plotscript' can use these log files to 
%       generate figures of the same kind used in the paper. 
%%%%

clear;
clc;

N = 100; %Number of Runs

% Parameters
tau  = 1;
tauf = 1e8;
sig  = 10;
k    = 3; %# of centers
n    = 2;
q    = 2; %#constraints per center

%% load data set

US=load('top50citySets.mat');
C=[US.lon, US.lat];
R=(0.1/sqrt(3.14))*sqrt(US.SqMi);

CB=[-80 34; -80 38; -92 37; -90 40; -115 45; -110 40] ;            
RB=[2 3 4 3 4 4];

projfun=@(mat)projball2(mat,C,R);
proj_c =@(mat)proj(mat,CB,RB,q);
X0 = zeros(k,n); 

for ii = 1:N
  % Starting point
    for i = 1:k
        X0(i,:)=CB(i,:)+RB(i)*rand(1,2)/sqrt(2);
    end

    %% DCA 
    tic
    [X,N1,iterlogsdca]=constrainedDCA2V2(C,R,X0,proj_c,tau,sig,tauf,q);
    LogsDCA(N-ii+1).time = toc;
    LogsDCA(N-ii+1).logs = iterlogsdca;
    [c,I]   = cost(C,X,projfun);
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
    [c2,I2]  = cost(C,X2,projfun);
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
    [c3,I3]       = cost(C,X3,projfun);
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
function P = proj(X,CB,RB,q)
    [k,n] = size(X);
    Pb    = zeros(k*q,n);
    for j=1:k
        Pb(2*j-1:2*j,:) = projball2(X(j,:),CB(2*j-1:2*j,:),RB(2*j-1:2*j));
    end
    P = Pb;
end

function [FileName] = CreateUniqueFileName(FileName)
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