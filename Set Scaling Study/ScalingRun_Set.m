%%%%
% Scaling Run - Test DCA, BDCA and BDCA adaptive for a constrained set
% clustering problem on selected dimensions, with a selected number of
% randomly chosen ball sets of radius 0.1. This generates a variable 
% 'Logs' which contains detailed information on all of the solves.
% 
% Constraints are given by the intersections of two balls with centers given 
% by variable CB and radius given by variable RB. By default we have 4
% constraints centered each formed by the intersection of two ball of radius 
% 1. See Example 7.5 for details.
%
% The parameters for the script
%
% N          - Number of times to rerun each numpoints/dimension pair for a
%              different random initial point.
% Dimensions - Set of dimensions to test
% NumPoints  - Set of the number of points to test for each dimension
% savetofile - Whether to save output log file 'Logs' to a .mat file or not
%
% Output: 
%
% Logs       - Log file containing detailed information about all of the
%              runs. A plotting script can use this
%              log file to generate figures of the same kind used in the
%              paper. 
%%%%

clc
clear

N = 100; %Number of Running Times

Dimensions = [2,3,5,10]; %Set of Dimensions to test
NumPoints  = [50,100,250,500,1000,5000,10000,50000]; %Set of Points to test
       
tau  = 1;  % penalty term
sig  = 10; % penalty growth parameter
tauf = 1e8;
ep   = 1e-6;
k    = 4; % number of seperate constraints
q    = 2; % number of overlaps in the constraints

 
CB_main = [1 5 5 5 5 5 5 5 5 5 ; 2 6 5 5 5 5 5 5 5 5 ; 5 4 1 2 3 1 2 3 1 2; 4 4 1 2 3 1 2 3 1 2  ; 8 5 9 8 7 9 8 7 9 8;...
    8 4 9 8 7 9 8 7 9 8;9 8 1 6 9 1 6 9 1 6;8 8 1 6 9 1 6 9 1 6 ]; %centers of constraints, defined up to 10D
RB_main = [1 1 1 1 1 1 1 1]; % Radius of centers

savetofile = true; %Whether to save to file or not

Logs = cell(length(Dimensions),length(NumPoints));

for dimiter = 1:length(Dimensions)
   for ptiter = 1:length(NumPoints)
       dim    = Dimensions(dimiter); 
       numpts = NumPoints(ptiter); 

       A = unifrnd(0,10,numpts,dim); % Data Centers
       R = 0.1*ones(numpts,1);         % Data Set Radius

       CB = CB_main(:,1:dim);
       RB = RB_main;

       projfun=@(mat)projball2(mat,A,R);
       proj_c =@(mat)proj(mat,CB,RB,q);

       X0       = zeros(k,dim);
       for runiter = 1:N          
            j = 1;
            for i = 1:2:2*k
                X0(j,:) = CB(i,:)+RB(i)*rand(1,dim)/sqrt(2);
                j       = j+1; 
            end
            %% DCA
            tic
            [X,N1,iterlogsdca]=constrainedDCA2V2(A,R,X0,proj_c,tau,sig,tauf,q);
            LogsDCA(N-runiter+1).time = toc;
            LogsDCA(N-runiter+1).logs = iterlogsdca;

            [c,~] = cost(A,X,projfun);
            dcaiter = 0;
            for i = 1:N1
                dcaiter = dcaiter + iterlogsdca(i).dcaiter;
            end
            LogsDCA(N-runiter+1).totaliter = dcaiter;
            LogsDCA(N-runiter+1).outeriter = N1;
            LogsDCA(N-runiter+1).cost      = c;

            %% BDCA -- Adaptive
            tic
            [X2,outeriter,iterlogsbdca]=constrainedBDCA2V2(A,R,X0,proj_c,tau,sig,tauf,q,true);
            LogsBDCA_adapt(N-runiter+1).time = toc;
            LogsBDCA_adapt(N-runiter+1).logs = iterlogsbdca;
            [c2,~] = cost(A,X2,projfun);
            bdcaiter = 0;
            for i = 1:outeriter
                bdcaiter = bdcaiter + iterlogsbdca(i).dcaiter;
            end
            LogsBDCA_adapt(N-runiter+1).totaliter = bdcaiter; 
            LogsBDCA_adapt(N-runiter+1).outeriter = outeriter; 
            LogsBDCA_adapt(N-runiter+1).cost      = c2;
            
            %% BDCA -- Non Adaptive
            tic
            [X2,outeriter,iterlogsbdca]=constrainedBDCA2V2(A,R,X0,proj_c,tau,sig,tauf,q,false);
            LogsBDCA(N-runiter+1).time = toc;
            LogsBDCA(N-runiter+1).logs = iterlogsbdca;
            [c2,~] = cost(A,X2,projfun);
            bdcaiter = 0;
            for i = 1:outeriter
                bdcaiter = bdcaiter + iterlogsbdca(i).dcaiter;
            end
            LogsBDCA(N-runiter+1).totaliter = bdcaiter; 
            LogsBDCA(N-runiter+1).outeriter = outeriter; 
            LogsBDCA(N-runiter+1).cost      = c2;
       end
       Log                  = struct;
       Log.dim              = dim; 
       Log.numpoints        = numpts;
       Log.numRuns          = N;
       Log.DCALogs          = LogsDCA; 
       Log.BDCAadaptiveLogs = LogsBDCA_adapt;
       Log.BDCALogs         = LogsBDCA; 
       
       Logs{dimiter,ptiter} = Log;
    end
end

if(savetofile)
    filename = CreateUniqueFileName('Profiling/SetResults');
    save(filename,'Logs');
end
%% Projection function
function P=proj(X,CB,RB,q)
[k,n]=size(X);
Pb = zeros(k*q,n);
for j=1:k
    Pb(2*j-1:2*j,:) = projball2(X(j,:),CB(2*j-1:2*j,:),RB(2*j-1:2*j));
end
P=Pb;
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

