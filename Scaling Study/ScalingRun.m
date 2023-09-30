%%%%
% Scaling Run - Test DCA, BDCA and BDCA adaptive for a constrained
% clustering problem on selected dimensions, with a selected number of
% randomly chosen points. This generates a variable 'Logs' which contains 
% detailed information on all of the solves.
% 
% Constraints are given by balls with centers given by variable C and
% radius given by variable r. By default we have 3 constraints centered at
% C_1 = [1,5....1,5], C_2 = [6,4,6,4.....,6,4], C_3 = [8,8,....,8]
% with radius 1. See Example 7.2 
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

N          = 100; %Number of Running Times
Dimensions = [2,3,5,10,20]; % Set of dimensions to test
NumPoints  = [50,100,500,1000,5000,10000,50000];     % Set of points to test
savetofile = true; %Whether to save to file or not     

% parameters for the solver 
tau  = 1;  % projection penalty
sig  = 10; % penalty growth parameter
tauf = 1e8;

%preallocate Log File
Logs = cell(length(Dimensions),length(NumPoints));

for dimiter = 1:length(Dimensions)
   for ptiter = 1:length(NumPoints)
       dim    = Dimensions(dimiter); 
       numpts = NumPoints(ptiter); 
       A = unifrnd(0,10,numpts,dim);
       
       [m,n] = size(A);
       C     = [ones(1,dim)*3-2*(-1).^(0:dim-1); ones(1,dim)*5+ (-1).^(0:dim-1); ones(1,dim)*8]; %Constraint Centers
       [k,~] = size(C);
        
       q       = 1; 
       r       = ones(1,k);  %Constraint Radius
       projfun = @(mat)proj(mat,C,r);
       X       = zeros(k,n);
       for runiter = 1:N          
            for i =1:k
                X(i,:) = C(i,:)+r(i)*rand(1,n)/sqrt(2);
            end
            X0=X;
            %% DCA
            tic
            [X,N1,iterlogsdca]        = constrainedDCA2V2(A,X0,projfun,tau,sig,tauf,q,[]);
            LogsDCA(N-runiter+1).time = toc;
            LogsDCA(N-runiter+1).logs = iterlogsdca;
            c = cost_opt(A,X);
            dcaiter = 0;
            for i = 1:N1
                dcaiter = dcaiter + iterlogsdca(i).dcaiter;
            end
            LogsDCA(N-runiter+1).totaliter = dcaiter;
            LogsDCA(N-runiter+1).outeriter = N1;
            LogsDCA(N-runiter+1).cost      = c;
            %% BDCA -- Adaptive
            tic
            [X2,outeriter,iterlogsbdca]=constrainedBDCA2V2(A,X0,projfun,tau,sig,tauf,q,true,[]);
            LogsBDCA_adapt(N-runiter+1).time = toc;
            LogsBDCA_adapt(N-runiter+1).logs      = iterlogsbdca;
            c2 =cost_opt(A,X2);
            bdcaiter = 0;
            for i = 1:outeriter
                bdcaiter = bdcaiter + iterlogsbdca(i).dcaiter;
            end
            LogsBDCA_adapt(N-runiter+1).totaliter = bdcaiter; 
            LogsBDCA_adapt(N-runiter+1).outeriter = outeriter; 
            LogsBDCA_adapt(N-runiter+1).cost      = c2;
            
            %% BDCA -- Non Adaptive
            tic
            [X2,outeriter,iterlogsbdca]=constrainedBDCA2V2(A,X0,projfun,tau,sig,tauf,q,false,[]);
            LogsBDCA(N-runiter+1).time = toc;
            LogsBDCA(N-runiter+1).logs      = iterlogsbdca;
            c2 = cost_opt(A,X2);
            bdcaiter = 0;
            for i = 1:outeriter
                bdcaiter = bdcaiter + iterlogsbdca(i).dcaiter;
            end
            LogsBDCA(N-runiter+1).totaliter = bdcaiter; 
            LogsBDCA(N-runiter+1).outeriter = outeriter; 
            LogsBDCA(N-runiter+1).cost      = c2;
       end
       Log = struct;
       Log.dim       = dim; 
       Log.numpoints = numpts;
       Log.DCALogs   = LogsDCA; 
       Log.BDCAadaptiveLogs = LogsBDCA_adapt;
       Log.BDCALogs  = LogsBDCA; 
       
       Logs{dimiter,ptiter} = Log;
    end
end

if(savetofile)
    filename = CreateUniqueFileName('profile/Results');
    save(filename,'Logs');
end

%% Projection function
function P=proj(X,C,r)
[k,n]=size(X);
P = zeros(k,n);
for i = 1:k
    P(i,:) = projball2(X(i,:),C(i,:),r(i)); 
end
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

