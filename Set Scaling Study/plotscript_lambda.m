%%%%
% Plotscript_Lambda - Generate figures of the lambda values
% throughout a BDCA run from the 'Logs' files resulting from the 
% 'Scaling Run' script. This is how those figures for the paper were
% generated. 
% 
% The parameters for the script: 
%
% Change the loaded file to select which 'Logs' file to generate figures
% for
% 
% select which cell(s) in the 'Logs' files to generate plots for
%%%%

clc
clear

load('profile/paper.mat'); % which saved 'Logs' file to plot

plotall = false; % If true generate a lambda plot for every numpoints/dim pair in the 'Logs' file

% dimension point pairs in the form [dim_1, pt_1; dim_2, pt_2; ... ] 
% to plot lambda for
dim_pt_pair  = [1,4; 1,8];

% particular run to generate the plot from 
run = 1;

numdim  = size(Logs,1);
numpts  = size(Logs,2);
numruns = length([Logs{1,1}.DCALogs.time]);

timeratio  = zeros(numruns,numpts);
points     = zeros(numpts,1); 
dimemsions = zeros(numdim,1);

for i = 1:numpts
    points(i) = Logs{1,i}.numpoints;
end

for i = 1:numdim
    dimemsions(i)   = Logs{i,1}.dim; 
end

if (plotall) % override dm_pt_pair with entire set of combinations
    dim_pt_pair = [repelem(dimemsions,numpts), repmat(points,numdim,1)];
end

%%
markersize = 8;
marker     = '.';

linemark = ':'; 

% Create the lambda plots
for i = 1:size(dim_pt_pair,1)
    dimm = dim_pt_pair(i,1); 
    pts  = dim_pt_pair(i,2);
    % BDCA 
    lambda    = vertcat(Logs{dimm,pts}.BDCALogs(run).logs.lambdaiter);
    taushifts = cumsum([Logs{dimm,pts}.BDCALogs(run).logs.dcaiter]);

    figure;
    plot(lambda,marker,"MarkerSize",markersize);
    hold on 
    xline(taushifts(1:end-1),linemark);
    hold off
    axis padded
    xlabel('Total Iterations','Interpreter','latex');
    ylabel('$\lambda$','Interpreter','latex');
    title("BDCA $\lambda$ values for Dimension "+ num2str(dimemsions(dimm))+ " and " + num2str(points(pts)) + " Points" ,'Interpreter','latex');
    fontsize(14,'points');
    %BDCA adaptive
    lambda = vertcat(Logs{dimm,pts}.BDCAadaptiveLogs(run).logs.lambdaiter);
    taushifts = cumsum([Logs{dimm,pts}.BDCAadaptiveLogs(run).logs.dcaiter]);
    
    figure;
    plot(lambda,marker,"MarkerSize",markersize);
    hold on 
    xline(taushifts(1:end-1),linemark);
    hold off
    axis padded
    xlabel('Total Iterations','Interpreter','latex');
    ylabel('$\lambda$','Interpreter','latex');
    title("BDCA Adaptive $\lambda$ values for Dimension "+ num2str(dimemsions(dimm))+ " and " + num2str(points(pts)) + " Points" ,'Interpreter','latex');
    fontsize(14,'points');
end