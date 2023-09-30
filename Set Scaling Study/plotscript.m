clc
clear
% close all

load('profile/paper.mat');

numdim  = size(Logs,1);
numpts  = size(Logs,2);
numruns = length([Logs{1,1}.DCALogs.time]);

timeratio = zeros(numruns,numpts);
points    = zeros(numpts,1); 
dimemsions      = zeros(numdim,1);

for i = 1:numpts
    points(i) = Logs{1,i}.numpoints;
end

for i = 1:numdim
    dimemsions(i)   = Logs{i,1}.dim; 
end


%%Plots ignoring variance and only plotting the average runtiem/iteration
%%count data that can be expected 
markersize = 20;
marker     = '.';

%% DCA/BDCA

avg_time_ratio = zeros(numdim,numpts); 
avg_iter_ratio = zeros(numdim,numpts);
legendstr      = cell(numdim,1);
for j = 1:numdim
   for i = 1:numpts
        avg_iter_ratio(j,i) = sum([Logs{j,i}.DCALogs.totaliter])./sum([Logs{j,i}.BDCALogs.totaliter]);
        avg_time_ratio(j,i) = sum([Logs{j,i}.DCALogs.time])./sum([Logs{j,i}.BDCALogs.time]);
   end
   legendstr{j} = "dim = " + num2str(dimemsions(j)); 
end

figure;
semilogx(points',avg_iter_ratio(1,:),marker,'MarkerSize',markersize);
hold on
for j = 2:numdim
    semilogx(points',avg_iter_ratio(j,:),marker,'MarkerSize',markersize);
end
yline(1,'k--');
hold off
xticks(points);
ylim([0,inf]);
ylim("padded");
xlabel('Number of Points');
ylabel('Average Iteration Ratio');
title("DCA/BDCA");
legend(legendstr,'Location','northwest');
fontsize(12,'points');
   
figure;
semilogx(points',avg_time_ratio(1,:),marker,'MarkerSize',markersize);
hold on
for j = 2:numdim
    semilogx(points',avg_time_ratio(j,:),marker,'MarkerSize',markersize);
end
yline(1,'k--');
hold off
xticks(points);
ylim([0,inf]);
ylim("padded");
xlabel('Number of Points');
ylabel('Average Time Ratio');
title("DCA/BDCA");
legend(legendstr,'Location','northwest');
fontsize(12,'points');

   

%% DCA/Adaptive BDCA

avg_time_ratio = zeros(numdim,numpts); 
avg_iter_ratio = zeros(numdim,numpts);
legendstr      = cell(numdim,1);
for j = 1:numdim
   for i = 1:numpts
        avg_iter_ratio(j,i) = sum([Logs{j,i}.DCALogs.totaliter])./sum([Logs{j,i}.BDCAadaptiveLogs.totaliter]);
        avg_time_ratio(j,i) = sum([Logs{j,i}.DCALogs.time])./sum([Logs{j,i}.BDCAadaptiveLogs.time]);
   end
   legendstr{j} = "dim = " + num2str(dimemsions(j)); 
end

figure;
semilogx(points',avg_iter_ratio(1,:),marker,'MarkerSize',markersize);
hold on
for j = 2:numdim
    semilogx(points',avg_iter_ratio(j,:),marker,'MarkerSize',markersize);
end
yline(1,'k--');
hold off
xticks(points);
% limsys = ylim; 
% ylim([0,limsys(2)]);
ylim([0,inf]);
ylim("padded");
xlabel('Number of Points');
ylabel('Average Iteration Ratio');
title("DCA/BDCA Adaptive");
legend(legendstr,'Location','northwest');
fontsize(12,'points');

figure;
semilogx(points',avg_time_ratio(1,:),marker,'MarkerSize',markersize);
hold on
for j = 2:numdim
    semilogx(points',avg_time_ratio(j,:),marker,'MarkerSize',markersize);
end
yline(1,'k--');
hold off
xticks(points); 
ylim([0,inf]);
ylim("padded");
xlabel('Number of Points');
ylabel('Average Time Ratio');
title("DCA/BDCA Adaptive");
legend(legendstr,'Location','northwest');
fontsize(12,'points');