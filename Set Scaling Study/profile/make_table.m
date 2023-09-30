clc
clear
format short

load('paper.mat');

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

mean_dca = zeros(numdim,numpts);
mean_bdca   = zeros(numdim,numpts);
mean_adapt  = zeros(numdim,numpts);

for j = 1:numdim
   for i = 1:numpts
        mean_dca(j,i)   = mean([Logs{j,i}.DCALogs.time]);
        mean_bdca(j,i)  = mean([Logs{j,i}.BDCALogs.time]);
        mean_adapt(j,i) =mean([Logs{j,i}.BDCAadaptiveLogs.time]);
   end  
end

std_dca = zeros(numdim,numpts);
std_bdca   = zeros(numdim,numpts);
std_adapt  = zeros(numdim,numpts);
for j = 1:numdim
   for i = 1:numpts
        std_dca(j,i)   = std([Logs{j,i}.DCALogs.time],0,2);
        std_bdca(j,i)  = std([Logs{j,i}.BDCALogs.time],0,2);
        std_adapt(j,i) =std([Logs{j,i}.BDCAadaptiveLogs.time],0,2);
   end  
end

DCA = zeros(numdim*2,numpts);
BDCA = zeros(numdim*2,numpts);
BDCA_adapt = zeros(numdim*2,numpts);
for j = 1:numdim
    DCA(2*j-1:2*j,:)=[mean_dca(j,:);std_dca(j,:)];    
end

for j = 1:numdim
    BDCA(2*j-1:2*j,:)=[mean_bdca(j,:);std_bdca(j,:)];    
end

for j = 1:numdim
    BDCA_adapt(2*j-1:2*j,:)=[mean_adapt(j,:);std_adapt(j,:)];    
end
