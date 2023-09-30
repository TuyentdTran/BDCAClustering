%% README 
%
%% Overview 
% Here we include the code used to generate the data and figures used in
% the paper. There are folders corresponding to the different examples in
% the paper, with the main functions and scripts documented. We will now
% list each folder and the relevant scripts inside that should be run.
% 
% The logging data used for the examples in the paper is always stored in a
% profile folder with name paper.mat. So for example: 
%
% *Scaling Study/profile/paper.mat*
%
% contains the data for Example 7.2 in the paper. 
%
% Also, all folders will have a script named plotscript, that generates
% performance figures for that examples log files. 
%
%% 50 Cities Set Clusters 
% Used to generate the results of Example 7.3 in the paper.
% The relevant scripts to run are as follows: 
%  
% # *Testing50Cities_Nruns.m*
% Generates Log files for the performance data in the profile folder.
% # *Test_50UScities.m*
% Basic timing results (not used in the paper) and generates a figure 
% of the solution. 
% # *plotscript.m*
% Plot performance figures for the generated logfiles.
%
%% Cities Set clusters
% Used to generate the results of Example 7.4 in the paper.
% The relevant scripts to run are as follows: 
%  
% # *Testing_USCities_N.m*
% Generates Log files for the performance data in the profile folder.
% # *Testing_UScities.m*
% Basic timing results (not used in the paper) and generates a figure 
% of the solution. 
% # *plotscript.m*
% Plot performance figures for the generated logfiles.
%
%% 76 city
% Used to generate the results of Example 7.1 in the paper.
% The relevant scripts to run are as follows: 
%  
% # *result.m*
% Generates Log files for the performance data in the profile folder.
% # *new_76.m*
% Basic timing results (not used in the paper) and generates a figure 
% of the solution. 
% # *plotscript.m*
% Plot performance figures for the generated logfiles.
%
%% Scaling Study
% Used to generate the results of Example 7.2 in the paper.
% The relevant scripts to run are as follows: 
%  
% # *ScalingRun.m*
% Generates Log files for the performance data in the profile folder.
% # *plotscript.m*
% Plot performance figures for the generated logfiles.
%
%% Set Scaling Study
% Used to generate the results of Example 7.5 in the paper.
% The relevant scripts to run are as follows: 
%  
% # *ScalingRun_Set.m*
% Generates Log files for the performance data in the profile folder.
%
% # *plotscript.m*
% Plot performance figures for the generated logfiles.
%