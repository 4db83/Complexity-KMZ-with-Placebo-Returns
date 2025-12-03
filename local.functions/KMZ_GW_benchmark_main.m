clear; clc
% print start time
% disp(datetime('now'))
addpath(genpath('./KMZ-local.functions/'))           % set path to db functions

%**************************************************************************
% Parameters Setting
%**************************************************************************

% training window list
trnwin_list = [12, 60, 120];

% Standardization = True
stdize  = 0;
% Demeaning = False
demean  = 1;
% save the results
SAVEON  = 1;
% PATH
save_path = './individual-files_stdize-0/';
if ~exist(save_path); mkdir(save_path); end

%**************************************************************************
% Get GW benchmark
%**************************************************************************

for trnwin = trnwin_list
    GW_benchmark_function(trnwin,stdize,demean,save_path,saveon)
end