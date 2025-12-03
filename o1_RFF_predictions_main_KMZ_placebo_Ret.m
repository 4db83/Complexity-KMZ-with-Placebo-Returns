clc; clear; tic;
warning('off', 'MATLAB:print:FigureTooLargeForPage'); warning('off', 'MATLAB:unreachableCode'); %#ok 
% --------------------------------------------------------------------------------------------------
set(groot,'defaultAxesXTickLabelRotationMode','manual');set(groot,'defaultLineLineWidth',2);set(groot,'defaultAxesFontSize',14)% set(groot,'defaultAxesFontName','Times New Roman')
% PATH 2 TOOLBOX: WIN (\) AND MAC/UNIX (/) % addpath(genpath('./utility.Functions'))                   % local path to utility.Functions
path_2_KMZ_local_functions = '../complexity KMZ - local.functions';
if (exist( path_2_KMZ_local_functions, 'dir')==7)
  addpath( genpath(path_2_KMZ_local_functions) ) % set path to db functions
else
  addpath(genpath('./local.functions')) % set path to db functions
end
% % To restore to defaultpath, call: restoredefaultpath %
% --------------------------------------------------------------------------------------------------
% CALL: get_all_db_toolbox_function_calls.m FROM code directory, collects needed db_toolbox functions
% --------------------------------------------------------------------------------------------------
% PRINT START TIME
disp(datetime('now'))
% NUMBER OF WORKERS. Set OUTSIDE OF LOOP  % delete(myCluster.Jobs)
start_parpool_with

% --------------------------------------------------------------------------------------------------
% FIX THE Y PLACEBO DATA SEED FOR ALL SIMULATIONS AT ONE VALUE. 
% --------------------------------------------------------------------------------------------------
% NUMBER OF SIMULATIONS: default is 1e3;
nSim = 1e3;
% Use value larger than 1000, so not to overlap with w weights RNDs 1:1000.
% placebo_seed = 1001;
% placebo_seed = 1111;
placebo_seed = 1234;

%**************************************************************************
% OUTPUT DIRECTORY NAME
%**************************************************************************
RFF_output_name = strcat( 'RFF_output_placebo_Ret_', num2str(placebo_seed), '/'); 
OUTPUT_DIR      = set_dir( strcat( './_nsims_source_', RFF_output_name ) );

%**************************************************************************
% Predictions for 1000 simulations with random seeds from 1 to 1000
% Note: Parallelization or HPC can be used to expedite the for loop
%**************************************************************************
% TRAINING WINDOW LIST
% trnwin_list = [ 12 60 120 ];
trnwin_list = [ 12 60 120 ];

% define P and z grids
Plist   = [2:2:20 24:12:(9*11) 10*(10:10:90) 1e3:1e3:3e3 4e3:2e3:12e3];
% lamlist = [1e-8 1e-6 1e-4 0.1 kron(10.^(-1:4), 2:2:10)];
% orginal SHRINKAGE PARAMETERS LAMBDA (Z) 
lamlist = 10.^(-3:1:3);

% gamma in Random Fourier Features
gamma = 2;
% Standardization the dependent variable Y
stdize_Y = 1;
% DEMEANING = FALSE
for demean = [0 1]
  % RUN FIRST RFFs now 
  for trnwin = trnwin_list
    % FIRST: RUN RFF ridge regressions: this takes time
    file_list = file_names_to_get(OUTPUT_DIR, trnwin, stdize_Y, demean, nSim);
    nFiles = length(file_list);
    fprintf("Computing: %.5g/%.5g.\n", [nFiles nSim])
    sep('*')
    % for ( ii = 1:nFiles )      
    parfor ( ii = 1:nFiles )
      seed_randn = file_list(ii);
      KMZ_tryrff_v2_function_for_each_sim_placebo_Ret(gamma, trnwin, seed_randn, stdize_Y, demean, OUTPUT_DIR, Plist, lamlist, placebo_seed)
    end
    % SECOND: RUN the Get GW benchmarks: this always uses lamlist = [0 10.^([-3:1:3])] and is quick 
    KMZ_GW_benchmark_function_placebo_Ret(trnwin, stdize_Y, demean, OUTPUT_DIR, placebo_seed);
  end
end




% --------------------------------------------------------------------------------------------------
% function to get existing file names to only compute those that don't already exist. Convenient
% when the cluster breaks down or runs out of time. 
% --------------------------------------------------------------------------------------------------
function list_of_seed_names = file_names_to_get(OUTPUT_DIR, trnwin, stdize_Y, demean, nSim)
% make function handle to easier create strings
fsp = @(trnwin,stdize_Y,demean) ( strcat('maxP-12000-trnwin-', num2str(trnwin), ...
        '-gamma-2-stdize-', num2str(stdize_Y), '-demean-', num2str(demean), '-v2') );
% string of file path
save_path       = strcat( OUTPUT_DIR, fsp(trnwin,stdize_Y,demean), '.' );
sep('*'); fprintf( save_path + " "); 
files_tmp       = dir([save_path]);
files_listed    = vertcat(char(files_tmp.name));
% GET FILE NAMES NUMBER ONLY AND CONVERT TO DOUBLE
existing_seeds  = sort(str2double(strrep(strrep(cellstr(files_listed),'iSim',''),'.mat','')));
seed_list_full  = (1:nSim)';
file_list_diff  = setdiff(seed_list_full, existing_seeds)';
% output
list_of_seed_names = file_list_diff;



end 







































% EOF