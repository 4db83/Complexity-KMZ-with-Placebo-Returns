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
% SET TO 1 TO DELETE THEM IF STORAGE IS AN ISSUE
DELETE_INDIVIDUAL_FILES   = 1;
% --------------------------------------------------------------------------------------------------
make_combined_files       = 1; % set to one if not already combined. 
make_performance_measures = 1;

% --------------------------------------------------------------------------------------------------
% FIX THE Y PLACEBO DATA SEED FOR ALL SIMULATIONS AT ONE VALUE. 
% --------------------------------------------------------------------------------------------------
% Use values larger than 1000, to not overlap with w weights RNDs 1:1000.
seed_set = [ 1001:1:2010 ];

for ii = 1:length(seed_set)
placebo_seed = seed_set(ii);
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


% **************************************************************************************************
% PATH to individual data files where the 1000 sims are stored from o1_RFF_predictions_main_KMZ
% **************************************************************************************************
RFF_output_name = strcat( 'RFF_output_placebo_Ret_', num2str(placebo_seed), '/');
NSIMS_SOURCE  = set_dir( strcat( './_nsims_source_', RFF_output_name ) );
COMBINED_PATH = set_dir( strcat( './_combined_RFF_', RFF_output_name ) );
GW_OUTPUT_DIR = NSIMS_SOURCE;


% MAIN Looping through to get all the individual files for all Trnwin etd.
stdize_Y  = 1;
% for demean = [ 0 1 ]
% for trnwin = [ 12 60 120 ]  
for demean = 0
for trnwin = 12 

%**************************************************************************
% Choices of parameters
%**************************************************************************
gamma   = 2; if gamma == 0.5; gamma_str = '0pt5'; else; gamma_str = num2str(gamma);end
% SAVE THE RESULTS
saveon  = 0;
% PORTFOLIO EVALUATION: FULL SAMPLE (this is not important for aggregration of output)
% --------------------------------------------------------------------------------------------------
subbeg = 1940;
% subbeg = 1931; % subbeg = 1935;
subend = 2020;
% --------------------------------------------------------------------------------------------------
% RESTRICT OUT-OF-SAMPLE EVALUATION PERIOD: 1940 onwards ∀ Tᵢₛ ∈ [12 60 120].
% Use subsample set to zero
subsamp = 0;
% max number of Random Fourier Features (RFFs)
maxP = 12000;
% the number of simulations is 1000
nSim = 1000; 
% SAVING STRING
para_str = strcat('maxP-', num2str(maxP), '-trnwin-', num2str(trnwin), '-gamma-', num2str(gamma), ...
                  '-stdize-', num2str(stdize_Y), '-demean-', num2str(demean), '-v2');

%**************************************************************************
% LOAD iSim1.mat file AND GW PARAMETERS only
%**************************************************************************
iSim1_file_path = strcat(NSIMS_SOURCE, para_str, '/iSim1.mat');
load( iSim1_file_path, 'T','nP','nL','Y', 'Plist', 'dates','lamlist' )
Y_Sim1 = Y;

% NOTE THAT THIS PATH DOES NOT CONTAIN THE FULL LIST OF FILES AFTER COMBINING THE DATA ... 
files_listing = dir( strcat(NSIMS_SOURCE, para_str, '/*.mat') );

% load benchmark GW(2008) "kitchen sink" regression results
if make_performance_measures
  disp('Loading Goyal Welch Benchmark...')
  load([GW_OUTPUT_DIR 'gybench-trnwin-' num2str(trnwin) '-stdize-' num2str(stdize_Y)  '-demean-' num2str(demean) '.mat'])
  disp('Done: Loading Goyal Welch Benchmark ...')
end 

%**************************************************************************
% Collect results of 1000 simulations
%**************************************************************************
combine_out_file_save = ['trnwin-' num2str(trnwin) '-gamma-' num2str(gamma) '-stdize-' num2str(stdize_Y) '-demean-' num2str(demean) '-Yprd_timing'];
comb_filename = strcat(COMBINED_PATH, combine_out_file_save, '.mat'); 

% --------------------------------------------------------------------------------------------------
% LOOP THROUGH THE INDIVIDUAL SIMULATIONS RESULTS AND COMBINE ALL RESULTS INTO 4-D ARRAYS
% --------------------------------------------------------------------------------------------------
if ~make_combined_files
    % if the combined result exists
    disp('Loading existing combined 4D data files, ... this about 66 seconds... ')
    tic; load(comb_filename); toc
else
    % if the combined result doesn't exist
    disp('--------------------------------------------------------------------------------------\n');
    disp('Reading in separate individual files, ... ')
    tic
    Yprd_collect    = nan(T,nP,nL,nSim);        % predicted Y
    Bnrm_collect    = nan(T,nP,nL,nSim);        % beta norm  
    % for (s = 1:nSim)
    for (s = 1:length(files_listing))
        % tic
        if mod(s, 100) == 0 ; disp( num2str(s, ' %4d') ); end 
        % load data
        filename = strcat(NSIMS_SOURCE, para_str , '/', files_listing(s).name);
        % tmp_load = load(filename, 'Yprd', 'Bnrm');
        tmp_load = load(filename);
        % check that the Y sim values are all the same
        if any(tmp_load.Y ~= Y_Sim1) 
          sep; disp("Y from different sims not matching --> results invalid"); sep
          head([Y_Sim1' tmp_load.Y'], 10); sep
          break;   % BREAK --> stop the loop
        end

        Yprd_collect(:,:,:,s)   = tmp_load.Yprd;
        Bnrm_collect(:,:,:,s)   = tmp_load.Bnrm;
        
        % make new y_sim1 for newly loaded data to check recursively
        Y_Sim1 = tmp_load.Y;
    end
    toc
    clearvars Yprd Bnrm
    Yprd = Yprd_collect;
    Bnrm = Bnrm_collect; % sum(Beta_hat.^2); ie., (∥β∥₂)²​ = ∑βᵢ² or β′β​

    Bnrmbar     = nanmean(Bnrm,4);              % nanmedian(Bnrm,4);
    Bnrmbar     = squeeze(nanmean(Bnrmbar,1));  % squeeze(nanmedian(Bnrmbar,1));
    timing      = Yprd.*Y';
    % NOT NEEDED
    % pihat       = [squeeze(nanmean(Yprd(:,nP,nL,:),4)) squeeze(nanmean(Yprd(:,1,1,:),4))]; % this is for P=12000 P=2 or 0? whatever the first entry is
    
    clearvars Yprd_collect Bnrm_collect
    % ADD -nocompression TAG, much faster, only a litte larger: takes 10 secs... 
    disp(['Saving data to combined data to ' comb_filename] )
    disp(' ... without compression this takes some time ... ')
    
    % SAVING Yprd, timing, 1092x46x7x1000 matrices 
    % tic; save(comb_filename, 'Yprd', 'timing', 'Bnrmbar', '-v7.3', '-nocompression'); toc
    % TO SAVE SPACE, DO NOT SAVE timing, JUST COMPUTE IT AFTERWARDS AS: timing = Yprd.*Y'
    tic; save(comb_filename, 'Yprd', 'Y', 'Bnrmbar', '-v7.3', '-nocompression'); toc % saves little space ... but
end
disp('--------------------------------------------------------------------------------------\n');
disp('Data Loaded')

% ADJUST THE SIZES CORRECTLY
% nP = length(Plist); nL = length(lamlist); T  = length(dates);

% COMPUTE TIMING IF NOT SAVED PREVIOUSLY
timing = Yprd.*Y'; 

% --------------------------------------------------------------------------------------------------
% DELETE INDIVIDUAL FILES (ONLY IF NEEDED FOR SPACE REASONS)
% --------------------------------------------------------------------------------------------------
if DELETE_INDIVIDUAL_FILES == 1
   disp("Deleting Sim Files")
  if length(files_listing) > 1
    for s = 2:nSim % keep the first iSim1.mat;
      filename = strcat(NSIMS_SOURCE, para_str , '/', files_listing(s).name);
      delete(filename);
%      disp("Deleting" + filename)
    end
  end
end
% --------------------------------------------------------------------------------------------------

% PERCENTILE LIST
pctlist     = [1 2.5 5 25 50 75 95 97.5 99];
% Generate exhibits
nPct = length(pctlist);
nSub = length(subbeg);
% leave this at 1
ss = 1; 

% Evaluation period
% LOC_EVAL_PERIOD = find(dates>=subbeg(ss)*100 & dates<=(subend(ss)+1)*100);
Iss   = (dates>=subbeg(ss)*100 & dates<=(subend(ss)+1)*100);
Inan  = ~isnan(Y' + timing(:,1,1,1)); 
Isel  = Iss & Inan; % Isel  = logical(Iss.*Inan);

% first subsample start
first_ss = find(Isel,1,'first');
% find first Non_nan value
firstNONnan = find(Inan,1,'first');
% GENERIC OUTPUT NAME
NAME_SUFFIX = ['trnwin-' num2str(trnwin) '-gamma-' gamma_str '-stdize-' num2str(stdize_Y) '-demean-' num2str(demean) '-' num2str(subbeg(ss)) '-' num2str(subend(ss))];
performance_filename = strcat(COMBINED_PATH, NAME_SUFFIX, '_KMZ_port_eval.mat');
ss_dates = char(num2str(dates(Iss)));
% FIRST NON-NAN EVALUATION PERIODS DATE
date0 = num2str( dates(firstNONnan) ) ;

dates_strng = [ss_dates(1,1:end-2) ':' ss_dates(1,end-1:end) ' - ' ss_dates(end,1:end-2) ':' ss_dates(end,end-1:end)];
if str2double(date0(1:4)) > str2double(dates_strng(1:4)); dates_strng(1:4) = date0(1:4); end

sep;  % print Sample evaluation information
fprintf('Data Start: %d. ', subbeg)
fprintf('Training window size T=%d. ', trnwin );
fprintf('Effective Portfolio Evaluation period: %s. \n', dates_strng); 
sep
% ADJUST DATES TO BE OVER THE NON-NAN PERIOD
fprintf('First non-nan Value at: %d.\n', firstNONnan)
fprintf('  Sub-sample starts at: %d.\n', first_ss)
fprintf('Effective Evaluation Period starts at: %d.\n', first_ss)

% **************************************************************************************************
% COMPUTE PORTFOLIO PERFORMANCE MEASURES
% SNR in Tibshiran 2022 etc, is SNR = (∥β∥₂)²/Var(uhat)
% **************************************************************************************************

% --------------------------------------------------------------------------------------------------
% COMPUTE/GENERATE PERFORMANCE MEASUREMENTS 
% --------------------------------------------------------------------------------------------------
if ~make_performance_measures
  % % OLD if file exists: if isfile(performance_filename)
  disp('Loading data KMZ portfolio performance file, ... ')
  tic; load(performance_filename); toc;
else
  disp('--------------------------------------------------------------------------------------\n');
  disp('Computing portfolio performance measures from separate simulations, ... ')
  % Performance initialization
  ER          = nan(nP,nL);
  SR          = nan(nP,nL);
  Vol         = nan(nP,nL);
  IR          = nan(nP,nL);
  IR_tstat    = nan(nP,nL);
  alpha       = nan(nP,nL);
  R2          = nan(nP,nL);
  IRstd       = nan(nP,nL);
  % already trimmed to correct size of subbeg
  Ytmp        = Y(:,Isel)';
  II          = ones(length(Ytmp),1);
  % add extra measures
  maxLoss     = nan(nP,nL);
  Skew        = nan(nP,nL);
  SR_tstat    = nan(nP,nL);
  IR_vs_GWz1000 = nan(nP,nL);
  IR_tstat_vs_GWz1000 = nan(nP,nL);
  Bias        = nan(nP,nL);
  MSE         = nan(nP,nL);
  for P=1:nP
    for L=1:nL
      ERtmp         = nan(nSim,1);
      Voltmp        = nan(nSim,1);
      IRtmp         = nan(nSim,1);
      IRttmp        = nan(nSim,1);
      alphatmp      = nan(nSim,1);
      R2tmp         = nan(nSim,1);
      IRstdtmp      = nan(nSim,1);
      % already trimmed to correct size of subbeg
      timtmp        = squeeze(timing(Isel,P,L,:));
      Yprdtmp       = squeeze(Yprd(Isel,P,L,:));
      uhat_tmp      = Ytmp - Yprdtmp;
      % add extra measures
      maxLosstmp    = nan(nSim,1);
      Skewtmp       = nan(nSim,1);
      SR_tstattmp   = nan(nSim,1);
      IR_vs_GWz1000_tmp       = nan(nSim,1);
      IR_tstat_vs_GWz1000_tmp = nan(nSim,1);
      % add Bias and MSE for each sim seperately
      Biastmp       = nan(nSim,1);
      MSEtmp        = nan(nSim,1);
      % do we need this?
      % loc         = find(~isnan(Yprdtmp(:,1)+Ytmp));
      % for (i=1:nSim) % regstats(y, X, model)
      parfor (i=1:nSim) % regstats(y, X, model)
        stats       = regstats(timtmp(:,i),Ytmp,'linear',{'tstat','r'});   % takes care of nans
        SRtmp(i)    = sharpe(timtmp(:,i),0);                               % takes care of nans
        ERtmp(i)    = nanmean(timtmp(:,i));
        Voltmp(i)   = nanstd(timtmp(:,i));
        IRtmp(i)    = stats.tstat.beta(1)/nanstd(stats.r);
        IRttmp(i)   = stats.tstat.t(1);
        alphatmp(i) = stats.tstat.beta(1);
        IRstdtmp(i) = nanstd(stats.r);
        R2tmp(i)    = 1-nanvar( Yprdtmp(:,i) - Ytmp ) / nanvar(Ytmp);
        % add extra measures
        maxLosstmp(i)  = -nanmin(timtmp(:,i));
        Skewtmp(i)     = skewness(timtmp(:,i));
        stats_SRtmp    = regstats(timtmp(:,i), II, 1, {'tstat'}); % add a one to not include a constant term in the regressions
        SR_tstattmp(i) = stats_SRtmp.tstat.t(1);
        % now add IR vs GW(z=1000)
        % NOT NEEDED stats_vs_GWz1000            = regstats( timtmp(loc,i), timing_gy(loc,end), 'linear',{'tstat','rsquare','r'});
        OLS_vs_GWz1000              = regstats( timtmp(:,i), timing_gy(Isel,end), 'linear',{'tstat','rsquare','r'});
        IR_vs_GWz1000_tmp(i)        = sharpe( OLS_vs_GWz1000.tstat.beta(1) + OLS_vs_GWz1000.r, 0);
        IR_tstat_vs_GWz1000_tmp(i)  = OLS_vs_GWz1000.tstat.t(1);
        Biastmp(i)  = mean(uhat_tmp(:,i));
        MSEtmp(i)   = mean(uhat_tmp(:,i).^2);
      end
      % Average them now over the 1000 Nsims
      SR(P,L)       = nanmean(SRtmp);
      ER(P,L)       = nanmean(ERtmp);
      Vol(P,L)      = nanmean(Voltmp);
      IR(P,L)       = nanmean(IRtmp);
      IR_tstat(P,L) = nanmean(IRttmp);
      alpha(P,L)    = nanmean(alphatmp);
      R2(P,L)       = nanmean(R2tmp);
      IRstd(P,L)    = nanmean(IRstdtmp);
      % add extra measures
      maxLoss(P,L)  = nanmean(maxLosstmp);
      Skew(P,L)     = nanmean(Skewtmp);
      SR_tstat(P,L) = nanmean(SR_tstattmp);
      IR_vs_GWz1000(P,L)       = nanmean(IR_vs_GWz1000_tmp);
      IR_tstat_vs_GWz1000(P,L) = nanmean(IR_tstat_vs_GWz1000_tmp);
      Bias(P,L)     = nanmean(Biastmp);
      MSE(P,L)      = nanmean(MSEtmp);
    end
    disp(['p=' num2str(P)])
  end
  toc
  % SAVE THE DATA
  % --------------------------------------------------------------------------------------------
  save( performance_filename, 'Isel','Iss','Inan','SR','ER','Vol','IR','IR_tstat','alpha','R2','IRstd', ...
    'maxLoss', 'Skew', 'SR_tstat', 'IR_vs_GWz1000', 'IR_tstat_vs_GWz1000', 'Bias', 'MSE', ...
    '-v7.3', '-nocompression');
  % --------------------------------------------------------------------------------------------
end

% loop though all files up to here
end; end

fprintf(' Finished reading single files and constructing portfolio measues... \n')

% % CLOSE THE LOG FILE FOR READING THE INDIVIDUAL FILES
% diary off

end








 
% %%
% clc  
% p36 = squeeze(Yprd);
% p36_mu = mean(p36,3);

% hold on;
% plot(Y','Color',.8*ones(1,3))
% plot(p36_mu)
% hold off































% EOF
