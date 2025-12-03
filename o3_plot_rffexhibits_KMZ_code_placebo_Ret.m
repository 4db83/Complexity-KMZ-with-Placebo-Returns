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
% start_parpool_with

PRNT_PDF  = 1;
PRNT_XLS  = 1;
PLOT_db   = 0;

Monitor_Pos = 1;              % change to where plots need to be displayed

% --------------------------------------------------------------------------------------------------
% FIX THE Y PLACEBO DATA SEED FOR ALL SIMULATIONS AT ONE VALUE. 
% --------------------------------------------------------------------------------------------------
% Use value larger than 1000, because w weights use 1:1000.
placebo_seed = 1001;
placebo_seed = 1111;

% **************************************************************************************************
% PATH to individual data files where the 1000 sims are stored from o1_RFF_predictions_main_KMZ
% **************************************************************************************************
RFF_output_name = strcat( 'RFF_output_placebo_Ret_', num2str(placebo_seed), '/');
NSIMS_SOURCE  = set_dir( strcat( './_nsims_source_', RFF_output_name ) );
COMBINED_PATH = set_dir( strcat( './_combined_RFF_', RFF_output_name ) );
GW_OUTPUT_DIR = NSIMS_SOURCE;
OOS_EVAL_OUTPUT_DIR = set_dir( strcat( './_oos_eval_results_KMZ_', RFF_output_name ));

% MAIN Looping through to get all the individual files for all Trnwin etd.
stdize_Y  = 1;
% for demean = [ 0 1 ]
% for trnwin = [ 12 60 120 ]  
demean = 0;
trnwin = 12;

% **************************************************************************************************
% Choices of parameters
% **************************************************************************************************
gamma   = 2; if gamma == 0.5; gamma_str = '0pt5'; else; gamma_str = num2str(gamma);end
% SAVE THE RESULTS
saveon  = 0;
% PORTFOLIO EVALUATION: FULL SAMPLE (this is not important for aggregration of output)
% --------------------------------------------------------------------------------------------------
subbeg = 1940;
% subbeg = 1931; % subbeg = 1935; % subbeg = 1940;
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

% **************************************************************************************************
% LOAD GW PARAMETERS AND BENCHMARK iSim1.mat file only
% **************************************************************************************************
iSim1_file_path = strcat(NSIMS_SOURCE, para_str, '/iSim1.mat');
% load( iSim1_file_path, 'T','nP','nL','Y', 'Plist', 'dates','lamlist' )
load( iSim1_file_path, 'Y', 'Plist', 'dates','lamlist' )
% nP = length(Plist); nL = length(lamlist); T  = length(dates);
iSim1_all = load( iSim1_file_path );

% NOTE THAT THIS PATH DOES NOT CONTAIN THE FULL LIST OF FILES AFTER COMBINING THE DATA ... 
files_listing = dir( strcat(NSIMS_SOURCE, para_str, '/*.mat') );

% load benchmark GW(2008) "kitchen sink" regression
disp('Loading Goyal Welch Benchmark...')
load([GW_OUTPUT_DIR 'gybench-trnwin-' num2str(trnwin) '-stdize-' num2str(stdize_Y)  '-demean-' num2str(demean) '.mat'])
disp('Done: Loading Goyal Welch Benchmark ...')

% **************************************************************************************************
% Collect results of 1000 simulations
% **************************************************************************************************
combine_out_file_save = ['trnwin-' num2str(trnwin) '-gamma-' num2str(gamma) '-stdize-' num2str(stdize_Y) , ...
                        '-demean-' num2str(demean) '-Yprd_timing'];
comb_filename = strcat(COMBINED_PATH, combine_out_file_save, '.mat'); 

% --------------------------------------------------------------------------------------------------
% LOAD COMBINED DATA 
% --------------------------------------------------------------------------------------------------
disp('Loading existing combined 4D data files, ... this takes about 11 seconds ... ')
tic; load(comb_filename); toc
disp('--------------------------------------------------------------------------------------');
disp('Data Loaded')

%% REMOVE THE COLUMN WITH Z = 0 AND USE THE Z = 1E-8 FOR NUMERICAL STABILITY. 
% lamlist(2)    = [];Yprd(:,:,2,:) = [];Bnrmbar(:,2)  = [];
nP = length(Plist); nL = length(lamlist); T  = length(dates);
L1000 = find(lamlist==1000);

% COMPUTE TIMING IF NOT SAVED PREVIOUSLY
timing = Yprd.*Y'; 

% leave this at 1
ss = 1; 
% Evaluation period
% LOC_EVAL_PERIOD = find(dates>=subbeg(ss)*100 & dates<=(subend(ss)+1)*100);
Iss   = (dates>=subbeg(ss)*100 & dates<=(subend(ss)+1)*100);
Inan  = ~isnan(Y' + timing(:,1,1,1)); 
Isel  = Iss & Inan; % Isel  = logical(Iss.*Inan);

% first subsample start
first_ss    = find(Isel,1,'first');
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

% START LOGGING THE OUTPUT IN DIARY
diary_filename = strcat(OOS_EVAL_OUTPUT_DIR, 'read_single_files_diary_', NAME_SUFFIX, '.txt');
if exist(diary_filename, 'file'); delete(diary_filename); end
diary(diary_filename)

% **************************************************************************************************
% COMPUTE PORTFOLIO PERFORMANCE MEASURES
% SNR in Tibshiran 2022 etc, is SNR = (∥β∥₂)²/Var(uhat)
% **************************************************************************************************

% --------------------------------------------------------------------------------------------------
% COMPUTE/GENERATE PERFORMANCE MEASUREMENTS 
% --------------------------------------------------------------------------------------------------
disp('Loading data KMZ portfolio performance file, ... this takes about 1 seconds ... ')
load(performance_filename);
  
% % ************************************************************************************************
% if ~make_performance_measures
%   % % OLD if file exists: if isfile(performance_filename)
%   disp('Loading data KMZ portfolio performance file, ... ')
%   tic; load(performance_filename); toc;
% else
%   disp('--------------------------------------------------------------------------------------\n');
%   disp('Computing portfolio performance measures from separate simulations, ... ')
%   % Performance initialization
%   Bias        = nan(nP,nL);
%   ER          = nan(nP,nL);
%   SR          = nan(nP,nL);
%   Vol         = nan(nP,nL);
%   IR          = nan(nP,nL);
%   IR_tstat    = nan(nP,nL);
%   alpha       = nan(nP,nL);
%   R2          = nan(nP,nL);
%   IRstd       = nan(nP,nL);
%   % already trimmed to correct size of subbeg
%   Ytmp        = Y(:,Isel)';
%   II          = ones(length(Ytmp),1);
%   % add extra measures
%   maxLoss     = nan(nP,nL);
%   Skew        = nan(nP,nL);
%   SR_tstat    = nan(nP,nL);
%   IR_vs_GWz1000 = nan(nP,nL);
%   IR_tstat_vs_GWz1000 = nan(nP,nL);
%   for P=1:nP
%     for L=1:nL
%       ERtmp         = nan(nSim,1);
%       Voltmp        = nan(nSim,1);
%       IRtmp         = nan(nSim,1);
%       IRttmp        = nan(nSim,1);
%       alphatmp      = nan(nSim,1);
%       R2tmp         = nan(nSim,1);
%       IRstdtmp      = nan(nSim,1);
%       % already trimmed to correct size of subbeg
%       timtmp        = squeeze(timing(Isel,P,L,:));
%       Yprdtmp       = squeeze(Yprd(Isel,P,L,:));
%       uhat_tmp      = Ytmp - Yprdtmp;
%       % add extra measures
%       maxLosstmp    = nan(nSim,1);
%       Skewtmp       = nan(nSim,1);
%       SR_tstattmp   = nan(nSim,1);
%       IR_vs_GWz1000_tmp       = nan(nSim,1);
%       IR_tstat_vs_GWz1000_tmp = nan(nSim,1);
%       % add Bias and MSE for each sim seperately
%       Biastmp       = nan(nSim,1);
%       MSEtmp        = nan(nSim,1);
%       % do we need this?
%       % loc         = find(~isnan(Yprdtmp(:,1)+Ytmp));
%       for (i=1:nSim) % regstats(y, X, model)
%         stats       = regstats(timtmp(:,i),Ytmp,'linear',{'tstat','r'});   % takes care of nans
%         SRtmp(i)    = sharpe(timtmp(:,i),0);                               % takes care of nans
%         ERtmp(i)    = nanmean(timtmp(:,i));
%         Voltmp(i)   = nanstd(timtmp(:,i));
%         IRtmp(i)    = stats.tstat.beta(1)/nanstd(stats.r);
%         IRttmp(i)   = stats.tstat.t(1);
%         alphatmp(i) = stats.tstat.beta(1);
%         IRstdtmp(i) = nanstd(stats.r);
%         R2tmp(i)    = 1-nanvar( Yprdtmp(:,i) - Ytmp ) / nanvar(Ytmp);
%         % add extra measures
%         maxLosstmp(i)  = -nanmin(timtmp(:,i));
%         Skewtmp(i)     = skewness(timtmp(:,i));
%         stats_SRtmp    = regstats(timtmp(:,i), II, 1, {'tstat'}); % add a one to not include a constant term in the regressions
%         SR_tstattmp(i) = stats_SRtmp.tstat.t(1);
%         % now add IR vs GW(z=1000)
%         % NOT NEEDED stats_vs_GWz1000            = regstats( timtmp(loc,i), timing_gy(loc,end), 'linear',{'tstat','rsquare','r'});
%         OLS_vs_GWz1000              = regstats( timtmp(:,i), timing_gy(Isel,end), 'linear',{'tstat','rsquare','r'});
%         IR_vs_GWz1000_tmp(i)        = sharpe( OLS_vs_GWz1000.tstat.beta(1) + OLS_vs_GWz1000.r, 0);
%         IR_tstat_vs_GWz1000_tmp(i)  = OLS_vs_GWz1000.tstat.t(1);
%         Biastmp(i)  = mean(uhat_tmp(:,i));
%         MSEtmp(i)   = mean(uhat_tmp(:,i).^2);
%       end
%       % Average them now over the 1000 Nsims
%       SR(P,L)       = nanmean(SRtmp);
%       ER(P,L)       = nanmean(ERtmp);
%       Vol(P,L)      = nanmean(Voltmp);
%       IR(P,L)       = nanmean(IRtmp);
%       IR_tstat(P,L) = nanmean(IRttmp);
%       alpha(P,L)    = nanmean(alphatmp);
%       R2(P,L)       = nanmean(R2tmp);
%       IRstd(P,L)    = nanmean(IRstdtmp);
%       % add extra measures
%       maxLoss(P,L)  = nanmean(maxLosstmp);
%       Skew(P,L)     = nanmean(Skewtmp);
%       SR_tstat(P,L) = nanmean(SR_tstattmp);
%       IR_vs_GWz1000(P,L)       = nanmean(IR_vs_GWz1000_tmp);
%       IR_tstat_vs_GWz1000(P,L) = nanmean(IR_tstat_vs_GWz1000_tmp);
%       Bias(P,L)     = nanmean(Biastmp);
%       MSE(P,L)      = nanmean(MSEtmp);
%     end
%     disp(['p=' num2str(P)])
%   end
%   toc
%   % SAVE THE DATA
%   % --------------------------------------------------------------------------------------------
%   save( performance_filename, 'Isel','Iss','Inan','SR','ER','Vol','IR','IR_tstat','alpha','R2','IRstd', ...
%     'maxLoss', 'Skew', 'SR_tstat', 'IR_vs_GWz1000', 'IR_tstat_vs_GWz1000', 'Bias', 'MSE', ...
%     '-v7.3', '-nocompression');
%   % --------------------------------------------------------------------------------------------
% end
% % ************************************************************************************************

% --------------------------------------------------------------------------------------------
sep('='); disp("Data-Info.: " + para_str(12:end-3))
% --------------------------------------------------------------------------------------------

% ==================================================================================================
% COMPUTE PORTFOLIO PERFORMANCE MEASURES BY FIRST AGGREGATING RFF FORECASTS AND PLOT
% ==================================================================================================
dates_oos = dates(Isel,:);
mean_timing = nanmean(timing(Isel,:,:,:),4);
mean_Yprd   = nanmean(Yprd(  Isel,:,:,:),4);
Ytmp        =            Y(  Isel)';
II  = ones(length(mean_timing),1);
% [~,nP,nL]  = size(mean_timing);
SR0     = nan(nP,nL);
ER0     = nan(nP,nL);
Vol0    = nan(nP,nL);
IR0     = nan(nP,nL);
alpha0  = nan(nP,nL);
R20     = nan(nP,nL);
Skew0   = nan(nP,nL);
SR0_tstat = nan(nP,nL);
IR0_tstat = nan(nP,nL);
IR0_vs_GWz1000 = nan(nP,nL);
IR0_tstat_vs_GWz1000 = nan(nP,nL);
% add uhat and Bias and MSE(uhat)
uhat_0  = Ytmp - mean_Yprd;
Bias0   = squeeze(nanmean(uhat_0));
MSE0    = squeeze(nanmean(uhat_0.^2));

for jj = 1:nL
  SR0(:,jj)     = sharpe(   mean_timing(:,:,jj), 0)';
  ER0(:,jj)     = nanmean(  mean_timing(:,:,jj)   )';
  Vol0(:,jj)    = nanstd(   mean_timing(:,:,jj)   )';
  Skew0(:,jj)   = skewness( mean_timing(:,:,jj)   )';
  maxLoss0(:,jj)= -nanmin(  mean_timing(:,:,jj)   )';
  for ii = 1:nP % regstats(y, X, model)
    % SR tstat
    OLS_II_ii         = regstats(mean_timing(:,ii,jj), II, 1, {'tstat'}); % add a one to not include a constant term in the regressions
    SR0_tstat(ii,jj)  = OLS_II_ii.tstat.t(1);
    % IR Market
    OLS_xRet_ii = regstats(mean_timing(:,ii,jj), Ytmp,'linear',{'tstat','r'});
    % CAPM: R(ii) = α + β*R(M) + ε. Using α/σ(ε) = E[R(ii)-β*R(M)]/σ[R(ii)-β*R(M)] =
    % E[α+ε]/σ[α+ε] ≈ Sharpe(α+ε), because Sharpe function uses T not T-1 in DF of variance computation.
    % IR0(ii,jj)      = OLS_xRet_ii.tstat.beta(1)/nanstd(OLS_xRet_ii.r);  % Using α/σ(ε)
    IR0(ii,jj)      = sharpe( OLS_xRet_ii.tstat.beta(1) + OLS_xRet_ii.r,0 );  % using Sharpe function as in KMZ
    IR0_tstat(ii,jj)= OLS_xRet_ii.tstat.t(1);
    alpha0(ii,jj)   = OLS_xRet_ii.tstat.beta(1);
    % IR vs GW(z=1000)
    OLS_vs_GWz1000              = regstats(mean_timing(:,ii,jj), timing_gy(Isel,end), 'linear',{'tstat','rsquare','r'});
    IR0_vs_GWz1000(ii,jj)       = sharpe( OLS_vs_GWz1000.tstat.beta(1) + OLS_vs_GWz1000.r, 0);
    IR0_tstat_vs_GWz1000(ii,jj) = OLS_vs_GWz1000.tstat.t(1);
    % R2
    R20(ii,jj) = 1-nanvar(mean_Yprd(:,ii,jj)-Ytmp)/nanvar(Ytmp);
  end
end

%% MAKE A TABLE FROM ARRAYS AND FIND MAX VALUES FOR EACH ONE
decml  = 6; PRINT_COMP = 0;
zNames = strrep(cellstr(strcat('z=',char(num2str(lamlist')))),' ','');
% zNames = strrep(zNames,'z=1e-08','z=0');
pNames = strrep(cellstr(strcat('P=',char(num2str(Plist')))),' ','');

% MAKE A TABLE
% compute beta Norm (average)
tNorm_beta    = array2table( Bnrmbar      ,'VariableNames',zNames,'RowNames', pNames);
% kmz aggregation
tbl_SR_kmz    = array2table(sqrt(12)*SR   ,'VariableNames',zNames,'RowNames', pNames);
tbl_ER_kmz    = array2table( ER           ,'VariableNames',zNames,'RowNames', pNames);
tbl_Vol_kmz   = array2table( Vol          ,'VariableNames',zNames,'RowNames', pNames);
tbl_IR_kmz    = array2table(sqrt(12)*IR   ,'VariableNames',zNames,'RowNames', pNames);
tbl_IRt_kmz   = array2table( IR_tstat     ,'VariableNames',zNames,'RowNames', pNames);
tbl_alpha_kmz = array2table( alpha        ,'VariableNames',zNames,'RowNames', pNames);
tbl_R2_kmz    = array2table( R2           ,'VariableNames',zNames,'RowNames', pNames);
tbl_Bias_kmz  = array2table( Bias         ,'VariableNames',zNames,'RowNames', pNames);
tbl_MSE_kmz   = array2table( MSE          ,'VariableNames',zNames,'RowNames', pNames);
tbl_SNR_kmz   = array2table( Bnrmbar./(Vol.^2),'VariableNames',zNames,'RowNames', pNames);
% correct aggregation
tbl_SR_0      = array2table(sqrt(12)*SR0  ,'VariableNames',zNames,'RowNames', pNames);
tbl_ER_0      = array2table( ER0          ,'VariableNames',zNames,'RowNames', pNames);
tbl_Vol_0     = array2table( Vol0         ,'VariableNames',zNames,'RowNames', pNames);
tbl_IR_0      = array2table(sqrt(12)* IR0 ,'VariableNames',zNames,'RowNames', pNames);
tbl_IRt_0     = array2table( IR0_tstat    ,'VariableNames',zNames,'RowNames', pNames);
tbl_alpha_0   = array2table( alpha0       ,'VariableNames',zNames,'RowNames', pNames);
tbl_R2_0      = array2table( R20          ,'VariableNames',zNames,'RowNames', pNames);
tbl_Bias_0    = array2table( Bias0        ,'VariableNames',zNames,'RowNames', pNames);
tbl_MSE_0     = array2table( MSE0         ,'VariableNames',zNames,'RowNames', pNames);
tbl_SNR_0     = array2table( Bnrmbar./(Vol0.^2),'VariableNames',zNames,'RowNames', pNames);

% % remove the z=1e-08 and z=1e-06 entries --> not needed
% % Bnrmbar(:,1:2) = [];
% PM_name = {'SR';'ER';'Vol';'IR';'IRt';'alpha';'R2';'Bias';'MSE';'SNR'};
% for kii = 1:length(PM_name)
%   eval(['tbl_' char(PM_name(kii)) '_kmz(:,{''z=1e-08'';''z=1e-06''}) = [];'])
%   eval(['tbl_' char(PM_name(kii)) '_0(:,{''z=1e-08'';''z=1e-06''}) = [];'])
% end
% tNorm_beta(:,1:2) = [];
% nL = size(tbl_alpha_0,2); 
% zNames(1:2) = [];

% ================================================================================================
% print performance measures to screen
% ================================================================================================
% Sharpe Ratios (SRs)
[pSR_kmz, zSR_kmz, maxSR_kmz, pSR_0, zSR_0, maxSR_0] = function_print_max_values(tbl_SR_kmz,tbl_SR_0,"SR", Plist,trnwin,PRINT_COMP);
% ------------------------------------------------------------------------------------------------
% Information Ratios (IRs)
% ------------------------------------------------------------------------------------------------
[pIR_kmz, zIR_kmz, maxIR_kmz, pIR_0, zIR_0, maxIR_0] = function_print_max_values(tbl_IR_kmz,tbl_IR_0,"IR", Plist,trnwin,PRINT_COMP);
% alpha (LINEAR --> MUST BE SAME)
% ------------------------------------------------------------------------------------------------
[palpha_kmz, zalpha_kmz, maxalpha_kmz, palpha_0, zalpha_0, maxalpha_0] = function_print_max_values(tbl_alpha_kmz,tbl_alpha_0,"α", Plist,trnwin,PRINT_COMP);
% ------------------------------------------------------------------------------------------------
% alpha t-stat (called IR-tstat)
% ------------------------------------------------------------------------------------------------
[pIRt_kmz, zIRt_kmz, maxIRt_kmz, pIRt_0, zIRt_0, maxIRt_0] = function_print_max_values(tbl_IRt_kmz,tbl_IRt_0,"IR-tstat", Plist,trnwin,PRINT_COMP);
% ------------------------------------------------------------------------------------------------
% Expected Return (E[R]) (LINEAR --> MUST BE SAME)
[pER_kmz, zER_kmz, maxER_kmz, pER_0, zER_0, maxER_0] = function_print_max_values(tbl_ER_kmz,tbl_ER_0,"ER", Plist,trnwin,PRINT_COMP);
% ------------------------------------------------------------------------------------------------
%  R²
[pR2_kmz, zR2_kmz, maxR2_kmz, pR2_0, zR2_0, maxR2_0] = function_print_max_values(tbl_R2_kmz,tbl_R2_0,"R2", Plist,trnwin,PRINT_COMP);
% ------------------------------------------------------------------------------------------------
%  MSE
[pMSE_kmz, zMSE_kmz, maxMSE_kmz, pMSE_0, zMSE_0, maxMSE_0] = function_print_max_values(tbl_MSE_kmz,tbl_MSE_0,"MSE", Plist,trnwin,PRINT_COMP,1);
% ------------------------------------------------------------------------------------------------
sep('=')
fprintf('Finished reading single files and constructing portfolio measures ... \n'); sep('=')

%% PLOTS ALL ------------------------------------------------------------------------------------ 

Leg_all_pos = 3;
Leg_all_col = 1;
% Leg_SR_pos = 5;
% Leg_SR_col = 1;
% Leg_all_buffer = [55 -15]; 
Leg_all_buffer = [-10 -20]; 

% define color palat to use
clr_palette = clr((1:nL)*1,1);
clrs2 = [ 0.50 0.50 0.50 ;  
          0.00 0.70 0.70 ];   
clr_palette = make_cellarray([lines(7); clrs2]);

if nL == 7; clr_palette = make_cellarray(lines(7)) ; end 
if nL >  9; clr_palette = clr((1:nL)*1,1) ; end 

% COMPLEXITY VECTOR IS
c = (Plist/trnwin)'; 
% some plotting controls  
Cline = .7;
% for Xlim0 = [ 0 ] 
if trnwin == 12;  xos = 5;  end
if trnwin == 60;  xos = 1;  end
if trnwin == 120; xos =.5;  end
Xlim0 =  0;
XlimUp = 50;  
% linewidth(s) of the two lines
LW  = 2.5; MK   = 'none';
LW0 = 2.0; MK0  = 'none';
LType = '-.'; xLbl = '$c$'; 
DK = .7;
SF = 15;
SFb = -.24; % move subtitel 
XLb = -.05; % move x-label
RSP =   85; % rowspace for figures
FNs = SF; clf; 
% LEGEND FONT SIZE
leg_FNTs = FNs-3;
F1 = figure(1); clf; % F1.Color= .6*ones(1,3); %F1 = colordef(F1,'black') % FIG WITH DARK BACKGROUND
monitors = get(0,'MonitorPositions'); F1.OuterPosition = monitors(Monitor_Pos,:);
Figure_Title = ['Evaluation Measures. Tis = ' num2str(trnwin) '. OOS: ' dates_strng ...
                ' (Stdize_Y=' num2str(stdize_Y) ', Demean =' num2str(demean) ')' ];
set(F1,'Name',Figure_Title,'WindowState','maximized');

tiledlayout(5,2,'TileSpacing', 'compact', 'Padding', 'loose');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if trnwin == 12 ; vLTH1 = 1/3*1000; vLTH2 = 3/4*1000; end 
if trnwin == 60 ; vLTH1 = 1/3*200;  vLTH2 = 3/4*200; end 
if trnwin == 120; vLTH1 = 1/3*100;  vLTH2 = 3/4*100; end 
decml_kmz = "%.4g"; decml_0 = "%.4g"; 

nexttile(1) % SR
hold on;
  pp1 = plot(c, tbl_SR_kmz.Variables ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
  set(pp1, {'Color'}, clr_palette);  
  hVa_clrs = get(pp1, 'Color');
  if PLOT_db
    pp1 = plot(c, tbl_SR_0.Variables, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9); 
  end
hold off;
big_CLRS = reshape([pp1.Color],3,length(pp1))';
% DEFINE COLOR PALETTE TO BE USED
CLRS4Labels = cell2mat(hVa_clrs);
% FIX COLOR ORDER to be the same, just a bit darker in the dashed line plots
for LL = 1:length(hVa_clrs);  pp1(LL).Color = hVa_clrs{LL}*DK; end
XL = xlabel(xLbl,'interpreter','latex');set(XL,'Units','normalized'); set(XL,'Position',[.5 XLb]);
if Xlim0;  xlim([0 XlimUp]); end
box on; addgrid
setyticklabels(0:.1:.7, 1, FNs); % ylim([0 .65])
% if trnwin == 60; setxticklabels(0:40:200); end
if (trnwin == 60)&&(~Xlim0); setxticklabels(0:40:200); end
ax = gca; ylim([-.064 ax.YLim(2)+.00]); hline(0); 
ST = subtitle('(a) Sharpe Ratio (SR)','Interpreter','LaTex','FontSize', SF); set(ST,'Units','normalized'); set(ST,'Position',[.5 SFb]);
if Xlim0; line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', Cline); end
xlim([0-xos c(end)+xos])
% ADD A VERTICAL LINE AT MAX(SR)
if PLOT_db
xL1_0 = xline( c(pSR_0), 'Color', CLRS4Labels(zSR_0,:)*DK,'Label',  ... % 'HandleVisibility', 'off', ...
  strcat('max[SR] $=', num2str(maxSR_0,'%2.4f'),' ~[c = ', num2str(c(pSR_0(1)),decml_0),',', ...
  tbl_SR_0(pSR_0,zSR_0).Properties.VariableNames,'$]'), ...
  'LabelOrientation','horizontal','Alpha',1, 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2 );
end
xL1_kmz = xline( c(pSR_kmz), 'Color', CLRS4Labels(zSR_kmz,:),'Label',... % 'HandleVisibility', 'off', ...
  strcat('max[SR$_{\mathrm{KMZ}}$] $=', num2str(maxSR_kmz,'%2.4f'),' ~[c = ', num2str(c(pSR_kmz(1)),decml_kmz),',', ...
  tbl_SR_kmz(pSR_kmz,zSR_kmz).Properties.VariableNames,'$]'), ...
  'LabelOrientation','horizontal','Alpha',1, 'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment', 'right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2 );
% adjust labeling position
if PLOT_db
  if c(pSR_0(1))   > vLTH1;   xL1_0.LabelHorizontalAlignment = 'center';  end
  if c(pSR_0(1))   > vLTH2;   xL1_0.LabelHorizontalAlignment = 'left';    end
end
if c(pSR_kmz(1)) > vLTH1; xL1_kmz.LabelHorizontalAlignment = 'center';  end
if c(pSR_kmz(1)) > vLTH2; xL1_kmz.LabelHorizontalAlignment = 'left';    end
yline(max(maxSR_kmz,maxSR_0),'Color','m','LineWidth',.66,'LineStyle','-')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile(2) % ER (LINEAR Measure, some as ER0) KMZ don't scale/annualize E(R) by 12
hold on;
  pp1 = plot(c, tbl_ER_0.Variables, '-','linewidth',LW0, 'Marker', MK0, 'Markersize', 9); 
  set(pp1, {'Color'}, clr_palette); 
hold off;
% FIX COLOR ORDER to be the same, just a bit darker in the dashed line plots
for LL = 1:length(hVa_clrs); pp1(LL).Color = hVa_clrs{LL}; end
if Xlim0;  xlim([0 XlimUp]); end
box on; addgrid 
if stdize_Y == 0
  ylim([-.825e-5 .9e-4])
  set(gca,'FontName','Times New Roman','Fontsize',FNs,'YColor','k','XColor','k')
else
  setyticklabels([0:.02:.08], 2, FNs); % ylim([0 .085])
  % if trnwin == 60; setxticklabels(0:40:200); end
  if (trnwin == 60)&&(~Xlim0); setxticklabels(0:40:200); end
  ax = gca; ylim([-.0073 ax.YLim(2)+.00]);
end
if (trnwin == 60)&&(~Xlim0); setxticklabels(0:40:200); end
hline(0); 
ST = subtitle('(b) Expected Return (ER)','Interpreter','LaTex', 'FontSize', SF); set(ST,'Units','normalized'); set(ST,'Position',[.5 SFb]);
if Xlim0; line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', Cline); end
XL = xlabel(xLbl,'interpreter','LaTex'); set(XL,'Units','normalized'); set(XL,'Position',[.5 XLb]);
xlim([0-xos c(end)+xos])
% ADD A VERTICAL LINE AT MAX(ER)
xL2_0 = xline(c(pER_0),'Color', CLRS4Labels(zER_0,:),'LabelOrientation','horizontal','Alpha',1,'Label', ...
  strcat('max[ER] $=',num2str(maxER_0,'%2.4f'),' ~[c = ', num2str(c(pER_0),decml_0),',', ...
  tbl_ER_0(pER_0,zER_0).Properties.VariableNames,'$]'), ...
 'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2);  
% adjust labeling position
if c(pER_0(1))   > vLTH1;   xL2_0.LabelHorizontalAlignment = 'center';  end
if c(pER_0(1))   > vLTH2;   xL2_0.LabelHorizontalAlignment = 'left';    end
yline(max(maxER_kmz,maxER_0),'Color','m','LineWidth',.66,'LineStyle','-')
% print2pdf('here')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile(3) % IR
hold on;
  pp1 = plot(c, tbl_IR_kmz.Variables ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
  set(pp1, {'Color'}, clr_palette); 
  hVa_clrs = get(pp1, 'Color');
  if PLOT_db
    pp1 = plot(c, tbl_IR_0.Variables, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9); 
  end
hold off;
% FIX COLOR ORDER to be the same, just a bit darker in the dashed line plots
for LL = 1:length(hVa_clrs);  pp1(LL).Color = hVa_clrs{LL}*DK; end
% xlabel(xLbl,'interpreter','latex')
XL = xlabel(xLbl,'interpreter','latex'); set(XL,'Units','normalized'); set(XL,'Position',[.5 XLb]);
if Xlim0;  xlim([0 XlimUp]); end
box on; addgrid
setyticklabels(0:.1:.4, 1, FNs); % ylim([0 .35])
% if trnwin == 60; setxticklabels(0:40:200); end
if (trnwin == 60)&&(~Xlim0); setxticklabels(0:40:200); end
ax = gca; ylim([-.036 ax.YLim(2)+.0]); hline(0);
ST = subtitle('(c) Information Ratio (IR)','Interpreter','LaTex','FontSize', SF); set(ST,'Units','normalized'); set(ST,'Position',[.5 SFb]);
if Xlim0; line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', Cline); end
xlim([0-xos c(end)+xos])
% ADD A VERTICAL LINE AT MAX(IR)
if PLOT_db
xL3_0 = xline( c(pIR_0), 'Color', CLRS4Labels(zIR_0,:)*DK,'Label',  ... % 'HandleVisibility', 'off', ...
  strcat('max[IR] $=', num2str(maxIR_0,'%2.4f'),' ~[c = ', num2str(c(pIR_0(1)),decml_0),',', ...
  tbl_IR_0(pIR_0,zIR_0).Properties.VariableNames,'$]'), ...
  'LabelOrientation','horizontal','Alpha',1, 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2 );
end
xL3_kmz = xline( c(pIR_kmz), 'Color', CLRS4Labels(zIR_kmz,:),'Label',... % 'HandleVisibility', 'off', ...
  strcat('max[IR$_{\mathrm{KMZ}}$] $=', num2str(maxIR_kmz,'%2.4f'),' ~[c = ', num2str(c(pIR_kmz(1)),decml_kmz),',', ...
  tbl_IR_kmz(pIR_kmz,zIR_kmz).Properties.VariableNames,'$]'), ...
  'LabelOrientation','horizontal','Alpha',1, 'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment', 'right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2 );
% adjust labeling position
if PLOT_db
  if c(pIR_0(1))   > vLTH1;   xL3_0.LabelHorizontalAlignment = 'center';  end
  if c(pIR_0(1))   > vLTH2;   xL3_0.LabelHorizontalAlignment = 'left';    end
end
if c(pIR_kmz(1)) > vLTH1; xL3_kmz.LabelHorizontalAlignment = 'center';  end
if c(pIR_kmz(1)) > vLTH2; xL3_kmz.LabelHorizontalAlignment = 'left';    end
yline(max(maxIR_kmz,maxIR_0),'Color','m','LineWidth',.66,'LineStyle','-')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile(4) % alpha 
hold on;
  pp1 = plot(c, tbl_alpha_0.Variables, '-','linewidth',LW0, 'Marker', MK0, 'Markersize', 9); 
  set(pp1, {'Color'}, clr_palette); 
hold off;
% FIX COLOR ORDER to be the same, just a bit darker in the dashed line plots
for LL = 1:length(hVa_clrs); pp1(LL).Color = hVa_clrs{LL}; end
XL = xlabel(xLbl,'interpreter','LaTex');set(XL,'Units','normalized'); set(XL,'Position',[.5 XLb]);
if Xlim0;  xlim([0-10.5 XlimUp]); end
xlim([0-xos c(end)+xos])
box on; addgrid 
% setyticklabels([0:.0001:.0003], 4, FNs); 
ylim([-.812e-5 .9e-4])
if stdize_Y == 0 
  set(gca,'FontName','Times New Roman','Fontsize',FNs,'YColor','k','XColor','k')
else
  setyticklabels([0:.01:.06], 2, FNs); 
  ax = gca; ylim([-.0055 ax.YLim(2)+.0]);
end 
hline(0); 
if (trnwin == 60)&&(~Xlim0); setxticklabels(0:40:200); end
ST = subtitle('(d) Alpha ($\alpha$)','Interpreter','LaTex','FontSize', SF); set(ST,'Units','normalized'); set(ST,'Position',[.5 SFb]);
if Xlim0; line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', Cline); end
xlim([0-xos c(end)+xos])
% ADD A VERTICAL LINE AT MAX(ER)
xL4_0 = xline(c(palpha_0),'Color', CLRS4Labels(zalpha_0,:),'LabelOrientation','horizontal','Alpha',1,'Label', ...
  strcat('max[$\alpha$] $=',num2str(maxalpha_0,'%2.4f'),' ~[c = ', num2str(c(palpha_0),decml_0),',', ...
  tbl_alpha_0(palpha_0,zalpha_0).Properties.VariableNames,'$]'), ...
 'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2);  
% adjust labeling position
if c(palpha_0(1))   > vLTH1;   xL4_0.LabelHorizontalAlignment = 'center';  end
if c(palpha_0(1))   > vLTH2;   xL4_0.LabelHorizontalAlignment = 'left';    end
yline(max(maxalpha_kmz,maxalpha_0),'Color','m','LineWidth',.66,'LineStyle','-')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile(5) % R2 % Two axis
hold on; 
  hVa = plot(c, tbl_R2_kmz.Variables ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
  set(hVa, {'Color'}, clr_palette); 
  hVa_clrs = get(hVa, 'Color');
  setyticklabels([-4:1:3], 0, FNs); % ylim([-3 .35])
  hline(0)
  yline(maxR2_kmz,'Color','m','LineWidth',.66,'LineStyle','-')
  if PLOT_db
yyaxis right
  hV0 =  plot(c, tbl_R2_0.Variables, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9); 
  set(hV0, {'Color'}, clr_palette)
  set(gca, 'YColor', 'k','Fontsize', FNs) % Set right y-axis color to black
  % FIX COLOR ORDER to be the same, just a bit darker in the dashed line plots
  for LL = 1:length(hVa_clrs); hV0(LL).Color = hVa_clrs{LL}*DK; end
  ax = gca; ax.YLim = [-.3 .11];
  % ylim([-.1 .02])  % R2ax = get(gca);  % R2_Ytkcs = R2ax.YTickLabel;
  hline(0); yline(maxR2_0,'Color','m','LineWidth',.66,'LineStyle','-');
hold off;
  end
box on; addgrid 
% % setyticklabels([-.2:.05:.0], 2, FNs); ylim([-.2 .025])
ST = subtitle('(e) $R^2$','Interpreter','LaTex','FontSize', SF); set(ST,'Units','normalized'); set(ST,'Position',[.5 SFb]);
if Xlim0; line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', Cline); end
XL = xlabel(xLbl,'interpreter','latex');set(XL,'Units','normalized'); set(XL,'Position',[.5 XLb]);
xlim([0-xos c(end)+xos])
if (trnwin == 60)&&(~Xlim0); setxticklabels(0:40:200); end
% ADD A VERTICAL LINE AT MAX(SR)
if PLOT_db
xL7_0 = xline( c(pR2_0), 'Color', CLRS4Labels(zR2_0,:)*DK,'Label',  ... % 'HandleVisibility', 'off', ...
  strcat('max[$R^2] =', num2str(maxR2_0,'%2.4f'),' ~[c = ', num2str(c(pR2_0(1)),decml_0),',', ...
  tbl_R2_0(pR2_0,zR2_0).Properties.VariableNames,'$] (RHS)'), ...
  'LabelOrientation','horizontal','Alpha',1, 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2 );
end
xL7_kmz = xline( c(pR2_kmz), 'Color', CLRS4Labels(zR2_kmz,:),'Label',... % 'HandleVisibility', 'off', ...
  strcat('max[$R^2_{\mathrm{KMZ}}] =', num2str(maxR2_kmz,'%2.4f'),' ~[c = ', num2str(c(pR2_kmz(1)),decml_kmz),',', ...
  tbl_R2_kmz(pR2_kmz,zR2_kmz).Properties.VariableNames,'$] (LHS)'), ...
  'LabelOrientation','horizontal','Alpha',1, 'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment', 'right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2 );
% adjust labeling position
if PLOT_db
  if c(pR2_0(1))   > vLTH1;   xL7_0.LabelHorizontalAlignment = 'center';  end
  if c(pR2_0(1))   > vLTH2;   xL7_0.LabelHorizontalAlignment = 'left';    end
end
if c(pR2_kmz(1)) > vLTH1; xL7_kmz.LabelHorizontalAlignment = 'center';  end
if c(pR2_kmz(1)) > vLTH2; xL7_kmz.LabelHorizontalAlignment = 'left';    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile(6) % |beta|
pp1 = plot(c, tNorm_beta.Variables ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
  set(pp1, {'Color'}, clr_palette); 
  hVa_clrs = get(pp1, 'Color');
XL = xlabel(xLbl,'interpreter','latex'); set(XL,'Units','normalized'); set(XL,'Position',[.5 XLb]);
if Xlim0;  xlim([0 XlimUp]); end
xlim([0-xos c(end)+xos])
box on; addgrid 
setyticklabels([0:.01:.05], 2, FNs); 
if stdize_Y == 1; setyticklabels([0:1:3], 0, FNs); end
if trnwin == 60; setxticklabels(0:40:200); end
ST = subtitle('(f) $\| \hat\beta \|$','Interpreter','LaTex','FontSize', SF); set(ST,'Units','normalized'); set(ST,'Position',[.5 SFb]);
xlim([0-xos c(end)+xos])
% --------------------------------------------------------------------------------------------------
% ADD LEGEND
% --------------------------------------------------------------------------------------------------
% leg_name = [cellstr(strcat(strcat('$z=',num2str((lamlist)')),'$'));'$c=1$']; 
leg_name = [cellstr(strcat('$',zNames,'$'))];
if stdize_Y == 0
  lh  = addlegend([pp1], leg_name,3,FNs-3,[],2/3,[-0 -18]);
else
  % lh  = addlegend([LegHndl], leg_name,3,FNs-3,[],2/3,[-0 -13]);
  lh =    addlegend([pp1], leg_name,Leg_all_pos,leg_FNTs,[],.44,Leg_all_buffer, 2, Leg_all_col);  
  lh.LineWidth = 1/2; % adjust line width of box around legend
end   % LLW = findobj(lh,'Type','Line'); set(LLW,'LineWidth',4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile(7) % alpha t-stat (IR t-stat)
hold on;
  pp1 = plot(c, tbl_IRt_kmz.Variables ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
  set(pp1, {'Color'}, clr_palette); 
  hVa_clrs = get(pp1, 'Color');
  if PLOT_db
    pp1 = plot(c, tbl_IRt_0.Variables, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9);  
  end
hold off;
% FIX COLOR ORDER to be the same, just a bit darker in the dashed line plots
for LL = 1:length(hVa_clrs); pp1(LL).Color = hVa_clrs{LL}*DK; end
XL = xlabel(xLbl,'interpreter','latex');set(XL,'Units','normalized'); set(XL,'Position',[.5 XLb]);
if Xlim0;  xlim([0 XlimUp]); end
box on; addgrid 
setyticklabels([0:1:4], 0, FNs);  % ylim([0 3.5])
% if trnwin == 60; setxticklabels(0:40:200); end
if (trnwin == 60)&&(~Xlim0); setxticklabels(0:40:200); end
ax = gca; ylim([-.38 ax.YLim(2)+.0]); hline(0); 
ST = subtitle('(g) Alpha $t-$statistic','Interpreter','LaTex','FontSize', SF); set(ST,'Units','normalized'); set(ST,'Position',[.5 SFb]);
if Xlim0; line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', Cline); end
xlim([0-xos c(end)+xos])
% ADD A VERTICAL LINE AT MAX(SR)
if PLOT_db
xL5_0 = xline( c(pIRt_0), 'Color', CLRS4Labels(zIRt_0,:)*DK,'Label',  ... % 'HandleVisibility', 'off', ...
  strcat('max[$\alpha(t-$stat)$]=', num2str(maxIRt_0,'%2.4f'),' ~[c = ', num2str(c(pIRt_0(1)),decml_0),',', ...
  tbl_IRt_0(pIRt_0,zIRt_0).Properties.VariableNames,'$]'), ...
  'LabelOrientation','horizontal','Alpha',1, 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2 );
end
xL5_kmz = xline( c(pIRt_kmz), 'Color', CLRS4Labels(zIRt_kmz,:),'Label',... % 'HandleVisibility', 'off', ...
  strcat('max[$\alpha(t-$stat$_{\mathrm{KMZ}})$)$]=', num2str(maxIRt_kmz,'%2.4f'),' ~[c = ', num2str(c(pIRt_kmz(1)),decml_kmz),',', ...
  tbl_IRt_kmz(pIRt_kmz,zIRt_kmz).Properties.VariableNames,'$]'), ...
  'LabelOrientation','horizontal','Alpha',1, 'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment', 'right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2 );
% adjust labeling position
if PLOT_db
  if c(pIRt_0(1))   > vLTH1;   xL5_0.LabelHorizontalAlignment = 'center';  end
  if c(pIRt_0(1))   > vLTH2;   xL5_0.LabelHorizontalAlignment = 'left';    end
end
if c(pIRt_kmz(1)) > vLTH1; xL5_kmz.LabelHorizontalAlignment = 'center';  end
if c(pIRt_kmz(1)) > vLTH2; xL5_kmz.LabelHorizontalAlignment = 'left';    end
yline(max(maxSR_kmz,maxSR_0),'Color','m','LineWidth',.66,'LineStyle','-')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile(8) % Volatility % KMZ don't scale/annualize vola by sqrt(12)
hold on;  % Two axis
pp1 =  plot(c, tbl_Vol_kmz.Variables ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9);
  set(pp1, {'Color'}, clr_palette); 
  hVa_clrs = get(pp1, 'Color');
  % setyticklabels([], 3, FNs); % 
  if stdize_Y == 0
    YTcks = [-0:.01:.02]; yticks(YTcks); 
    ylim(minmax(YTcks));    
    % ax = gca; ax.YTickLabel(end-1) = {'0.010'};
  else
    setyticklabels([0:2:10], 0, FNs);
  end
  % set(gca, 'YTick', YTcks); % YTcklabs = get(gca,'YTickLabel'); YTcklabs(1) = {'-1.0'}
if PLOT_db  
yyaxis right
  hV0 =  plot(c, tbl_Vol_0.Variables, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9); 
  % FIX COLOR ORDER to be the same, just a bit darker in the dashed line plots
  for LL = 1:length(hVa_clrs); hV0(LL).Color = hVa_clrs{LL}*DK; end
  if stdize_Y == 0
    ax = gca; ax.YLim = [-6e-4 6e-4];
    if trnwin == 12 && demean == 1
      ax = gca; ax.YLim = [4e-4 10e-4];
    end
  else
    ax = gca; ax.YLim = [-.8 .7];
  end  
hold off;
end
set(gca,'FontName','Times New Roman','Fontsize',FNs,'YColor','k','XColor','k')
% set(gca,'YTickLabel',YTcklabs);
hline(0)
% xlabel(xLbl,'interpreter','latex')
XL = xlabel(xLbl,'interpreter','latex');set(XL,'Units','normalized'); set(XL,'Position',[.5 XLb]);
if (trnwin == 60)&&(~Xlim0); setxticklabels(0:40:200); end
if Xlim0;  xlim([0 XlimUp]); end
xlim([0-xos c(end)+xos])
box on; addgrid 
% setyticklabels([0:1:5], 0, FNs); % ylim([0 .5])
% subtitle('Volatility','Interpreter','Latex')
ST = subtitle('(h) Volatility','Interpreter','LaTex','FontSize', SF); set(ST,'Units','normalized'); set(ST,'Position',[.5 SFb]);
if Xlim0; line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', Cline); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ADD PLOT ROWSPACE FOR SUBTITLES TO BE VISIBLE
addrowspace(RSP)
% pdf_name = strcat(para_str, '_eval_plots');
pdf_name = strcat('Tis=', num2str(trnwin), '_eval_all_plots_',  num2str(subbeg), '-', num2str(subend) );
if Xlim0;     pdf_name = strcat(pdf_name , '_short_Xaxis=',     num2str(XlimUp)); end

pdf_name_eval_all = strcat(OOS_EVAL_OUTPUT_DIR, NAME_SUFFIX,'_eval_all');
if PRNT_PDF;  plot2pdf(pdf_name_eval_all); end
% ---------------------------------------------------------------------------------------------

%% PLOTS FOR SR ONLY 
Leg_SR_pos = 5;
Leg_SR_col = 1;
% Leg_all_buffer = [55 -15]; 
% ---------------------------------------------------------------------------------------------
FNs = 16;  % FNs = 16;
% define performanc emeasure to evaluate
% PFM = 'SR'; PFM = 'ER'; 
% PFMs = {'SR','ER'}
% PFMs = {'ER'};
% PFMs = {'IR'};
% PFMs = {'alpha'};
PFMs = {'SR'};

ff = 1:length(PFMs);
PFM = char(PFMs(ff)); 
% PFM = 'MSE';
% if strcmp(PFM,'ER'); 

if strcmp(PFM,'ER') 
  YTcks_PFM = [-.00:.01:.04]; YTcks_decmls = 2;   osw = .90;
elseif strcmp(PFM,'MSE') 
  YTcks_PFM = [0:1:3];        YTcks_decmls = 0;   osw = .93;
elseif strcmp(PFM,'IR')  % IR use these 
  YTcks_PFM = [-0:.1:.4];      YTcks_decmls = 2;   osw = .90;
elseif strcmp(PFM,'alpha')  % IR use these 
  YTcks_PFM = [-0:.01:.03];    YTcks_decmls = 2;   osw = .90;
else % SR use these 
  YTcks_PFM = [-1:.1:.7];      YTcks_decmls = 1;   osw = .94;
end

eval( strcat('pPFM_0   = p',      PFM,'_0;'  ) );
eval( strcat('pPFM_kmz = p',      PFM,'_kmz;') );
eval( strcat('zPFM_0   = z',      PFM,'_0;'  ) );
eval( strcat('zPFM_kmz = z',      PFM,'_kmz;') );
eval( strcat('tbl_PFM_0 = tbl_',  PFM,'_0;'  ) );
eval( strcat('tbl_PFM_kmz = tbl_',PFM,'_kmz;') );
eval( strcat('maxPFM_0 = max',    PFM,'_0;'  ) );

F2 = figure(2); clf; % F1.Color= .6*ones(1,3); %F1 = colordef(F1,'black') % FIG WITH DARK BACKGROUND
monitors = get(0,'MonitorPositions'); F2.OuterPosition = monitors(Monitor_Pos,:);
Figure_Title_2 = ['Sharpe Ratio (' PFM '). Tis = ' num2str(trnwin) '. OOS: ' dates_strng ... 
                  ' (Stdize_Y=' num2str(stdize_Y) ', Demean =' num2str(demean) ')' ];
set(F2,'Name',Figure_Title_2,'WindowState','maximized');

hold on;
  pla = plot(c, tbl_PFM_kmz.Variables ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9);  
  set(pla, {'Color'}, clr_palette);
  hVa_clrs = get(pla, 'Color');
  if ~(strcmp(PFM,'ER') || strcmp(PFM,'alpha'))
    pl0 = plot(c, tbl_PFM_0.Variables   , LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9);  
    set(pl0, {'Color'}, clr_palette);
    % FIX COLOR ORDER to be the same, just a bit darker in the dashed line plots
    for LL = 1:length(hVa_clrs); pl0(LL).Color = hVa_clrs{LL}*DK; end
  end
hold off;

setplotdims([.10 1/5 .82 .77], 3/4); 
xlabel(xLbl,'interpreter','latex')
% leg_name = [cellstr(strcat(strcat('$z=',num2str((lamlist)')),'$'));'$c=1$'];  
% leg_name = strrep(leg_name,'0.6000000000000001','0.6');
if trnwin == 12;   xos = .1; end
if trnwin == 60;   xos = .1; end
if trnwin == 120;  xos = .1; end
xlim([0-xos c(end)+xos])
box on; addgrid
setyticklabels(YTcks_PFM, YTcks_decmls, FNs); % ylim([0 .52]); 
ylim([-.02  max(YTcks_PFM)+.00]); 
% ylim([-.002  max(YTcks_PFM)+.00]); % for alpha
hline(0); tickshrink(.4)
set(gca,'TickLength',[0.02 0.02])

% ADD VERTICAL LINE AT C = 1;
xl = xline(1,'Color', .5*ones(1,3),'LineWidth',2/3,'LineStyle','--');
xl.Alpha = 1;  % Set alpha to 1 (fully opaque)
% xl = xline(8);
yline(maxPFM_0,'Color','m','LineWidth',.66,'LineStyle','-')

% ADD A VERTICAL LINE AT MAX (SR) of correct strategy
decml_fmt   = "%.3g"; 
label_name1 = "max(" + PFM + ") $=" + num2str(maxPFM_0, decml_fmt ) + "$"; 
label_name2 = "[$c=" + num2str(c(pPFM_0), decml_fmt ) + ",~" + tbl_SR_0(pPFM_0,zPFM_0).Properties.VariableNames + "$]";
label_name  = label_name1 + "~" + label_name2;

% max position of labels
PSt = max(ylim)*osw;
if c(pPFM_0) < 8 
  xLa = xline( c(pPFM_0), 'Color', big_CLRS(zPFM_0,:)*DK,'LabelOrientation','horizontal','Alpha',1, ...% 'HandleVisibility', 'off', ...
  'Label', label_name1, 'LabelHorizontalAlignment','right','LineWidth',1.5,'Interpreter','latex','FontSize',FNs-2);
  text(c(pPFM_0)+.15, PSt, 0, label_name2 , 'Color', big_CLRS(zPFM_0,:)*DK, 'Interpreter','latex','FontSize',FNs-3)
  % text(c(pPFM_0)+.15, xLa.Parent.Position(2)-.035, 0, label_name2 , 'Color', big_CLRS(zPFM_0,:)*DK, 'Interpreter','latex','FontSize',FNs-3)
else
  xLa = xline( c(pPFM_0), 'Color', big_CLRS(zPFM_0,:)*DK,'LabelOrientation','horizontal','Alpha',1, ...% 'HandleVisibility', 'off', ...
  'Label', label_name1, 'LabelHorizontalAlignment', 'left','LineWidth',1.5,'Interpreter','latex','FontSize',FNs-2);
  text(c(pPFM_0)-3.85, PSt, 0, label_name2 , 'Color', big_CLRS(zPFM_0,:)*DK, 'Interpreter','latex','FontSize',FNs-3)
  % text(c(pPFM_0)-3.03, xLa.Parent.Position(2)-.035, 0, label_name2 , 'Color', big_CLRS(zPFM_0,:)*DK, 'Interpreter','latex','FontSize',FNs-3)
end
tickshrink(.2);

% YLINE AT MAXSR DO NOT PLOT IT HERE 
% PSt = xLa.Parent.Position(2)+.007;
if (trnwin == 12 && demean == 0)
  text(3, PSt,        0, label_name1 , 'Color', big_CLRS(zPFM_0,:)*DK, 'Interpreter','latex','FontSize',FNs-3)
  text(3, PSt -.035,  0, label_name2 , 'Color', big_CLRS(zPFM_0,:)*DK, 'Interpreter','latex','FontSize',FNs-3)
  % text(3, PSt*osw,  0, label_name2 , 'Color', big_CLRS(zPFM_0,:)*DK, 'Interpreter','latex','FontSize',FNs-3)
end
% subtitle('SR'); % title('SR','FontWeight','normal')
% L0 = line(1*ones(1,2),[-1e9,1e9],'Linewidth',1,'LineStyle','-','color','k'); %'HandleVisibility', 'off');

switch trnwin 
  case 12 
    % BRK = breakxaxis([10.35 989.6], .02); 
    % BRK = breakxaxis([50.35 979.6], .02); 
    % BRK = breakxaxis([104.35 6690.6], .02); 
    BRK = breakxaxis([14.35 990.6], .02); 
    % BRK = breakxaxis([201 430], .02); 
    % set(BRK.leftAxes ,'TickLength',tcklength/2*ones(1,2))
    % set(BRK.rightAxes,'TickLength',tcklength/1*ones(1,2))  
    % end
  case  60 
    BRK = breakxaxis([14.35   190.6], .02); 
    % set(BRK.leftAxes ,'TickLength',tcklength*ones(1,2))
    % set(BRK.rightAxes,'TickLength',tcklength/1*ones(1,2))  
    % end
  case 120
    BRK = breakxaxis([14.35   90.6], .02); 
    % set(BRK.leftAxes ,'TickLength',tcklength*ones(1,2))
    % set(BRK.rightAxes,'TickLength',tcklength/1*ones(1,2))  
end
xlim([0-xos c(end)+xos])
% ADJUST TICKLENGTH
tcklength = 0.01;
% get(BRK.leftAxes ,'TickLength');
% get(BRK.rightAxes,'TickLength')
set(BRK.leftAxes ,'TickLength', tcklength*ones(1,2)*.6 ); get(BRK.leftAxes ,'TickLength');
set(BRK.rightAxes,'TickLength', tcklength*ones(1,2)*.9 ); get(BRK.rightAxes,'TickLength');

% --------------------------------------------------------------------------------------------------
% ADD LEGEND
% --------------------------------------------------------------------------------------------------
%lh = addlegend([pla],leg_name, 5, FNs-3,[],.55);
lh  = addlegend([pla],leg_name, Leg_SR_pos, FNs-4.0,[],.44,[0 0],1.5, Leg_SR_col);
LLW = findobj(lh,'Type','Line'); set(LLW,'LineWidth',2);
% addlegend([LegHndl], leg_name,3,FNs-3,[],2/3,[-0 -31]);

% PRINT TO PDF FILE
% --------------------------------------------------------------------------------------------------
% pdf_name_SR = strcat(OUTPUT_DIR, NAME_SUFFIX,'_SRs_narrow');
pdf_name_PFM = strcat(OOS_EVAL_OUTPUT_DIR, NAME_SUFFIX, '_', PFM, 's');
if PRNT_PDF;  plot2pdf(pdf_name_PFM); end
% tickshrink(.1)
% **************************************************************************************************  


%% Generate Table 1: Comparison with Welch and Goyal (2008) and mkt
% **************************************************************************************************
out = [];
% SUB-SAMPLE SET SUBSAMPLE % note this, only sets the subsample to all non-nan values so gives
% different OOS for T={12,60,120}
I_KMZ = Isel; % Isel or Inan;
% I_KMZ = Inan;
% SubSmpl = find(~isnan(Y'+Yprd_gy(:,1))); 
beg_date_str = num2str(dates(find(I_KMZ,1,'first')));
end_date_str = num2str(dates(end));
% --------------------------------------------------------------------------------------------------
date_range_T1 = strcat(beg_date_str(1:4), '- ',end_date_str(1:4));
disp("Date range for Table I: " + date_range_T1)
% --------------------------------------------------------------------------------------------------
% GW(z=0)
% --------------------------------------------------------------------------------------------------
gy_z0      = timing_gy(I_KMZ,1);
Yprd_gy_z0 =   Yprd_gy(I_KMZ,1);
II_z0      = ones(length(gy_z0),1);
r2tmp_z0   = 1-nanvar( Y(I_KMZ)'-Yprd_gy_z0 ) / nanvar(Y(I_KMZ)');
           % regstats(y, X, model)
stats_z0   = regstats(gy_z0, II_z0,1,'tstat');
% GW(z=0 )vs mkt
stats_z0_vsMKT   = regstats(gy_z0, Y(I_KMZ)','linear',{'tstat','rsquare','r'}); 
tmp1res_z0_vsMKT = stats_z0_vsMKT.tstat.beta(1) + stats_z0_vsMKT.r;
out_z0 = [out ; 100*r2tmp_z0 sqrt(12)*sharpe(gy_z0,0) stats_z0.tstat.t(1) sqrt(12)*sharpe(tmp1res_z0_vsMKT,0) stats_z0_vsMKT.tstat.t(1) nan(1,2) -min(gy_z0) skewness(gy_z0)];
% --------------------------------------------------------------------------------------------------
% GW(z=1000)
% --------------------------------------------------------------------------------------------------
gy_z1000      = timing_gy(I_KMZ,end);  
Yprd_gy_z1000 =   Yprd_gy(I_KMZ,end);
r2tmp_z1000   = 1-nanvar( Y(I_KMZ)'-Yprd_gy_z1000 ) / nanvar(Y(I_KMZ)');
stats_z1000   = regstats(gy_z1000, II_z0,1,'tstat');
% GW(z=1000) vs mkt
stats_z1000_vsMKT   = regstats(gy_z1000, Y(I_KMZ)','linear',{'tstat','rsquare','r'}); 
tmp1res_z1000_vsMKT = stats_z1000_vsMKT.tstat.beta(1) + stats_z1000_vsMKT.r;
% out     = [out ; 100*r2tmp sqrt(12)*sharpe(tmp1,0) stats.tstat.t(1) nan(1,4) min(tmp1) skewness(tmp1)];
out_GW_z1000 = [out_z0 ; 100*r2tmp_z1000 sqrt(12)*sharpe(gy_z1000,0) stats_z1000.tstat.t(1) sqrt(12)*sharpe(tmp1res_z1000_vsMKT,0) stats_z1000_vsMKT.tstat.t(1) nan(1,2) -min(gy_z1000) skewness(gy_z1000)];

% --------------------------------------------------------------------------------------------------
% RFF(z=1000) High complexity
% --------------------------------------------------------------------------------------------------
RFF_timing = nanmean(timing(I_KMZ,end,L1000,:),4);  % R_pi,t  --> average prediction over 1000 RFFS
RFF_Yprd   = nanmean(  Yprd(I_KMZ,end,L1000,:),4);  % pi,t    --> average prediction over 1000 RFFS
r2_RFF     = 1 - nanvar( Y(I_KMZ)'-RFF_Yprd ) / nanvar(Y(I_KMZ)');
stats_RFF  = regstats(RFF_timing, II_z0,1,'tstat');
% RFF(z = 1000) High complexity vs mkt IR
stats_RFF_vsMKT = regstats(RFF_timing, Y(I_KMZ)','linear',{'tstat','rsquare','r'}); 
res_RFF_vsMKT   = stats_RFF_vsMKT.tstat.beta(1) + stats_RFF_vsMKT.r;
% RFF(z = 1000) High complexity vs GW(z=1000)
stats_RFF_vsGWz1000 = regstats(RFF_timing, gy_z1000, 'linear',{'tstat','rsquare','r'}); 
res_RFF_vsGWz1000   = stats_RFF_vsGWz1000.tstat.beta(1) + stats_RFF_vsGWz1000.r;
% --------------------------------------------------------------------------------------------------

% MAKE OUTPUT TABLES LAST ROW FOR 
% --------------------------------------------------------------------------------------------------
% what KMZ do NOTE SHARPE FUNCTION USES MLE VARIANCE that is why sharpe(res_RFF_vsMKT,0) ~= nanmean(res_RFF_vsMKT)./nanstd(res_RFF_vsMKT)
out_RFF   = [ 100*r2_RFF sqrt(12)*sharpe(RFF_timing,0) stats_RFF.tstat.t(1) ...
              sqrt(12)*sharpe(res_RFF_vsMKT,0)         stats_RFF_vsMKT.tstat.t(1) ...
              sqrt(12)*sharpe(res_RFF_vsGWz1000,0)     stats_RFF_vsGWz1000.tstat.t(1) ...
              -min(RFF_timing) skewness(RFF_timing)];

% using my tables above based on average forecasts
out_RFF0  = [ 100*R20(end,end) sqrt(12)*SR0(end,end) SR0_tstat(end,end), ...
              sqrt(12)*IR0(end,end)                  IR0_tstat(end,end)   , ...
              sqrt(12)*IR0_vs_GWz1000(end,end)       IR0_tstat_vs_GWz1000(end,end)   , ...
              maxLoss0(end,end) Skew0(end,end) ];

sep;disp(['Tis = ' num2str(trnwin) ': Their mixed approach. OOS = ' date_range_T1]);sep

% NOW ADD WHAT IS SHOWN IN VoC PLOTS: PARFOR I=1:NSIM FOR RFF(Z=1000)
out_RFF_VoC_plots = [100*R2(end,L1000) sqrt(12)*SR(end,L1000) SR_tstat(end,L1000), ...
                sqrt(12)*IR(end,L1000)                IR_tstat(end,L1000)   , ...
                sqrt(12)*IR_vs_GWz1000(end,L1000)     IR_tstat_vs_GWz1000(end,L1000)   , ...
                maxLoss(end,L1000) Skew(end,L1000)];

out_RFF_all   = [ out_GW_z1000; out_RFF; out_RFF_VoC_plots ];

lst([ out_GW_z1000; out_RFF; nan(1,size(out_RFF,2)); out_RFF_VoC_plots ], '%2.4f')
sep

% MAKE ONE BIG STATS OUTPUT FILE COMARING THE AGGREGATION OUTPUT [~,nC,nZ]  = size(mean_timing);
rowNan = NaN(1,2*nL+3); colNan = NaN(nP,1);
P_allout = [ [NaN(1,2) lamlist NaN lamlist];    rowNan;
  [ Plist' c sqrt(12)*SR colNan sqrt(12)*SR0];  rowNan;
  [ Plist' c sqrt(12)*IR colNan sqrt(12)*IR0];  rowNan;
  [ Plist' c alpha       colNan alpha0];        rowNan;
  [ Plist' c IR_tstat    colNan IR0_tstat];     rowNan;
  [ Plist' c R2          colNan R20];           rowNan;
  [ Plist' c ER          colNan ER0];           rowNan;
  [ Plist' c Vol         colNan Vol0];          rowNan;
  [ Plist' c Bias        colNan Bias0];         rowNan;
  [ Plist' c MSE         colNan MSE0];          rowNan;
  ];

% PRINT TO EXCEL
if PRNT_XLS == 1
  % PRINT SMALL STATS TO EXCEL 
  xls_Table_1_name = strcat(OOS_EVAL_OUTPUT_DIR,NAME_SUFFIX,'-KMZ-Table-I.xlsx');
  writematrix(out_RFF_all,  xls_Table_1_name,'Range','B2');
  % cell names
  xls_col_names = {'R2','SR','t-stat.','IR','t-stat.','IR vs. GW(z=1000)','t-stat.','Max-Loss','Skew'};
  xls_row_names = {'Linear (z=0)';'Linear (z=1000)';'Nonlinear (z=1000)';'Nonlinear VoC Plots'};
  % writecell
  writecell(xls_col_names, xls_Table_1_name,'Range','B1');
  writecell(xls_row_names, xls_Table_1_name,'Range','A2');

  xlsoutName = strcat(OOS_EVAL_OUTPUT_DIR,NAME_SUFFIX,'-KMZ-Table-I-VoC_Full.xlsx');

  % which columns to put the lables in the excel file
  Col_KMZ = number2excel_column(2+4);
  Col_0   = number2excel_column(nL+1+2+4); 

  % PRINT BIG STATS TO EXCEL
  writematrix(P_allout,           xlsoutName,'Range','A1');
  writematrix('SR(KMZ)',          xlsoutName,'Range',Col_KMZ + num2str(0*(nP+1)+2) );
  writematrix('IR(KMZ)',          xlsoutName,'Range',Col_KMZ + num2str(1*(nP+1)+2) );
  writematrix('Alpha(KMZ)',       xlsoutName,'Range',Col_KMZ + num2str(2*(nP+1)+2) );
  writematrix('Alpha-tstat(KMZ)', xlsoutName,'Range',Col_KMZ + num2str(3*(nP+1)+2) );
  writematrix('R2(KMZ)',          xlsoutName,'Range',Col_KMZ + num2str(4*(nP+1)+2) );
  writematrix('ER(KMZ)',          xlsoutName,'Range',Col_KMZ + num2str(5*(nP+1)+2) );
  writematrix('Vol(KMZ)',         xlsoutName,'Range',Col_KMZ + num2str(6*(nP+1)+2) );
  writematrix('Bias(KMZ)',        xlsoutName,'Range',Col_KMZ + num2str(7*(nP+1)+2) );
  writematrix('MSE(KMZ)',         xlsoutName,'Range',Col_KMZ + num2str(8*(nP+1)+2) );

  writematrix('P',                xlsoutName,'Range','A2');
  writematrix('z',                xlsoutName,'Range','B1');
  writematrix('c',                xlsoutName,'Range','B2');
  writematrix('SR',               xlsoutName,'Range',Col_0 + num2str(0*(nP+1)+2) );
  writematrix('IR',               xlsoutName,'Range',Col_0 + num2str(1*(nP+1)+2) );
  writematrix('Alpha',            xlsoutName,'Range',Col_0 + num2str(2*(nP+1)+2) );
  writematrix('Alpha-tstat',      xlsoutName,'Range',Col_0 + num2str(3*(nP+1)+2) );
  writematrix('R2',               xlsoutName,'Range',Col_0 + num2str(4*(nP+1)+2) );
  writematrix('ER',               xlsoutName,'Range',Col_0 + num2str(5*(nP+1)+2) );
  writematrix('Vol',              xlsoutName,'Range',Col_0 + num2str(6*(nP+1)+2) );
  writematrix('Bias',             xlsoutName,'Range',Col_0 + num2str(7*(nP+1)+2) );
  writematrix('MSE',              xlsoutName,'Range',Col_0 + num2str(8*(nP+1)+2) );

  fprintf(' Finished writing to excel ... \n')
end

diary off


 
% -----------------------------------------------------------------------------------------------
% loop though all files up to here

% end; end

fprintf(" ALL DONE \n")

%% Extra stuff
% **************************************************************************************************








% EOF