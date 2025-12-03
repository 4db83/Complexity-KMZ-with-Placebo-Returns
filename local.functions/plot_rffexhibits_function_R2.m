% This file is in Step2_RFFexhibits. FIRST RUN PREDICTIONS_MAIN FROM (STEP 1)
% \EmpiricalAnalysis\Step1_Predictions
clear; clc; 
set(groot,'defaultLineLineWidth',2 ); set(groot,'defaultAxesXTickLabelRotationMode','manual') 
set(groot,'defaultAxesFontSize' ,14); set(groot,'defaultAxesFontName','Times New Roman')
set(groot,'DefaultAxesTitleFontWeight', 'normal')
% set(F1,'Name','GW regressors','WindowState','maximized','Position',[1268 1  1268  2178]);
addpath(genpath('D:/matlab.tools/db.toolbox/db'))% add path to db functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma   = 2;
% trnwin  = [120];
% for trnwin  = [12 60 120]
for trnwin  = [12]
stdize  = 1;
ADD_zGrid = 0;
PRNT_PDF  = 0;
COPY_FILE = 0;

% DIARY LOG FILE OF OUTPUT
log_filename0 = strcat('Log File Tis=', num2str(trnwin),'. Check');
if ADD_zGrid; log_filename0 = [log_filename0, ' Finer zGrid ']; end
diary off;      % Turn off diary if it's currently on
%**************************************************************************
% Choices of parameters
%**************************************************************************
if gamma == 0.5; gamma_str = '0pt5'; else; gamma_str = num2str(gamma);end

% Demeaning = False
demean  = 0;
% save the results
saveon  = 0;
% Use subsample set to zero
subsamp = 0;
% max number of Random Fourier Features (RFFs)
maxP = 12000;
% the number of simulations is 1000
nSim = 1000; 

% saving string
para_str = strcat('maxP-', num2str(maxP), '-trnwin-', num2str(trnwin), '-gamma-', num2str(gamma), '-stdize-', num2str(stdize), '-demean-', num2str(demean), '-v2');

% Local path
RFF_PATH = 'tryrff_v2_SeparateSims';       % default, same as KMZ: lamlist = 
% RFF_PATH = 'tryrff_v2_SeparateSims_zGrid';
% is_zGrid =  contains(RFF_PATH,'zGrid'); % check which data series is used. 
if ADD_zGrid; RFF_PATH = strcat(RFF_PATH, '_zGrid'); end

% gybench_datadir = '../Step1_Predictions/tryrff_v2_SeparateSims/';
gybench_datadir = strcat('../Step1_Predictions/', RFF_PATH, '/');
figdir          = strcat('./RFF_Empirical_figures/', para_str,'/');
% datadir         = '../Step1_Predictions/tryrff_v2_SeparateSims/';
datadir         = strcat('../Step1_Predictions/', RFF_PATH, '/');
save_path       = strcat(datadir, para_str);
combined_data_save_path = './combined_data/';

if ADD_zGrid; combined_data_save_path = './combined_data_zGrid/'; end

% BUILD THE OUTPUT FOLDERS
% mkdir(figdir);
if ~exist(combined_data_save_path); mkdir(combined_data_save_path); end

%**************************************************************************
% Load parameters and benchmark
%**************************************************************************
% filename = strcat([save_path '/iSim1.mat']);
% load the baseline Sim1 data which inludes more variables
% disp('Loading Sim1 ...')
load( strcat([save_path '/iSim1.mat']), 'T','nP','nL','Y', 'Plist', 'log_lamlist','dates','lamlist' );
% disp('Done: Loading Sim1 ...')
% note that this path does not contain the full list of files after combining the data ... 
files_listing   = dir([save_path '/*.mat']);

%**************************************************************************
% Collect results of 1000 simulations
%**************************************************************************
% load the benchmark of Welch and Goyal (2008) "kitchen sink" regression
disp('Loading Goyal Welch Benchmark...')
load([gybench_datadir 'gybench-trnwin-' num2str(trnwin) '-stdize-' num2str(stdize)  '-demean-' num2str(demean) '.mat'])
disp('Done: Loading Goyal Welch Benchmark ...')

% Y_B_file_save = ['trnwin-' num2str(trnwin) '-gamma-' num2str(gamma) '-stdize-' num2str(stdize) '-demean-' num2str(demean) '-Y-B'];
% filename = strcat(combined_data_save_path, Y_B_file_save, '.mat');

combine_out_file_save = ['trnwin-' num2str(trnwin) '-gamma-' num2str(gamma) '-stdize-' num2str(stdize) '-demean-' num2str(demean) '-Yprd_timing'];
comb_filename = strcat(combined_data_save_path, combine_out_file_save, '.mat'); 

% loop through the individual simulations results and combine all results into 4-D arrays
if isfile(comb_filename)
    % if the combined result exists
    if ADD_zGrid; disp('Loading existing combined 4D data files, ... this takes about 60 seconds ... ');
    else          disp('Loading existing combined 4D data files, ... this takes about 10 seconds ... ')
    end
    tic; load(comb_filename); toc
else
    % if the combined result doesn't exist
    disp('--------------------------------------------------------------------------------------\n');
    disp('Reading in seperate individual files, ... this takes about 25 seconds ... ')
    tic
    Yprd_collect    = nan(T,nP,nL,nSim); % predicted Y
    Bnrm_collect    = nan(T,nP,nL,nSim); % beta norm  
    for s = 1:nSim
        % tic
        if mod(s, 100) == 0 ; disp( num2str(s, ' %4d') ); end 
        % load data
        filename = strcat([save_path '/' files_listing(s).name]);
        load(filename, 'Yprd', 'Bnrm');

        Yprd_collect(:,:,:,s)   = Yprd;
        Bnrm_collect(:,:,:,s)   = Bnrm;
    end
    toc
    clearvars Yprd Bnrm
    Yprd = Yprd_collect;
    Bnrm = Bnrm_collect;

    Bnrmbar     = nanmean(Bnrm,4);  % nanmedian(Bnrm,4);
    Bnrmbar     = squeeze(nanmean(Bnrmbar,1)); % squeeze(nanmedian(Bnrmbar,1));
    timing      = Yprd.*Y';
    pihat       = [squeeze(nanmean(Yprd(:,nP,nL,:),4)) squeeze(nanmean(Yprd(:,1,1,:),4))]; % this is for P=12000 P=2 or 0? whatever the first entry is
    
    % clearvars 
    clearvars Yprd_collect Bnrm_collect
    % ADD -nocompression TAG, much faster, only a litte larger: takes 10 secs... 
    disp(['Saving data to combined data to ' comb_filename] )
    disp(' ... this takes about 10 seconds ... ')
    % SAVING Yprd, timing, 1092x46x7x1000 matrics 
    tic; save(comb_filename, 'Yprd', 'timing', 'Bnrmbar', '-v7.3', '-nocompression'); toc
end
disp('Data Loaded')
% PERCENTILE LIST
pctlist     = [1 2.5 5 25 50 75 95 97.5 99];
% Generate exhibits
%**************************************************************************
% Portfolio Evaluation: Full Sample
%**************************************************************************
if subsamp==1
    subbeg = [1926 1926 1975];
    subend = [2020 1974 2020];
else
    subbeg = 1926;
    subend = 2020;
end
nPct = length(pctlist);
nSub = length(subbeg);
% RESTRICT OUT-OF-SAMPLE EVALUATION PERIOD: 1940 onwards ∀ Tᵢₛ ∈ [12 60 120].
% subbeg = 1940;  %  
% subbeg = 1975;

%% for subbeg = [1926 1940 1975]
% MAIN LOOP of KMZ performance evaluation for each simulated series ( and Plist, lamlis values)
% for ss=1:nSub
ss = 1; 
% Evaluation period
locev       = find(dates>=subbeg(ss)*100 & dates<=(subend(ss)+1)*100);
% Suffix
suffix      = ['trnwin-' num2str(trnwin) '-gamma-' num2str(gamma) '-stdize-' num2str(stdize) '-demean-' num2str(demean) '-' num2str(subbeg(ss)) '-' num2str(subend(ss))];
performance_filename = strcat(combined_data_save_path, suffix, '_KMZ_port_eval.mat');
ss_dates = char(num2str(dates(locev)));
% first non-nan evaluation periods date
date0 = num2str(dates(find_first_non_nan(timing(:,1,1,1))));
dates_strng = [ss_dates(1,1:end-2) ':' ss_dates(1,end-1:end) ' - ' ss_dates(end,1:end-2) ':' ss_dates(end,end-1:end)];
if str2double(date0(1:4)) > str2double(dates_strng(1:4)); dates_strng(1:4) = date0(1:4); end

% start logging the output DIARY
log_filename = [log_filename0, ' ', strrep(strrep(dates_strng,':','.'),' ',''), '.txt']; 
if exist(log_filename, 'file'); delete(log_filename); end
diary(log_filename)

sep;  fprintf('Portfolio Evaluation period: %s. ', dates_strng) 
      fprintf('Sample size Tis =  %d\n', trnwin ); sep

% start new parpool
start_parpool_with;

% GENERATE PERFORMANCE MEASUREMENTS ---------------------------------------------------------------
  if isfile(performance_filename)
  % if 0
    disp('loading data KMZ_port_eval file, ... this takes about 1 seconds ... ')
    tic; load(performance_filename); toc;
  else
   disp('--------------------------------------------------------------------------------------\n');
   disp('Computing portfolio performance measures from separate simulations, ... this takes about 248 seconds... ')
  % Performance initialization
    ER          = nan(nP,nL);
    SR          = nan(nP,nL);
    Vol         = nan(nP,nL);
    IR          = nan(nP,nL);
    IRt         = nan(nP,nL);
    alpha       = nan(nP,nL);
    R2          = nan(nP,nL);
    IRstd       = nan(nP,nL);
    Ytmp        = Y(locev)';
    II          = ones(length(Ytmp),1);
    % add extra measures
    maxLoss     = nan(nP,nL);
    Skew        = nan(nP,nL);
    SR_tstat    = nan(nP,nL);
    IR_vs_GWz1000 = nan(nP,nL);
    IR_tstat_vs_GWz1000 = nan(nP,nL);
    for p=1:nP
      for l=1:nL
        ERtmp         = nan(nSim,1);
        Voltmp        = nan(nSim,1);
        IRtmp         = nan(nSim,1);
        IRttmp        = nan(nSim,1);
        alphatmp      = nan(nSim,1);
        R2tmp         = nan(nSim,1);
        IRstdtmp      = nan(nSim,1);
        timtmp        = squeeze(timing(locev,p,l,:));
        Yprdtmp       = squeeze(Yprd(locev,p,l,:));
        % add extra measures
        maxLosstmp    = nan(nSim,1);
        Skewtmp       = nan(nSim,1);
        SR_tstattmp   = nan(nSim,1);
        IR_vs_GWz1000_tmp       = nan(nSim,1);
        IR_tstat_vs_GWz1000_tmp = nan(nSim,1);
        parfor i=1:nSim
          stats       = regstats(timtmp(:,i),Ytmp,'linear',{'tstat','r'});
          SRtmp(i)    = sharpe(timtmp(:,i),0);
          ERtmp(i)    = nanmean(timtmp(:,i));
          Voltmp(i)   = nanstd(timtmp(:,i));
          IRtmp(i)    = stats.tstat.beta(1)/nanstd(stats.r);
          IRttmp(i)   = stats.tstat.t(1);
          alphatmp(i) = stats.tstat.beta(1);
          IRstdtmp(i) = nanstd(stats.r);
          loc         = find(~isnan(Yprdtmp(:,i)+Ytmp));
          R2tmp(i)    = 1-var(Yprdtmp(loc,i)-Ytmp(loc),'omitnan')/var(Ytmp(loc),'omitnan');   
          % add extra measures
          maxLosstmp(i)  = -nanmin(timtmp(:,i));
          Skewtmp(i)     = skewness(timtmp(:,i));
          stats_SRtmp    = regstats(timtmp(:,i), II, 1, {'tstat'}); % add a one to not include a constant term in the regressions
          SR_tstattmp(i) = stats_SRtmp.tstat.t(1);
          % now add IR vs GW(z=
          stats_vs_GWz1000 = regstats(timtmp(loc,i),timing_gy(loc,end), 'linear',{'tstat','rsquare','r'}); 
          IR_vs_GWz1000_tmp(i) = sharpe( stats_vs_GWz1000.tstat.beta(1) + stats_vs_GWz1000.r, 0);
          IR_tstat_vs_GWz1000_tmp(i) = stats_vs_GWz1000.tstat.t(1);
        end
        % Average them now over the 1000 Nsims
        SR(p,l)       = nanmean(SRtmp);
        ER(p,l)       = nanmean(ERtmp);
        Vol(p,l)      = nanmean(Voltmp);
        IR(p,l)       = nanmean(IRtmp);
        IRt(p,l)      = nanmean(IRttmp);
        alpha(p,l)    = nanmean(alphatmp);
        R2(p,l)       = nanmean(R2tmp);
        IRstd(p,l)    = nanmean(IRstdtmp);
        % add extra measures
        maxLoss(p,l)  = nanmean(maxLosstmp);
        Skew(p,l)     = nanmean(Skewtmp);
        SR_tstat(p,l) = nanmean(SR_tstattmp);
        IR_vs_GWz1000       = nanmean(IR_vs_GWz1000_tmp);
        IR_tstat_vs_GWz1000 = nanmean(IR_tstat_vs_GWz1000_tmp);
      end
      disp(['p=' num2str(p)])
    end
    toc
    % SAVE THE DATA
    % --------------------------------------------------------------------------------------------
    save( performance_filename, 'locev','SR','ER','Vol','IR','IRt','alpha','R2','IRstd', ...
        'maxLoss', 'Skew', 'SR_tstat', 'IR_vs_GWz1000', 'IR_tstat_vs_GWz1000', ... 
        '-v7.3', '-nocompression');
    % --------------------------------------------------------------------------------------------
    % % 'SRpct','ERpct','volpct','IRpct','IRtpct','alphapct','R2pct',...
    % % 'Bnrmbar','timing', 'pihat', 'Yprd', '-v7.3', '-nocompression');
  end
% end  

%% COMPUTE PORTFOLIO PERFORMANCE MEASURES AS OVER MEANS OF RFF FORECASTS
  % ss_dates = char(num2str(dates(locev)));
  % % [ss_dates(1,1:end-2) ':' ss_dates(1,end-1:end)] - [ss_dates(end,1:end-2) ':' ss_dates(end,end-1:end)]
  % dates_strng = [ss_dates(1,1:end-2) ':' ss_dates(1,end-1:end) ' - ' ss_dates(end,1:end-2) ':' ss_dates(end,end-1:end)];
  % fprintf(' Portfolio Evaluation period: %s. \n', dates_strng ) 
  dates_oos = dates(locev,:);
  mean_timing = mean(timing(locev,:,:,:),4);
  mean_Yprd   = mean(Yprd(  locev,:,:,:),4);
  Ytmp        =         Y(  locev)';
  [~,nC,n3]  = size(mean_timing); 
  SR0     = nan(nC,n3);  
  ER0     = nan(nC,n3); 
  Vol0    = nan(nC,n3); 
  IR0     = nan(nC,n3); 
  IRt0    = nan(nC,n3); 
  alpha0  = nan(nC,n3); 
  R20     = nan(nC,n3); 

  for jj = 1:n3
    SR0(:,jj)     = sharpe( mean_timing(:,:,jj), 0)';
    ER0(:,jj)     = nanmean(mean_timing(:,:,jj)   )';
    Vol0(:,jj)    = nanstd( mean_timing(:,:,jj)   )';
    for ii = 1:nC
      stats_ii    = regstats(mean_timing(:,ii,jj),Ytmp,'linear',{'tstat','r'});
      IR0(ii,jj)  = stats_ii.tstat.beta(1)/nanstd(stats_ii.r);
      IRt0(ii,jj) = stats_ii.tstat.t(1);
      alpha0(ii,jj) = stats_ii.tstat.beta(1);
      R20(ii,jj)    = 1-var(mean_Yprd(:,ii,jj)-Ytmp,'omitnan')/var(Ytmp,'omitnan');
      % R2tmp(i)      = 1-var( Yprdtmp(loc,i)-Ytmp(loc) ,'omitnan')/var(Ytmp(loc),'omitnan'); 
    end
  end

% MAKE A TABLE FROM ARRAYS AND FIND MAX VALUES FOR EACH ONE
  zVarNames = strrep(cellstr(strcat('z=',char(num2str(lamlist')))),' ','');
  zVarNames = strrep(zVarNames,'z=0.6000000000000001','z=0.6');
  pNames    = strrep(cellstr(strcat('P=',char(num2str(Plist')))),' ','');

  if ADD_zGrid 
    % this adds the KMZ suffix to the zVarnmaes, to distinguish KMZ from these
    % zVarNames(1:7) = strcat(zVarNames(1:7),'(KMZ)');
    % REMOVE THE KMZ ENTRIES FOR 0,1,10,100, AND 1000 to avoid double plotting them
    N0 = [1:3 9:length(lamlist)]; 
    zVarNames = zVarNames(N0); 
    SR0    = SR0   (:,N0);      SR    = SR   (:,N0);
    ER0    = ER0   (:,N0);      ER    = ER   (:,N0);
    Vol0   = Vol0  (:,N0);      Vol   = Vol  (:,N0);
    IR0    = IR0   (:,N0);      IR    = IR   (:,N0);
    IRt0   = IRt0  (:,N0);      IRt   = IRt  (:,N0);
    alpha0 = alpha0(:,N0);      alpha = alpha(:,N0);
    R20    = R20   (:,N0);      R2    = R2   (:,N0);
  end
  % MAKE A TABLE 
  SR_m      = array2table(sqrt(12)*SR0  ,'VariableNames',zVarNames,'RowNames', pNames);
  ER_m      = array2table( ER0          ,'VariableNames',zVarNames,'RowNames', pNames);
  Vol_m     = array2table( Vol0         ,'VariableNames',zVarNames,'RowNames', pNames);
  IR_m      = array2table(sqrt(12)* IR0 ,'VariableNames',zVarNames,'RowNames', pNames);
  IRt_m     = array2table( IRt0         ,'VariableNames',zVarNames,'RowNames', pNames);
  alpha_m   = array2table( alpha0       ,'VariableNames',zVarNames,'RowNames', pNames);
  R2_m      = array2table( R20          ,'VariableNames',zVarNames,'RowNames', pNames);
  
  SR_kmz    = array2table(sqrt(12)*SR   ,'VariableNames',zVarNames,'RowNames', pNames);
  ER_kmz    = array2table( ER           ,'VariableNames',zVarNames,'RowNames', pNames);
  Vol_kmz   = array2table( Vol          ,'VariableNames',zVarNames,'RowNames', pNames);
  IR_kmz    = array2table(sqrt(12)*IR   ,'VariableNames',zVarNames,'RowNames', pNames);
  IRt_kmz   = array2table( IRt          ,'VariableNames',zVarNames,'RowNames', pNames);
  alpha_kmz = array2table( alpha        ,'VariableNames',zVarNames,'RowNames', pNames);
  R2_kmz    = array2table( R2           ,'VariableNames',zVarNames,'RowNames', pNames);
  
  %% 
  KMZ_strng_SR = '';
  clc
  % Sharpe Ratios (SRs)
  [pSR, zSR] =find(max(max(SR)) ==SR);
  [pSR0,zSR0]=find(max(max(SR0))==SR0);
  maxSR   = max(max(SR_kmz.Variables));
  maxSR0  = max(max(SR_m.Variables));
  % max max
  mmaxSR = sqrt(12)*max(max(max(SR)),max(max(SR0)));
  if max(max(SR)) > max(max(SR0))
    pSR_ = pSR; zSR_ = zSR; KMZ_strng_SR = strcat(' (KMZ)',KMZ_strng_SR);
  else;  pSR_ = pSR0; zSR_ = zSR0; end
  sep
  fprintf('Max SR(KMZ) = %4.6f. \n', sqrt(12)*max(max(SR))  )
  fprintf('Max SR(0)   = %4.6f. \n', sqrt(12)*max(max(SR0)) )
  % print_table(SR_m(pSR0,zSR0),6,1,' Max Sharpe Ratios (SRs)')
  % SR_m2 = table_max(SR_m,'SR(0)');
  % print_table(sortrows(SR_m2,'SR(0)','descend'), 6, 1,[],3);
  % PRINT: Sharpe Ratios (SRs)
% ------------------------------------------------------------------------------------------------    
  % nd = 60;
  % fprintf(strcat(repmat('-',1, nd),'\n'))
  % fprintf('             z        SR(KMZ)   SR0 \n')
  % fprintf(strcat(repmat('-',1, nd),'\n'))
  % % [sqrt(12)*max(SR)' sqrt(12)*max(SR0)'] Takes it row-wise
  % fprintf(' %18.3f   %4.6f    %4.6f \n', [lamlist(N0); sqrt(12)*max(SR(:,N0)); sqrt(12)*max(SR0)])
  % fprintf(strcat(repmat('-',1, nd),'\n'))
% ------------------------------------------------------------------------------------------------  
  % % MAX IR 
  [pIR, zIR] =find(max(max(IR))==IR);
  [pIR0,zIR0]=find(max(max(IR0))==IR0);
  mmaxIR = sqrt(12)*max(max(max(IR)),max(max(IR0)));
  if max(max(IR)) > max(max(IR0)); pIR_ = pIR; zIR_ = zIR; else; pIR_ = pIR0; zIR_ = zIR0; end
  sep
  fprintf('Max IR(KMZ) = %4.6f. \n', sqrt(12)*max(max(IR))  )
  fprintf('Max IR(0)   = %4.6f. \n', sqrt(12)*max(max(IR0)) )
  % print_table(IR_m(pIR0,zIR0),6,1,'Max IR');
  % IR_m2 = table_max(IR_m,'IR');
  % print_table(sortrows(IR_m2,'IR','descend'), 6 ,1,[],3)
% ------------------------------------------------------------------------------------------------
  % max alpha 
  mmax_alpha = max( max(max(alpha)), max(max(alpha0)) );
  [pa0,za0]   = find(max(max(alpha0))==alpha0);
  sep
  fprintf('Max alpha(KMZ) = %4.6f. \n', max(max(alpha))  )
  fprintf('Max alpha(0)   = %4.6f. \n', max(max(alpha0)) )
  % print_table(alpha_m(pa0,za0),6,1,'Max Alpha');
  % alpha_m2 = table_max(alpha_m,'alpha');
  % print_table(sortrows(alpha_m2,'alpha','descend'), 6 ,1,[],3)
% ------------------------------------------------------------------------------------------------
  % max alpha t-stat IRt 
  mmax_IRt = max( max(max(IRt)), max(max(IRt0)) );
  [pat0,zat0]   = find(max(max(IRt0))==IRt0);
  sep
  fprintf('Max alpha-tstat(KMZ) = %4.6f. \n', max(max(IRt))  )
  fprintf('Max alpha-tstat(0)   = %4.6f. \n', max(max(IRt0)) )
  % print_table(IRt_m(pat0,zat0),6,1,'Max Alpha t-stat');
  % IRt_m2 = table_max(IRt_m,'alpha t-stat');
  % print_table(sortrows(IRt_m2,'alpha t-stat','descend'), 6 ,1,[],3)
% --------------------------------------------------------------------------------------------------  
  % max ER
  mmax_ER = max( max(max(ER)), max(max(ER0)) );
  [pER0,zER0]   = find(max(max(ER0))==ER0);
  sep
  fprintf('Max ER(KMZ) = %4.6f. \n', max(max(ER))  )
  fprintf('Max ER(0)   = %4.6f. \n', max(max(ER0)) )
  % print_table(ER_m(pER0,zER0),6,1,'Max E(R)');
  % ER_m2 = table_max(ER_m,'E(R)');
  % print_table(sortrows(ER_m2,'E(R)','descend'), 6 ,1,[],3)
% ------------------------------------------------------------------------------------------------
  % max R2
  KMZ_strng_R2 = '';
  mmax_R2 = max( max(max(R2)), max(max(R20)) );
  [pR2,zR2]     = find(max(max(R2))==R2);
  [pR20,zR20]   = find(max(max(R20))==R20);
  if max(max(R2)) > max(max(R20)); pR2_ = pR2; zR2_ = zR2; KMZ_strng_R2 = strcat(' (KMZ)',KMZ_strng_R2); else; pR2_ = pR20; zR2_ = zR20; end
  sep
  fprintf('Max R2(KMZ) = %4.6f. \n', max(max(R2))  )
  fprintf('Max R2(0)   = %4.6f. \n', max(max(R20)) )
  if max(max(R2)) > max(max(R20));  print_table(R2_kmz(pR2,zR2),6,1,'KMZ''s Max R²');
  else;                             print_table(R2_m(pR20,zR20),6,1,'Max R²'); end

% PLOTS 
CLRS = [0         0.4470    0.7410;
        0.8500    0.3250    0.0980;
        0.9290    0.6940    0.1250;
        0.4940    0.1840    0.5560;
        0.4660    0.6740    0.1880;
        0.3010    0.7450    0.9330;
        0.6350    0.0780    0.1840];
 
% some plotting controls  
for Xlim0 = [0] 
XlimUp = 25;  
PLOT_ALL = 1;
PLOT_INDIVIDUAL = 1;
% linewidth(s) of the two lines
LW  = 2.5; MK   = 'none';
LW0 = 2.0; MK0  = 'none';
LType = '-.'; c = (Plist/trnwin)'; xLbl = '$c$'; 
DK = .7;

if PLOT_ALL 
  FNs = 16; 
  clf; F1 = figure(1); set(F1, 'Name',['Evaluation Measures. Tis = ' num2str(trnwin) '. OOS: ' dates_strng], 'WindowState','maximized','Position',[1268 1  1268  2178]);
  tiledlayout(5,2,'TileSpacing', 'compact', 'Padding', 'loose');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nexttile(1) % 1 SR
  hold on;
  if ~ADD_zGrid 
    pp1 = plot(Plist/trnwin,sqrt(12)*SR ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
    hVa_clrs = get(pp1, 'Color');
  end
    pp1 = plot(Plist/trnwin,sqrt(12)*SR0, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9); 
  hold off;
  big_CLRS = reshape([pp1.Color],3,length(pp1))';
  for LL = 1:length(hVa_clrs)
    pp1(LL).Color = hVa_clrs{LL}*DK;
  end
  xlabel(xLbl,'interpreter','latex')
  if Xlim0;  xlim([0 XlimUp]); end
  box on; addgrid
  setyticklabels(0:.1:.5, 2, FNs)
  subtitle('SR')
  if Xlim0; line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', 1.5); end
  % ADD A VERTICAL LINE AT MAX(SR0)
  if ~ADD_zGrid
    xL1 = xline(c(pSR_), 'Color', CLRS(zSR_,:)*DK, ... % 'HandleVisibility', 'off', ...
      'Label', strcat('max(SR) $=', num2str(mmaxSR,'%2.4f'),' ~[c = ', num2str(c(pSR0(1))),',', SR_m(pSR0,zSR0).Properties.VariableNames,'$]'), 'LabelOrientation','horizontal','Alpha',1, ...
      'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment', 'center','LineWidth',1,'Interpreter','latex','FontSize',FNs-2);
    if trnwin == 120; xL1.LabelHorizontalAlignment = 'right'; end
    % if trnwin == 12; xL1.LabelHorizontalAlignment = 'center'; end      
  else 
    % strcat('max(SR) $=', num2str(mmaxSR,4),', ~ c = ', num2str(c(pSR0(1))), '$.')
    xL1 = xline(c(pSR_), 'Color', big_CLRS(zSR_,:), ... %'HandleVisibility', 'off', ...
      'Label', strcat('max(SR) $=', num2str(mmaxSR,'%2.4f'),' ~[c = ', num2str(c(pSR0(1))),',', SR_m(pSR0,zSR0).Properties.VariableNames,'$]'), 'LabelOrientation','horizontal','Alpha',1, ...
      'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment', 'right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2);
      if trnwin == 12; xL1.LabelHorizontalAlignment = 'left'; end
      % xL1.LabelHorizontalAlignment = 'right';
      yline(mmaxSR,'Color',0*ones(1,3),'LineWidth',1,'LineStyle','-')
  end
  yline(mmaxSR,'Color','m','LineWidth',.66,'LineStyle','-')
  if trnwin == 60; xL1.LabelHorizontalAlignment = 'right'; end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  nexttile(3) % 3 alpha 
  hold on;
  if ~ADD_zGrid
     % plot(Plist/trnwin,alpha , '- ','linewidth', LW , 'Marker', MK , 'Markersize', 9); 
  end
     plot(Plist/trnwin,alpha0, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9); 
  hold off;
  xlabel(xLbl,'interpreter','latex')
  if Xlim0;  xlim([0 XlimUp]); end
  box on; addgrid
  setyticklabels([0:.01:.03], 2, FNs); %ylim([0 .025])
  subtitle('Alpha')
  % line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', 1.5);
  if Xlim0; line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', 1.5); end
  % ADD A VERTICAL LINE AT 
  if ~ADD_zGrid
  xL2 = xline(c(pa0),'Color',CLRS(za0,:),'LabelOrientation','horizontal','Alpha',1, ...
   'Label',strcat('$\max(\alpha) =',num2str(mmax_alpha,'%2.4f'),' ~[c = ', num2str(c(pa0),'%2.1f'),',', alpha_m(pa0,za0).Properties.VariableNames,'$]'), ...
   'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2);  
  if trnwin == 12; xL2.LabelHorizontalAlignment = 'left'; end
  if trnwin == 60; xL2.LabelHorizontalAlignment = 'right'; end
  else
  xL2 = xline(c(pa0),'Color',big_CLRS(za0,:),'LabelOrientation','horizontal','Alpha',1, ...
   'Label',strcat('$\max(\alpha) =',num2str(mmax_alpha,'%2.4f'),' ~[c = ', num2str(c(pa0),'%2.1f'),',', alpha_m(pa0,za0).Properties.VariableNames,'$]'), ...
   'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2);  
    if trnwin == 12; xL2.LabelHorizontalAlignment = 'left'; end
  end
  yline(mmax_alpha,'Color','m','LineWidth',.66,'LineStyle','-')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  nexttile(2) % 2 IR
  hold on;
  if ~ADD_zGrid 
    % plot(Plist/trnwin,sqrt(12)*IR ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
    pp1 = plot(Plist/trnwin,sqrt(12)*IR ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
    hVa_clrs = get(pp1, 'Color');
  end
    % plot(Plist/trnwin,sqrt(12)*IR0, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9); 
    pp1 = plot(Plist/trnwin,sqrt(12)*IR0, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9); 
  hold off;
  for LL = 1:length(hVa_clrs)
    pp1(LL).Color = hVa_clrs{LL}*DK;
  end
  xlabel(xLbl,'interpreter','latex')
  if Xlim0;  xlim([0 XlimUp]); end
  box on; addgrid
  setyticklabels(0:.1:.4, 1, FNs); ylim([0 .33])
  subtitle('Information Ratio')
  % line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', 1.5);
  if Xlim0; line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', 1.5); end
  % ADD A VERTICAL LINE AT MAX(SR0)
  if ~ADD_zGrid
    xL3 = xline(c(pIR_), 'Color', CLRS(zIR_,:)*DK, 'LabelOrientation','horizontal','Alpha',1, ...
      'Label', strcat('max(IR) $=', num2str(mmaxIR,'%2.4f'),' ~[c = ', num2str(c(pIR0(1)),'%2.1f'),',', IR_m(pIR0,zIR0).Properties.VariableNames,'$]'), ...
      'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment', 'right', 'LineWidth',1,'Interpreter','latex','FontSize',FNs-2);  
    if trnwin == 120; xL3.LabelHorizontalAlignment = 'right'; end
    if trnwin == 12;  xL3.LabelHorizontalAlignment = 'left';  end
  else
    xL3 = xline(c(pIR_(1)),'Color',big_CLRS(zIR_(1),:),'LabelOrientation','horizontal','Alpha',1, ...
      'Label',strcat('max(IR) $=',num2str(mmaxIR,'%2.4f'),' ~[c = ', num2str(c(pIR0(1)),'%2.1f'),',', IR_m(pIR0,zIR0).Properties.VariableNames,'$]'), ...
      'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2);  
    if trnwin == 12; xL3.LabelHorizontalAlignment = 'left'; end
    % xL2.LabelHorizontalAlignment = 'left';    
    % yline(mmaxIR,'Color',0*ones(1,3),'LineWidth',1,'LineStyle','-')
  end
  yline(mmaxIR,'Color','m','LineWidth',.66,'LineStyle','-')
  % if trnwin == 60; xL3.LabelHorizontalAlignment = 'right'; end
  % xlim([0 2])
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  nexttile(4) % 4 alpha t-stat
  hold on;
  if ~ADD_zGrid
    pp1 = plot(Plist/trnwin,IRt ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9);  
    hVa_clrs = get(pp1, 'Color');
  end
    pp1 = plot(Plist/trnwin,IRt0, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9);  
  hold off;
  for LL = 1:length(hVa_clrs)
    pp1(LL).Color = hVa_clrs{LL}*DK;
  end
  xlabel(xLbl,'interpreter','latex')
  if Xlim0;  xlim([0 XlimUp]); end
  box on; addgrid 
  setyticklabels([0:1:3], 0, FNs); 
  subtitle('Alpha $t$-statistic','Interpreter','Latex')
  % line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', 1.5);
  if Xlim0; line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', 1.5); end
  % ADD A VERTICAL LINE AT 
  if ~ADD_zGrid
  xL4 = xline(c(pat0),'Color',CLRS(zat0,:)*DK,'LabelOrientation','horizontal','Alpha',1, ...
   'Label',strcat('max(alpha $t-$stat.) $ =',num2str(mmax_IRt,'%2.4f'),' ~[c = ', num2str(c(pat0),'%2.1f'),',', alpha_m(pat0,zat0).Properties.VariableNames,'$]'), ...
   'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2);  
  if trnwin == 12; xL4.LabelHorizontalAlignment = 'left'; end
  else
  xL4 = xline(c(pat0),'Color',big_CLRS(zat0,:),'LabelOrientation','horizontal','Alpha',1, ...
   'Label',strcat('max(alpha $t-$stat.) $ =',num2str(mmax_IRt,'%2.4f'),' ~[c = ', num2str(c(pat0),'%2.1f'),',', alpha_m(pat0,zat0).Properties.VariableNames,'$]'), ...
   'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2);  
    if trnwin == 12; xL4.LabelHorizontalAlignment = 'left'; end
    if trnwin == 12; xL4.LabelVerticalAlignment = 'bottom'; end
  end
  yline(mmax_IRt,'Color','m','LineWidth',.66,'LineStyle','-')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  nexttile(5) % 5 R2
  % hold on;
  % if ~ADD_zGrid
  %   pl0 = plot(Plist/trnwin,R2  ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
  % end
  %         plot(Plist/trnwin,R20 , LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9); 
  % hold off; hline(0)

  hold on; % scale/annualize vola by 12, KMZ don't do this, so leave it at 1
  % if ~ADD_zGrid
  % yyaxis left 
  hVa =  plot(Plist/trnwin, R2 ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
    % setyticklabels([0:1:5], 0, FNs); % ylim([0 .5])
    setyticklabels([-.5:.1:.2], 1, FNs); ylim([-.5 .25])
    hVa_clrs = get(hVa, 'Color');
    hline(0)
  yyaxis right
  hV0 =  plot(Plist/trnwin, R20, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9); 
    % ylim([-.05 .05])
    set(gca, 'YColor', 'k') % Set right y-axis color to black
    % setyticklabels([0:.1:3], 0, FNs); % ylim([0 .5])
    for LL = 1:length(hVa_clrs)
      hV0(LL).Color = hVa_clrs{LL}*DK;
    end
    ylim([-.04 .02])
    hline(0)
  hold off;

  xlabel(xLbl,'interpreter','latex')
  if Xlim0;  xlim([0 XlimUp]); end
  box on; addgrid 
  % setyticklabels([-.2:.05:.0], 2, FNs); ylim([-.2 .025])
  subtitle('$R^2$','Interpreter','Latex')
  if Xlim0; line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', 1.5); end
  % ADD A VERTICAL LINE AT 
  if ~ADD_zGrid
  R2_CLRS = CLRS(zSR_,:)*DK; 
  if trnwin == 12; R2_CLRS = CLRS(zSR_,:); end
  xL5 = xline(c(pR2_),'Color',R2_CLRS,'LabelOrientation','horizontal','Alpha',1, ...
   'Label',strcat('$\max(R^2) =',num2str(mmax_R2,'%2.4f'),' ~[c = ', num2str(c(pR2_),'%2.0f'),',', alpha_m(pR2_,zR2_).Properties.VariableNames,'$]',KMZ_strng_R2), ...
   'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','LineWidth',1,'Interpreter','latex','FontSize',FNs-2);  
    if ismember(trnwin, [60, 120]); xL5.LabelHorizontalAlignment = 'right'; end
    % if trnwin == 60; xL5.LabelHorizontalAlignment = 'right'; end
  else
  xL5 = xline(c(pR2_),'Color',big_CLRS(zR2_,:),'LabelOrientation','horizontal','Alpha',1, ...
   'Label',strcat('$\max(R^2) =',num2str(mmax_R2,'%2.4f'),' ~[c = ', num2str(c(pR2_),'%2.0f'),',', alpha_m(pR2_,zR2_).Properties.VariableNames,'$]',KMZ_strng_R2), ...
   'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2);  
    if trnwin == 12; xL5.LabelHorizontalAlignment = 'left'; end
  end
  yline(mmax_R2,'Color','m','LineWidth',.66,'LineStyle','-')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  nexttile(7) % 7 ER (LINEAR Measure, some as ER0)
  hold on;
  if ~ADD_zGrid
    % plot(Plist/trnwin,ER ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
  end
    plot(Plist/trnwin,ER0, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9); 
  hold off;
  xlabel(xLbl,'interpreter','latex')
  if Xlim0;  xlim([0 XlimUp]); end
  box on; addgrid 
  setyticklabels([0:.01:.04], 2, FNs); % ylim([0 .035])
  subtitle('Expected Return','Interpreter','Latex')
  if Xlim0; line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', 1.5); end
  % ADD A VERTICAL LINE AT 
  if ~ADD_zGrid
  xL7 = xline(c(pER0),'Color',CLRS(zER0,:),'LabelOrientation','horizontal','Alpha',1, ...
   'Label',strcat('$\max(E(R_t)) =',num2str(mmax_ER,'%2.4f'),' ~[c = ', num2str(c(pER0),'%2.1f'),',', alpha_m(pER0,zER0).Properties.VariableNames,'$]'), ...
   'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2);  
  if trnwin == 12;  xL7.LabelHorizontalAlignment = 'left';  end
  if ismember(trnwin, [60, 120]);  xL7.LabelVerticalAlignment = 'top';  end
  else
  xL7 = xline(c(pER0),'Color',big_CLRS(zER0,:),'LabelOrientation','horizontal','Alpha',1, ...
   'Label',strcat('$\max(E(R_t)) =',num2str(mmax_ER,'%2.4f'),' ~[c = ', num2str(c(pER0),'%2.1f'),',', alpha_m(pER0,zER0).Properties.VariableNames,'$]'), ...
   'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','LineWidth',1,'Interpreter','latex','FontSize',FNs-2);  
    if trnwin == 12; xL7.LabelHorizontalAlignment = 'left'; end
  end
  yline(mmax_ER,'Color','m','LineWidth',.66,'LineStyle','-')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  nexttile(6) % 6 |beta|
  pla = plot(Plist/trnwin, Bnrmbar ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
  % pl0 = plot(Plist/trnwin, Bnrmbar ,  '--','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
  xlabel(xLbl,'interpreter','latex')
  if Xlim0;  xlim([0 XlimUp]); end
  box on; addgrid 
  setyticklabels([0:1:3], 0, FNs);
  subtitle('$\| \hat\beta \|$','Interpreter','Latex')
  % ----------------------------------------------------------------------------------------------
  if ~ADD_zGrid 
    leg_name = [cellstr(strcat(strcat('$z=',num2str((lamlist)')),'$'));'$c=1$'];  
    % leg_name = [cellstr(strcat(strcat('$\log_{10}(z)=',num2str(log10(lamlist)')),'$'));'$c=1$'];  
    if Xlim0
      % L0 = line(1*ones(1,2),[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', 1.5); %'HandleVisibility', 'off');
      L0 = line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', 1.5); 
      % leg_name = [cellstr(strcat(strcat('$\log_{10}(z)=',num2str(log10(lamlist)')),'$'));'$c=1$'];  
      % leg_name = [cellstr(strcat(strcat('$z=',num2str((lamlist)')),'$'));'$c=1$'];  
      lh = addlegend([pla; L0],leg_name,3,FNs-3);
    else 
      % leg_name = [cellstr(strcat(strcat('$\log_{10}(z)=',num2str(log10(lamlist)')),'$'))];  
      lh = addlegend([pla],leg_name,3,FNs-3,[],1);
      % lh = addlegend([pl0],leg_name,3,FNs-7,[],.5,[],4);
    end
    LLW = findobj(lh,'Type','Line');set(LLW,'LineWidth',2);
  end
  % ----------------------------------------------------------------------------------------------
  
  nexttile(8) % 8 Volatility
  hold on; % scale/annualize vola by 12, KMZ don't do this, so leave it at 1
  % if ~ADD_zGrid
  % yyaxis left 
  hVa =  plot(Plist/trnwin, Vol ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
    setyticklabels([0:1:6], 0, FNs); % ylim([0 5.5])
    if trnwin ==  12; setyticklabels([0:1:5], 0, FNs); ylim([0 5.5]); end
    hVa_clrs = get(hVa, 'Color');
  yyaxis right
  hV0 =  plot(Plist/trnwin, Vol0, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9); 

  ylim([-.1 .40]);
  if trnwin == 120 ; ylim([-.1 .30]); end
  if trnwin ==  60 ; ylim([-.1 .40]); end
    set(gca, 'YColor', 'k') % Set right y-axis color to black
    % setyticklabels([0:.1:3], 0, FNs); % ylim([0 .5])
    for LL = 1:length(hVa_clrs)
      hV0(LL).Color = hVa_clrs{LL}*DK;
    end
    hline(0)
  hold off;
  xlabel(xLbl,'interpreter','latex')
  if Xlim0;  xlim([0 XlimUp]); end
  box on; addgrid 
  % setyticklabels([0:1:5], 0, FNs); % ylim([0 .5])
  subtitle('Volatility','Interpreter','Latex')
  if Xlim0; line([1,1],[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', 1.5); end

  % pdf_name = strcat(para_str, '_eval_plots');
  pdf_name = strcat('Tis=', num2str(trnwin), '_eval_all_plots_', num2str(subbeg), '-', num2str(subend) );
  if Xlim0;     pdf_name = strcat(pdf_name , '_short_Xaxis'); end
  if ADD_zGrid; pdf_name = strcat(pdf_name , '_zGrid'); end
  % if PRNT_PDF; plot2pdf(pdf_name); end
  % plot2pdf(strcat(pdf_name,'_2nd'))
end
end
% copyfile('*.pdf', 'D:\_research\_current\complexity\graphics');
%%
% --------------------------------------------------------------------------------------------------
% MAKE SEPERATE PLOTS FOR SR OR IR ONLY
% --------------------------------------------------------------------------------------------------
if PLOT_INDIVIDUAL
for TYPE_ = [0] % if 0 ---> db average first
FNs = 16; 
F2 = figure(2); clf; % tiledlayout(5,2,'TileSpacing', 'compact', 'Padding', 'loose');
set(F2,'Name',['Sharpe Ratio (SR). Tis = ' num2str(trnwin) '. OOS: ' dates_strng], 'WindowState','maximized','Position',[1268 1  1268  2178]);
  hold on;
  % if ~ADD_zGrid 
  %   pl0 = plot(Plist/trnwin,sqrt(12)*SR ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9);  
  % end
  % if TYPE_ == 0
  % pla = plot(Plist/trnwin,sqrt(12)*SR ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9);  
  % pl0 = plot(Plist/trnwin,sqrt(12)*SR0, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9);  
  pla = plot(Plist/trnwin,R2 ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9);  
  pl0 = plot(Plist/trnwin,R20, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9);  
  % else
  %   pla = plot(Plist/trnwin,sqrt(12)*SR ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9);  
  % end

  hVa_clrs = get(pl0, 'Color');
  for LL = 1:length(hVa_clrs)
    pl0(LL).Color = hVa_clrs{LL}*DK;
  end
  % end
  hold off;
  setplotdims([.10 2/3 .82 .17], 3/4); 
  xlabel(xLbl,'interpreter','latex')
  leg_name = [cellstr(strcat(strcat('$z=',num2str((lamlist)')),'$'));'$c=1$'];  
  leg_name = strrep(leg_name,'0.6000000000000001','0.6');
  box on; addgrid
  % setyticklabels(0:.1:.6, 1, FNs); % ylim([0 .52]); 
  setyticklabels(-1:.2:.2, 1, FNs); % ylim([0 .52]); 
  tickshrink(.5)
  % set(gca,'TickLength',[0.02 0.02])
  
  % 0 suffix --> db or correct way, average first then Evaluate.
  % ADD A VERTICAL LINE AT MAX (SR) 
  if TYPE_ == 0
  xLa = xline( c(pSR0), 'Color', big_CLRS(zSR0,:)*DK,'LabelOrientation','horizontal','Alpha',1, ...% 'HandleVisibility', 'off', ...
    'Label', strcat('max(SR) $=', num2str(maxSR0,4),'~[c = ', num2str(c(pSR0)),',',alpha_m(pSR0,zSR0).Properties.VariableNames,'$]'), ...
    'LabelHorizontalAlignment', 'right','LineWidth',1.5,'Interpreter','latex','FontSize',FNs-2);
  else
  xLa = xline( c(pSR), 'Color', big_CLRS(zSR,:),'LabelOrientation','horizontal','Alpha',1, ...% 'HandleVisibility', 'off', ...
    'Label', strcat('max(SR) $=', num2str(maxSR,4),'~[c = ', num2str(c(pSR)),',',alpha_m(pSR,zSR).Properties.VariableNames,'$]'), ...
    'LabelHorizontalAlignment', 'right','LineWidth',1.5,'Interpreter','latex','FontSize',FNs-2);
  if trnwin == 60  ; xLa.LabelHorizontalAlignment = 'left';   end
  end

  if trnwin == 120 ; xLa.LabelHorizontalAlignment = 'left' ;  end
  if trnwin == 12  ; xLa.LabelHorizontalAlignment = 'right';  end
  % xLa.LabelHorizontalAlignment = 'right';
  yline(mmaxSR,'Color','m','LineWidth',.66,'LineStyle','-')
  
  % subtitle('SR'); % title('SR','FontWeight','normal')
  L0 = line(1*ones(1,2)      ,[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', 1.5); %'HandleVisibility', 'off');
  
  % if ADD_zGrid; xlim([0 1000.5]); end
  % adjust ticklenght
  tcklength = 0.015;
  % set(BRK.rightAxes,'TickLength',tcklength/2*ones(1,2))  
  %   if trnwin ==  12; breakxaxis([71 250], .02); end 
  %   if trnwin ==  60; breakxaxis([9  190], .02); end
  %   if trnwin == 120; breakxaxis([7   90], .02); end
  % else
  
  % if ~ADD_zGrid
    if trnwin ==  12 
      BRK = breakxaxis([10 990], .02); 
      % BRK = breakxaxis([151 480], .02); 
      % BRK = breakxaxis([201 430], .02); 
      set(BRK.leftAxes ,'TickLength',tcklength*ones(1,2))
      set(BRK.rightAxes,'TickLength',tcklength/2*ones(1,2))  
    end
  % end

  if trnwin ==  60 
    BRK = breakxaxis([11  190], .02); 
    set(BRK.leftAxes ,'TickLength',tcklength*ones(1,2))
    set(BRK.rightAxes,'TickLength',tcklength/1*ones(1,2))  
  end

    
  if trnwin == 120
    BRK = breakxaxis([7   90], .02); 
    set(BRK.leftAxes ,'TickLength',tcklength*ones(1,2))
    set(BRK.rightAxes,'TickLength',tcklength/1.5*ones(1,2))  
  end

  % add legend
  if trnwin == 60 
    % if ~ADD_zGrid; lh = addlegend([pl0; L0],leg_name, 5, FNs-3); end
    if ~ADD_zGrid; lh = addlegend([pla],leg_name, 5, FNs-3); end
  end
  % addlegend([pl0; L0], leg_name, 5, 10)

  % print to pdf file
  % PRNT_PDF = 1
  % pdf_name_ind_SR = strcat('Tis=', num2str(trnwin), '_eval_SR_plots_', num2str(subbeg), '-', num2str(subend), '_Type_', num2str(TYPE_) );
  pdf_name_ind_SR = strcat('Tis=', num2str(trnwin), '_eval_SR_plots_', num2str(subbeg), '-', num2str(subend), '_Type_Combined'  );
  if ADD_zGrid; pdf_name_ind_SR = strcat(pdf_name_ind_SR,'_zGrid'); end
  if PRNT_PDF;  plot2pdf(pdf_name_ind_SR); end

% %% not needed for now, 
% F3 = figure(3); clf; % tiledlayout(5,2,'TileSpacing', 'compact', 'Padding', 'loose');
% set(F3, 'Name',['Information Ratio (IR). OOS: ' dates_strng], 'WindowState','maximized','Position',[1268 1  1268  2178]);
%  % IR
%   hold on;
%   if ADD_zGrid
%     pl0 = plot(Plist/trnwin,sqrt(12)*IR ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
%   else
%     pl0 = plot(Plist/trnwin,sqrt(12)*IR ,  '- ','linewidth',LW , 'Marker', MK , 'Markersize', 9); 
%     plot(Plist/trnwin,sqrt(12)*IR0, LType,'linewidth',LW0, 'Marker', MK0, 'Markersize', 9); 
%   end
%   hold off;
%   xlabel(xLbl,'interpreter','latex')
%   leg_name = [cellstr(strcat(strcat('$\log_{10}(z)=',num2str(log10(lamlist)')),'$'))];  
%   box on; addgrid
%   setyticklabels(0:.1:.6, 1, FNs); 
%   L0 = line(1*ones(1,2)      ,[-1e9,1e9],'LineStyle','-','color', 'k','Linewidth', 1.5);%  'HandleVisibility', 'off');
%   setplotdims([.10 2/3 .82 .22])
%   % ADD A VERTICAL LINE AT MAX (IR) % uncmment the visbility bit to remove the double entry of
%   xL = xline(c(pIR0(1)), 'Color', big_CLRS(zIR0(1),:), ... %'HandleVisibility', 'off', ...
%     'Label', strcat('max(IR) $=', num2str(mmaxIR,4),', ~ c = ', num2str(c(pIR0(1))), '$.'), 'LabelOrientation','horizontal','Alpha',1, ...
%     'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','LineWidth',1.5,'Interpreter','latex','FontSize',FNs-2);
%   yline(mmaxIR,'Color','m','LineWidth',.66,'LineStyle','-')
%   if trnwin == 12; breakxaxis([91  800], .02);  end
%   if trnwin == 60; breakxaxis([9   190], .02);  end
%   constantLines  = findall(gca, 'Type', 'ConstantLine');
%   pdf_name_ind_IR = strcat('Tis=', num2str(trnwin), '_eval_IR_plots_', num2str(subbeg), '-', num2str(subend) );
%   if ADD_zGrid; pdf_name_ind_IR = strcat(pdf_name_ind_IR,'_zGrid'); end
%   if PRNT_PDF; plot2pdf(pdf_name_ind_IR); end
end

if COPY_FILE 
  % copy all *.pdf files to D:\_research\_current\complexity\graphics
  copyfile('*.pdf', 'D:\_research\_current\complexity\graphics');
  % diary off; 
  copyfile('*.txt', 'D:\_research\_current\complexity\graphics');
end


%% make my own Table 1 tables 
% clc
%**************************************************************************
% Generate Table 1: Comparison with Welch and Goyal (2008) and mkt
%**************************************************************************
out = [];
% sub-sample set subsample 
SubSmpl = find(~isnan(Y'+Yprd_gy(:,1))); 
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% GW(z=0)
% --------------------------------------------------------------------------------------------------
gy_z0      = timing_gy(:,1);
Yprd_gy_z0 =   Yprd_gy(:,1);
r2tmp_z0   = 1-nanvar(Y(SubSmpl)'-Yprd_gy_z0(SubSmpl))/nanvar(Y(SubSmpl)');
stats_z0   = regstats(gy_z0,ones(size(Y')),1,'tstat');
% GW(z=0 )vs mkt
stats_z0_vsMKT   = regstats(gy_z0,Y','linear',{'tstat','rsquare','r'}); 
tmp1res_z0_vsMKT = stats_z0_vsMKT.tstat.beta(1) + stats_z0_vsMKT.r;
out_z0 = [out ; 100*r2tmp_z0 sqrt(12)*sharpe(gy_z0,0) stats_z0.tstat.t(1) sqrt(12)*sharpe(tmp1res_z0_vsMKT,0) stats_z0_vsMKT.tstat.t(1) nan(1,2) -min(gy_z0) skewness(gy_z0)];
% --------------------------------------------------------------------------------------------------
% GW z = 1000 
% --------------------------------------------------------------------------------------------------
gy_z1000      = timing_gy(:,end);  
Yprd_gy_z1000 =   Yprd_gy(:,end);
r2tmp_z1000   = 1-nanvar(Y(SubSmpl)'-Yprd_gy_z1000(SubSmpl))/nanvar(Y(SubSmpl)');
stats_z1000   = regstats(gy_z1000,ones(size(Y')),1,'tstat');
% GW(z=1000) vs mkt
stats_z1000_vsMKT   = regstats(gy_z1000,Y','linear',{'tstat','rsquare','r'}); 
tmp1res_z1000_vsMKT = stats_z1000_vsMKT.tstat.beta(1) + stats_z1000_vsMKT.r;
% out     = [out ; 100*r2tmp sqrt(12)*sharpe(tmp1,0) stats.tstat.t(1) nan(1,4) min(tmp1) skewness(tmp1)];
out_z1000 = [out_z0 ; 100*r2tmp_z1000 sqrt(12)*sharpe(gy_z1000,0) stats_z1000.tstat.t(1) sqrt(12)*sharpe(tmp1res_z1000_vsMKT,0) stats_z1000_vsMKT.tstat.t(1) nan(1,2) -min(gy_z1000) skewness(gy_z1000)];

% --------------------------------------------------------------------------------------------------
% RFF(z = 1000) High complexity
% --------------------------------------------------------------------------------------------------
RFF_timing = nanmean(timing(:,end,end,:),4);  % R_pi,t  --> average prediction over 1000 RFFS
RFF_Yprd   = nanmean(  Yprd(:,end,end,:),4);  % pi,t    --> average prediction over 1000 RFFS
r2_RFF     = 1 - nanvar( Y(SubSmpl)'-RFF_Yprd(SubSmpl) ) / nanvar(Y(SubSmpl)');
stats_RFF  = regstats(RFF_timing, ones(size(Y')),1,'tstat');
% RFF(z = 1000) High complexity vs mkt IR
stats_RFF_vsMKT = regstats(RFF_timing,Y','linear',{'tstat','rsquare','r'}); 
res_RFF_vsMKT   = stats_RFF_vsMKT.tstat.beta(1) + stats_RFF_vsMKT.r;
% RFF(z = 1000) High complexity vs GW(z=1000)
stats_RFF_vsGWz1000 = regstats(RFF_timing,timing_gy(:,end),'linear',{'tstat','rsquare','r'}); 
res_RFF_vsGWz1000   = stats_RFF_vsGWz1000.tstat.beta(1) + stats_RFF_vsGWz1000.r;
out_RFF   = [out_z1000 ; 100*r2_RFF sqrt(12)*sharpe(RFF_timing,0) stats_RFF.tstat.t(1) ...
                         sqrt(12)*sharpe(res_RFF_vsMKT,0)         stats_RFF_vsMKT.tstat.t(1) ...
                         sqrt(12)*sharpe(res_RFF_vsGWz1000,0)     stats_RFF_vsGWz1000.tstat.t(1) ...
                         -min(RFF_timing) skewness(RFF_timing)];

sep;disp(['Tis = ' num2str(trnwin) ': Their mixed approach']);sep

% disp([sqrt(12)*sharpe(tmp1res,0) stats.tstat.t(1)])
% Now Add what is plotted from: parfor i=1:nSim for RFF(z=1000)
out_RFF_plots = [100*R2(end,end) sqrt(12)*SR(end,end) SR_tstat(end,end), ...
                  sqrt(12)*IR(end,end)                IRt(end,end)   , ...
                  sqrt(12)*IR_vs_GWz1000(end,end)     IR_tstat_vs_GWz1000(end,end)   , ...
                  maxLoss(end,end) Skew(end,end), ...
                  ];

out_RFF_all   = [out_RFF; out_RFF_plots ];
lst(out_RFF_all, 4)
sep
% print to xls if needed
suffix = ['trnwin-' num2str(trnwin) '-gamma-' gamma_str '-stdize-' num2str(stdize) '-demean-' num2str(demean) '-' num2str(subbeg(ss)) '-' num2str(subend(ss))];
writematrix(out_RFF_all,[suffix '-out.xlsx'],'Range','B2');

% %% compare to pure loop and pure average approach
% SRt_stat    = regstats(mean_Ret_pi(:,end,end), ones(size(mean_Ret_pi(:,end,end))) ,1,'tstat')
% Table1_rff  =  [R2(end,end)  sqrt(12)*SR(end,end)  SRt_stat.tstat.t(1) sqrt(12)*IR(end,end) IRt(end,end) alpha(end,end)] 
% 
% Table1_rff0 = [100*R20(end,end), ... 
%               sqrt(12)*SR0(end,end) 
%               sqrt(12)*IR0(end,end) 
%               IRt0(end,end) 
%               alpha0(end,end)
%               ] 
diary off




end 

end




























% end % functon end

