% ************************************************************************************************
if ~make_performance_measures
  % % OLD if file exists: if isfile(performance_filename)
  disp('Loading data KMZ portfolio performance file, ... ')
  tic; load(performance_filename); toc;
else
  disp('--------------------------------------------------------------------------------------\n');
  disp('Computing portfolio performance measures from separate simulations, ... ')
  % Performance initialization
  Bias        = nan(nP,nL);
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
      for (i=1:nSim) % regstats(y, X, model)
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
% ************************************************************************************************
