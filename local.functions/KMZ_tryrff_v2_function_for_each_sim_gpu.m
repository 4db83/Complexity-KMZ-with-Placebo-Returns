function [] = KMZ_tryrff_v2_function_for_each_sim_gpu(gamma, trnwin, iSim, stdize, demean, OUT_DIR_NAME)

%**************************************************************************
% Setup
%**************************************************************************
nSim = 1;
maxP = 12000;
trainfrq = 1;
% OUT_DIR_NAME = './temp1/';
Plist = [2:2:22 24:12:(9*12) 12*(10:10:100) 1500 2e3:2e3:12e3];
log_lamlist = -3:1:3;
lamlist = 10.^(log_lamlist);
nL = length(lamlist);
nP = length(Plist);
para_str = strcat('maxP-', num2str(maxP), '-trnwin-', num2str(trnwin), ...
    '-gamma-', num2str(gamma), '-stdize-', num2str(stdize), '-demean-', num2str(demean), '-v2');
save_path = strcat(OUT_DIR_NAME, para_str);
if ~exist(save_path, 'dir'); mkdir(save_path); end

%**************************************************************************
% Load Data
%**************************************************************************
load KMZ_GYdata.mat

X = [X lagmatrix(Y,1)];
X = volstdbwd(X,[]); % still on CPU
Y2 = 0;
for j = 1:12
    Y2 = Y2 + lagmatrix(Y.^2, j);
end
Y2 = Y2 / 12;

if stdize == 1    
    Y = Y ./ sqrt(Y2);
end

Y = Y(37:end);
X = X(37:end,:);
dates = dates(37:end,:);
T = length(Y);
X = X';
Y = Y';
d = size(X,1);

%**************************************************************************
% GPU Conversion
%**************************************************************************
X = gpuArray(X);
Y = gpuArray(Y);

%**************************************************************************
% Initialize Outputs
%**************************************************************************
Yprd = nan(T, nP, nL, nSim, 'gpuArray');
Bnrm = nan(T, nP, nL, nSim, 'gpuArray');

s = iSim;
disp(s);
rng(s);
W = gpuArray(randn(max(Plist), d));

for p = 1:nP
    P = floor(Plist(p)/2);
    wtmp = W(1:P, :);
    Z = [cos(gamma * wtmp * X); sin(gamma * wtmp * X)];
    Yprdtmp = nan(T, nL, 'gpuArray');
    Bnrmtmp = nan(T, nL, 'gpuArray');

    for t = trnwin+1:T
        trnloc = (t - trnwin):t-1;
        Ztrn = Z(:, trnloc);
        Ytrn = Y(trnloc);
        Ztst = Z(:, t);

        if demean == 1
            Ymn = mean(Ytrn, 'omitnan');
            Zmn = mean(Ztrn, 2, 'omitnan');
        else
            Ymn = 0;
            Zmn = 0;
        end

        Ytrn = Ytrn - Ymn;
        Ztrn = Ztrn - Zmn;
        Ztst = Ztst - Zmn;
        Zstd = std(Ztrn, 0, 2, 'omitnan');
        Ztrn = Ztrn ./ Zstd;
        Ztst = Ztst ./ Zstd;

        if t == trnwin+1 || mod(t-trnwin-1, trainfrq) == 0
            if P <= trnwin
                B = ridgesvd(gather(Ytrn'), gather(Ztrn'), lamlist * trnwin); % gather to CPU for SVD
            else
                B = get_beta(gather(Ytrn'), gather(Ztrn'), lamlist); % gather to CPU
            end
            B = gpuArray(B); % push back to GPU
        end

        Yprdtmp(t, :) = B' * Ztst + Ymn;
        Bnrmtmp(t, :) = sum(B.^2);
    end

    Yprd(:, p, :, 1) = Yprdtmp;
    Bnrm(:, p, :, 1) = Bnrmtmp;
end

%**************************************************************************
% Save
%**************************************************************************
Yprd = gather(Yprd);
Bnrm = gather(Bnrm);
if iSim == 1
    save([save_path '/iSim' num2str(iSim) '.mat']);
else 
    save([save_path '/iSim' num2str(iSim) '.mat'], 'Yprd', 'Bnrm');
end

end
