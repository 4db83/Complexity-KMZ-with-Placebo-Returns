function B = ridgesvd_gpu(Y, X, lambda)
% GPU-compatible version of ridge regression via SVD
% Inputs:
%   Y       Tx1 (gpuArray or CPU)
%   X       TxP (gpuArray or CPU)
%   lambda  Lx1 (CPU)
% Output:
%   B       PxL (gpuArray)

% Check for NaNs
if any(isnan(X(:))) || any(isnan(Y))
    error('missing data')
end

% Convert to GPU if not already
if ~isa(X, 'gpuArray')
    X = gpuArray(X);
end
if ~isa(Y, 'gpuArray')
    Y = gpuArray(Y);
end

[T, P] = size(X);
L = length(lambda);

% SVD on GPU
[U, S, V] = svd(X, 'econ');  % econ is faster for T < P
D = diag(S);                 % Singular values (on GPU)

% Handle complement zeros
if T >= P
    compl = gpuArray.zeros(P, T - P);
else
    compl = gpuArray.zeros(P - T, T);
end

B = gpuArray.nan(P, L);

% Ridge solution for each lambda
for l = 1:L
    if T >= P
        B(:, l) = V * [diag(D ./ (D.^2 + lambda(l))), compl] * (U' * Y);
    else
        B(:, l) = V * [diag(D ./ (D.^2 + lambda(l))); compl] * (U' * Y);
    end
end

end
