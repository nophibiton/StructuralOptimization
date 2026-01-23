function [K, dK] = covHamming(hyp, x, z, ii)
% covHamming - GPML-compatible Hamming-distance kernel
%
% k(x,z) = sf2 * exp(-hdist(x,z)/ell)
% hyp = [ log(ell)
%         log(sf) ]

if nargin < 2
    K = '2'; return;  % number of hyperparameters
end

ell = exp(hyp(1));
sf2 = exp(2*hyp(2));

if nargin < 3 || isempty(z)
    z = x;
end
xeqz = isempty(z);
dg = ischar(z) && strcmp(z, 'diag');

% ---- Diagonal mode ----
if dg
    K = sf2 * ones(size(x,1),1);
    if nargout > 1, dK = zeros(size(K)); end
    return;
end

% ---- Compute pairwise Hamming distances safely ----
n = size(x,1);
if xeqz
    hdist = zeros(n);
    for a = 1:n
        xa = x(a,:);
        for b = a:n
            d = sum(xa ~= x(b,:));
            hdist(a,b) = d;
            hdist(b,a) = d;
        end
    end
else
    m = size(z,1);
    hdist = zeros(n,m);
    for a = 1:n
        xa = x(a,:);
        for b = 1:m
            hdist(a,b) = sum(xa ~= z(b,:));
        end
    end
end

% ---- Covariance matrix ----
K = sf2 * exp(-hdist / ell);

% ---- Derivative mode ----
if nargin == 4 && ~isempty(ii)
    switch ii
        case 1  % ∂K/∂log(ell)
            dK = K .* (hdist / ell);
        case 2  % ∂K/∂log(sf)
            dK = 2 * K;
        otherwise
            error('covHamming: invalid hyperparameter index');
    end
else
    dK = [];
end
end