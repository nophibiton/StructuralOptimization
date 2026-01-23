function [X_new,EI_x] = EI(surrogateModel,X,y)

% unpack parameters of the surrogate model
hyp      = surrogateModel.hyp;
meanfunc = surrogateModel.meanfunc;
covfunc  = surrogateModel.covfunc;
likfunc  = surrogateModel.likfunc;

% sampling of input values/other sampling techniques can be used
X_samples = (0:0.001:15)';

[mu, s2] = gp(hyp, @infGaussLik, ...
                      meanfunc, covfunc, likfunc, X, y, X_samples);

% calculate max expected improvement
f_star = min(y);         % best so far
sigma  = sqrt(s2);
Z      = (f_star-mu)./sigma;
EI_x   = (f_star-mu).*normcdf(Z) + sigma.*normpdf(Z);

% return point/X_new with max expected improvement
[~, idx] = max(EI_x);
X_new    = X_samples(idx);

end