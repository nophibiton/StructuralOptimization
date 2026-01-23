function [modelparams] = gp_model(X, y)


% create Kriging model

% define parameters
meanfunc = @meanZero;
covfunc  = @covSEiso; 
likfunc  = @likGauss;

% initialize
hyp.cov  = [log(0.5); log(1.0)];
hyp.mean = [];
hyp.lik  = log(1e-6);


hyp = minimize(hyp, @gp, -50, @infGaussLik, meanfunc, ...
                            covfunc, likfunc, X, y);

% define parameters of the surrogate model
modelparams.hyp       = hyp;
modelparams.meanfunc  = meanfunc;
modelparams.covfunc   = covfunc;
modelparams.likfunc   = likfunc;

% [mu, s2] = gp(hyp, @infGaussLik, ...
%                       meanfunc, covfunc, likfunc, X, y, X);


end