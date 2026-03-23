% Performs reliability-based design optimization using decoupled approach
% The example uses Sequential optimization and reliability assessment
% (SORA) described in [1] and solves Example 11.3.5 in [2].
%
% Reference:
%
% [1] Du, X., & Chen, W. (2004). Sequential Optimization and Reliability 
% Assessment Method for Efficient Probabilistic Design. Journal of 
% Mechanical Design, 126(2), 225–233. https://doi.org/10.1115/1.1649968
%
% [2] Melchers, R. E., & Beck, A. T. (2018). Structural Reliability 
% Analysis and Prediction (2nd ed.). John Wiley & Sons.
%
clc,clear; format compact; format shortG;

options = optimoptions('fmincon','display','none','algorithm','sqp');

targetbeta = 4.0;

objparams   = struct;

lsffunc   = @g_fun;

lsfparams = struct;

% define probability data
probdata       = struct;
probdata.Xdist = {    'normal',     'normal',  'normal'};
probdata.Xmu   = [          40,           50,      1000];
probdata.Xstd  = [           5,           2.5,      200];
probdata.R     = eye(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get distribution parameters
distparams          = get_dist_params(probdata);
probdata.distparams = distparams;

Ro   = modify_corr(probdata.R,probdata);
Lo   = (chol(Ro))';
iLo  = inv(Lo);

hist            = struct;
hist.fobj_val   = [];
max_iter        = 20;
eps1            = 10^-4;

fobj_val_old = inf;

k     = 1;
ximpp_prev  = probdata.Xmu';
si          = zeros(size(ximpp_prev));
dvstar_prev = 1;
while true
    % (1) Deterministic optimization
    % perform deterministic optimization
    fun = @(d) calc_obj(d,objparams); % objective function

    % modify probabilistic constraint into a deterministic
    det_g_fun = @(d)eval_lsf(ximpp_prev, d, false, lsffunc, lsfparams,distparams);
    nonlcon   = @(d)deal(-1*det_g_fun(d),[]); % constraint function
    
    % perform optimization
    d0            = dvstar_prev;
    [dvstar,fobj_val] = fmincon(fun,d0,[],[],[],[],0,5,nonlcon,options);

    % (2) Reliability optimization
    % perform inverse reliability analysis
    ximpp = invform_hmv(ximpp_prev,dvstar,targetbeta,probdata,Lo,iLo,lsffunc,lsfparams,distparams);

    % check convergence
    if abs(fobj_val_old - fobj_val) < eps1
        exit_message = 'Objective function value converges.';
        break;
    end
    if k == max_iter
        exit_message = 'Maximum iteration is observed.';
        break;
    end

    % store previous iteration results
    ximpp_prev   = ximpp;
    dvstar_prev  = dvstar;
    fobj_val_old = fobj_val;

    k = k+1;
end
dvstar
disp(exit_message)



