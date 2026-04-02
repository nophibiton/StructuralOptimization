% Performs reliability-based design optimization using decoupled approach.
% The example uses Sequential optimization and reliability assessment
% (SORA) described in [1]. The example solves Example 7.2 Short column in [2].
%
% Reference:
%
% [1] Du, X., & Chen, W. (2004). Sequential Optimization and Reliability 
% Assessment Method for Efficient Probabilistic Design. Journal of 
% Mechanical Design, 126(2), 225–233. https://doi.org/10.1115/1.1649968
%
% [2] Aoues, Y., & Chateauneuf, A. (2009). Benchmark study of numerical 
% methods for reliability-based design optimization. Structural and 
% Multidisciplinary Optimization, 41(2). https://doi.org/10.1007/S00158-009-0412-2
%
clc,clear; format compact; format shortG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Short column Normal case %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i_cov = [0.05,0.10,0.15]

% define objective function parameters
objparams   = struct;
calc_obj    = @objfun_column_lsf;

% define constraint function parameters
targetbeta = 3.0;
lsffunc    = @short_column_lsf;
lsfparams.ind_d_det = [];
lsfparams.ind_d_ran = 1:2; % index of the random design variablve
ind_d_ran_x         = 5:6; % index of the random variable for the design

% define optimization options
options = optimoptions('fmincon','display','none','algorithm','sqp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hist            = struct;
hist.dvstar     = [];
hist.fobj_val   = [];
max_iter        = 20;
eps1            = 10^-4;

fobj_val_old = inf;

dvstar_prev = [2,2];

    mu_h = dvstar_prev(1);
    mu_b = dvstar_prev(2);
    cov_bh = i_cov;
    
    % define base units
    kN  = 1;
    m   = 1;
    kPa = kN/m^2;
    MPa = 1000*kPa;

    % define probability data
    probdata       = struct;
    probdata.Xdist = { 'normal','normal','normal','normal','normal','normal'};
    probdata.Xmu   = [ 2500*kN, 250*kN*m,125*kN*m, 40*MPa, mu_h, mu_b];
    probdata.Xstd  = [ 0.2*2500*kN, 0.3*250*kN*m, 0.3*125*kN*m, 0.1*40*MPa, cov_bh*mu_h,cov_bh*mu_b];
    probdata.R     = eye(6);
    % get distribution parameters
    distparams          = get_dist_params(probdata);
    probdata.distparams = distparams;

Ro   = modify_corr(probdata.R,probdata);
Lo   = (chol(Ro))';
iLo  = inv(Lo);

ximpp_prev  = probdata.Xmu';
k           = 1;

while true
    % (1) Deterministic optimization
    % perform deterministic optimization
    fun = @(d) calc_obj(d,objparams); % objective function

    % modify probabilistic constraint into a deterministic
    if ~isempty(lsfparams.ind_d_ran)
        lsfparams.si = dvstar_prev(lsfparams.ind_d_ran) - ximpp_prev(ind_d_ran_x)';
    end
    lsfparams.isDeterministic = true;
    det_g_fun = @(d)eval_lsf(ximpp_prev,d,false, lsffunc, lsfparams,distparams);
    nonlcon   = @(d)deal([-1*det_g_fun(d);0.5*d(1)-d(2);d(2)-2*d(1)],[]); % constraint function
    
    % perform optimization
    d0                = dvstar_prev;
    [dvstar,fobj_val] = fmincon(fun,d0,[],[],[],[],[0,0],[15,15],nonlcon,options);

    % (2) Reliability assessment
    % perform inverse reliability analysis
    mu_h           = dvstar(1);
    mu_b           = dvstar(2);
    probdata.Xmu   = [ 2500*kN, 250*kN*m,125*kN*m, 40*MPa, mu_h, mu_b];
    probdata.Xstd  = [ 0.2*2500*kN, 0.3*250*kN*m, 0.3*125*kN*m, 0.1*40*MPa, cov_bh*mu_h,cov_bh*mu_b];
    distparams          = get_dist_params(probdata);
    probdata.distparams = distparams;
    Ro   = modify_corr(probdata.R,probdata);
    Lo   = (chol(Ro))';
    iLo  = inv(Lo);

    lsfparams.isDeterministic = false;
    ximpp = invform_hmv(ximpp_prev,dvstar,targetbeta,...
                           probdata,Lo,iLo,lsffunc,lsfparams,distparams);

    hist.fobj_val = [hist.fobj_val;fobj_val];
    hist.dvstar   = [hist.dvstar; dvstar];

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


fprintf('\n');
disp(['Short column results, Normal case (COV=' num2str(i_cov) '):'])
disp(['Optimum: ' num2str(dvstar)]);
disp(['fstar: ' num2str(fobj_val)]);


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Short column Non-normal case %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i_cov = [0.05,0.10,0.15]

% define objective function parameters
objparams   = struct;
calc_obj    = @objfun_column_lsf;

% define constraint function parameters
targetbeta = 3.0;
lsffunc    = @short_column_lsf;
lsfparams.ind_d_det = [];
lsfparams.ind_d_ran = 1:2; % index of the random design variablve
ind_d_ran_x         = 5:6;  % index of the random variable for the design

% define optimization options
options = optimoptions('fmincon','display','none','algorithm','sqp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hist            = struct;
hist.fobj_val   = [];
max_iter        = 20;
eps1            = 10^-4;

fobj_val_old = inf;

dvstar_prev = [2,2];

    mu_h = dvstar_prev(1);
    mu_b = dvstar_prev(2);
    cov_bh = i_cov;
    
    % define base units
    kN  = 1;
    m   = 1;
    kPa = kN/m^2;
    MPa = 1000*kPa;

    % define probability data
    probdata       = struct;
    probdata.Xdist = { 'Gumbel','Gumbel','Gumbel','weibull','Lognormal','Lognormal'};
    probdata.Xmu   = [ 2500*kN, 250*kN*m,125*kN*m, 40*MPa, mu_h, mu_b];
    probdata.Xstd  = [ 0.2*2500*kN, 0.3*250*kN*m, 0.3*125*kN*m, 0.1*40*MPa, cov_bh*mu_h,cov_bh*mu_b];
    probdata.R     = eye(6);
    % get distribution parameters
    distparams          = get_dist_params(probdata);
    probdata.distparams = distparams;
    Ro   = modify_corr(probdata.R,probdata);
    Lo   = (chol(Ro))';
    iLo  = inv(Lo);


ximpp_prev  = probdata.Xmu';
k           = 1;

while true
    % (1) Deterministic optimization
    % perform deterministic optimization
    fun = @(d) calc_obj(d,objparams); % objective function

    % modify probabilistic constraint into a deterministic
    if ~isempty(lsfparams.ind_d_ran)
        lsfparams.si = dvstar_prev(lsfparams.ind_d_ran) - ximpp_prev(ind_d_ran_x)';
    end
    lsfparams.isDeterministic = true;
    det_g_fun = @(d)eval_lsf(ximpp_prev,d,false, lsffunc, lsfparams,distparams);
    nonlcon   = @(d)deal([-1*det_g_fun(d);0.5*d(1)-d(2);d(2)-2*d(1)],[]); % constraint function
    
    % perform optimization
    d0                = dvstar_prev;
    [dvstar,fobj_val] = fmincon(fun,d0,[],[],[],[],[0,0],[15,15],nonlcon,options);

    % (2) Reliability assessment
    % perform inverse reliability analysis
    mu_h = dvstar(1);
    mu_b = dvstar(2);
    probdata.Xmu   = [ 2500*kN, 250*kN*m,125*kN*m, 40*MPa, mu_h, mu_b];
    probdata.Xstd  = [ 0.2*2500*kN, 0.3*250*kN*m, 0.3*125*kN*m, 0.1*40*MPa, cov_bh*mu_h,cov_bh*mu_b];
    distparams          = get_dist_params(probdata);
    probdata.distparams = distparams;
    Ro   = modify_corr(probdata.R,probdata);
    Lo   = (chol(Ro))';
    iLo  = inv(Lo);

    lsfparams.isDeterministic = false;
    ximpp = invform_hmv(ximpp_prev,dvstar,targetbeta,...
                           probdata,Lo,iLo,lsffunc,lsfparams,distparams);

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

fprintf('\n');
disp(['Short column results, Nonnormal case (COV=' num2str(i_cov) '):'])
disp(['Optimum: ' num2str(dvstar)]);
disp(['fstar: ' num2str(fobj_val)]);

end
