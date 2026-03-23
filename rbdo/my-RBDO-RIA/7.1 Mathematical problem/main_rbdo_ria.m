% Performs reliability-based design optimization using double-loop
% approach. The reliability index approach (RIA) with iHLRF algorithm
% in [1] for the reliability analysis is applied. The example
% solves Example 7.1 Mathematical problem in [2].
%
% Reference:
%
% [1] Der Kiureghian, A. (2005). First- and Second-Order Reliability 
% Methods. In Engineering Design Reliability Handbook. CRC press.
%
% [2] Aoues, Y., & Chateauneuf, A. (2009). Benchmark study of numerical 
% methods for reliability-based design optimization. Structural and 
% Multidisciplinary Optimization, 41(2). https://doi.org/10.1007/S00158-009-0412-2
%
clc,clear; format compact; format shortG;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear limit state %
%%%%%%%%%%%%%%%%%%%%%%%%%

% define probability data
probdata       = struct;
probdata.Xdist = {    'normal',     'normal'};
probdata.Xmu   = [         5.0,          3.0];
probdata.Xstd  = [     0.3*5.0,      0.3*3.0];
probdata.R     = eye(2);

% define objective function parameters
objparams   = struct;
calc_obj    = @objfun_nonlinear_lsf;

% define constraint function parameters
constparams = struct;
constparams.targetbeta = -norminv(0.01);
constparams.lsffunc    = @nonlinear_lsf;
constparams.probdata   = probdata;

fun     = @(x) calc_obj(x,objparams);     % objective function
nonlcon = @(x) calc_const(x,constparams); % constraint function

A   = [];   b = [];    % coefficients of inequality constraints
Aeq = []; beq = [];    % coefficients of equality constraints

x0  = [2,2];                % intial guess
lb  = [0,0]; ub  = [15,15]; % upper and lower bounds

options = optimoptions('fmincon','display','none','algorithm','sqp');
[x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

disp('Nonlinear limit state results: ')
disp(['Optimum: ' num2str(x)]);
disp(['fstar: ' num2str(fval)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Highly nonlinear limit state %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define probability data
probdata       = struct;
probdata.Xdist = {    'normal',     'normal'};
probdata.Xmu   = [         5.0,          3.0];
probdata.Xstd  = [     0.3*5.0,      0.3*3.0];
probdata.R     = eye(2);

% define objective function parameters
objparams   = struct;
calc_obj    = @objfun_nonlinear_lsf;

% define constraint function parameters
constparams = struct;
constparams.targetbeta = -norminv(0.01);
constparams.lsffunc    = @highly_nonlinear_lsf;
constparams.probdata   = probdata;

fun     = @(x) calc_obj(x,objparams);     % objective function
nonlcon = @(x) calc_const(x,constparams); % constraint function

A   = [];   b = [];    % coefficients of inequality constraints
Aeq = []; beq = [];    % coefficients of equality constraints

x0  = [2,2];                % intial guess
lb  = [0,0]; ub  = [15,15]; % upper and lower bounds

options = optimoptions('fmincon','display','none','algorithm','sqp');
[x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

disp('Highly nonlinear limit state results: ')
disp(['Optimum: ' num2str(x)]);
disp(['fstar: ' num2str(fval)]);