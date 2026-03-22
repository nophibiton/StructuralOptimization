% Performs reliability-based design optimization using double-loop
% approach. The performance measure approach with hybrid-mean value method
% in [1] for the inverse reliability analysis is applied. The example
% solves Example 7.2 Short column in [2].
%
% Reference:
%
% [1] Youn, B. D., Choi, K. K., & Park, Y. H. (2003). Hybrid analysis 
% method for reliability-based design optimization. Journal of Mechanical 
% Design, Transactions of the ASME, 125(2). https://doi.org/10.1115/1.1561042
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

% define probability data
probdata       = struct;

% define objective function parameters
objparams   = struct;
calc_obj    = @objfun_column_lsf;

% define constraint function parameters
constparams = struct;
constparams.targetbeta = 3.0;
constparams.cov        = i_cov;
constparams.isNormal   = true; % Normal case
constparams.lsffunc    = @short_column_lsf;
constparams.probdata   = probdata;

fun     = @(x) calc_obj(x,objparams);     % objective function
nonlcon = @(x) calc_const(x,constparams); % constraint function

A   = [];   b = [];    % coefficients of inequality constraints
Aeq = []; beq = [];    % coefficients of equality constraints

x0  = [0.5,0.5];                % intial guess
lb  = [0.1,0.1]; ub  = []; % upper and lower bounds

options = optimoptions('fmincon','display','none','algorithm','sqp');
[x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

fprintf('\n');
disp(['Short column results, Normal case (COV=' num2str(constparams.cov) '):'])
disp(['Optimum: ' num2str(x)]);
disp(['fstar: ' num2str(fval)]);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Short column Non-normal case %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i_cov = [0.05,0.10,0.15]

% define probability data
probdata       = struct;

% define objective function parameters
objparams   = struct;
calc_obj    = @objfun_column_lsf;

% define constraint function parameters
constparams = struct;
constparams.targetbeta = 3.0;
constparams.cov        = i_cov;
constparams.isNormal   = false; % Non-normal case
constparams.lsffunc    = @short_column_lsf;
constparams.probdata   = probdata;

fun     = @(x) calc_obj(x,objparams);     % objective function
nonlcon = @(x) calc_const(x,constparams); % constraint function

A   = [];   b = [];    % coefficients of inequality constraints
Aeq = []; beq = [];    % coefficients of equality constraints

x0  = [0.5,0.5];                % intial guess
lb  = [0.1,0.1]; ub  = []; % upper and lower bounds

options = optimoptions('fmincon','display','none','algorithm','sqp');
[x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

fprintf('\n');
disp(['Short column results, Nonnormal case (COV=' num2str(constparams.cov) '):'])
disp(['Optimum: ' num2str(x)]);
disp(['fstar: ' num2str(fval)]);

end
