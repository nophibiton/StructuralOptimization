% Performs reliability-based design optimization using single-loop
% approach. This single loop approach implements the study in [1]. 
% The example solves Example 4.1 Mathematical problem in [1].
%
% Reference:
%
% [1] Liang, J., Mourelatos, Z. P., & Tu, J. (2004). A Single-Loop Method 
% for Reliability-Based Design Optimization. 419–430. 
% https://doi.org/10.1115/DETC2004-57255
%

clc,clear; format compact; format shortG;
global xk_prev xk_buffer last_dk_evaluated hist
hist = struct;
hist.dk = [];

xk_prev=[];
xk_buffer=[];
last_dk_evaluated=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.1 A Mathematical Example %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define probability data
probdata       = struct;

% define objective function parameters
objparams   = struct;
calc_obj    = @objfun_mathematical_lsf;

% define constraint function parameters
constparams = struct;
constparams.targetbeta = 3.0; % for all constraints, j=1,2,3
constparams.lsffunc{1}    = @confun_g1_lsf;
constparams.lsffunc{2}    = @confun_g2_lsf;
constparams.lsffunc{3}    = @confun_g3_lsf;
constparams.probdata   = probdata;

fun     = @(x)calc_obj(x,objparams);     % objective function
nonlcon = @(x)calc_const(x,constparams); % constraint function

A   = [];   b = [];         % coefficients of inequality constraints
Aeq = []; beq = [];         % coefficients of equality constraints

x0  = [5,5];                % intial guess
lb  = [0,0]; ub  = [10,10]; % upper and lower bounds

options = optimoptions('fmincon','OutputFcn',@outfun, ...
                                     'display','none','algorithm','sqp');
[x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

fprintf('\n4.1 A mathematical example results:\n')
fmt = ['xstar = [%.4f' repmat(' %.4f', 1, numel(x)-1) '].\n'];
fprintf(fmt,x);
fprintf('fstar = %.4f.\n',fval);
