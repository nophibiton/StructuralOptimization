% Performs reliability-based design optimization using double-loop
% approach. The reliability index approach (RIA) with iHLRF algorithm
% in [1] for the reliability analysis is applied. The example
% solves Example 11.3.5 in [2].
%
% Reference:
%
% [1] Der Kiureghian, A. (2005). First- and Second-Order Reliability 
% Methods. In Engineering Design Reliability Handbook. CRC press.
%
% [2] Melchers, R. E., & Beck, A. T. (2018). Structural Reliability 
% Analysis and Prediction (2nd ed.). John Wiley & Sons.
%
clc,clear; format compact; format shortG;

% define probability data
probdata       = struct;
probdata.Xdist = {    'normal',     'normal',  'normal'};
probdata.Xmu   = [          40,           50,      1000];
probdata.Xstd  = [           5,           2.5,      200];
probdata.R     = eye(3);

% define objective function parameters
objparams   = struct;

% define constraint function parameters
constparams = struct;
constparams.targetbeta = 4.0;
constparams.lsffunc    = @g_fun;
constparams.probdata   = probdata;

fun     = @(x) calc_obj(x,objparams);     % objective function
nonlcon = @(x) calc_const(x,constparams); % constraint function

A   = [];   b = [];    % coefficients of inequality constraints
Aeq = []; beq = [];    % coefficients of equality constraints

x0  = 1;               % intial guess
lb  = 0; ub  = 5;      % upper and lower bounds

options = optimoptions('fmincon','display','none','algorithm','sqp');
[x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

disp(['Optimum: ' num2str(x)]);