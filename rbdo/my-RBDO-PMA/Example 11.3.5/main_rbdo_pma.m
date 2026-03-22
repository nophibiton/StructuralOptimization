% Performs reliability-based design optimization using double-loop
% approach. The performance measure approach (PMA) with hybrid-mean value method
% in [1] for the inverse reliability analysis is applied. The example
% solves Example 11.3.5 in [2].
%
% Reference:
%
% [1] Youn, B. D., Choi, K. K., & Park, Y. H. (2003). Hybrid analysis 
% method for reliability-based design optimization. Journal of Mechanical 
% Design, Transactions of the ASME, 125(2). https://doi.org/10.1115/1.1561042
%
% [2] Melchers, R. E., & Beck, A. T. (2018). Structural Reliability 
% Analysis and Prediction (2nd ed.). John Wiley & Sons.
%
clc,clear; format compact; format shortG;

objparams   = struct;
constparams = struct;
constparams.targetbeta = 4.0;

fun     = @(x) calc_obj(x,objparams);     % objective function
nonlcon = @(x) calc_const(x,constparams); % constraint function

A   = [];   b = [];    % coefficients of inequality constraints
Aeq = []; beq = [];    % coefficients of equality constraints

x0  = 1;               % intial guess
lb  = 0; ub  = 5;      % upper and lower bounds

options = optimoptions('fmincon','display','none','algorithm','sqp');
[x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

disp(['Optimum: ' num2str(x)]);