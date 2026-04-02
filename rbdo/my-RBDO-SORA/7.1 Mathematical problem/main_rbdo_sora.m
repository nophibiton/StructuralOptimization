% Performs reliability-based design optimization using decoupled approach.
% The example uses Sequential optimization and reliability assessment
% (SORA) described in [1]. The example solves Example 7.1 Mathematical 
% problem in [2].
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

%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear limit state %
%%%%%%%%%%%%%%%%%%%%%%%%%   

% define objective function parameters
objparams   = struct;
calc_obj    = @objfun_nonlinear_lsf;

% define constraint function parameters
lsffunc    = @nonlinear_lsf;
lsfparams  = struct;
targetbeta = -norminv(0.01);

% define probability data
probdata       = struct;
probdata.Xdist = {    'normal',     'normal'};
probdata.Xmu   = [         5.0,          3.0];
probdata.Xstd  = [     0.3*5.0,      0.3*3.0];
probdata.R     = eye(2);

% define optimization options
options = optimoptions('fmincon','display','none','algorithm','sqp');

% plots
set(0,'defaultAxesFontSize',12)
% set(0,'defaultTextFontName','Times New Roman')
% set(0,'defaultAxesFontName','Times New Roman')
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get distribution parameters
distparams          = get_dist_params(probdata);
probdata.distparams = distparams;

Ro   = modify_corr(probdata.R,probdata);
Lo   = (chol(Ro))';
iLo  = inv(Lo);

hist            = struct;
hist.dvstar     = [];
hist.fobj_val   = [];
max_iter        = 20;
eps1            = 10^-4;

fobj_val_old = inf;

k           = 1;
dv0         = [2,4];
ximpp_prev  = probdata.Xmu';
dvstar_prev = dv0;

ind_d_ran_x         = [];
lsfparams.ind_d_det = 1:2;
lsfparams.ind_d_ran = [];

while true
    % (1) Deterministic optimization
    % perform deterministic optimization
    fun = @(d) calc_obj(d,objparams); % objective function

    % modify probabilistic constraint into a deterministic
    if ~isempty(lsfparams.ind_d_ran)
        lsfparams.si = dvstar_prev(ind_d_ran) - ximpp_prev(ind_d_ran_x);
    end
    det_g_fun = @(d)eval_lsf(ximpp_prev,d,false, lsffunc, lsfparams,distparams);
    nonlcon   = @(d)deal(-1*det_g_fun(d),[]); % constraint function
    
    % perform optimization
    d0                = dvstar_prev;
    [dvstar,fobj_val] = fmincon(fun,d0,[],[],[],[],[0,0],[15,15],nonlcon,options);

    % (2) Reliability assessment
    % perform inverse reliability analysis
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
set(0,'defaultAxesFontSize',12)
% set(0,'defaultTextFontName','Times New Roman')
% set(0,'defaultAxesFontName','Times New Roman')
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')

f1 = figure;
set(gcf,'units','inches','position',[1,1,6,3.3])
tiledlayout(1,2);

ax1 = nexttile;
xx1 = linspace(1,10,100);
xx2 = linspace(1,10,100);
[X1, X2] = meshgrid(xx1, xx2);
for i=1:size(X1,1)
    for j=1:size(X1,2)
        xx = [X1(i,j),X2(i,j)];  
        fobj_cont(i,j) = fun(xx);
        [gxx(i,j), ~] = nonlcon(xx);
    end
end

% plot objective function
contour(X1, X2, fobj_cont,10, '-','LineColor', [0.8 0.8 0.8], ...
                                              'LineWidth', 0.8); hold on;
lgdstr{1} = '';

% plot constraint function
contour(X1, X2, gxx, [0, 0], 'k-','LineWidth', 2); hold on;
lgdstr{2} = '$g_1(x,d)$';

% plot starting point
scatter(dv0(1),dv0(2),30,Marker="o", ...
    MarkerEdgeColor='k',MarkerFaceColor='r'); hold on;
lgdstr{3} = '$d_0$ - starting pt.';

% plot iteration
plot([dv0(1);hist.dvstar(1:end,1)],[dv0(2);hist.dvstar(1:end,1)],'k-',LineWidth=1.0, ...
    Marker='o',MarkerFaceColor='k', MarkerSize=2); hold on;
lgdstr{4} = '';

% plot optimum
scatter(hist.dvstar(end,1),hist.dvstar(end,2),100,Marker="pentagram", ...
    MarkerEdgeColor='k',MarkerFaceColor='r'); hold on;
lgdstr{5} = strcat('$x^* = $', num2str(dvstar(1),'%.3f'), ',',num2str(dvstar(2),'%.3f') );

xlim([1,10]);
ylim([1,10]);

xlabel('$x_1$',Interpreter='latex');
ylabel('$x_2$',Interpreter='latex');
legend(lgdstr,Interpreter="latex", Location="southoutside");

title('Nonlinear LSF')
box on;


disp('Nonlinear limit state results: ')
disp(['Optimum: ' num2str(dvstar)]);
disp(['fstar: ' num2str(fobj_val)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Highly nonlinear limit state %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define objective function parameters
objparams  = struct;
calc_obj   = @objfun_nonlinear_lsf;

% define constraint function parameters
targetbeta = -norminv(0.01);
lsffunc    = @highly_nonlinear_lsf;

% define probability data
probdata       = struct;
probdata.Xdist = {    'normal',     'normal'};
probdata.Xmu   = [         5.0,          3.0];
probdata.Xstd  = [     0.3*5.0,      0.3*3.0];
probdata.R     = eye(2);

% define optimization options
options = optimoptions('fmincon','display','none','algorithm','sqp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get distribution parameters
distparams          = get_dist_params(probdata);
probdata.distparams = distparams;

Ro   = modify_corr(probdata.R,probdata);
Lo   = (chol(Ro))';
iLo  = inv(Lo);

hist            = struct;
hist.dvstar     = [];
hist.fobj_val   = [];
max_iter        = 20;
eps1            = 10^-4;

fobj_val_old = inf;

k           = 1;
dv0         = [3,2];
ximpp_prev  = probdata.Xmu';
dvstar_prev = dv0;

ind_d_ran_x         = [];
lsfparams.ind_d_det = 1:2;
lsfparams.ind_d_ran = [];

while true
    % (1) Deterministic optimization
    % perform deterministic optimization
    fun = @(d) calc_obj(d,objparams); % objective function

    % modify probabilistic constraint into a deterministic
    if ~isempty(lsfparams.ind_d_ran)
        lsfparams.si = dvstar_prev(ind_d_ran) - ximpp_prev(ind_d_ran_x);
    end
    det_g_fun = @(d)eval_lsf(ximpp_prev, d, false, lsffunc, lsfparams,distparams);
    nonlcon   = @(d)deal(-1*det_g_fun(d),[]); % constraint function
    
    % perform optimization
    d0            = dvstar_prev;
    [dvstar,fobj_val] = fmincon(fun,d0,[],[],[],[],[0,0],[15,15],nonlcon,options);

    % (2) Reliability assessment
    % perform inverse reliability analysis
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots

ax2 = nexttile;

xx1 = linspace(0,5,100);
xx2 = linspace(0,5,100);
[X1, X2] = meshgrid(xx1, xx2);
for i=1:size(X1,1)
    for j=1:size(X1,2)
        xx = [X1(i,j),X2(i,j)];  
        fobj_cont(i,j) = fun(xx);
        [gxx(i,j), ~] = nonlcon(xx);
    end
end

% plot objective function
contour(X1, X2, fobj_cont,10, '-','LineColor', [0.8 0.8 0.8], ...
                                              'LineWidth', 0.8); hold on;
lgdstr{1} = '';

% plot constraint function
contour(X1, X2, gxx, [0, 0], 'k-','LineWidth', 2); hold on;
lgdstr{2} = '$g_1(x,d)$';

% plot starting point
scatter(dv0(1),dv0(2),30,Marker="o", ...
    MarkerEdgeColor='k',MarkerFaceColor='r'); hold on;
lgdstr{3} = '$d_0$ - starting pt.';

% plot iteration
plot([dv0(1);hist.dvstar(1:end,1)],[dv0(2);hist.dvstar(1:end,1)],'k-',LineWidth=1.0, ...
    Marker='o',MarkerFaceColor='k', MarkerSize=2); hold on;
lgdstr{4} = '';

% plot optimum
scatter(hist.dvstar(end,1),hist.dvstar(end,2),100,Marker="pentagram", ...
    MarkerEdgeColor='k',MarkerFaceColor='r'); hold on;
lgdstr{5} = strcat('$x^* = $', num2str(dvstar(1),'%.3f'), ',',num2str(dvstar(2),'%.3f') );

xlim([0,5]);
ylim([0,5]);

xlabel('$x_1$',Interpreter='latex');
ylabel('$x_2$',Interpreter='latex');
legend(lgdstr,Interpreter="latex", Location="southoutside");

title('Highly nonlinear LSF')
box on;

exportgraphics(gcf,'figure.pdf','ContentType','vector');

disp('Highly nonlinear limit state results: ')
disp(['Optimum: ' num2str(dvstar)]);
disp(['fstar: ' num2str(fobj_val)]);