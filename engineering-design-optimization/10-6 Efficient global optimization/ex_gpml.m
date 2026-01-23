clc, clear; format compact; format shortG;

addpath(genpath(strcat(pwd,'\gpml'     )) );

set(0,'defaultAxesFontSize',12)
set(0,'defaultTextFontName','Cambria Math')
set(0,'defaultAxesFontName','Cambria Math')

% exact function
f   = @(x) exp(-0.1*x) .* sin(x); 

% sample points
xi  = [0.5, 2, 2.5, 9, 10]';
yi  = f(xi);

% create Kriging model
meanfunc = @meanZero;
covfunc  = @covSEiso; 
likfunc  = @likGauss;

hyp.cov  = [log(0.5); log(1.0)];
hyp.mean = [];
hyp.lik  = log(1e-6);

hyp = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, ...
                            covfunc, likfunc, xi, yi);
[mu_tr, s2_tr] = gp(hyp, @infGaussLik, ...
                      meanfunc, covfunc, likfunc, xi, yi, xi);

ytrain = yi;
% Performance metrics
err = ytrain - mu_tr;

RMSE = sqrt(mean(err.^2));
MAE  = mean(abs(err));

SS_res = sum(err.^2);
SS_tot = sum((ytrain - mean(ytrain)).^2);
R2 = 1 - SS_res/SS_tot;

fprintf('Training metrics:\n');
fprintf('RMSE = %.4f\n', RMSE);
fprintf('MAE  = %.4f\n', MAE);
fprintf('R2   = %.4f\n', R2);


figure;
set(gcf,'units','inches','position',[1,1,4.0,3.0]);



% exact function
xx     = 0:0.01:15;
yexact = f(xx);
plot(xx,yexact,'k-','LineWidth',1.0); hold on;
lgdstr{1} = 'Exact';

[yy_pred, s2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc,xi, yi, xx');
mu    = yy_pred;


plot(xx,yy_pred,'color',[0.53 0.81 0.92],'LineWidth',1.5); hold on;
lgdstr{2} = 'Kriging';


% sample points
s = scatter(xi,yi,'ro','MarkerFaceColor','r'); hold on;
lgdstr{3} = 'Sample points';


xx = xx(:);
mu = mu(:);
s2 = s2(:);
sigma = sqrt(s2);

fill([xx; flipud(xx)], ...
     [mu + 1*sigma; flipud(mu - 1*sigma)], ...
     [0.53 0.81 0.92], ...   % sky blue
     'FaceAlpha', 0.3, ...
     'EdgeColor', 'none');


ylim([-1.5,1.5])
xlim([0,12])

box on;


rmpath(genpath(strcat(pwd,'\gpml'     ))   );
