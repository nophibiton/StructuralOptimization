clc, clear; format compact; format shortG;

addpath(genpath(strcat(pwd,'\gpml'     )) );
warning off

%
% Performs Efficient global optimization. Solves example 10.10 (page 417) of the
% book of Martins and Ning (2021)
%
% Joaquim R. R. A. Martins and Andrew Ning. Engineering Design Optimization.
% Cambridge University Press, 2021. ISBN: 9781108833417.
% 

% define functions
f         = @(x) exp(-0.1*x) .* sin(x);
surrogate = @gp_model;
infill    = @EI;

% convergence criteria
max_iter = 10;
tau      = 10^-6;

% initial sample points
xi  = [0.5, 2, 2.5, 9, 10]';

% evaluate sample points
fi  = f(xi);

% best point so far
[fstar,idx_star] = min(fi);
xstar            = xi(idx_star);

k = 1;

hist(k).X     = xi;
hist(k).y     = fi;
hist(k).xstar = xstar;
hist(k).fstar = fstar;


while true

    % construct surrogate model
    surrogateModel = surrogate(xi,fi);

    % maximize expected improvement
    [xk,EI_x] = infill(surrogateModel,xi,fi);

    hist(k).EI_x  = EI_x;
    hist(k).surrogateModel = surrogateModel;

    % evaluate true function at predicted optimum
    fk = f(xk);

    % add new point
    xi = [xi; xk];
    fi = [fi;fk];

    % best point so far
    [fstar,idx_star] = min(fi);
    xstar            = xi(idx_star);


    k = k+1;

    hist(k).X     = xi;
    hist(k).y     = fi;
    hist(k).xstar = xstar;
    hist(k).fstar = fstar;

    if (k >= max_iter) || (max(EI_x) < tau), break; end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plots

for k = [1,2,3,7]

    set(0,'defaultAxesFontSize',12)
    set(0,'defaultTextFontName','Cambria Math')
    set(0,'defaultAxesFontName','Cambria Math')

    figure;
    set(gcf,'units','inches','position',[1,1,4.0,3.0]);

    % input values
    x      = 0:0.001:15;

    % (1) plot exact function
    yexact = f(x);
    plot(x,yexact,'k-','LineWidth',1.0); hold on;

    % (2) plot predicted by surrogate model
    surrogateModel = hist(k).surrogateModel;

    % unpack parameters of the surrogate model
    hyp            = surrogateModel.hyp;
    meanfunc       = surrogateModel.meanfunc;
    covfunc        = surrogateModel.covfunc;
    likfunc        = surrogateModel.likfunc;

    [mu, s2] = gp(hyp, @infGaussLik, meanfunc, covfunc, ...
                                             likfunc,hist(k).X, hist(k).y, x');
    ypred = mu;
    plot(x,ypred,'color',[0.53 0.81 0.92],'LineWidth',1.5); hold on;

    x  = x(:); mu = mu(:); s2 = s2(:); sigma = sqrt(s2);
    fill([x; flipud(x)], ...
         [mu + 1*sigma; flipud(mu - 1*sigma)], ...
         [0.53 0.81 0.92], ...   % sky blue
         'FaceAlpha', 0.3, ...
         'EdgeColor', 'none');


    % (3) plot sample points
    s = scatter(hist(1).X,hist(1).y,40,'rp','MarkerFaceColor','r'); hold on;
    x_added = hist(k).X(~ismember(hist(k).X, hist(1).X));
    y_added = hist(k).y(~ismember(hist(k).y, hist(1).y));
    s = scatter(x_added,y_added,20,'ko','MarkerFaceColor','k'); hold on;
    legend({'Exact','GP Model','','Initial points','Added points'});


    ylim([-1.5,2.5])
    xlim([0,12])
    xlabel('x')
    ylabel('f(x)')
    box on;

    title(strcat('k = ', num2str(k)))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Animation

set(0,'defaultAxesFontSize',12)
set(0,'defaultTextFontName','Cambria Math')
set(0,'defaultAxesFontName','Cambria Math')

figure;
set(gcf,'units','inches','position',[1,1,4.0,3.0],'Color','w');

scale = 1;

% (1) plot exact function
yexact = f(x);
exact_plot = plot(x,yexact,'k-','LineWidth',1.3);hold on;

% (2) plot initial points
initial_points = scatter(hist(1).X,hist(1).y,40,'rp','MarkerFaceColor','r'); 

% (3) plot GP
gp_plot = plot(NaN,NaN,'color',[0.53 0.81 0.92],'LineWidth',1.5);
gp_fill = fill(NaN, ...
             NaN, ...
             [0.53 0.81 0.92], ...   % sky blue
             'FaceAlpha', 0.3, ...
             'EdgeColor', 'none');

% (4) plot added points
s = scatter(NaN,NaN,20,'ko','MarkerFaceColor','k'); 

legend({'Exact','Initial points','GP Model','','Added points'},Box="off");

ylim([-1.5,2.5])
xlim([0,15])
xlabel('x')
ylabel('f(x)')

uistack([initial_points,s,exact_plot],'top');

for k = 1:8

    % input values
    x      = 0:0.001:15;

    % (2) plot predicted by surrogate model
    surrogateModel = hist(k).surrogateModel;

    % unpack parameters of the surrogate model
    hyp            = surrogateModel.hyp;
    meanfunc       = surrogateModel.meanfunc;
    covfunc        = surrogateModel.covfunc;
    likfunc        = surrogateModel.likfunc;

    [mu, s2] = gp(hyp, @infGaussLik, meanfunc, covfunc, ...
                                             likfunc,hist(k).X, hist(k).y, x');
    ypred = mu;
    gp_plot.XData = x;
    gp_plot.YData = ypred;

    x  = x(:); mu = mu(:); s2 = s2(:); sigma = sqrt(s2);
    gp_fill.XData = [x; flipud(x)];
    gp_fill.YData = [mu + 1*sigma; flipud(mu - 1*sigma)];

    % (3) plot sample points
    x_added = hist(k).X(~ismember(hist(k).X, hist(1).X));
    y_added = hist(k).y(~ismember(hist(k).y, hist(1).y));
    s.XData = x_added;
    s.YData = y_added;

    title(strcat('k = ', num2str(k)))

    drawnow;

    frame    = getframe(gcf);
    img      = frame2im(frame);
    imgSmall = imresize(img, scale, 'nearest');   % downsample
    img      = imresize(imgSmall, size(img(:,:,1)), 'nearest');  % upsample
    [A,map] = rgb2ind(img,256);

    if k == 1
        imwrite(A,map,'kriging_example.gif','gif','LoopCount',inf,'DelayTime',0.7);
    else
        imwrite(A,map,'kriging_example.gif','gif','WriteMode','append','DelayTime',0.7);
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rmpath(genpath(strcat(pwd,'\gpml'     ))   );