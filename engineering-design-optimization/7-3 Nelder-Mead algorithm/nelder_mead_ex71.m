clc, clear; format compact; format longG;
%
% Performs Nelder-Mead algorithm. Solves example 7.1 (page 291) of the
% book of Martins and Ning (2021)
%
% Joaquim R. R. A. Martins and Andrew Ning. Engineering Design Optimization.
% Cambridge University Press, 2021. ISBN: 9781108833417.
% 

% define Nelder-Mead input parameters
x0           = [-1,2];
params.tau_x = 10^-3;
params.tau_f = 10^-3; 

% define parameters of objective function
func.params = {};
func.fobj   = @bean;

[xstar, fstar, hist] = nelder_mead(x0, func, params)


%%%Plot
x1 = linspace(-2.5,2.5);
x2 = linspace(-1,3);
[X1,X2] = meshgrid(x1,x2);

for ii=1:size(X1,1)
    for jj=1:size(X2,2)
        Z(ii,jj) = bean([X1(ii,jj),X2(ii,jj) ],params);
    end
end


figure;
set(gcf, 'units','inches','position',[1,1,3.5,3.5]);

contour(X1,X2,Z,15); hold on;

for ii=1:length(hist)
    V = hist(ii).x;
    patch('Faces',[1 2 3], ...
          'Vertices',V, ...
          'FaceColor','none', ...
          'EdgeColor','r', ...
          'EdgeAlpha',0.2, ...
          'LineWidth',1.5);
    hold on;
end

xlabel('x_1'); ylabel('x_2');

xlim([-2.5,2.5]);
ylim([-1,3]);


exportgraphics(gcf,'figure.pdf','ContentType','vector');

%%%

figure;
set(gcf, 'units','inches','position',[1,1,3.5,3.5]);

contour(X1,X2,Z,15); hold on;
scatter(xstar(1),xstar(2),Marker="pentagram", ...
        Color='r',MarkerFaceColor='r',MarkerEdgeColor='r'); hold on;

V = hist(1).x;
h = patch('Faces',[1 2 3], ...
      'Vertices',V, ...
      'FaceColor','none', ...
      'EdgeColor','r', ...
      'EdgeAlpha',0.2, ...
      'LineWidth',1.5);

for k = 1:length(hist)

    V = hist(k).x;
    set(h,'Vertices',V);

    title(sprintf('Step %d', k));

    drawnow;

    frame = getframe(gcf);
    img = frame2im(frame);
    [A,map] = rgb2ind(img,256);

    if k == 1
        imwrite(A,map,'contour.gif','gif','LoopCount',inf,'DelayTime',0.5);
    else
        imwrite(A,map,'contour.gif','gif','WriteMode','append','DelayTime',0.5);
    end
end

