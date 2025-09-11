clc, clear; format compact; format longG;
%
% Performs BFGS. Solves Bean Function (page 580) of the
% book of Martins and Ning (2021)
%
% Joaquim R. R. A. Martins and Andrew Ning. Engineering Design Optimization.
% Cambridge University Press, 2021. ISBN: 9781108833417.
% 

% initial condition and parameters
x0 = [-1;2];
tau = 10^-6;

% define parameters of objective function and its derivative
params = {};
func.params = {};
func.myfun = @bean;
func.myderiv = @beanderiv;

% parameters of bracketing/strong wolfe line search algorithm
mu1 = 10^-4;
mu2 = 0.7;
sigma = 2;
max_iter = 5;
ls.params{1} = mu1;
ls.params{2} = mu2;
ls.params{3} = sigma;
ls.params{4} = max_iter;
ls.func = @line_search_SW;

[x_star, f_star,x_hist] = BFGS(x0, tau, func, ls);

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
contour(X1,X2,Z); hold on;

plot(x_hist(1,:),x_hist(2,:),Marker="o",Color='r',MarkerFaceColor='r')

xlabel('x_1')
ylabel('x_2')

%%%


