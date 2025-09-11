clc, clear; format compact; format longG;
%
% Performs BFGS. Solves example 4.10 (page 112) of the
% book of Martins and Ning (2021)
%
% Joaquim R. R. A. Martins and Andrew Ning. Engineering Design Optimization.
% Cambridge University Press, 2021. ISBN: 9781108833417.
% 

% initial condition and parameters
x0 = [10;4];
tau = 10^-6;

% define parameters of objective function and its derivative
beta = 5.0;
func.params.beta = beta;
func.myfun = @myfun;
func.myderiv = @myderiv;

% parameters of bracketing/strong wolfe line search algorithm
mu1 = 10^-4;
mu2 = 0.9;
sigma = 2;
max_iter = 10;
ls.params{1} = mu1;
ls.params{2} = mu2;
ls.params{3} = sigma;
ls.params{4} = max_iter;
ls.func = @line_search_SW;

[x_star, f_star,x_hist] = BFGS(x0, tau, func, ls);

%%%Plot
x1 = linspace(-5,12.5);
x2 = linspace(-7.5,7.5);
[X1,X2] = meshgrid(x1,x2);

params.beta = beta;

for ii=1:size(X1,1)
    for jj=1:size(X2,2)
        Z(ii,jj) = myfun([X1(ii,jj),X2(ii,jj) ],params);
    end
end

figure;
contour(X1,X2,Z); hold on;

plot(x_hist(1,:),x_hist(2,:),Marker="o",Color='r',MarkerFaceColor='r')

xlabel('x_1')
ylabel('x_2')

%%%


