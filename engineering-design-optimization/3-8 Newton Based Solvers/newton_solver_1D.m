clc, clear; format compact; format longG;
%
% Performs Newton's method for a 1 dimensional case. Solves 
% example 3.7 (page 68) of the book of Martins and Ning (2021)
%
% Joaquim R. R. A. Martins and Andrew Ning. Engineering Design Optimization.
% Cambridge University Press, 2021. ISBN: 9781108833417.
% 

u(1) = 1.5;
k = 1;

tol = 10^-3;

while true

    u(k+1) = u(k) - myfun(u(k) ) / myderiv(u(k) ); % Equation 3.24
    
    if abs(u(k+1)-u(k)) <= tol, break; end

    k = k+1;
end


figure;

x = 0:0.01:1.5;
y = myfun(x); 
r = myfun(u);

plot(x,y); hold on;

plot(u,r,LineStyle='none', Marker='o',Color='r',MarkerFaceColor='g');

yline(0);


ylim([-5,16])

xlabel('u'); ylabel('r');
box on;


