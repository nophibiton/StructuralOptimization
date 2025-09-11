clc, clear; format compact; format longG;

%
% Performs Newton's method for a 1 dimensional case using SECANT Method. Solves 
% example 3.7 (page 68) of the book of Martins and Ning (2021)
%
% Joaquim R. R. A. Martins and Andrew Ning. Engineering Design Optimization.
% Cambridge University Press, 2021. ISBN: 9781108833417.
% 
u(1) = 1.5;
k = 1;

tol = 10^-3;
delta = 10^-5;

while true

    r_k1 = myfun(u(k) + delta );
    r_k = myfun(u(k) );
    du = (r_k1 - r_k)/ (delta);
    u(k+1) = u(k) - r_k / du;    % Equation 3.26
    
    if abs(u(k+1)-u(k)) <= tol, break; end

    k = k+1;
end

x = 0:0.01:1.5;
y = myfun(x); 
r = myfun(u);

plot(x,y); hold on;

plot(u,r,LineStyle='none', Marker='o',Color='r',MarkerFaceColor='g');

yline(0);



xlabel('u');
ylabel('r');

ylim([-5,16])

box on;


