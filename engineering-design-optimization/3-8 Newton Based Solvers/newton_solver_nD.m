clc, clear; format compact; format longG;

%
% Performs Newton's method for a n dimensional case. Solves 
% example 3.8 (page 70) of the book of Martins and Ning (2021)
%
% Joaquim R. R. A. Martins and Andrew Ning. Engineering Design Optimization.
% Cambridge University Press, 2021. ISBN: 9781108833417.
% 

u = [2;3] ;
k = 1;

tol = 10^-3;

while true

    J = myjacobian(u(:,k) );
    r_k = myfuns(u(:,k) );
    
    du = linsolve(J,-r_k');  % Solves Equation 3.30

    u(:,k+1) = u(:,k) + du;  % Equation 3.31

    if norm(u(:,k+1) - u(:,k) ) <= tol, break; end

    k = k+1;
end



