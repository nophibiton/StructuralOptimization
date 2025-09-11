clc, clear; format compact; format longG;

% Performs backtracking line search. Solves example 4.8 (page 101) of the
% book of Martins and Ning (2021)
%
% Joaquim R. R. A. Martins and Andrew Ning. Engineering Design Optimization.
% Cambridge University Press, 2021. ISBN: 9781108833417.
% 
%


xk = [-1.25; 1.25];
pk = [4; 0.75];

rho = 0.7;
mu1 = 10^-4;
alpha_init = 1.2;

params{1} = mu1;
params{2} = rho;

alpha = backtrack(xk, pk, @myfun, @myderiv, alpha_init, params);

function alpha = backtrack(xk, pk, myfun, myderiv, alpha_init, params)
%
% Implements backtracking line search algorith. The following steps follows
% Algorithm 4.2 of the Martins and Ning (2021).
%
%
    mu1 = params{1};
    rho = params{2};

    alpha = alpha_init;
    phi0 = myfun(xk);
    dphi = myderiv(xk)' * pk;

    while true 

        alpha = alpha*rho;
        phi   = myfun (xk + alpha*pk);
        if phi <= phi0 + mu1*alpha*dphi, break; end

    end    

end