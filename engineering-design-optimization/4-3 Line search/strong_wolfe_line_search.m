clc, clear; format compact; format longG;
%
% Performs backtracking line search. Solves example 4.9 (page 107) of the
% book of Martins and Ning (2021)
%
% Joaquim R. R. A. Martins and Andrew Ning. Engineering Design Optimization.
% Cambridge University Press, 2021. ISBN: 9781108833417.
% 

xk = [-1.25; 1.25];
pk = [4; 0.75];

mu1 = 10^-4;
mu2 = 0.9;
sigma = 2;

max_iter = 10;

fun_params = {};

params{1} = mu1;
params{2} = mu2;
params{3} = sigma;
params{4} = max_iter;
params{5} = fun_params;

alpha_init = 1.2;

aalpha = 0.01:0.01:1.2;

for ii=1:length(aalpha), pphi(ii) = myfun(xk + aalpha(ii)*pk); end
plot(aalpha,pphi);

[alpha_star] = line_search_SW(xk, pk, @myfun, @myderiv, alpha_init, params);

function [alpha_star] = line_search_SW(xk, pk, myfun, myderiv, alpha_init, params)
%
% Implements bracketing strategy for line search algorith. The following steps follows
% Algorithm 4.3 (page 104) of the Martins and Ning (2021).
%
% Strong wolfe (SW) conditions
% (SW1) phi(alpha) <= phi(0) + mu1 * alpha * dphi(0)  <- sufficient decrease
% (SW2) |dphi(alpha)| <= mu2 * |dphi(0)|  <- sufficient curvature
%

    mu1 = params{1};
    mu2 = params{2};
    sigma = params{3};
    %max_iter = params{4};
    fun_params = params{5};
    

    alpha1 = 0;
    alpha2 = alpha_init;

    phi0 = myfun(xk,fun_params);
    dphi0 = myderiv(xk,fun_params)' * pk;

    first = true;
    % BRACKETTING phase: find (a) INTERVAL that contains 
    % strong wolfe conditions or (b) a POINT that satisfies it
    while true

        phi2 = myfun (xk + alpha2*pk, fun_params);
        phi1 = myfun (xk + alpha1*pk, fun_params);

        if (phi2 > phi0 + mu1*alpha2*dphi0) || (~first && (phi2 >= phi1))
            % bracketing: SW1 is not met but possible min in interval can 
            % found
            alpha_star = pinpoint(alpha1, alpha2, xk, pk, @myfun, @myderiv, phi0, dphi0,params);
            break;
        end

        dphi2 = myderiv(xk + alpha2*pk,fun_params)' * pk;
        if abs(dphi2) <= -mu2*dphi0
            % SW1 and SW2 condtions are met
            alpha_star = alpha2;
            break;
        elseif dphi2 >= 0
            % positive slope but SW2 condition is not satisfied
            alpha_star = pinpoint(alpha2, alpha1, xk, pk, @myfun, @myderiv, phi0, dphi0,params);
            break;
        else
            % increase step length
            alpha1 = alpha2;
            alpha2 = sigma*alpha2;
        end

        first = false;
    end


end

function [alpha_star] = pinpoint(alpha_l, alpha_h, xk, pk, myfun, myderiv, phi0, dphi0,params)
%
% Implements pinpointing strategy for line search algorith. The following steps follows
% Algorithm 4.4 (page 106) of the Martins and Ning (2021).
% 
    mu1 = params{1};
    mu2 = params{2};
    %sigma = params{3};
    max_iter = params{4};
    fun_params = params{5};

    phi_l = myfun(xk + alpha_l*pk,fun_params);

    k = 0;
    while true
        
        % bisection metod or you can use interpolation (quad, cubic)
        alpha_p = 0.5*(alpha_l+alpha_h); % uses bisection method
        
        phi_p = myfun(xk + alpha_p*pk,fun_params);
        if phi_p > phi0+mu1*alpha_p*dphi0 || phi_p > phi_l 
            alpha_h = alpha_p;
        else
            dphi_p = myderiv(xk + alpha_p*pk,fun_params)' * pk;
            if abs(dphi_p) <= -mu2*dphi0
                % SW1 and SW2 are satisfied
                alpha_star = alpha_p;
                break;
            elseif dphi_p*(alpha_h-alpha_l) >=0
                alpha_h = alpha_l;
            end
            alpha_l = alpha_p;
        end

        if k==max_iter, alpha_star = alpha_p; break; end

        k = k+1;
    end


end