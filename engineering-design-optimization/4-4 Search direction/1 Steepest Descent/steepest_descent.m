function [x_star, f_star, x_hist] = steepest_descent(x0, tau, func, ls)
%
% Implements steepest descent. The following steps follows
% Algorithm 4.5 (page 111) of the Martins and Ning (2021).
%
% Joaquim R. R. A. Martins and Andrew Ning. Engineering Design Optimization.
% Cambridge University Press, 2021. ISBN: 9781108833417.
% 

    params = func.params;
    myfun = func.myfun;
    myderiv = func.myderiv;

    line_search = ls.func;
    ls.params{length(ls.params)+1} = params;

    k = 0;
    xk = x0;

    max_iter = 500;
    
    x_hist = xk;
    alpha_hist = [];
    df_hist = [];
    pk_hist = [];

    while true
        df = myderiv(xk,params);
        pk = -df/norm(df);

        df_hist = [df_hist,df];
        pk_hist = [pk_hist,pk];

        if k==0, alpha_init = 1;
        else, alpha_init = alpha_hist(end)*((df_hist(:,end-1)'*pk_hist(:,end-1))/(df_hist(:,end)'*pk_hist(:,end)) ); end

        alpha_k = line_search(xk, pk, myfun, myderiv, alpha_init, ls.params);
        alpha_hist = [alpha_hist,alpha_k];

        xk = xk + alpha_k*pk;
        x_hist = [x_hist, xk];

        k = k+1;

        if k==max_iter || norm(df,"inf") < tau
            x_star = xk;
            f_star = myfun(x_star,params);
            break; 
        end

    end


end