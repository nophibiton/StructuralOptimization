function [x_star, f_star, x_hist] = BFGS(x0, tau, func, ls)
%
% Implements BFGS. The following steps follows
% Algorithm 4.7 (page 129) of the Martins and Ning (2021).
%
% Joaquim R. R. A. Martins and Andrew Ning. Engineering Design Optimization.
% Cambridge University Press, 2021. ISBN: 9781108833417.
% 
    params = func.params;
    myfun = func.myfun;
    myderiv = func.myderiv;

    line_search = ls.func;
    ls.params{length(ls.params)+1} = params;

    max_iter = 50;
    xk = x0;
    

    x_hist = xk;
    alpha_hist = [];
    df_hist = [];
    pk_hist = [];


    k = 0;
    I = eye(length(x0));
    while true
        df = myderiv(xk,params);
        df_hist = [df_hist,df];

        if k==0
            H = (1/norm(df))*I;
        else
            s = x_hist(:,end) - x_hist(:,end-1);
            y = df_hist(:,end) - df_hist(:,end-1);
            sigma = 1 / (y'*s);
            H = (I - sigma*s*y') * H * (I - sigma*y*s') + sigma*s*s';
        end
        
        pk = -H*df;

        alpha_init = 1;

        alpha_k = line_search(xk, pk, myfun, myderiv, alpha_init, ls.params);
        alpha_hist = [alpha_hist,alpha_k];

        xk = xk + alpha_k*pk;
        x_hist = [x_hist, xk];

        k = k + 1;

        if k==max_iter || norm(df,"inf") < tau
            x_star = xk;
            f_star = myfun(x_star,params);
            break; 
        end 
        

    end


end