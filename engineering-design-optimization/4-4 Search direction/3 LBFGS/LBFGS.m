function [x_star, f_star, x_hist] = LBFGS(x0, tau, func, ls)
%
% Implements L-BFGS. The following steps follows
% Algorithm 4.8 (page 132) of the Martins and Ning (2021).
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


    n_var = length(xk);
    m = 3; % LBFGS parameter

    alfa_arr = zeros(1,m);
    sigma_arr = zeros(1,m);
    s_arr = zeros(n_var,m);
    y_arr = zeros(n_var,m);


    df = myderiv(x0,params);
    df_hist = [df_hist,df];

    k = 0;
    I = eye(length(x0));
    while true


        if k==0
            H = (1/norm(df))*I;

            pk = -H*df;

        else

            d = df;
            km = max(1,m-k+1);

            for i = m : -1 : km
                alfa_arr(1,i) = sigma_arr(1,i) * s_arr(:,i)' * d;
                d = d - alfa_arr(1,i) * y_arr(:,i);
            end        

            H_init = ((s_arr(:, km)'*y_arr(:,km)) / ((y_arr(:,km)'*y_arr(:,km))) ) * I;
            d = H_init*d;

            for i = km : m
                beta = sigma_arr(1, i) * y_arr(:,i)' * d;
                d = d + s_arr(:, i) * (alfa_arr(1, i) - beta);
            end

            pk = -d;

        end

        alpha_init = 1;

        alpha_k = line_search(xk, pk, myfun, myderiv, alpha_init, ls.params);
        alpha_hist = [alpha_hist,alpha_k];

        xk = xk + alpha_k*pk;
        x_hist = [x_hist, xk];

        df = myderiv(xk,params);
        df_hist = [df_hist,df];

        s = alpha_k*pk;
        y = df_hist(:,end) - df_hist(:,end-1);        
        sigma = 1 / (y'*s);

        %  Move each sigma, s, and y entry up one spot in their respective arrays.
        for i = 1 : m-1
            sigma_arr(1, i) = sigma_arr(1, i+1);
            s_arr(:, i) = s_arr(:, i+1);
            y_arr(:, i) = y_arr(:, i+1);
        end
        %  Set the the last entry in the rho, s, and y arrays to the most recent
        %  sigma, s, and y values, respectively.
        sigma_arr(1, m) = sigma;
        s_arr(:, m) = s;
        y_arr(:, m) = y;

        k = k + 1;

        if k==max_iter || norm(df,"inf") < tau
            x_star = xk;
            f_star = myfun(x_star,params);
            break; 
        end 
        

    end


end