function [g_val,dg_val] = g_fun(x,d,params)

    % unpack
    %d      = params.d;

    sigmay = x(1);
    W      = x(2);
    M      = x(3);
    
    g_val = d*sigmay*W - M;


    %%% derivative of g(X,params) if available
    dg_val = [];

end