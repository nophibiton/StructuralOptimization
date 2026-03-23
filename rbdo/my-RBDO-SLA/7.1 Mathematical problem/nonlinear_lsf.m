function [g_val,dg_val] = nonlinear_lsf(x,params)

    % unpack
    d      = params.d;
    d1     = d(1);
    d2     = d(2);

    X1     = x(1);
    X2     = x(2);
    
    g_val = (1/5)*d1*d2*X2^2 - X1;


    %%% derivative of g(X,params) if available
    dg_val = [];

end