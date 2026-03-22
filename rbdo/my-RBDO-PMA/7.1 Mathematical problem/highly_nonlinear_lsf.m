function [g_val,dg_val] = highly_nonlinear_lsf(x,params)

    % unpack
    d      = params.d;
    d1     = d(1);
    d2     = d(2);

    X1     = x(1);
    X2     = x(2);
    
    g_val = d1*d2*X2 - log(X1);


    %%% derivative of g(X,params) if available
    dg_val = [];

end