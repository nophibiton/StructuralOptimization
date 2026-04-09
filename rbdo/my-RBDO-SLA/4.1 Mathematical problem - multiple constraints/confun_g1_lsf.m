function [g_val,dg_val] = confun_g1_lsf(x,params)

    % unpack
    X1       = x(1);
    X2       = x(2);
    
    % calculate constraints
    g_val    = ((X2*X1^2)/20) - 1;

    %%% derivative of g(X,params) if available
    dg_val = [];

end