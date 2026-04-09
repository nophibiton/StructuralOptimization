function [g_val,dg_val] = confun_g2_lsf(x,params)

    % unpack
    X1       = x(1);
    X2       = x(2);
    
    % calculate constraints
    g_val    = ((X1+X2-5)^2/30) + ((X1-X2-12)^2/120) - 1;

    %%% derivative of g(X,params) if available
    dg_val = [];

end