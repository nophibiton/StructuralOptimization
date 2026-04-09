function [g_val,dg_val] = confun_mathematical_lsf(x,params)

    % unpack
    X1       = x(1);
    X2       = x(2);
    
    % calculate constraints
    g_val    = zeros(1,3);
    g_val(1) = ((X2*X1^2)/20) - 1;
    g_val(2) = ((X1+X2-5)^2/30) + ((X1-X2-12)^2/120) - 1;
    g_val(3) = (80/(X1^2+8*X2+5)) - 1;


    %%% derivative of g(X,params) if available
    dg_val = [];

end