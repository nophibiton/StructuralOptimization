function [f_val,df_val] = objfun_mathematical_lsf(d,params)

    % unpack
    d1     = d(1);
    d2     = d(2);
    
    f_val  = d1+d2;

    %%% derivative of objective function f(), if available
    df_val = [];

end