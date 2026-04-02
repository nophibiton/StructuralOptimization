function [g_val,dg_val] = nonlinear_lsf(x,d,params)

    % unpack
    ind_d_det = params.ind_d_det;
    ind_d_ran = params.ind_d_ran;
    
    if ~isempty(ind_d_ran)
        d_ran = d(ind_d_ran) - params.si;
    end

    d_det = d(ind_d_det);

    %d      = params.d;
    d1     = d_det(1);
    d2     = d_det(2);

    X1     = x(1);
    X2     = x(2);
    
    g_val  = (1/5)*d1*d2*X2^2 - X1;


    %%% derivative of g(X,params) if available
    dg_val = [];

end