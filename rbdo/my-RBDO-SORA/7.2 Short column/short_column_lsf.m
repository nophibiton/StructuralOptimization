function [g_val,dg_val] = short_column_lsf(x,d,params)

    % unpack
    isDeterministic = params.isDeterministic;
    ind_d_det = params.ind_d_det;
    ind_d_ran = params.ind_d_ran;    

    if isDeterministic
        if ~isempty(ind_d_ran)
            d_ran = d(ind_d_ran) - params.si;
        end
        d_det  = d(ind_d_det);

        h      = d_ran(1);
        b      = d_ran(2);        
    else
        h      = x(5);
        b      = x(6);
    end

    F      = x(1);
    M1     = x(2);
    M2     = x(3);
    fy     = x(4);    
    
    g_val = 1 - ( (4*M1)/(fy*b*h^2) ) - ( (4*M2)/(fy*h*b^2) ) - ( (F^2)/(fy*b*h)^2 );


    %%% derivative of g(X,params) if available
    dg_val = [];

end