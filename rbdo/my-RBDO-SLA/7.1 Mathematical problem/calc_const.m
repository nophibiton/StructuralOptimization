function [c,ceq] = calc_const(dv,params)

global xk_prev xk_buffer last_x_evaluated

targetbeta = params.targetbeta;
probdata   = params.probdata;
lsffunc    = params.lsffunc;

% define limit state parameters/decision variables
lsfparams.d     = dv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get distribution parameters
distparams          = get_dist_params(probdata);
probdata.distparams = distparams;

Ro   = mod_corr(probdata.R,probdata);
Lo   = (chol(Ro))';
iLo  = inv(Lo);

if isempty(xk_prev)

    % set starting point
    x    = probdata.Xmu(:);
    u    = trans_x2u(x,probdata,iLo);  % transform x to u
        
    % compute limit state function and gradients
    [~,dgx] = eval_lsf(x,true,lsffunc,lsfparams,distparams);
    Jux     = calc_jacobian(x,u,probdata,Lo,iLo);
    Jxu     = inv(Jux);    
    dGu     = (dgx'*Jxu)';
    
    xk = probdata.Xmu + (-dGu/norm(dGu))'.*probdata.Xstd*targetbeta;
    xk_prev = xk(:); % initialize
end

% only recompute if dk changed
if isempty(last_x_evaluated) || norm(dv - last_x_evaluated) > 1e-12

    x    = xk_prev;
    u    = trans_x2u(x,probdata,iLo);  % transform x to u
    
    % compute limit state function and gradients
    [~,dgx] = eval_lsf(x,true,lsffunc,lsfparams,distparams);
    Jux     = calc_jacobian(x,u,probdata,Lo,iLo);
    Jxu     = inv(Jux);
    dGu     = (dgx'*Jxu)';
        
    xk = probdata.Xmu + (-dGu/norm(dGu))'.*probdata.Xstd*targetbeta;

    xk_buffer        = xk(:); % store it in buffer
    last_x_evaluated = dv;    % update tracker
else
    % reuse previous computation
    xk = xk_buffer;
end

isGrad = false;
gxk    = eval_lsf(xk,isGrad,lsffunc,lsfparams,distparams);

gpstar = gxk;

% inequality nonlinear constraint
c = -gpstar;

% equality nonlinear constraint
ceq = [];

end