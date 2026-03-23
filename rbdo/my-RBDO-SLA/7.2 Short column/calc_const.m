function [c,ceq] = calc_const(dv,params)

global xk_prev xk_buffer last_dk_evaluated

targetbeta = params.targetbeta;
cov_bh     = params.cov;
isNormal   = params.isNormal;
% probdata   = params.probdata;
lsffunc    = params.lsffunc;

% define limit state parameters/decision variables
lsfparams.d     = dv;
mu_h            = dv(1);
mu_b            = dv(2);

% define base units
kN  = 1;
m   = 1;
kPa = kN/m^2;
MPa = 1000*kPa;

% define probability data
probdata       = struct;
if isNormal
    probdata.Xdist = { 'normal','normal','normal','normal','normal','normal'};
else
    probdata.Xdist = { 'Gumbel','Gumbel','Gumbel','weibull','Lognormal','Lognormal'};
end
probdata.Xmu   = [ 2500*kN, 250*kN*m,125*kN*m, 40*MPa,mu_h, mu_b];
probdata.Xstd  = [ 0.2*2500*kN, 0.3*250*kN*m, 0.3*125*kN*m, 0.1*40*MPa, cov_bh*mu_h,cov_bh*mu_b];
probdata.R     = eye(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get distribution parameters
distparams          = get_dist_params(probdata);
probdata.distparams = distparams;

Ro   = modify_corr(probdata.R,probdata);
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
if isempty(last_dk_evaluated) || norm(dv - last_dk_evaluated) > 1e-12

    x    = xk_prev;
    u    = trans_x2u(x,probdata,iLo);  % transform x to u
    
    % compute limit state function and gradients
    [~,dgx] = eval_lsf(x,true,lsffunc,lsfparams,distparams);
    Jux     = calc_jacobian(x,u,probdata,Lo,iLo);
    Jxu     = inv(Jux);
    dGu     = (dgx'*Jxu)';
        
    xk = probdata.Xmu + (-dGu/norm(dGu))'.*probdata.Xstd*targetbeta;

    xk_buffer         = xk(:); % store it in buffer
    last_dk_evaluated = dv;    % update tracker
else
    % reuse previous computation
    xk = xk_buffer;
end
% dv
% xk_prev
% xk_buffer
% last_dk_evaluated
% xk

isGrad = false;
gxk    = eval_lsf(xk,isGrad,lsffunc,lsfparams,distparams);

gpstar = gxk;

% inequality nonlinear constraint
c = [-gpstar;... 
    0.5*dv(1)-dv(2);... 
    dv(2)-2*dv(1)];

% equality nonlinear constraint
ceq = [];

end