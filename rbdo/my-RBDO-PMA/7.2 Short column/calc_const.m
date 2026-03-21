function [c,ceq] = calc_const(dv,params)

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

Ro   = mod_corr(probdata.R,probdata);
Lo   = (chol(Ro))';
iLo  = inv(Lo);

% set starting point
x    = probdata.Xmu(:);
u    = trans_x2u(x,probdata,iLo);  % transform x to u

% set iteration number
i    = 1;

% set tolerances
eps1     = 10^-3;
eps2     = 10^-3;
max_iter = 20;

hist = struct;
hist.x = [];
hist.u = [];

Guold      = inf;
uold       = inf;

while true

    % transform u to x
    x   = trans_u2x(u,probdata, Lo);

    % compute jacobians at x = xi
    Jux = calc_jacobian(x,u,probdata,Lo,iLo);
    Jxu = inv(Jux);

    % compute limit state function and gradients
    isGrad   = true;
    [gx,dgx] = eval_lsf(x,isGrad,lsffunc,lsfparams,distparams);
    Gu       = gx;
    dGu      = (dgx'*Jxu)';

    % compute alpha
    alpha = -dGu/norm(dGu);

    hist.x = [hist.x,x];
    hist.u = [hist.u,u];

    % check convergence
    if norm(u-uold) < eps1 && abs(Gu-Guold) < eps2
        exit_message = 'MPTP found.';
        break;
    end

    if i == max_iter
        exit_message = 'Maximum iteration is observed.';
        break;
    end

    Guold = Gu;
    uold  = u;   

    % determine new iteration point
    % unew  = targetbeta*alpha; % Advanced mean value method
    if i==1
        unew       = targetbeta*alpha;
        alpha_old2 = alpha;
    elseif i==2
        unew       = targetbeta*alpha;
        alpha_old1 = alpha;
    else
        zeta = dot(alpha-alpha_old1,alpha_old1-alpha_old2);
        if zeta>0
            alpha_comb = alpha;
        else
            alpha_comb = ( alpha + alpha_old1 + alpha_old2 )...
                            /norm( alpha + alpha_old1 + alpha_old2);
        end
        unew       = targetbeta*alpha_comb; % Hybrid mean value method
        alpha_old2 = alpha_old1;
        alpha_old1 = alpha;
    end

    u     = unew;
    i     = i + 1;

end

% alpha = -dGu/norm(dGu)
% ustar = u
% xstar = trans_u2x(ustar,probdata, Lo)
gpstar = gx;

% inequality nonlinear constraint
c = [-gpstar;... 
    0.5*dv(1)-dv(2);... 
    dv(2)-2*dv(1)]; 
% c = -gpstar;

% equality nonlinear constraint
ceq = [];

end