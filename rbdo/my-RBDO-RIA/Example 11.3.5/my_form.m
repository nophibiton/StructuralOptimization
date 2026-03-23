clc,clear; format compact; format shortG;


lsfparams       = struct;
lsfparams.d     = 1;
lsffunc         = @g_fun;

Xdist = {    'normal',     'normal',  'normal'};
Xmu   = [          40,           50,      1000];
Xstd  = [           5,           2.5,      200];

R     = eye(3);

probdata       = struct;
probdata.Xdist = Xdist;
probdata.Xmu   = Xmu;
probdata.Xstd  = Xstd;
probdata.R     = R;

% get distribution parameters
distparams          = get_dist_params(probdata);
probdata.distparams = distparams;

Ro   = mod_corr(R,probdata);
Lo   = (chol(Ro))';
iLo  = inv(Lo);

% set starting point
x    = Xmu(:);
u    = trans_x2u(x,probdata,iLo);  % transform x to u

% set iteration number
i    = 0;

% set tolerances
eps1     = 10^-3;
eps2     = 10^-3;
max_iter = 20;

hist = struct;
hist.x = [];
hist.u = [];

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

    % set scale parameter
    if i == 0
        Go = Gu;
    end

    % compute alpha
    alpha = -dGu/norm(dGu);

    hist.x = [hist.x,x];
    hist.u = [hist.u,u];

    % check convergence
    if abs(Gu/Go) < eps1 && norm(u- alpha'*u*alpha) < eps2
        exit_message = 'MPP found.';
        break;
    end
    if i == max_iter
        exit_message = 'Maximum iteration is observed.';
        break;
    end

    % determine search direction
    d = (Gu/norm(dGu) + alpha'*u)*alpha-u;

    % determine step size
    %lambda = 1.0;
    lambda = calc_step_size(u,Gu,dGu,d, lsffunc,lsfparams,distparams,probdata,Lo);

    % determine new iteration point
    unew = u + lambda*d;
    u    = unew;
    i    = i + 1;

end

alpha = -dGu/norm(dGu)
ustar = u
xstar = trans_u2x(ustar,probdata, Lo)
beta  = alpha'*ustar
