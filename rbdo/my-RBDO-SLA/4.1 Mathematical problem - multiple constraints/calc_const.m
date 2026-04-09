function [c,ceq] = calc_const(dv,params)

global xk_prev xk_buffer last_dk_evaluated

% unpack parameters
targetbeta      = params.targetbeta;
lsffunc         = params.lsffunc;

% define limit state parameters/decision variables
lsfparams.d     = dv;
mu_x1           = dv(1);
mu_x2           = dv(2);

% define probability data
probdata        = struct;
probdata.Xdist  = { 'normal','normal'};
probdata.Xmu    = [ mu_x1, mu_x2];
probdata.Xstd   = [   0.3,   0.3];
probdata.R      = eye(2);

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
        
    for j=1:numel(lsffunc)

        % compute limit state function and gradients
        [~,dgx] = eval_lsf(x(:),true,lsffunc{j},lsfparams,distparams);
        Jux     = calc_jacobian(x(:),u,probdata,Lo,iLo);
        Jxu     = inv(Jux);    
        dGu     = (dgx'*Jxu)';

        xk(j,:) = probdata.Xmu + (-dGu/norm(dGu))'.*probdata.Xstd*targetbeta;
        xk_prev(j,:) = xk(j,:); % initialize
    end

end

% only recompute if mu_k or dk changed (see flowchat of Fig.5 in Ref[1])
if isempty(last_dk_evaluated) || norm(dv - last_dk_evaluated) > 1e-12

    for j=1:numel(lsffunc)
        x    = xk_prev(j,:);
        u    = trans_x2u(x,probdata,iLo);  % transform x to u
        
        % compute limit state function and gradients
        [~,dgx] = eval_lsf(x(:),true,lsffunc{j},lsfparams,distparams);
        Jux     = calc_jacobian(x(:),u,probdata,Lo,iLo);
        Jxu     = inv(Jux);
        dGu     = (dgx'*Jxu)';
            
        xk(j,:) = probdata.Xmu + (-dGu/norm(dGu))'.*probdata.Xstd*targetbeta;
        xk_buffer(j,:) = xk(j,:); % store it in buffer
    end

    last_dk_evaluated = dv;    % update tracker
else
    % reuse previous computation
    xk = xk_buffer;
end

% evaluate lsf at xk
gpstar = zeros(numel(lsffunc),1);
for j=1:numel(lsffunc)
    isGrad = false;
    x      = xk(j,:);
    gxk    = eval_lsf(x(:),isGrad,lsffunc{j},lsfparams,distparams);
    
    gpstar(j) = gxk;
end

% inequality nonlinear constraint
c = -gpstar;

% equality nonlinear constraint
ceq = [];

end