function [xstar, fstar, hist] = nelder_mead(x0, func, params)
%
% Implements Nelder-Mead algorithm. The following steps follows
% Algorithm 7.1 (page 290) of the Martins and Ning (2021).
%
% Joaquim R. R. A. Martins and Andrew Ning. Engineering Design Optimization.
% Cambridge University Press, 2021. ISBN: 9781108833417.
% 

xstar = inf;
fstar = inf;
hist  = struct;

% unpack parameters
tau_x     = params.tau_x;
tau_f     = params.tau_f;

objparams = func.params;
calc_obj  = func.fobj;

% create simplex with edge length l
n = length(x0);
X = zeros(n+1,n);
l = 1.0; % edge length
for i=1:n
    s = zeros(1,n);
    for j=1:n
        if j==i
            s(j) = (1/(n*sqrt(2)))*(sqrt(n+1)-1) + (l/sqrt(2));
        else
            s(j) = (1/(n*sqrt(2)))*(sqrt(n+1)-1);
        end
    end

    X(i,:) = x0 + s;
end
X(end,:) = x0; 

% evaluate objective function of the current simplex f(xi)
fobj = zeros(n+1,1);
for i=1:n+1
    fobj(i) = calc_obj(X(i,:),objparams);
end

k = 1;

[~,idx] = min(fobj);
xstar       = X(idx,:);
hist(k).xstar = xstar;
hist(k).x     = X;

while true

    % order from lowest(best) to highest f(xi)
    [fobj, idx] = sort(fobj);
    X           = X(idx,:);

    % calculate centroid excluding the worst point x_n+1
    xc = (1/n)* sum( X(1:n,:) );

    % perfrom reflection (alpha=1)
    xr = xc + (xc - X(n+1,:));

    % evaluate reflected point f(xr)
    f_xr = calc_obj(xr, objparams);

    if f_xr < fobj(1)
            % perform expansion (alpha=2)
            xe = xc + 2*(xc - X(n+1,:));
    
            % evaluate expanded point f(xe)
            f_xe = calc_obj(xe, objparams);
            if f_xe < fobj(1)
                % accept expansion point
                X(n+1,:)  = xe; 
                fobj(n+1) = f_xe;
            else
                % accept reflected point
                X(n+1,:)  = xr; 
                fobj(n+1) = f_xr;
            end

    elseif f_xr < fobj(n)      % is the reflected point better than the 2nd worst
            % accept reflected point
            X(n+1,:)  = xr;     
            fobj(n+1) = f_xr;
    else
            if f_xr > fobj(n+1) % is the reflected point worse than the worst
                % perform inside contraction
                xic = xc - 0.5*(xc - X(n+1,:));

                % evaluate expanded point f(xic)
                f_xic = calc_obj(xic, params);
                if f_xic < fobj(n+1)
                    % accept inside contraction
                    X(n+1,:)  = xic; 
                    fobj(n+1) = f_xic;
                else
                    % shrink (gamma=0.5)
                    for j=1:n
                        X(j+1,:)  = X(1,:) + 0.5*(X(j+1,:)-X(1,:));
                        fobj(j+1) = calc_obj(X(j+1,:), objparams);
                    end
                end

            else

                % perform outside contraction
                xoc = xc + 0.5*(xc - X(n+1,:));

                % evaluate expanded point f(xoc)
                f_xoc = calc_obj(xoc, params);
                if f_xoc < f_xr
                    % accept outside contraction
                    X(n+1,:)  = xoc; 
                    fobj(n+1) = f_xoc;
                else
                    % shrink (gamma=0.5)
                    for j=1:n
                        X(j+1,:)  = X(1,:) + 0.5*(X(j+1,:)-X(1,:));
                        fobj(j+1) = calc_obj(X(j+1,:), objparams);
                    end
                end

            end
                    
    end

    k = k+1;
    [fstar,idx]   = min(fobj);
    xstar         = X(idx,:);
    hist(k).xstar = xstar;
    hist(k).x     = X;


    % check convergence
    delta_x = sum( vecnorm(X(1:n,:) - X(n,:), 2, 2) );
    delta_f = sqrt( (sum((fobj - mean(fobj)).^2 ))/(n+1));
    if (delta_x < tau_x) && (delta_f < tau_f), break; end

end

