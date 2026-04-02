function [gx,dgx] = eval_lsf(x,d, isGrad, lsffunc,lsfparams,distparams)


func = lsffunc;

% evaluate limit state function at x
gx   = func(x,d,lsfparams);

% calculate gradient using FFD
dgx = zeros(size(x));
if isGrad
    for i = 1:length(x)
        delta  = distparams(i,2)/200;
        dx     = x;
        dx(i)  = dx(i) + delta;
        gdx    = func(dx,d,lsfparams);
        dgx(i) = (gdx-gx)/delta;
    end
end

end