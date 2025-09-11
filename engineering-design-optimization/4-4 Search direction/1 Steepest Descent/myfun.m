function y = myfun(x,params)

beta = params.beta;

y = x(1).^2 + beta*x(2).^2;

end