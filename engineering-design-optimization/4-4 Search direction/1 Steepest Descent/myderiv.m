function g = myderiv(x,params)


beta = params.beta;

g(1,1) = 2*x(1);
g(2,1) = beta*2*x(2);

end