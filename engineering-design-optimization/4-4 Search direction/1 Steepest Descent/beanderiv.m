function g = beanderiv(x, params)

g(1,1) = 2*x(1) - 2*x(1)*(2*x(2) - x(1)^2) -2;
g(2,1) = -2 + 6*x(2) - 2*x(1)^2;

end