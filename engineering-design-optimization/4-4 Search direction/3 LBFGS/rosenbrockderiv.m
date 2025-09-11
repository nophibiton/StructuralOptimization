function g = rosenbrockderiv(x, params)

g(1,1) = 2*x(1) - 400*x(1)*(x(2)-x(1)^2 ) - 2;
g(2,1) = 200*x(2) - 200*x(1)^2;

end