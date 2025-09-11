function grad = myderiv(x,params)

grad(1,1) = 0.6*x(1)^5 - 6*x(1)^3 + 10*x(1) + 0.5*x(2);
grad(2,1) = 0.4*x(2)^3 + 6*x(2) - 9 + 0.5*x(1);


end