function y = myfun(x,params)

y = 0.1*x(1)^6 - 1.5*x(1)^4 + 5*x(1)^2 + 0.1*x(2)^4 + 3*x(2)^2 ...
            - 9*x(2) + 0.5*x(1)*x(2);

end