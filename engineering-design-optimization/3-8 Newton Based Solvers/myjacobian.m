function [J] = myjacobian(u)

J(1,2) = 1;
J(2,2) = 1;

J(1,1) = 1/u(1)^2 ;
J(2,1) = -1/(2*sqrt(u(1)) );

end