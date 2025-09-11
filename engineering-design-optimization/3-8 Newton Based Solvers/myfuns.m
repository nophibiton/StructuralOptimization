function [r] = myfuns(u)

r(1) = u(2) - 1./u(1);
r(2) = u(2) - sqrt(u(1) );

end