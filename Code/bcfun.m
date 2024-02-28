function [pa,qa,pb,qb] = bcfun(~,ua,~,ub,~,u3)
% Boundary conditions - single traveling front:
% u(x -> -inf,t) = u3 ,  u(x -> +inf) = 0
pa = ua - u3;
qa = 0;
pb = ub;
qb = 0;