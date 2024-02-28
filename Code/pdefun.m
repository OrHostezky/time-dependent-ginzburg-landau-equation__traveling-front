function [c,s,f] = pdefun(~,~,u,dudx,u3)
% Defining the 1D-TDGLE PDE: c*dudt = dudx + f(u)
c = 1;
s = dudx;
f = -u*(u-1)*(u-u3); % Reaction function f(u)