function [c,xi0,xfront,tau] = fit_front(u_tx,xmesh,tspan,t0,u3)
% (1) Computes x values of the traveling front xfront (where u = u3/2) for
% each time point in tspan, (2) fits for the velocity c and spatial shift
% xi0 between fit and origin at t = 0, such that u(xi = xi0) = u3/2 from
% time t0 > tspan(1) onwards, and (3) computes the time tau (relative to
% tspan(1)) of convergence to fit. If no convergence is detected during
% tspan, tau = -1
%%% NOTE: It is convenient to use tspan such that tspan(1) = 0

% Compute xfront
xfront = nan(1,length(tspan));
for i = 1:length(xfront)
    u_ti = u_tx(i,:);
    xfront(i) = xmesh(abs(u_ti-u3/2) == min(abs(u_ti-u3/2)));
end

% Fit xfront for c and xi0
idx = find(tspan >= t0,1); % Ignore times before t0
c = polyfit(tspan(idx:end),xfront(idx:end),1);
xi0 = c(2); % Spatial shift at t = 0
c = c(1);

% Compute convergence time
dxi = xfront - c*tspan - xi0; % Traveling coordinate shift xi - xi0
tau = tspan(find(abs(dxi) <= 10*min(abs(dxi)),1)) - tspan(1);
if isempty(tau), tau = -1; end % No convergence
