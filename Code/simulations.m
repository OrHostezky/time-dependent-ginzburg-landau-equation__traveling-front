%%% Or Hostezky
%%% 1D-TDGLE with a polynomial reaction function f(u) = -u(u - 1)(u - u3)
clear; clc;
xmesh = -100:0.02:100; % Spatial points for evaluation

% Different initial conditions
u0x_1 = @(x,u3) u3/pi * (pi/2 - atan(x));
u0x_2 = @(x,u3) u3/2 * (1 - sign(x).*min(abs(x),1));
u0x_3 = @(x,u3) u3/2 * (1 - erf(x));
u0x_4 = @(x,u3) u3/2 * (1 - erf(30*x));

% Analytical asymptotic waveform
wave_fun = @(xmesh,u3) u3 ./ (1 + exp(u3/sqrt(2)*xmesh)); 

%%% NOTE: Under the above assignments, each of the simulation blocks can
%%% run independently (1st block can run without)


%% 1st block: Plot reaction function f(u)
u3 = 2.5; % Reaction parameter
u = 0:0.01:u3;
f = -u.*(u-1).*(u-u3);

figure 
plot(u,f,'DisplayName','$f(u)$')
hold on
plot([0 1 u3],[0 0 0],'o','DisplayName','Fixed points')
plot([0 u3],[0 0],'-k','Linewidth',1.5,'HandleVisibility','off')
xlim([0 u3])
xlabel('$u$','Interpreter','latex','FontSize',12)
ylabel('$f(u)$','Interpreter','latex','FontSize',12)
legend('Location','northwest','Interpreter','latex','FontSize',9)

saveas(gcf,'../Data_&_Plots/Figure_1.fig')
saveas(gcf,'../Data_&_Plots/Figure_1.png')


%% 2nd block: Specific case - waveform + evaluating c (c > 0)
u3 = 2 + sqrt(2);  % Reaction parameter
tspan = 0:0.02:2; % Temporal span
u0x = u0x_4;       % Initial condition

% Solve PDE
u_tx = pdepe(0,@(x,t,u,dudx) pdefun(x,t,u,dudx,u3),@(x) u0x(x,u3), ...
    @(a,ua,b,ub,t) bcfun(a,ua,b,ub,t,u3),xmesh,tspan);

% Compute x values of the front, fit for c and xi0
t0 = 2; % Ignore transient in fit
[c,xi0,xfront,~] = fit_front(u_tx,xmesh,tspan,t0,u3);

% Plot
tidx = 1:10:101; % Sample every 10 time-steps
cmap = parula(length(tidx));
figure

% Traveling wavefront
subplot(2,2,[1 2])
for i = 1:length(tidx)
    plot(xmesh,u_tx(tidx(i),:),'.','Color',cmap(i,:), ...
        'DisplayName',['$t = ' num2str(tspan(tidx(i))) '$'])
    hold on
end
plot(xfront(tidx),u3/2*ones(1,length(tidx)),'ok','DisplayName','$u_t(x) = u_3/2$')
box on
xlim([-1 4])
ylim([0 u3])
xlabel('$x$','Interpreter','latex','FontSize',12)
ylabel('$u_t(x)$','Interpreter','latex','FontSize',12)
legend('Interpreter','latex','FontSize',9)

% xfront vs time
subplot(2,2,3)
plot(tspan,xfront,'.','DisplayName','Solution')
hold on
plot(tspan,c*tspan+xi0,'-k','DisplayName','Linear fit')
xlim([0 5])
ylim([0 5*c+xi0])
xlabel('$t$','Interpreter','latex','FontSize',12)
ylabel('$x_{front}(t)$','Interpreter','latex','FontSize',12)
legend('Location','northwest','Interpreter','latex','FontSize',9)

% Waveform at different times
subplot(2,2,4)
for i = 1:length(tidx)
    plot(xmesh-c*tspan(tidx(i))-xi0,u_tx(tidx(i),:),'.','Color',cmap(i,:), ...
        'HandleVisibility','off')
    hold on
end
plot(-2.5:0.05:2.5,wave_fun(-2.5:0.05:2.5,u3),'--k','DisplayName','Analytical')
xlim([-2.5 2.5])
ylim([0 u3])
xlabel('$\xi_t - \xi_0$','Interpreter','latex','FontSize',12)
ylabel('$u(\xi_t - \xi_0)$','Interpreter','latex','FontSize',12)
legend('Interpreter','latex','FontSize',9)

% Save
saveas(gcf,'../Data_&_Plots/Figure_2.fig')
saveas(gcf,'../Data_&_Plots/Figure_2.png')
save('../Data_&_Plots/Figure_2.mat')


%% 3rd block: Invariance of solution at t -> inf to initial condition
tspan = 0:0.5:5;                 % Temporal span
u0x = {u0x_1 u0x_2 u0x_3 u0x_4}; % Initial conditions

% Solve PDE for each initial condition, compute xi0 values
u_tx = nan(length(tspan),length(xmesh),4);
xi0 = nan(1,4);
for i = 1:4
    u_tx(:,:,i) = pdepe(0,@(x,t,u,dudx) pdefun(x,t,u,dudx,u3),@(x) u0x{i}(x,u3), ...
        @(a,ua,b,ub,t) bcfun(a,ua,b,ub,t,u3),xmesh,tspan);
    t0 = 2; % Ignore transient in fit
    [c,xi0(i),~,~] = fit_front(u_tx(:,:,i),xmesh,tspan,t0,u3);
end

% Plot
cmap = parula(length(u0x));
figure

% Initial conditions
subplot(2,1,1)
plot(xmesh,u_tx(1,:,1),'.','Color',cmap(1,:))
hold on
plot(xmesh,u_tx(1,:,2),'.','Color',cmap(2,:))
plot(xmesh,u_tx(1,:,3),'.','Color',cmap(3,:))
plot(-5:5e-4:5,u0x{4}(-5:5e-4:5,u3),'.','Color',cmap(4,:))
plot(-5:0.05:5,wave_fun(-5:0.05:5,u3),'--k')
xlim([-5 5])
ylim([0 u3])
xlabel('$x$','Interpreter','latex','FontSize',12)
ylabel('$u_0(x)$','Interpreter','latex','FontSize',12)

% Solutions at t = 5
subplot(2,1,2)
plot(xmesh-c*tspan(end)-xi0(1),u_tx(length(tspan),:,1),'.','Color', ...
    cmap(1,:),'DisplayName','$u_0(x) = u_3/\pi \big(\pi/2 - \textrm{tan}^{-1}(x)\big)$')
hold on
plot(xmesh-c*tspan(end)-xi0(2),u_tx(length(tspan),:,2),'.','Color', ...
    cmap(2,:),'DisplayName','$u_0(x) = u_3/2 \big(1 - \textrm{sign}(x)\textrm{min}\{1,|x|\} \big)$')
plot(xmesh-c*tspan(end)-xi0(3),u_tx(length(tspan),:,3),'.','Color', ...
    cmap(3,:),'DisplayName','$u_0(x) = u_3/2 \big(1 - \textrm{erf}(x)\big)$')
plot(xmesh-c*tspan(end)-xi0(4),u_tx(length(tspan),:,4),'.','Color', ...
    cmap(4,:),'DisplayName','$u_0(x) = u_3/2 \big(1 - \textrm{erf}(30x)\big)$')
plot(-2.5:0.05:2.5,wave_fun(-2.5:0.05:2.5,u3),'--k','DisplayName', ...
    'Analytical asymptotic wavefront')
xlim([-2.5 2.5])
ylim([0 u3])
xlabel('$\xi_{t = 5} - \xi_0$','Interpreter','latex','FontSize',12)
ylabel('$u(\xi_{t = 5} - \xi_0)$','Interpreter','latex','FontSize',12)
legend('Interpreter','latex','FontSize',9)

% Save
saveas(gcf,'../Data_&_Plots/Figure_3.fig')
saveas(gcf,'../Data_&_Plots/Figure_3.png')
save('../Data_&_Plots/Figure_3.mat')


%% 4th block:  Verification of the 'area rule' and convergence
u3 = 1.005:0.005:3;                     % Reaction parameter
tspan1 = 0:0.01:20; tspan2 = 0:0.01:10; % Temporal spans (near / far u3 = 1)
t0_1 = 10; t0_2 = 5;                    % Corresponding times from which to fit
idx = round(length(u3)/4);              % Boundary between computation segments

% Initial condition (avoid broadcast u0x)
u0x = cell(1, length(u3) - length(idx));
u0x(:) = {u0x_3};

% Solve PDE and compute c and tau at each u3 value
c = nan(1,length(u3));
tau_smpl = nan(1,length(u3));
parfor i = 1:idx
    u_tx = pdepe(0,@(x,t,u,dudx) pdefun(x,t,u,dudx,u3(i)),@(x) u0x{i}(x,u3(i)), ...
        @(a,ua,b,ub,t) bcfun(a,ua,b,ub,t,u3(i)),xmesh,tspan1);
    [c(i),~,~,tau] = fit_front(u_tx,xmesh,tspan1,t0_1,u3(i));
    tau_smpl(i) = tau;
end
parfor i = idx+1:length(u3)
    u_tx = pdepe(0,@(x,t,u,dudx) pdefun(x,t,u,dudx,u3(i)),@(x) u0x{i}(x,u3(i)), ...
        @(a,ua,b,ub,t) bcfun(a,ua,b,ub,t,u3(i)),xmesh,tspan2);
    [c(i),~,~,tau] = fit_front(u_tx,xmesh,tspan2,t0_2,u3(i));
    tau_smpl(i) = tau;
end

% Average tau in batches (several u3 values per point)
npts = 20;                                           % No. of batches
du3 = (max(u3) - min(u3)) / npts;                    % Batch size
u3pts = u3(1) + du3/2 : du3 : u3(end) - du3/2;       % Batch locations
idx = 1;                                             % Batch no.
tau = nan(1,npts); dtau = tau;                       % Convergence time and error
w = exp(-8*(((1:npts) - (npts+1)/2) / (npts-1)).^2); % Weights

for i = 1 : round(length(u3)/npts) : round(length(u3) - npts + 1)
    tau_batch = tau_smpl(i : i + round(length(u3)/npts) - 1); % Batch
    conv = tau_batch>=0; % Converged cases
    tau(idx) = sum(w(conv).*tau_batch(conv)) / round(length(u3)/npts); % Weighted mean
    dtau(idx) = std(tau_batch(conv),w(conv)); % Weighted STD
    idx = idx + 1;
end

% Plot
figure

% Velocities
subplot(1,2,1)
plot(u3,c,'.','Color','#4DBEEE','DisplayName','Numerical')
hold on
plot(u3,(u3-2)/sqrt(2),'--k','DisplayName','Analytical')
xlim([1 3])
ylim([-inf inf])
xlabel('$u_3$','Interpreter','latex','FontSize',12)
ylabel('$c(u_3)$','Interpreter','latex','FontSize',12)
legend('Location','northwest','Interpreter','latex','FontSize',9)

% Convergence times
subplot(1,2,2)
errorbar(u3pts,tau,dtau,'.')
xlim([1 3])
ylim([0 inf])
xlabel('$u_3$','Interpreter','latex','FontSize',12)
ylabel('$\tau(u_3)$','Interpreter','latex','FontSize',12)

% Save
saveas(gcf,'../Data_&_Plots/Figure_4.fig')
saveas(gcf,'../Data_&_Plots/Figure_4.png')
save('../Data_&_Plots/Figure_4.mat')

