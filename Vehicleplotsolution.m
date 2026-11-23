clear; clc; close all;

%% Parameters
U   = 28;
Cf  = -1.232e5;
Cr  = -1.042e5;
m   = 1.67e3;
Iz  = 2.1e3;
a   = 0.99;
b   = 1.7;

A = [0 1 0 0;
     0 (Cf+Cr)/(m*U) -(Cf+Cr)/m (a*Cf-b*Cr)/(m*U);
     0 0 0 1;
     0 (a*Cf-b*Cr)/(Iz*U) -(a*Cf-b*Cr)/Iz ...
       (a^2*Cf+b^2*Cr)/(Iz*U)];

B = [0;
     -Cf/m;
     0;
     -a*Cf/Iz];

Bd = [0; 
      (a*Cf-b*Cr)/m-U^2;
      0;
      (a^2*Cf+b^2*Cr)/Iz];

%% Initial condition
y0 = [4; 1; 3; 2];

%% Simulation
Tspan = [0 10];

[t,y] = ode45(@(t,y) vehicle_dyn(t,y,A,B),Tspan,y0);

%% Recompute control for plotting
u = zeros(length(t),1);

for k=1:length(t)
    u(k) = feedback_law(y(k,:)');
end

%% Plots
figure;

subplot(2,2,1)
plot(t,y(:,1),'LineWidth',2.5)
grid on
xlabel('$t$','Interpreter','latex','FontSize',22)
ylabel('$e$','Interpreter','latex','FontSize',22)
set(gca,'FontSize',18)

subplot(2,2,2)
plot(t,y(:,2),'LineWidth',2.5)
grid on
xlabel('$t$','Interpreter','latex','FontSize',22)
ylabel('$\dot{e}$','Interpreter','latex','FontSize',22)
set(gca,'FontSize',18)

subplot(2,2,3)
plot(t,y(:,3),'LineWidth',2.5)
grid on
xlabel('$t$','Interpreter','latex','FontSize',22)
ylabel('$e_{\theta}$','Interpreter','latex','FontSize',22)
set(gca,'FontSize',18)

subplot(2,2,4)
plot(t,y(:,4),'LineWidth',2.5)
grid on
xlabel('$t$','Interpreter','latex','FontSize',22)
ylabel('$\dot{e}_{\theta}$','Interpreter','latex','FontSize',22)
set(gca,'FontSize',18)
%figure
%plot(t,u,'LineWidth',1.5)
%grid on
%xlabel('t (s)')
%ylabel('u(t)')
%title('Steering control')

%% ===== Local functions =====

function dy = vehicle_dyn(~,y,A,B)

    u = feedback_law(y);

    dy = A*y + B*u;

end

function u = feedback_law(y)

    % Gain K
    K = -[0.1166, 0.6243, 0.6490, 0.5674];

    % Produit Ky
    Ky = K * y;

    % Saturation sat_i(x) = max(-1, min(1, i*x))
    sat = @(x, i) max(-1, min(1, i*x));

    % Inner integral
    f1 = @(nu) sat((-nu/6) * Ky, 1/3);
    I1 = integral(f1, 0, 1, 'ArrayValued', true);

    % Middle integral
    f2 = @(tau) tau .* sat(tau .* I1, 2/3);
    I2 = integral(f2, 0, 1, 'ArrayValued', true);

    % Outer integral
    f3 = @(r) sat(r .* I2, 7/6);
    I3 = integral(f3, 0, 1, 'ArrayValued', true);

    % Control input
    u = -54 * I3;

end


% axes (IMPORTANT : applique aux ticks aussi)
%set(gca, ...
%    'FontSize',30, ...
 %   'FontName','Times')

% labels (LaTeX + taille plus grande)
%xlabel('$t$', ...
%    'Interpreter','latex', ...
 %   'FontSize',30, ...
 %   'FontName','Times')

%ylabel('$x(t)$', ...
%    'Interpreter','latex', ...
 %   'FontSize',30, ...
  %  'FontName','Times')

% légende (LaTeX + taille cohérente)
%hleg = legend({'$x_1$','$x_2$'}, ...
 %   'Interpreter','latex', ...
  %  'FontSize',30, ...
  %  'Location','best');

% légende (LaTeX + taille cohérente)
%hleg = legend({'continuous NN control','linear control'}, ...
 %   'FontSize',30, ...
 %   'Location','best');

% export EPS (comme dans le script)
set(gcf,'PaperPositionMode','auto')
print('-depsc2','-loose','figure_states.eps')