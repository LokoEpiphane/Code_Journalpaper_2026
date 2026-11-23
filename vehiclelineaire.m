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

%% Feedback gain
K = -0.11*[0.1166, 0.6243, 0.6490, 0.5674];

%% Initial condition
y0 = [4; 1; 3; 2];

%% Simulation
Tspan = [0 10];

[t,y] = ode45(@(t,y) vehicle_dyn(t,y,A,B,K),Tspan,y0);

%% Control history
u = zeros(length(t),1);

for k = 1:length(t)
    u(k) = K*y(k,:)';
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

%% Control plot
%figure
%plot(t,u,'LineWidth',1.5)
%grid on
%xlabel('t (s)')
%ylabel('u(t)')
%title('Control input')

%% Closed-loop poles
Acl = A + B*K;

disp('Closed-loop eigenvalues:')
disp(eig(Acl))

%% Local function
function dy = vehicle_dyn(~,y,A,B,K)

u = K*y;

dy = A*y + B*u;

end

