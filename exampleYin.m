clear; clc; close all;

% Première matrice Q
Q1 = [0.1129 0.1131 0.0772 0.1108;
      0.1131 0.1149 0.0374 0.1113;
      0.0772 0.0374 1.6263 0.1139;
      0.1108 0.1113 0.1139 0.1189];

% Deuxième matrice Q
Q2 = [0.0021 0.0043 0.0133 0.0028;
      0.0043 0.0280 -0.0723 0.0110;
      0.0133 -0.0723 4.2856 0.1928;
      0.0028 0.0110 0.1928 0.0380];

% Sous-matrices correspondant à (x1,x2,0,0)
A1 = Q1(1:2,1:2);
A2 = Q2(1:2,1:2);

% Paramètre
t = linspace(0,2*pi,1000);

%% Ellipse 1
[V1,D1] = eig(A1);
Y1 = [cos(t)./sqrt(D1(1,1));
      sin(t)./sqrt(D1(2,2))];
X1 = V1*Y1;

%% Ellipse 2
[V2,D2] = eig(A2);
Y2 = [cos(t)./sqrt(D2(1,1));
      sin(t)./sqrt(D2(2,2))];
X2 = V2*Y2;

%% Tracé
figure; hold on;

plot(X1(1,:),X1(2,:), ...
    'LineWidth',2, ...
    'DisplayName','Q_1');

plot(X2(1,:),X2(2,:), ...
    'LineWidth',2, ...
    'DisplayName','Q_2');

axis equal;
grid on;
xlabel('x_1');
ylabel('x_2');
title('Sections (x_1,x_2,0,0) des ellipsoïdes x^TQx=1');
legend('Location','best');