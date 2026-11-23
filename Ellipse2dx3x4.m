clear; clc; close all;

% Première matrice
Q1 = [0.1129 0.1131 0.0772 0.1108;
      0.1131 0.1149 0.0374 0.1113;
      0.0772 0.0374 1.6263 0.1139;
      0.1108 0.1113 0.1139 0.1189];

% Deuxième matrice
Q2 = [0.0021 0.0087 0.0106 0.0072;
      0.0087 0.0498 -0.0408 0.0360;
      0.0106 -0.0408 2.7933 0.1015;
      0.0072 0.0360 0.1015 0.0468];

% Blocs correspondant à (0,0,x3,x4)
A1 = Q1(3:4,3:4);
A2 = Q2(3:4,3:4);

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

plot(X1(1,:),X1(2,:), 'LineWidth',2, ...
     'DisplayName','Q_1');

plot(X2(1,:),X2(2,:), 'LineWidth',2, ...
     'DisplayName','Q_2');

axis equal;
grid on;

xlabel('x_3');
ylabel('x_4');

title('Sections (0,0,x_3,x_4) des ellipsoïdes');
legend('Location','best');