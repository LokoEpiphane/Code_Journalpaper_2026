% Region of attraction analysis on the pendulum system with NN controller
clear

%% parameters
% Nominal speed of the vehicle travels at.
U = 28; % m/s
% Model
% Front cornering stiffness for one wheel.
Ca1 = -61595; % unit: Newtons/rad
% Rear cornering stiffness for one wheel.
Ca3 = -52095; % unit: Newtons/rad

% Front cornering stiffness for two wheels.
Caf = Ca1*2; % unit: Newtons/rad
% Rear cornering stiffness for two wheels.
Car = Ca3*2; % unit: Newtons/rad

% Vehicle mass
m = 1670; % kg
% Moment of inertia
Iz = 2100; % kg/m^2

% Distance from vehicle CG to front axle
a = 0.99; % m
% Distance from vehicle CG to rear axle
b = 1.7; % m

g = 9.81;

% sampling period
dt = 0.02;

% Continuous-time state space matrices
% States are lateral displacement(e) and heading angle error(deltaPsi) and 
% their derivatives.
% Inputs are front wheel angle and curvature of the road.
Ac = [0 1 0 0; ...
      0, (Caf+Car)/(m*U), -(Caf+Car)/m, (a*Caf-b*Car)/(m*U); ...
      0 0 0 1; ...
      0, (a*Caf-b*Car)/(Iz*U), -(a*Caf-b*Car)/Iz, (a^2*Caf+b^2*Car)/(Iz*U)];
B1c = [0; 
      -Caf/m;
      0; ...
      -a*Caf/Iz];
% discribes how the curvature enters the system. 
B2c = [0; 
      (a*Caf-b*Car)/m-U^2;
      0; ...
      (a^2*Caf+b^2*Car)/Iz];

%% x^+ = AG*x + BG1*q + BG2*usat
AG = Ac*dt + eye(4);

% describes how q enters the system
BG1 = B1c*dt;
% describes how usat enters the system
BG2 = B1c*dt;
nG = size(AG, 1);
nu = 1;
nq = 1;

% the input to Delta: p = usat

%% load weights and biases of the NN controller     
fname = './vehicle_Wb_s32_tanh/Wb_s32_tanh/';
load([fname 'W1.csv'])
load([fname 'W2.csv'])
load([fname 'W3.csv'])
n0 = nG;
n1 = size(W1,1);
n2 = size(W2,1);
n3 = size(W3,1);
nphi = n1+n2+n3;

b1 = zeros(n1,1);
b2 = zeros(n2,1);
b3 = zeros(n3,1);

%% bounds for the inputs to the nonlinearity
xeq = zeros(4,1);
v1eq = W1*xeq + b1;
w1eq = tanh(v1eq);
v2eq = W2*w1eq + b2;
w2eq = tanh(v2eq);
v3eq = W3*w2eq + b3; % This is also u_*
% usat = sat(v3) = sat(u)

deltav1 = 0.6; 
v1up = deltav1 + v1eq;
v1lb = v1eq - deltav1;

% Note that tanh(v1) is sector-bounded by [alpha1,beta] = [tanh(r1)/r1,1]
alpha1sec = min((tanh(v1up)-tanh(v1eq))./(v1up-v1eq), (tanh(v1eq)-tanh(v1lb))./(v1eq-v1lb));
syms r
grad_tanh = jacobian(tanh(r),r);
alpha1off = min(double(subs(grad_tanh,r,v1lb)), double(subs(grad_tanh,r,v1up)));
beta = 1;
w1up = tanh(v1up);
w1lb = tanh(v1lb);

v2up = W2*1/2*(w1up+w1lb) + b2 + abs(W2)*1/2*abs(w1up-w1lb);
v2lb = W2*1/2*(w1up+w1lb) + b2 - abs(W2)*1/2*abs(w1up-w1lb);
alpha2off = [];
for i = 1:n2
    alpha2off = [alpha2off; min(double(subs(grad_tanh,r,v2lb(i))), double(subs(grad_tanh,r,v2up(i))))];
end
alpha2sec = min((tanh(v2up)-tanh(v2eq))./(v2up-v2eq), (tanh(v2eq)-tanh(v2lb))./(v2eq-v2lb));
w2up = tanh(v2up);
w2lb = tanh(v2lb);

v3up = W3*1/2*(w2up+w2lb) + b3 + abs(W3)*1/2*abs(w2up-w2lb);
v3lb = W3*1/2*(w2up+w2lb) + b3 - abs(W3)*1/2*abs(w2up-w2lb);
umax = 30/180*pi;
alpha3sec = min((sat(v3up,umax,-umax)-sat(v3eq,umax,-umax))./(v3up-v3eq), (sat(v3eq,umax,-umax)-sat(v3lb,umax,-umax))./(v3eq-v3lb));
if v3up > umax || v3lb < -umax
    alpha3off = 0;
else
    alpha3off = 1;
end

Alphasec = blkdiag(diag(alpha1sec),diag(alpha2sec),alpha3sec);
Alphaoff = blkdiag(diag(alpha1off),diag(alpha2off),alpha3off);
Beta = beta*eye(nphi);

%% multiplier for LTI uncertainty

delta_norm = 0.1;

Apsi = -eye(2);
npsi = size(Apsi,1);
Apsi = Apsi*dt + eye(npsi);
Bpsi1 = [1; 0];
Bpsi1 = Bpsi1*dt;
Bpsi2 = [0; 1];
Bpsi2 = Bpsi2*dt;
Cpsi = [0,1,0,0;0,0,0,1]';
Dpsi1 = [1,0,0,0]';
Dpsi2 = [0,0,1,0]';

%% filter for the off-by-one IQC 
% xo^+ = Ao*xo + Bov*v + Bow*w
% zo^+ = Co*xi + Dov*v + Dow*w
% denote v as the input to the activation fcn, 
% denote w as the output to the activation fcn.
Ao = zeros(nphi);
Bov = -Beta;
Bow = eye(nphi);
Co = [eye(nphi); zeros(nphi)];
Dov = [Beta; -Alphaoff];
Dow = [-eye(nphi); eye(nphi)];

% filter for sector IQC, notice that Dsecv = Dov, Dsecw = Dow
Csec = [zeros(nphi); zeros(nphi)];
Dsecv = [Beta; -Alphasec];
Dsecw = [-eye(nphi); eye(nphi)];

%% define X to gather the states x and xo, X = [x; xo; psi]
% X^+ = AX + B[u_sat;q;v;w]
A = blkdiag(AG, Ao, Apsi);
B = [BG2, BG1, zeros(nG,nphi), zeros(nG,nphi);...
     zeros(nphi,nu), zeros(nphi,nq), Bov, Bow;...
     Bpsi1, Bpsi2, zeros(npsi,nphi), zeros(npsi,nphi)];  
 
Rx = [eye(nG), zeros(nG, nphi+npsi+nphi+nq)];
Rxoff = [zeros(nphi,nG),eye(nphi),zeros(nphi,npsi+nphi+nq)];
Rpsi = [zeros(npsi,nG+nphi),eye(npsi),zeros(npsi,nphi+nq)];
Rusat = [zeros(nu,nG+nphi+npsi+n1+n2),eye(nu),zeros(nu,nq)];
Rq = [zeros(nq,nG+nphi+npsi+nphi),eye(nq)];
Rvphi = [blkdiag([W1,zeros(n1,nphi+npsi)], W2, [W3,zeros(n3,nu)]), zeros(nphi,nq)];
Rwphi = [zeros(nphi,nG+nphi+npsi),eye(nphi),zeros(nphi,nq)];

%% Convex Optimization - compute ROA (YALMIP)

P       = sdpvar(nG+nphi+npsi, nG+nphi+npsi, 'symmetric');
lamoff  = sdpvar(nphi,nphi,'diagonal');
lamsec  = sdpvar(nphi,nphi,'diagonal');
M11     = sdpvar(2,2,'symmetric');

Constraints = [];

%% PSD constraints
Constraints = [Constraints, P >= 1e-8*eye(nG+nphi+npsi)];
Constraints = [Constraints, M11 >= 0];
Constraints = [Constraints, diag(lamoff) >= 0];
Constraints = [Constraints, diag(lamsec) >= 0];

%% Uncertainty multiplier
M = blkdiag(delta_norm*M11, -1/delta_norm*M11);

%% Off-by-one IQC
Moff = [zeros(nphi), lamoff;
        lamoff, zeros(nphi)];

Rmido = [Co, Dov, Dow] * [Rxoff; Rvphi; Rwphi];

%% Sector IQC
Msec = [zeros(nphi), lamsec;
        lamsec, zeros(nphi)];

Rmidsec = [Csec, Dsecv, Dsecw] * [Rxoff; Rvphi; Rwphi];

Qphi = Rmido' * Moff * Rmido + ...
       Rmidsec' * Msec * Rmidsec;

%% IQC term
Qpi = [Rpsi; Rusat; Rq]' * ...
      [Cpsi, Dpsi1, Dpsi2]' * M * ...
      [Cpsi, Dpsi1, Dpsi2] * ...
      [Rpsi; Rusat; Rq];

%% Lyapunov difference
S = [A'*P*A - P, A'*P*B;
     B'*P*A,     B'*P*B];

RV = [Rx; Rxoff; Rpsi; Rusat; Rq; Rvphi; Rwphi];

QV = RV' * S * RV;

%% Main LMI
Constraints = [Constraints, ...
    QV + Qphi + Qpi <= 1e-8*eye(nG+nphi+npsi+nphi+nq)];

%% Inclusion constraints
for i = 1:n1

    LMI_i = [deltav1^2, W1(i,:), zeros(1,nphi+npsi);
             [W1(i,:), zeros(1,nphi+npsi)]', P];

    Constraints = [Constraints, LMI_i >= 0];

end

%% Objective
Objective = trace(P(1:nG,1:nG));

%% Solve with SDPT3
ops = sdpsettings('solver','sdpt3','verbose',1);
sol = optimize(Constraints, Objective, ops);

%% Save data
filename = ['vehicle_offby1andSector_LTIunc',datestr(now, 'dd-mmm-yyyy-HH-MM-PM'),'.mat'];
save(filename)

%%
if sol.problem == 0

    Pval = value(P);

    P11 = Pval(1:nG,1:nG);

    radii   = 1./sqrt(eig(P11));
    traceP  = trace(P11);
    volume  = 1/det(P11);

    disp('Optimization successful');
    disp('P = ');
    disp(Pval);

else
    disp('Optimization failed');
    disp(sol.info);

end


% pvar x1 x2 x3 x4
% x = [x1;x2;x3;x4];
% V = x'*P(1:nG,1:nG)*x;
% 
% subplot(1,2,1)
% domain1 = [-36, 36, -10, 10];
% pcontour(subs(V,[x3;x4],[0;0]),1,domain1,'b')
% xlim(domain1(1,1:2))
% ylim(domain1(1,3:4))
% hold on
% 
% subplot(1,2,2)
% domain2 = [-0.6, 0.6, -6, 6];
% pcontour(subs(V,[x1;x2],[0;0]),1,domain2,'b')
% hold on
% xlim(domain2(1,1:2))
% ylim(domain2(1,3:4))

Pval = value(P);

fprintf('Bloc 4x4 de P :\n');
disp(Pval(1:4,1:4));

fprintf('La matrice K :\n');
disp(W3*W2*W1);