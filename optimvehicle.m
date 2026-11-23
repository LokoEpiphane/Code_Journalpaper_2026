clc;
%% Données du problème
U = 28;
C1 = -1.2321e+3;
C2 = -1.0421e+3;
m = 1.6710;
Iz = 2.110;
a = 0.99;
b = 1.7;
A = [0 1 0 0;
0 (C1+C2)/(m*U) -(C1+C2)/(m) (a*C1-b*C2)/(m*U);
0 0 0 1;
0 (a*C1-b*C2)/(Iz*U) -(a*C1-b*C2)/(Iz) (a*a*C1+b*b*C2)/(Iz*U)];
B = [ 0;
-C1/m;
0;
-(a*C2)/Iz];
Bd=[ 0;
(a*C1-b*C2)/m-U*U;
0;
(a*a*C1+b*b*C2)/Iz];

alpha0 =1.1;
alpha1 = 0.6;
alpha2 = 1/6;
rho0 = 0.9;
rho1 = 1.4;
rho2 = 1/4;
K2 = -[0.1166 0.6243 0.6490 0.5674]; % vecteur ligne
%% Grilles de recherche
lambda1_values = linspace(310,320,20);
lambda2_values = linspace(299,310,20);
tau_grid = linspace(0,1,10);
lambda0_values = linspace(210,220,10);
%% stockage des résultats
best_trace = inf;
best_Q = [];
best_lambda1 = NaN;
best_lambda2 = NaN;
best_lambda0 = NaN;
best_gamma = inf;
%% boucle principale double
for lambda1 = lambda1_values
   for lambda2 = lambda2_values
       for lambda0 = lambda0_values

         % variable SDP
           Q = sdpvar(4,4,'symmetric');
         % blocs de la LMI
           M11 = (A + alpha0*alpha1*alpha2*B*K2)'*Q + Q*(A + alpha0*alpha1*alpha2*B*K2);
           M12 = Q*B + lambda0*rho0*alpha1*alpha2*K2';
           M13 = alpha0*Q*B+rho1*alpha2*lambda1*K2';
           M14 = alpha0*alpha1*Q*B + lambda2*rho2*K2';
           M21 = M12';
           M22 = -2*lambda0;
           M23 = lambda0*rho0;
           M24 = rho0*alpha1*lambda0;
           M31 = M13';
           M32 = M23';
           M33 = -2*lambda1;
           M34 = rho1*lambda1;

           M41 = M14';
           M42 = M24';
           M43 = M34';
           M44 = -2*lambda2;
           M = [M11 M12 M13 M14;
                M21 M22 M23 M24;
                M31 M32 M33 M34;
                M41 M42 M43 M44];
        % contraintes
          Constraints = [];
          Constraints = [Constraints, Q >= 1e-5*eye(4)];
          Constraints = [Constraints, M <= -1e-8*eye(7)];
        % contraintes en tau
            for tau = tau_grid
                a_tau = -(tau/3)*K2;
                LMI_tau = [1 a_tau;
                            a_tau' Q];
                Constraints = [Constraints, LMI_tau >= 1e-8*eye(5)];
            end
        % objectif
          Objective = trace(Q);
        % résolution
          options = sdpsettings('solver','sdpt3','verbose',0);
          sol = optimize(Constraints,Objective,options);
          if sol.problem == 0
             Q_val = value(Q);
             tr = trace(Q_val);
             if tr < best_trace
                 best_trace = tr;
                 best_Q = Q_val;
                 best_lambda1 = lambda1;
                 best_lambda2 = lambda2;
                 best_lambda0 = lambda0;
             end
%if gamma < best_gamma
% best_gamma = gamma;
%end
         end
      end
  end
end
d=det(inv(best_Q));
%% résultats
disp('Meilleur résultat trouvé:')
disp(['lambda1 = ', num2str(best_lambda1)])
disp(['lambda2 = ', num2str(best_lambda2)])
disp(['lambda0 = ', num2str(best_lambda0)])
disp(['trace(Q) = ', num2str(best_trace)])
disp('Q = ')
disp(best_Q)
disp('Le determinant de linverse de Q=')
disp(d)