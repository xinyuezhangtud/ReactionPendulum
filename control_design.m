clear all
close all
clc

load('reaction_pendulum.mat')
%% LQR

Ts = 0.05;
G_disc = ss(Phi, Gamma, C, D, Ts);
G_cont = ss(A,B,C,D);

% poles_disc = pole(G_disc);

poles_cont = pole(G_cont);

q1 = 100;
q2 = 100; 
q3 = 100;

Q =[q1 0 0;
    0 q2 0;
    0 0 q3];
R = 10000;
K_lqr = dlqr(Phi,Gamma, Q,R);

poles_cont_5 = 10*poles_cont;

p1_obsv = exp(poles_cont_5(1)*Ts);
p2_obsv = exp(-poles_cont_5(2)*Ts);
p3_obsv = exp(poles_cont_5(3)*Ts);

p_obsv = [p1_obsv,p2_obsv,p3_obsv]';


L = place(Phi', C', p_obsv)';
eig(Phi - L*C);

L_2 = [0 0; 
     L(2,1), L(2,2);
     0 0];

eigPhi_LC = eig(Phi - L*C);
eigPhi = eig(Phi);
eigPhi_LC_x2 = eig(Phi - L_2*C);


Phi_obsv = [Phi - Gamma*K_lqr Gamma*K_lqr;
            zeros(3,3) Phi-L*C];

% Phi_obsv = Phi - Gamma*K_lqr;
Gamma_obsv = [Gamma;
    zeros(size(Gamma))];

C_obsv = [1 0 0; 0 0 1];
% C_obsv = eye(6);

% sysObserver = ss(Phi_obsv, Gamma_obsv, C_obsv,[], Ts);
cl_sys_2 = ss(Phi-L_2*C, Gamma,C_obsv,[],Ts);
cl_sys = ss(Phi-L*C, Gamma,C_obsv,[],Ts);


%%
t = 0:Ts:5;
u = zeros(size(t));
[yObs, tOut, xObs] = lsim(cl_sys, u, t, [0.1 0.0 0.0]');
% [y, tOut, x] = lsim(G_disc, u, t, [0.1 0.0 0.0]');
[y2, tOut, x2] = lsim(cl_sys_2, u, t, [0.1 0.0 0.0]');
figure(10)
stairs(y2)
hold on
stairs(yObs)
legend('y2_1','y2_2','y_1','y_2')
% 
% figure(1)
% stairs(t, xObs(:,1))
% hold on
% stairs(t, x(:,1))
% stairs(t, x2(:,1))
% 
% title('theta')
% legend('full OBS', 'No OBS', ' Obsx2')
% 
% figure(2)
% stairs(t, xObs(:,2))
% hold on
% stairs(t, x(:,2))
% stairs(t, x2(:,2))
% 
% title('theta dot')
% legend('full OBS', 'No OBS', ' Obsx2')
% 
% figure(3)
% stairs(t, xObs(:,3))
% hold on
% stairs(t, x(:,3))
% stairs(t, x2(:,3))
% 
% title('phi dot')
% legend('full OBS', 'No OBS', ' Obsx2')
% 
% % figure(2)
% % stairs(t, y(:,4))
% % hold on
% % stairs(t, y(:,5))
% % stairs(t, y(:,6))
% title('e')