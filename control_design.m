clear all
close all
clc

load('reaction_pendulum.mat')
%% LQR
h=0.05;
Ts = 0.05;
G_disc = ss(Phi, Gamma, C, D, Ts);
G_cont = ss(A,B,C,D);

% poles_disc = pole(G_disc);

poles_cont = pole(G_cont);

q1 = 1e1;
q2 = 1e1; 
q3 = 1e-4;
% q1 = 4.5;
% q2 = 3.2;
% q3 = 0.000031;

Q =[q1 0 0;
    0 q2 0;
    0 0 q3];
R = 1e2;
% R=100;
K_lqr = dlqr(Phi,Gamma, Q,R);
disp(K_lqr);
poles_cont_5 = 8*poles_cont;
poles_disc = pole(G_disc);

p1_obsv = exp(-poles_cont_5(1)*Ts);
p2_obsv = exp(poles_cont_5(2)*Ts);
p3_obsv = exp(poles_cont_5(3)*Ts);

p_obsv = [p1_obsv, p2_obsv, p3_obsv]';

% L = place(Phi', C', p_obsv)';
p_obsv_disc = eig(Phi - Gamma*K_lqr)*0.5;
L = place(Phi', C', p_obsv_disc)';
% 
Phi_obsv = [Phi - Gamma*K_lqr Gamma*K_lqr;
            zeros(3,3) Phi-L*C];

% Phi_obsv = Phi - Gamma*K_lqr;
Gamma_obsv = [Gamma;
    zeros(size(Gamma))];

C_obsv = [1 0 0; 0 0 1];
C_obsv = eye(6);

sysObserver = ss(Phi_obsv, Gamma_obsv, C_obsv,[], Ts);
sys_cl = ss(Phi-Gamma*K_lqr, Gamma, G_disc.C, G_disc.D, Ts);
t = 0:Ts:50;
u = zeros(size(t));
y = lsim(sysObserver, u, t, [0.1 0.02 0.0 0.1 0.1 0.1]');
y = lsim(sys_cl, u, t, [0.01 0 0]');
figure(1)
stairs(t, y(:,1))
hold on
stairs(t, y(:,2))
% stairs(t, y(:,3))
title('states')
legend('theta', 'phi dot')
% figure(2)
% stairs(t, y(:,4))
% hold on
% stairs(t, y(:,5))
% stairs(t, y(:,6))
% title('e')