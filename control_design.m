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

q1 = 1;
q2 = 10; 
q3 = 100;

Q =[q1 0 0;
    0 q2 0;
    0 0 q3];
R = 10000;
K_lqr = dlqr(Phi,Gamma, Q,R);

poles_cont_5 = 5*poles_cont;

p1_obsv = exp(-poles_cont_5(1)*Ts);
p2_obsv = exp(poles_cont_5(2)*Ts);
p3_obsv = exp(poles_cont_5(3)*Ts);

p_obsv = [p1_obsv, p2_obsv, p3_obsv]';

L = place(Phi', C', p_obsv)';
eig(Phi - L*C);
% 
Phi_obsv = [Phi - Gamma*K_lqr Gamma*K_lqr;
            zeros(3,3) Phi-L*C];

% Phi_obsv = Phi - Gamma*K_lqr;
Gamma_obsv = [Gamma;
    zeros(size(Gamma))];

C_obsv = [1 0 0; 0 0 1];
C_obsv = eye(6);

sysObserver = ss(Phi_obsv, Gamma_obsv, C_obsv,[], Ts);

t = 0:Ts:5;
u = zeros(size(t));
y = lsim(sysObserver, u, t, [0.1 0.02 0.0 0.1 0.1 0.1]');

figure(1)
stairs(t, y(:,1))
hold on
stairs(t, y(:,2))
stairs(t, y(:,3))
title('states')
legend('theta', 'theta dot', 'phi dot')
figure(2)
stairs(t, y(:,4))
hold on
stairs(t, y(:,5))
stairs(t, y(:,6))
title('e')