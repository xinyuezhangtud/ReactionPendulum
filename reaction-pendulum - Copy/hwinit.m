% This script is executed every time that an experiment is initialized
%     clear all
%     close all
%     clc
% Define the sampling time here:
h = 0.05;

a21 = -9.9834;
a22 = -9.9834;
a24 = 0.0108;
a41 = 11.4724;
a42 = -0.0378;
a44 = -0.9255;
b2  = -3.0625;
b4  = 348.4161;
Ts = 0.05;

A = [0 1 0; 
     a21 a22 a24; 
     a41 a42 a44];
B = [0; b2; b4];
C = [1 0 0;
    0 0 1];
D = [0;0];

% Discretize system using 'zoh'
sys_CT = ss(A,B,C,D);
sys_DT = c2d(sys_CT, Ts, 'zoh');

Ts = 0.05;
G_disc = ss(sys_DT.A, sys_DT.B, sys_DT.C, sys_DT.D, Ts);
G_cont = ss(A,B,C,D);

% poles_disc = pole(G_disc);

poles_cont = pole(G_cont);

q1 = 1e1;
q2 = 1; 
q3 = 0.07;

Q =[q1 0 0;
    0 q2 0;
    0 0 q3];
R = 0.1;
K = dlqr(sys_DT.A,sys_DT.B, Q,R)

poles_cont_5 = 0.4*poles_cont;

p1_obsv = exp(poles_cont_5(1)*10*Ts);
p2_obsv = exp(poles_cont_5(2)*10*Ts);
p3_obsv = exp(poles_cont_5(3)*10*Ts);

p_obsv = [p1_obsv, p2_obsv,p3_obsv]'
% p_obsv =eig(sys_DT.A-sys_DT.B*K)
% p_obsv1 = [p_obsv(1);p_obsv(2);p_obsv(3)*100]
% p_obsv = [0.51 0.4 0.7]
L = place(sys_DT.A', sys_DT.C', p_obsv1)';

% Phi_obsv = [Phi - Gamma*K_lqr Gamma*K_lqr;
%             zeros(3,3) Phi-L*C];
% % Phi_obsv = Phi - Gamma*K_lqr;
% Gamma_obsv = [Gamma;
%     zeros(size(Gamma))];

% % sysObserver = ss(Phi_obsv, Gamma_obsv, C_obsv,[], Ts);
% cl_sys_2 = ss(Phi-L_2*C, Gamma,C_obsv,[],Ts);

cl_sys = ss(sys_DT.A-sys_DT.B*K, sys_DT.B,sys_DT.C,[],Ts);
OBS  = ss(sys_DT.A-L*sys_DT.C, sys_DT.B,sys_DT.C,[],Ts);

t = 0:Ts:5;
u = zeros(size(t));
[y, tOut, xObs] = lsim(cl_sys, u, t, [0.01 0.0 0.0]');
[yObs, tOut, xObs] = lsim(OBS, u, t, [0.1 0.02 0.01]');

% figure(10)
% stairs(yObs)
% hold on
% % stairs(yObs)
% legend('y2_1','y2_2','y_1','y_2')

% figure(2)
% plot(theta_obs)
% hold on
% plot(theta_meas)
% legend('observer','measured')