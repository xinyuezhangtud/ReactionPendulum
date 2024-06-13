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
%% downward equilibrium
% q1 = 1e4;
% q2 = 1e2; 
% q3 = 1e-3;
% 
% Q =[q1 0 0;
%     0 q2 0;
%     0 0 q3];
% R = 1e2;
%% upwards
% q1 = 1e5;
% q2 = 1e3; 
% q3 = 1e-6;
% Q =[q1 0 0;
%     0 q2 0;
%     0 0 q3];
% R = 1e1;

%% upwards
q1 =4.5;
q2 = 3.50; 
q3 = 0.0071;
Q =[q1 0 0;
    0 q2 0;
    0 0 q3];
R = 20;
%% 

K_lqr = dlqr(Phi,Gamma, Q,R);

% Pc = [-10, -4+2j, -4-2j];
% Pd = exp(Pc*Ts);
% K_lqr = place(Phi, Gamma, Pd);



disp(K_lqr);
poles_cont_5 = 10*poles_cont;
poles_disc = pole(G_disc);

p1_obsv = exp(-poles_cont_5(1)*Ts);
p2_obsv = exp(poles_cont_5(2)*Ts);
p3_obsv = exp(poles_cont_5(3)*Ts);

p_obsv = [p1_obsv, p2_obsv, p3_obsv]';

% L = place(Phi', C', p_obsv)';

p_obsv_disc = eig(Phi - Gamma*K_lqr)*1;
sys_cl = ss(Phi-Gamma*K_lqr, Gamma, G_disc.C, G_disc.D, Ts);
sys_cl_cont= d2c(sys_cl, 'zoh');


poles_cl_cont= pole(sys_cl_cont);

% poles_cl_cont_5 = 6*real(poles_cl_cont) +2*imag(poles_cl_cont)*i;

poles_cl_cont_5 = 2*(poles_cl_cont);

p1_obsv_dc = exp(poles_cl_cont_5(1)*Ts);
p2_obsv_dc = exp(poles_cl_cont_5(2)*Ts);
p3_obsv_dc = exp(poles_cl_cont_5(3)*Ts);

p_obsv_dc = [p1_obsv_dc;p2_obsv_dc;p3_obsv_dc];

L = place(Phi', C', p_obsv_dc)';
% 
Phi_obsv = [Phi - Gamma*K_lqr Gamma*K_lqr;
            zeros(3,3) Phi-L*C];

% Phi_obsv = Phi - Gamma*K_lqr;
Gamma_obsv = [Gamma;
    zeros(size(Gamma))];

C_obsv = eye(6);

sysObserver = ss(Phi_obsv, Gamma_obsv, C_obsv,[], Ts);
sys_cl = ss(Phi-Gamma*K_lqr, Gamma, G_disc.C, G_disc.D, Ts);
t = 0:Ts:50;
u = zeros(size(t));
y = lsim(sysObserver, u, t, [0.1 0.02 0.0 0.1 0.1 0.1]');
% y = lsim(sys_cl, u, t, [0.01 0 0]');
figure(1)
stairs(t, y(:,1))
hold on
stairs(t, y(:,2))
stairs(t, y(:,3))
title('states')
legend('theta', 'phi dot')
figure(2)
stairs(t, y(:,4))
hold on
stairs(t, y(:,5))
stairs(t, y(:,6))
title('e')