clear all
close all
clc

m_p     = 0.117;
m_d     = 0.019;
m_tot   = m_p + m_d;
l_cg    = 0.0987;
l       = 0.14298;
g       = 9.81;
k       = 0.06;
r       = 0.1;
% J_th    = (1/12) * m_p *l^2;
% J_phi   = (1/2) * m_d * r^2 + m_d * (l^2/4);
J_th    = 6.2533e-4;
J_phi   = 9.4559e-4

A = [0 1 0 0;
    -(m_tot*g*l_cg)/J_th 0 0 0;
    0 0 0 1;
    m_tot*g*l_cg/J_th 0 0 0];

B = [0;
    -k/J_th;
    0;
    k*(1/J_th + 1/J_phi)];

C = [1 0 0 0;
    0 0 1 0 ];

D = [0];

x_e = [pi 0 0 0];

disp('A');
disp(A);
disp('B')
disp(B);
disp('C')
disp(C);

disp('controllability');
disp(rank(ctrb(A, B)));

%% EXPORT A, B, C, D matrices
save reaction_pendulum.mat A B C D x_e

%% Discretize the system
Ts = 0.01;

G = ss(A,B,C,D);

G_dt = c2d(G, Ts, 'zoh');

[Phi, Gamma, C, D] = ssdata(G_dt);

disp('Phi');
disp(Phi);

disp('Gamma');
disp(Gamma);

disp('C');
disp(C);

disp('D');
disp(D);

disp('Controllability')
disp(rank(ctrb(Phi, Gamma)));

%% EXPORT Phi, Gamma, C, D matrices
save pendubot_dt.mat Phi Gamma C D x_e Ts

%% Simulation

x_0 = [1.2*pi 0 0 0];
t = 0:Ts:10;
u = 0*t;
y= lsim(G_dt, u, t, x_0);

figure(1)
plot(y)

figure(2)
step(G_dt)