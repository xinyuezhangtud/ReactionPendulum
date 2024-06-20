clear all
close all
clc
%% 
load('reaction_pendulum.mat')

h=0.05;
Ts = 0.05;

G_disc = ss(Phi, Gamma, C, D, Ts);
G_cont = ss(A,B,C,D);
poles_cont = pole(G_cont);
%% downward equilibrium

% Define weight matrices
q1 = 1e4;
q2 = 1e2; 
q3 = 1e-3;

%% upwards equilibrium

% Weight on the states
q1 = 4.0; 
q2 = 3.5; 
q3 = 0.0071; 

Q =[q1 0 0;
    0 q2 0;
    0 0 q3];

% Weight on the control action
R = 1e2;

%% Computing the LQR control gain and Observer gain

%Computing the LQR control gain 
K_lqr = dlqr(Phi,Gamma, Q,R);

% Computing poles of the closed loop system
p_obsv_disc = eig(Phi - Gamma*K_lqr);
sys_cl = ss(Phi-Gamma*K_lqr, Gamma, G_disc.C, G_disc.D, Ts);

% Making them continuous time
sys_cl_cont= d2c(sys_cl, 'zoh');

% Making them 6 times faster
poles_cl_cont= pole(sys_cl_cont);
poles_cl_cont_5 = 6*(poles_cl_cont);

p1_obsv_dc = exp(poles_cl_cont_5(1)*Ts);
p2_obsv_dc = exp(poles_cl_cont_5(2)*Ts);
p3_obsv_dc = exp(poles_cl_cont_5(3)*Ts);

p_obsv_dc = [p1_obsv_dc;p2_obsv_dc;p3_obsv_dc];

% Computing the Observer gain
L = place(Phi', C', p_obsv_dc)';


