clear all
close all
clc

%% System identification
m_p     = 0.117;
m_d     = 0.019;
m_tot   = m_p + m_d;
l_cg    = 0.0987;
l       = 0.14298;
g       = 9.81;
k       = 0.06;
r       = 0.1;
J_th    = (1/12) * m_p *l^2;
J_phi   = (1/2) * m_d * r^2 + m_d * (l^2/4);
J_th    = 6.2533e-4;
J_phi   = 9.4559e-4;
%%
a = m_tot*g*l_cg;
b = J_th;
c = J_phi;
d = k;
e = -0.62;

fun1 = @r_pendulum;
par = [a b c d e];
aux = {};

T = 0;

[A,B,C,D,K,x0] = r_pendulum(par, T)

m = idgrey('r_pendulum', par,'c',aux,T);
% 
No_input_data = load('/Users/mavi.ciabba/Desktop/VALE/Q4/Integration project/GIT_HUB/reaction-pendulum_new/no_input.mat')
y = [No_input_data.theta.Data No_input_data.phi_dot.Data];
u = No_input_data.i_m.Data;
data = iddata(y,u,0.05);
m_est = greyest(data,m);
