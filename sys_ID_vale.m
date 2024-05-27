close all
clear all
clc

load('/Users/mavi.ciabba/Desktop/VALE/Q4/Integration project/GIT_HUB/reaction-pendulum_new/chirp.mat')
% load('/Users/mavi.ciabba/Desktop/VALE/Q4/Integration project/GIT_HUB/reaction-pendulum_new/no_sin_01amp.mat')

input = u.Data(500:800);
y = [theta.Data(500:800) phi_dot.Data(500:800)]; 

% figure(1)
% plot(y)
%% 
% m_d     = 0.4;
% m       = 0.5;
% l       = 0.7;
% b_theta = 0.55;% best with 0.08
% b_phi   = 0.7 ;
% g       = 9.81;
% r       = 0.05;
% k       = 0.01;
m_d     = 0.4;
m       = 0.5;
l       = 0.7;
b_theta = 0.055;% best with 0.08
b_phi   = 0.07 ;
g       = 9.81;
r       = 0.05;
k       = 0.01;
J_phi   = (1/2) * m_d * r^2 + m_d * ((l^2)/4);
J_th    = (1/12) * (m) *l^2;
Ts = 0.05; 
l_cg = 0.1;

data = iddata(y, input, Ts);


theta0 = [m*g*l_cg; J_th; b_theta; b_phi; k; J_phi];
% theta0 = [4.18; 0.38; 0.09; 0.05; 1.2; 0.0034];

% [A,B,C,D]= reaction_wheel_pendulum(theta0)

sys= idgrey(@reaction_wheel_pendulum, theta0, 'c');
sys.Structure.Parameters(1).Minimum = 0;

opt = greyestOptions('InitialState','estimate','Display','on');

estimatedModel = greyest(data, sys, opt);

estimatedParams = getpvec(estimatedModel);
disp('Estimated Parameters:');
disp(estimatedParams);

compare(data,estimatedModel);

function [A, B, C, D] = reaction_wheel_pendulum(theta, Ts)
    mgl = theta(1);
    J_theta = theta(2);
    k = theta(5);
    J_phi = theta(6);
    b_theta = theta(3);
    b_phi = theta(4);

    A = [0 1 0 0; 
        -mgl/J_theta -b_theta/J_theta 0 b_phi/J_theta; 
        0 0 0 1;
        mgl/J_theta b_theta/J_theta 0 -b_phi/J_theta];

    B = [0; -k/J_theta; 0; (k/J_theta + k/J_phi)];
    C = [1 0 0 0; 0 0 0 1];
    D = [0;0];
end