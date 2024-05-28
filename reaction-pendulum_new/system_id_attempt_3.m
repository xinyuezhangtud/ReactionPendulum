close all
clear all
clc

load('no_sin_01amp.mat')
input2 = u.Data(400:900);
y2 = [theta.Data(400:900) phi_dot.Data(400:900)]; 
data2 = iddata(y2, input2, 0.05);


load('no_chirp_015amp.mat')
input = u.Data(400:900);
y = [theta.Data(400:900) phi_dot.Data(400:900)]; 


m_d     = 0.2;
m       = 0.311; %0.311
l       = 0.9163; %0.9163
b_phi   = 0.011; 
g       = 9.81;
r       = 0.17;
k       = 0.07;
J_phi   = (1/2) * m_d * r^2 + m_d * ((l^2)/4);
J_th    = (1/12) * (m) *l^2;
Ts = 0.05; 

data = iddata(y, input, Ts);


theta0 = [m*g*l; J_th; J_phi;k];
% theta0 = [m*g*l; J_th; b_theta; b_phi];
%% 

sys= idgrey(@reaction_wheel_pendulum, theta0, 'c');
sys.Structure.Parameters(1).Minimum = 0;

opt = greyestOptions('InitialState','estimate','Display','on');

estimatedModel = greyest(data2, sys, opt);

estimatedParams = getpvec(estimatedModel);
disp('Estimated Parameters:');
disp(estimatedParams);

figure(1)
compare(data,estimatedModel);

figure(2)
compare(data2,estimatedModel);

function [A, B, C, D] = reaction_wheel_pendulum(theta, Ts)
    mgl = theta(1);
    J_theta = theta(2);
    k = theta(4);
    J_phi = theta(3);

    A = [0 1 0; 
        -mgl/J_theta 0 0;
        mgl/J_theta 0 0];
    B = [0; -k/J_theta; k/J_theta + k/J_phi];

    C = [1 0 0; 0 0 1];
    D = [0;0];
end
