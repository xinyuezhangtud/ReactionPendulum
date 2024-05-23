close all
clear all
clc

load('no_chirp_015amp.mat')

input = u.Data(400:900);
y = [theta.Data(400:900) phi_dot.Data(400:900)]; 
 

m_d     = 0.16;
m       = 0.3572;
l       = 0.9;
b_theta = 0.08;
b_phi   = 0.1;
g       = 9.81;
r       = 0.09;
k       = 0.06;
J_phi   = (1/2) * m_d * r^2 + m_d * ((l^2)/4);
J_th    = (1/12) * (m) *l^2;
Ts = 0.05; 

data = iddata(y, input, Ts);

%% 
theta0 = [m*g*l; J_th; b_theta; b_phi; k; J_phi];
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


    A = [0 1 0 0; -mgl/J_theta -b_theta/J_theta 0 b_phi/J_theta; 0 0 0 1;
        mgl/J_theta b_theta/J_theta 0 -b_phi/J_theta];
    B = [0; -k/J_theta;0; k/J_theta + k/J_phi];
    % B = zeros(4,1);
    C = [1 0 0 0; 0 0 0 1];
    D = [0;0];
end
