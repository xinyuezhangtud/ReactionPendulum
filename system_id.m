close all
clear all
clc
load('chirp_new.mat')
x = 200;
% y = 100/0.05-1;
x1 = 400;
% load('no_chirp_015amp.mat')
load('no_chirp_02amp.mat')
input2 = u.Data(x1:end);
y2 = [theta.Data(x1:end) phi_dot.Data(x1:end)]; 
data2 = iddata(y2, input2, 0.05);
% load('no_sin_01amp.mat')
% load('no_chirp_015amp.mat')
% 
% load('no_sin_01amp.mat')
load('chirp_new_3.mat')
input = u.Data(x:end);
y = [theta.Data(x:end) phi_dot.Data(x:end)]; 


m_d     = 0.16;
m       = 0.311; %0.311
l       = 0.5; %0.9163
b_theta = 0.0834; %0.0834
b_phi   = 0.011; 
g       = 9.81;
r       = 0.09;
k       = 0.06;
J_phi   = (1/2) * m_d * r^2 + m_d * ((l^2)/4);
J_th    = (1/12) * (m) *l^2;
Ts = 0.05; 

data = iddata(y, input, Ts);


theta0 = [m*g*l; J_th; b_theta; b_phi; k; J_phi];
mgl=m*g*l;

% theta0 = [6.8969; 0.6429; 0.0681; 0.0057;1.9185;0.0058];
% theta0 = [6.8969; 3; 0.8; 0.01;1.9185;0.058];
%% 
a21=  -10.4142;
a22= -0.1793;
a24= b_phi/J_th;
a41=10.4142;
a42 =0.1793;
a44= -(b_phi/J_th) -b_phi/J_phi;
b2 = -k/J_th;
b4 = k/J_th + k/J_phi;

a21=  -10.4142;
a22= -0.1793;
a24= -0.0030;
a41=10.4142;
a42 =0.1793;
a44= -0.9433;
b2 = -3.0445;
b4 = 343.9576;

sys= idgrey(@reaction_wheel_pendulum, {a21,a22,a24,a42,a44,b2,b4}, 'c');
% sys.Structure.Parameters(1).Minimum = 0;

sys.Structure.Parameters(1).Free= true;
sys.Structure.Parameters(2).Free= true;
sys.Structure.Parameters(4).Free= true;
sys.Structure.Parameters(5).Free= true;

opt = greyestOptions('InitialState','estimate','Display','on');
estimatedModel = greyest(data, sys, opt);

estimatedParams = getpvec(estimatedModel);
disp('Estimated Parameters:');
disp(estimatedParams);

figure(1)
compare(data,estimatedModel);

figure(2)
compare(data2,estimatedModel);
%% 
[~, fit1, ~] = compare(data, estimatedModel);
%% 

function [A, B, C, D] = reaction_wheel_pendulum(a21,a22,a24,a42,a44,b2,b4,Ts)
    % mgl = theta(1);
    % J_theta = theta(2);
    % b_theta = theta(3);
    % b_phi = theta(4);
    % k = theta(5);
    % J_phi = theta(6);


    A = [0 1 0; 
        a21 a22 a24; 
        -a21 a42 a44];
    B = [0; b2; b4];
    % B = zeros(4,1);
    C = [1 0 0; 0 0 1];
    D = [0;0];
end
