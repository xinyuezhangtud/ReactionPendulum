close all
clear all
clc

load('no_sin_01amp.mat')

input = u.Data(100:700);
y = [theta.Data(100:700), phi_dot.Data(100:700)];

% Define physical parameters
m_d     = 0.16;
m       = 0.511;
l       = 0.7163;
b_theta = 0.08;
b_phi   = 0.08;
g       = 9.81;
r       = 0.09;
k       = 0.06;
J_phi   = (1/2) * m_d * r^2 + m_d * ((l^2)/4);
J_theta = (1/12) * (m) * l^2;
Ts      = 0.05;

% Create iddata object
data = iddata(y, input, Ts);
data.OutputName = {'Pendulum Position', 'Disk Velocity'};
data.OutputUnit = {'rad', 'rad/s'};
data.Tstart = 0;
data.TimeUnit = 's';
data.InputName = 'Current';

% Define estimation parameters
order = [2, 1, 4]; % [ny, nu, nx]
mgl = m*g*l;
parameters = {mgl; J_theta; b_theta; b_phi; k; J_phi};
initial_states = [y(1,1), 0, 0, y(2,1)]';
theta0 = {mgl; J_theta; b_theta; b_phi; k; J_phi};
theta0 = [6.8969; 0.6429; 0.0681; 0.0057;1.9185;0.0058];
T = 0;
% Create the idnlgrey model object
nonlinear_model = idnlgrey(@pendulum, order, theta0, initial_states, T);
setpar(nonlinear_model,'Fixed',{false false false false false false});
nonlinear_model = nlgreyest(data, nonlinear_model)%, 'Display', 'Full');
[nonlinear_b_est, dnonlinear_b_est] = getpvec(nonlinear_model, 'free');

% Compare the estimated model with data
compare(data, nonlinear_model);

% Function for the nonlinear grey-box model
function [dx, y] = pendulum(t, x, u, mgl, J_theta, b_theta, b_phi, k, J_phi, varargin)
    dx(1) = x(2);
    dx(2) = -(mgl*cos(x(1))/J_theta) - (b_theta/J_theta)*x(2) + (b_phi/J_theta)*x(4) - (k/J_theta) * u;
    dx(3) = x(4);
    dx(4) = (mgl*cos(x(1))/J_theta) + (b_theta/J_theta) * x(2) - (b_phi/J_theta)*x(4) - (b_phi/J_phi)*x(4) + ((k/J_theta) + (k/ J_phi))*u;

    % Output equations
    y = [x(1); x(4)];
end