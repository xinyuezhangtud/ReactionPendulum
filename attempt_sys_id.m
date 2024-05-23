close all
clear all
clc
%% 

% Load the data
load('no_sin_01amp.mat')
% load('no_input_02amp.mat')

% Prepare output data (assuming there is no input)
y = [theta.Data(80:end) phi_dot.Data(80:end)]; 

% Define constants
m_p = 0.170;
m_d = 0.03;
m = 0.355;
l = 0.9;
g = 9.81;
k = 0.06;
r = 0.05;
Ts = 0.05;

% Prepare the iddata object with no input
data = iddata(y, input, Ts);

% Initial guess for parameters [mgl, J_theta, k, J_phi]
theta0 = [1.4887; 0.0018; 1.5224; 0.0051];

% Create the grey-box model
sys = idgrey(@reaction_wheel_pendulum, theta0, 'c');

% Set parameter constraints
sys.Structure.Parameters(1).Minimum = 0;


% Set options for greyest to estimate initial state
opt = greyestOptions('InitialState', 'estimate', 'Display', 'on');

% Estimate the parameters
estimatedModel = greyest(data, sys, opt);

% Extract estimated parameters
estimatedParams = getpvec(estimatedModel);
disp('Estimated Parameters:');
disp(estimatedParams);

% Compare model output with data
compare(data, estimatedModel);

% Function defining the state-space model
function [A, B, C, D] = reaction_wheel_pendulum(theta, Ts)
    mgl = theta(1);
    J_theta = theta(2);
    k = theta(3);
    J_phi = theta(4);
    
    A = [0 1 0; -mgl/J_theta 0 0; mgl/J_theta 0 0];
    B = zeros(3,1);
    C = [1 0 0; 0 0 1];
    D = [0; 0];
end
