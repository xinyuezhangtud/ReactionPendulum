close all
clear all
clc

% Load the data
% load('no_sin_01amp.mat')
% load('no_input_02amp.mat')
load('no_input_fr.mat')
% input = u.Data(65:220);
% y = [theta.Data(65:220)]; 

input = u.Data(100:(30/0.05));
y = [theta.Data(100:(30/0.05))]; 


% Define constants
g = 9.81;

l = 0.09;
m = 1;
b = 0.1;

Ts = 0.05;

nx = 10; % Choose a number larger than the expected model order

% Prepare the iddata object with no input
data = iddata(y, input, Ts);
 

a_1 = -9.81/0.6;
a_2 = -0.012;

% Create the grey-box model
sys = idgrey(@pendulum_model, {a_1, a_2}, 'c');

% Set parameter constraints
% sys.Structure.Parameters(1).Minimum = 0;

% Set options for greyest to estimate initial state
opt = greyestOptions('InitialState', 'estimate', 'Display', 'on');
% opt = greyestOptions()
% Estimate the parameters
estimatedModel = greyest(data, sys, opt);

% Extract estimated parameters
estimatedParams = getpvec(estimatedModel);
disp('Estimated Parameters:');
disp(estimatedParams);

% Compare model output with data
compare(data, estimatedModel);


function [A, B, C, D] = pendulum_model(a_1,a_2, Ts)
    % g = theta(1);
    % g = 9.81;
    % l = theta(1);
    % m = theta(2);
    % b = theta(3);
    
    A = [0 1; a_1 a_2];
    B = [0; 0]; 
    C = [1 0]; 
    D = [0]; 
    % G_cont = ss(A,B,C,D);
    % G_disc = c2d(G_cont, 0.05, 'zoh');
    % A = G_disc.A;
    % B = G_disc.B;
end

