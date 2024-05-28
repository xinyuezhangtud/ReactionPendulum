% % Clear workspace and command window
% close all
% clear all
% clc
% 
% % Load datasets
% load('no_sin_01amp.mat')
% input2 = u.Data(400:900);
% y2 = [theta.Data(400:900) phi_dot.Data(400:900)]; 
% data2 = iddata(y2, input2, 0.05);
% 
% load('no_chirp_015amp.mat')
% input = u.Data(400:900);
% y = [theta.Data(400:900) phi_dot.Data(400:900)]; 
% data = iddata(y, input, 0.05);
% 
% % Physical parameters
% m_d     = 0.16;
% m       = 0.311;
% l       = 0.9163;
% b_theta = 0.0834;
% b_phi   = 0.011; 
% g       = 9.81;
% r       = 0.09;
% k       = 0.06;
% J_phi   = (1/2) * m_d * r^2 + m_d * ((l^2)/4);
% J_th    = (1/12) * (m) *l^2;
% 
% % Initial parameters for optimization
% theta0 = [m*g*l; J_th; b_theta; b_phi; k; J_phi];
% Ts = 0.05; 
% 
% % Set up grey box model
% sys = idgrey(@reaction_wheel_pendulum, theta0, 'c');
% sys.Structure.Parameters(1).Minimum = 0;
% 
% opt = greyestOptions('InitialState','estimate','Display','on');
% 
% % Optimization routine
% num_trials = 100; % Number of random trials
% best_params = theta0;
% best_score = Inf;
% 
% for i = 1:num_trials
%     % Randomize initial conditions
%     theta_init = theta0 .* (1 + 0.1 * randn(size(theta0)));
% 
%     % Estimate model
%     sys_init = idgrey(@reaction_wheel_pendulum, theta_init, 'c');
%     try
%         estimatedModel = greyest(data, sys_init, opt);
%         estimatedParams = getpvec(estimatedModel);
% 
%         % Validate on both datasets
%         [~, fit1, ~] = compare(data, estimatedModel);
%         [~, fit2, ~] = compare(data2, estimatedModel);
% 
%         % Combined score (you can adjust the weights as needed)
%         score = - (fit1 + fit2);
% 
%         if score < best_score
%             best_score = score;
%             best_params = estimatedParams;
%             disp('New best score found:');
%             disp(best_score);
%             disp('With parameters:');
%             disp(best_params);
%         end
%     catch
%         disp('Estimation failed for this set of initial conditions.');
%     end
% end
% 
% % Display best parameters
% disp('Best Parameters:');
% disp(best_params);
% 
% % Compare the best model with both datasets
% bestModel = idgrey(@reaction_wheel_pendulum, best_params, 'c');
% figure(1)
% compare(data, bestModel);
% figure(2)
% compare(data2, bestModel);
% 
% % Save best parameters for future use
% save('best_params.mat', 'best_params');
% 
% % Function definition
% function [A, B, C, D] = reaction_wheel_pendulum(theta, Ts)
%     mgl = theta(1);
%     J_theta = theta(2);
%     b_theta = theta(3);
%     b_phi = theta(4);
%     k = theta(5);
%     J_phi = theta(6);
% 
%     A = [0 1 0 0; 
%         -mgl/J_theta -b_theta/J_theta 0 b_phi/J_theta; 
%         0 0 0 1;
%         mgl/J_theta b_theta/J_theta 0 -b_phi/J_theta];
%     B = [0; -k/J_theta;0; k/J_theta + k/J_phi];
%     C = [1 0 0 0; 0 0 0 1];
%     D = [0;0];
% end

% Clear workspace and command window
close all
clear all
clc

% Load datasets
load('no_sin_01amp.mat')
input2 = u.Data(400:900);
y2 = [theta.Data(400:900) phi_dot.Data(400:900)]; 
data2 = iddata(y2, input2, 0.05);

load('no_chirp_015amp.mat')
input = u.Data(400:900);
y = [theta.Data(400:900) phi_dot.Data(400:900)]; 
data = iddata(y, input, 0.05);

% Initial physical parameters
m_d     = 0.16;
m       = 0.311;
l       = 0.9163;
b_theta = 0.0834;
b_phi   = 0.011; 
g       = 9.81;
r       = 0.09;
k       = 0.06;
J_phi   = (1/2) * m_d * r^2 + m_d * ((l^2)/4);
J_th    = (1/12) * (m) *l^2;

% Initial parameters for optimization
theta0 = [m*g*l; J_th; b_theta; b_phi; k; J_phi];
Ts = 0.05; 

% Set up grey box model
sys = idgrey(@reaction_wheel_pendulum, theta0, 'c');
sys.Structure.Parameters(1).Minimum = 0;

opt = greyestOptions('InitialState','estimate','Display','on');

% Optimization routine
num_trials = 100; % Number of random trials
best_params = theta0;
best_score = Inf;

for i = 1:num_trials
    % Randomize initial conditions with broader range
    theta_init = theta0 .* (1 + 0.5 * randn(size(theta0)));

    % Estimate model
    sys_init = idgrey(@reaction_wheel_pendulum, theta_init, 'c');
    try
        estimatedModel = greyest(data, sys_init, opt);
        estimatedParams = getpvec(estimatedModel);

        % Validate on both datasets
        [~, fit1, ~] = compare(data, estimatedModel);  % fit1 is the fit percentage for the first dataset
        [~, fit2, ~] = compare(data2, estimatedModel); % fit2 is the fit percentage for the second dataset

        % Define score based on fit percentages within desired range (75%-85%)
        fit_target = 80;
        fit_tolerance = 5;
        
        score1 = abs(fit_target - fit1) / fit_tolerance;
        score2 = abs(fit_target - fit2) / fit_tolerance;

        % Combined score (sum of normalized deviations)
        score = score1 + score2;

        % Check if this is the best score so far
        if score < best_score
            best_score = score;
            best_params = estimatedParams;
            disp('New best score found:');
            disp(best_score);
            disp('With parameters:');
            disp(best_params);
        end
    catch
        disp('Estimation failed for this set of initial conditions.');
    end
end

%% 

% Display best parameters
disp('Best Parameters:');
disp(best_params);

% Compare the best model with both datasets
bestModel = idgrey(@reaction_wheel_pendulum, best_params, 'c');
figure(1)
compare(data, bestModel);
figure(2)
compare(data2, bestModel);

% Save best parameters for future use
save('best_params.mat', 'best_params');
%% 

% Function definition
function [A, B, C, D] = reaction_wheel_pendulum(theta, Ts)
    mgl = theta(1);
    J_theta = theta(2);
    b_theta = theta(3);
    b_phi = theta(4);
    k = theta(5);
    J_phi = theta(6);

    A = [0 1 0 0; 
        -mgl/J_theta -b_theta/J_theta 0 b_phi/J_theta; 
        0 0 0 1;
        mgl/J_theta b_theta/J_theta 0 -b_phi/J_theta];
    B = [0; -k/J_theta;0; k/J_theta + k/J_phi];
    C = [1 0 0 0; 0 0 0 1];
    D = [0;0];
end

