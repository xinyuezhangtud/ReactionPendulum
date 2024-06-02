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
J_th    = (1/12) * (m) * l^2;

% Initial parameters for optimization
theta0 = [m*g*l; J_th; b_theta; b_phi; k; J_phi];
Ts = 0.05;

% Define search ranges for each parameter, exploring a wider range but ensuring positive values
param_ranges = {
    logspace(log10(0.1 * m * g * l), log10(10 * m * g * l), 5), ...  % mgl
    logspace(log10(0.1 * J_th), log10(10 * J_th), 5), ...           % J_theta
    logspace(log10(0.1 * b_theta), log10(10 * b_theta), 5), ...     % b_theta
    logspace(log10(0.1 * b_phi), log10(10 * b_phi), 5), ...         % b_phi
    logspace(log10(0.1 * k), log10(10 * k), 5), ...                 % k
    logspace(log10(0.1 * J_phi), log10(10 * J_phi), 5)              % J_phi
};

% Initialize best score
best_score = Inf;
best_params = theta0;

% Grid search over all parameter combinations
for mgl = param_ranges{1}
    for J_theta = param_ranges{2}
        for b_theta = param_ranges{3}
            for b_phi = param_ranges{4}
                for k = param_ranges{5}
                    for J_phi = param_ranges{6}
                        % Current parameter set
                        theta = [mgl, J_theta, b_theta, b_phi, k, J_phi];

                        % Compute score for current parameter set
                        score = modelFitScore(theta, data, data2, Ts);

                        % Update best parameters if current score is better
                        if score < best_score && score >= 0
                            best_score = score;
                            best_params = theta;
                            disp('New best score found:');
                            disp(best_score);
                            disp('With parameters:');
                            disp(best_params);
                        end
                    end
                end
            end
        end
    end
end

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

% Function definition
function score = modelFitScore(theta, data, data2, Ts)
    % Create grey-box model with given parameters
    sys = idgrey(@reaction_wheel_pendulum, theta, 'c');

    % Set stability threshold for the model estimation
    opt = greyestOptions('InitialState', 'estimate', 'Display', 'off');
    opt.StabilityThreshold.A = 1e-6;  % Custom stability threshold for A matrix

    try
        % Estimate model
        estimatedModel = greyest(data, sys, opt);
        
        % Validate on both datasets
        [~, fit1, ~] = compare(data, estimatedModel);
        [~, fit2, ~] = compare(data2, estimatedModel);

        % Define score based on fit percentages within desired range (70%-85%)
        fit_min = 70;
        fit_max = 85;

        % Penalize deviations from the desired range
        penalty = 0;
        penalty = penalty + sum(max(0, fit_min - fit1)) + sum(max(0, fit1 - fit_max));
        penalty = penalty + sum(max(0, fit_min - fit2)) + sum(max(0, fit2 - fit_max));

        % Combined score
        score = penalty + sum(abs(fit1 - (fit_min + fit_max) / 2)) + sum(abs(fit2 - (fit_min + fit_max) / 2));

    catch
        % If estimation fails, return a large score
        score = 1e6;
    end
end

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
