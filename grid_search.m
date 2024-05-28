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
theta0 = [6.8969; 0.6429; 0.0681; 0.0057;1.9185;0.0058];

Ts = 0.05;

% Define search ranges for each parameter using linspace for finer control
param_ranges = {
    linspace(0.01 * m * g * l, 1000 * m * g * l, 10), ...  % mgl
    linspace(0.01 * J_th, 1000 * J_th, 10), ...           % J_theta
    linspace(0.01 * b_theta, 1000 * b_theta, 10), ...     % b_theta
    linspace(0.01 * b_phi, 1000* b_phi, 10), ...         % b_phi
    linspace(0.01 * k, 1000 * k, 10), ...                 % k
    linspace(0.01 * J_phi, 1000 * J_phi, 10)              % J_phi
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
                        
                        % Print current parameter set being evaluated
                        % disp('Evaluating parameters:');
                        % disp(theta);

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
                            
                            % Check if all fits are above 70%
                            if allFitsAboveThreshold(best_params, data, data2, Ts, 70)
                                disp('Stopping criterion met. Exiting search.');
                                break;
                            end
                        end
                    end
                    if allFitsAboveThreshold(best_params, data, data2, Ts, 70)
                        break;
                    end
                end
                if allFitsAboveThreshold(best_params, data, data2, Ts, 70)
                    break;
                end
            end
            if allFitsAboveThreshold(best_params, data, data2, Ts, 70)
                break;
            end
        end
        if allFitsAboveThreshold(best_params, data, data2, Ts, 70)
            break;
        end
    end
    if allFitsAboveThreshold(best_params, data, data2, Ts, 70)
        break;
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
%% 

% Save best parameters for future use
save('best_params.mat', 'best_params');

% Function definition
function score = modelFitScore(theta, data, data2, Ts)
    % Create grey-box model with given parameters
    sys = idgrey(@reaction_wheel_pendulum, theta, 'c');
    
    % Create a stable predictor
    sys = init(sys);

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

        % Penalize deviations below the desired range
        penalty = 0;
        penalty = penalty + sum(max(0, fit_min - fit1)) + sum(max(0, fit_min - fit2));

        % Deviation from target fit percentage (77.5%)
        target_fit = (fit_min + 85) / 2;
        deviation = sum(abs(fit1 - target_fit)) + sum(abs(fit2 - target_fit));

        % Combined score
        score = penalty + deviation;

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
        mgl/J_theta b_theta/J_theta 0 -b_phi/J_theta-(b_phi/J_phi)];
    B = [0; -k/J_theta;0; k/J_theta + k/J_phi];
    C = [1 0 0 0; 0 0 0 1];
    D = [0;0];
end

function stop = allFitsAboveThreshold(params, data, data2, Ts, threshold)
    % Create grey-box model with given parameters
    sys = idgrey(@reaction_wheel_pendulum, params, 'c');
    
    % Create a stable predictor
    sys = init(sys);

    % Set stability threshold for the model estimation
    opt = greyestOptions('InitialState', 'estimate', 'Display', 'off');
    opt.StabilityThreshold.A = 1e-6;  % Custom stability threshold for A matrix

    try
        % Estimate model
        estimatedModel = greyest(data, sys, opt);
        
        % Validate on both datasets
        [~, fit1, ~] = compare(data, estimatedModel);
        [~, fit2, ~] = compare(data2, estimatedModel);

        % Check if all fits are above the threshold
        if all(fit1 >= threshold) && all(fit2 >= threshold)
            stop = true;
        else
            stop = false;
        end
    catch
        % If estimation fails, return false
        stop = false;
    end
end
