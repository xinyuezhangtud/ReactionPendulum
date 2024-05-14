% This script is executed every time that an experiment is initialized

% Define the sampling time here:
h = 0.05;

save('zero_rad') % upward position
save('pi_rad') % downward position
save('sixtythree_deg') % 63 from upward position (to the right) 
save('one_twentytwo_deg') % 122 degree from dwonward to the left, 58 deg from upward counter clockwise (our reference) 
save('deg_71') % 251 degree  clockwise from down, 71 deg from upward clockwise (our reference) 
save('deg_145') % 145 degree  clockwise from down, 35 deg from upward counter-clockwise (our reference) 

%% Data no input 
save('theta.mat') % output data of no input and initial displacement of angle
save('phi_dot') % output data of no input and initial displacement of angle
save('i_m') % output data of no input and initial displacement of angle
% Data no input 2
save('theta2') % output data of no input and initial displacement of angle
save('phi_dot2') % output data of no input and initial displacement of angle
save('i_m2') % output data of no input and initial displacement of angle

%% Data sine input amplitude 0.15
save('theta_sine_0.15.mat') % output data of no input and initial displacement of angle
save('phi_dot_sine_015') % output data of no input and initial displacement of angle
save('i_m_sine_015') % output data of no input and initial displacement of angle
save('u_sin_015') % output data of no input and initial displacement of angle

%% Data frequency sweep input amplitude 0.2
save('theta_sweep_1') % output data of no input and initial displacement of angle
save('phi_dot_sweep_1') % output data of no input and initial displacement of angle
save('i_m_sweep_1') % output data of no input and initial displacement of angle
save('u_sweep_1') % output data of no input and initial displacement of angle

%% Data frequency square wave input amplitude 0.15
save('theta_square') % output data of no input and initial displacement of angle
save('phi_dot_square') % output data of no input and initial displacement of angle
save('i_m_square') % output data of no input and initial displacement of angle
save('u_square') % output data of no input and initial displacement of angle
%% Data frequency step input gain 0.5
save('theta_step') % output data of no input and initial displacement of angle
save('phi_dot_step') % output data of no input and initial displacement of angle
save('i_m_step') % output data of no input and initial displacement of angle
save('u_step') % output data of no input and initial displacement of angle