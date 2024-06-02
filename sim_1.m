close all
clear all
clc


m_d     = 0.46;
m       = 0.7;
l       = 0.6163;
b_theta = 1;% best with 0.08
b_phi   = 2.9 ;
g       = 9.81;
r       = 0.05;
k       = 0.006;
J_phi   = (1/2) * m_d * r^2 + m_d * ((l^2)/4);
J_th    = (1/12) * (m) *l^2;
Ts = 0.05; 
l_cg = 0.1;
mgl = m*9.81*l_cg;
theta = [m*g*l; J_th; b_theta; b_phi; k; J_phi];


    mgl = theta(1);
    J_theta = theta(2);
    k = theta(5);
    J_phi = theta(6);
    b_theta = theta(3);
    b_phi = theta(4);
% 
    A = [-b_theta/J_theta, mgl/J_theta, b_phi/J_theta,                   0;
        1,                0,           0,                               0;
        b_theta/J_theta,   -mgl/J_theta -b_phi*((1/J_theta)+(1/J_phi)),  0;
        0,                0,            0,                              1;
    ];

%     A = [0 1 0 0; 
%         -mgl/J_theta -b_theta/J_theta 0 b_phi/J_theta; 
%         0 0 0 1;
%         mgl/J_theta b_theta/J_theta 0 -b_phi/J_theta];

    B = [-k/J_theta; 0; (k/J_theta + k/J_phi); 0];
    % B = zeros(4,1);
%     B = [0; -k/J_theta; 0; (k/J_theta + k/J_phi)]
    C = [0 0 1 0; 0 0 1 0];
    D = [0;0];

    sys = ss(A,B,C,D);
    t = 0:0.002:2;  % 201 points
u = zeros(1, size(t,2));
[y, tOut, x] = lsim(sys, u,t,[0.2, 0, 0 0])

figure(1)
plot(y)


%% 
% 
% x0 = [0.2,0,0,0]; %initial condition
% t_a = linspace(0,1,20); %time vector
% [t,x] = ode45(@(t,x)myfunction(t, x, u, mgl, J_theta, b_theta, b_phi, k, J_phi),t_a,x0);
% figure
% plot(t_a,x.','LineWidth',1.5);
% 
% 
% function [dx, y] = pendulum(t, x, u, mgl, J_theta, b_theta, b_phi, k, J_phi, varargin)
% State equations
% dx = zeros(4, 1);
% dx(1) = x(2);
% dx(2) = -(mgl*cos(x(1))/J_theta) - (b_theta/J_theta)*x(2) + (b_phi/J_theta)*x(4) - (k/J_theta) * u;
% dx(3) = x(4);
% dx(4) = (mgl*cos(x(1))/J_theta) + (b_theta/J_theta) * x(2) - (b_phi/J_theta)*x(4) + ((k/J_theta) + (k/ J_phi))*u;
% 
% Output equations
% y = [x(1); x(4)];
% % Output equation
% y = [x(1);
%     x(4)];
% 
% State equation
% dx = [x(2);
%     -(mgl*cos(x(1))/J_theta) -(b_theta/J_theta)*x(2) + (b_phi/J_theta)*x(4) -(k/J_theta)*u;
%     x(3);
%     (mgl*cos(x(1))/J_theta) + (b_theta/J_theta)*x(2) - (b_phi/J_theta)*x(4) + ((k/J_theta)+ (k/J_phi))*u;
%     ];
% end