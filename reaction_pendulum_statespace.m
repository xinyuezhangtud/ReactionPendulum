clear all
close all
clc

a21 = -9.9834;
a22 = -0.1314;
a24 = 0.0110;
a41 = 9.9834;
a42 = -0.0474;
a44 = -0.9340;
b2  = -3.0817;
b4  = 348.1863;

A = [0 1 0; 
     -a21 a22 a24; 
     -a41 a42 a44];

B = [0; b2; b4];

    % B = zeros(4,1);
C = [1 0 0; 0 0 1];
D = [0;0];


x_e = [pi/2; 0; 0 ;0];

disp('A');
disp(A);
disp('B')

disp(B);
disp('C')
disp(C);

disp('controllability');

disp(rank(ctrb(A, B)));

disp('observability')
disp(rank(obsv(A,C)));

%% EXPORT A, B, C, D matrices
save reaction_pendulum.mat A B C D x_e

%% Discretize the system
Ts = 0.05;

G = ss(A,B,C,D);

G_dt = c2d(G, Ts, 'zoh');

[Phi, Gamma, C, D] = ssdata(G_dt);

disp('Phi');
disp(Phi);

disp('Gamma');
disp(Gamma);

disp('C');
disp(C);

disp('D');
disp(D);

disp('Controllability')
disp(rank(ctrb(Phi, Gamma)));

%% EXPORT Phi, Gamma, C, D matrices
save reaction_pendulum.mat A B Phi Gamma C D x_e Ts

%% Simulation

% x_0 = [0.1 0 0 0]';
% t = 0:Ts:5;
% u = 0*t;
% y= lsim(G_dt, u, t, x_0);
% 
% figure(1)
% plot(y)
% legend('theta','phi')
