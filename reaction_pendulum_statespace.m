clear all
close all
clc

a21 = -9.9678;
a22 = -0.1416;
a24 =  0.0114;
a44 = -0.9971;
b2  = -3.0585;
b4  = 346.8628;

A = [0 1 0; 
     -a21 a22 a24; 
     0 0 a44];

B = [0; b2; b4];

    % B = zeros(4,1);
C = [1 0 0; 0 0 1];
D = [0;0];

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
save reaction_pendulum.mat A B C D

%% Discretize the system
h = 0.05;

G = ss(A,B,C,D);

G_dt = c2d(G, h, 'zoh');

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
save reaction_pendulum.mat A B Phi Gamma C D h

%% Simulation

% x_0 = [0.1 0 0 0]';
% t = 0:Ts:5;
% u = 0*t;
% y= lsim(G_dt, u, t, x_0);
% 
% figure(1)
% plot(y)
% legend('theta','phi')
