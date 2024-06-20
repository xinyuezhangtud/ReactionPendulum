clear all
close all
clc

%% INSPIRATION TAKEN FROM MPC CODE 
load('reaction_pendulum.mat')
load('InvSet_A.mat');
load('InvSet_b.mat')
%% 
a21 = -9.9678;
a22 = -0.1416;
a24 =  0.0114;
a44 = -0.9971;
b2  = -3.0585;
b4  = 346.8628;
Ts = 0.05;

A = [0 1 0; 
     -a21 a22 a24; 
     0 0 a44];
B = [0; b2; b4];
C = [1 0 0; 0 0 1];
D = [0;0];
% Discretize system using 'zoh'
Ts = 0.05;
sys_CT = ss(A,B,C,D);
sys_DT = c2d(sys_CT, Ts, 'zoh');
G_disc.A = sys_DT.A; 
G_disc.B = sys_DT.B;
G_disc.C = sys_DT.C;
G_disc.D = sys_DT.D;

weight.Q =[4 0 0;0 3.5 0; 0 0 0.0071];  
weight.R = 1e2;
[P, K, eigvals] = dare(G_disc.A, G_disc.B, weight.Q, weight.R, [], []);
weight.P = P;
T_sim = 500; % Simulation horizon and input constraints
umin = -4; 
umax = -umin;  
dim.nx = size(G_disc.A,1); % state dimension
dim.nu = size(G_disc.B,2); % input dimension
dim.ny = size(G_disc.C,1);
dim.N = 50;

% Define prediction model + cost
[S,T] = predmodgen(G_disc,dim);

% Cost definition
Q_bar = blkdiag(kron(eye(dim.N), weight.Q), weight.P);
H = S' * Q_bar * S + kron(eye(dim.N),weight.R);   
h = S' * Q_bar * T;

% MPC formulation
x = zeros(dim.nx,T_sim+1);
u_rec = zeros(dim.nu,T_sim);
G_disc.x0=[-0.2; 1.0; 0];
x(:,1) = G_disc.x0;

% Inequality Constraints
x1_constraint = 0.25; 
x2_constraint = 2; 
x3_constraint = 400; 

up =[x1_constraint; x2_constraint; x3_constraint];
lb = -up;
bounds = [up; -lb];
eps = repmat(bounds, dim.N+1,1);

f = [eye(3); -eye(3)];
F = kron(eye(dim.N+1), f);
eps_1 = F*T;
A_in = F*S;  

% Equality constraints (terminal)
Ae = S(end-dim.nx+1:end,:);
be = T(end-dim.nx+1:end,:);
% 
% % Compute terminal set using MPT3
% model = LTISystem(sys_DT);
% model.x.min = -[x1_constraint; x2_constraint; x3_constraint];
% model.x.max = [x1_constraint; x2_constraint; x3_constraint];
% model.u.max = umax;
% model.u.min = umin;
% model.x.penalty = QuadFunction(weight.Q);
% model.u.penalty = QuadFunction(weight.R);
% Tset = model.LQRSet;
% 
% system = LTISystem('A',[G_disc.A],'B', [G_disc.B],'C', G_disc.C,'D', G_disc.D);
% system.x.min = -[x1_constraint; x2_constraint; x3_constraint];
% system.x.max = [x1_constraint; x2_constraint; x3_constraint];
% system.u.min = umin;
% system.u.max = umax;
% InvSet = system.invariantSet()
% figure(4)
% InvSet.plot()
% xlabel('theta')
% ylabel('theta dot')
% zlabel('phi dot')

%% MPC
total_iterations = T_sim;
o = waitbar(0, 'Progress');


total_iterations = T_sim;

o = waitbar(0, 'Progress');
for k = 1 : T_sim
    Constraint = [];
    x_0 = x(:,k); 
%     u_horizon = sdpvar(dim.nu*dim.N,1);  % [5,1] where 5 is horizon, 1, input dim          
    % MPC using quadprog
    
    A_in_QP = [A_in; A_Sinv*Ae];
    b_in_QP = [(eps-eps_1*x_0); (b_Sinv- A_Sinv*be*x_0)];
    lb = umin*ones(dim.N,1);
    ub = umax*ones(dim.N,1);
    [u_horizon_QP,fval,exitflag,output]= quadprog(H, (h * x_0), A_in_QP, b_in_QP, [],[], lb, ub );
    exitflag
    u_rec(k)=u_horizon_QP(1);
    % Compute the state/output evolution
    x(:,k+1)=G_disc.A*x_0 + G_disc.B*u_rec(k);        

    waitbar(k /total_iterations , o, sprintf('Progress: %d%%', round(k / total_iterations * 100)));
end
close(o);


%% Plots

% time = (0:k) * 0.05;
% figure(10)
% stairs(time, x(1,:)),
% hold on
% % stairs(time,x(2,:)),
% % stairs(time, x(3,:)),
% xlabel('Time t [s]','Interpreter', 'latex'), 
% ylabel('$$\theta$$','Interpreter', 'latex'), 
% grid on;
% legend('$$N = 2$$','$$N = 5$$','$$N = 20$$','$$N = 50$$', 'Interpreter', 'latex');
% 
% figure(2)
% stairs(time, x(2,:)),
% hold on
% % stairs(time,x(2,:)),
% % stairs(time, x(3,:)),
% xlabel('Time t [s]','Interpreter', 'latex'), 
% ylabel('$$\dot{\theta}$$','Interpreter', 'latex'), 
% grid on;
% legend('$$N = 2$$','$$N = 5$$','$$N = 20$$','$$N = 50$$', 'Interpreter', 'latex');
% 
% figure(3)
% stairs(time, x(3,:)),
% hold on
% % stairs(time,x(2,:)),
% % stairs(time, x(3,:)),
% xlabel('Time t [s]','Interpreter', 'latex'), 
% ylabel('$$\dot{\phi}$$','Interpreter', 'latex'), 
% grid on;
% legend('$$N = 2$$','$$N = 5$$','$$N = 20$$','$$N = 50$$', 'Interpreter', 'latex');
% 
% time_lq = (0:k-1) * 0.05;
% figure(20)
% stairs(time_lq, u_rec)
% hold on
% xlabel('Time t [s]','Interpreter', 'latex'), 
% ylabel('Control input $$u$$','Interpreter', 'latex'), 
% legend('$$N = 2$$','$$N = 5$$','$$N = 20$$','$$N = 50$$', 'Interpreter', 'latex');
% grid on;


%% Define prediction model 

% Generation of prediction model -> Prediction matrix from initial state
function [S,T]=predmodgen(LTI,dim)

T=zeros(dim.nx*(dim.N+1),dim.nx); % nx = 4

for k=0:dim.N
    T(k*dim.nx+1:(k+1)*dim.nx,:)=LTI.A^k;
end

%Prediction matrix from input
S=zeros(dim.nx*(dim.N+1),dim.nu*(dim.N));
for k=1:dim.N
    for i=0:k-1
        S(k*dim.nx+1:(k+1)*dim.nx,i*dim.nu+1:(i+1)*dim.nu)=LTI.A^(k-1-i)*LTI.B;
    end
end

predmod.T=T;
predmod.S=S;

end