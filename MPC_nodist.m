clear all
close all
clc
%%
% % INSPIRATION TAKEN FROM MPC CODE 
% load('b_in_Tset.mat')
% load('A_in_Tset.mat')

%%
% INSPIRATION TAKEN FROM MPC CODE 
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

    % B = zeros(4,1);
C = [1 0 0; 0 0 1];
D = [0;0];

% Discretize system using 'zoh'
sys_CT = ss(A,B,C,D);
sys_DT = c2d(sys_CT, Ts, 'zoh');

G_disc.A = sys_DT.A; 
G_disc.B = sys_DT.B;
G_disc.C = sys_DT.C;
G_disc.D = sys_DT.D;


weight.Q =[5 0 0;0 4.5 0; 0 0 0.0071];   % weight on output
weight.R = 30;
[P, K, eigvals] = dare(G_disc.A, G_disc.B, weight.Q, weight.R, [], []);
weight.P = P;
% R = 10, N = 40 works, otherise terminal ste is a bit not converged
% Simulation horizon and input constraints
T_sim = 500;
umin = -4; % 0.1 A as rated current
umax = -umin;  
dim.nx = size(G_disc.A,1); % state dimension
dim.nu = size(G_disc.B,2); % input dimension
dim.ny = size(G_disc.C,1);
dim.N = 20;

% Define prediction model + cost
[S,T] = predmodgen(G_disc,dim);

% Cost definition
Q_bar = blkdiag(kron(eye(dim.N), weight.Q), weight.P);
H = S' * Q_bar * S + kron(eye(dim.N),weight.R);   
h = S' * Q_bar * T;

% MPC formulation
x = zeros(dim.nx,T_sim+1);
u_rec = zeros(dim.nu,T_sim);
G_disc.x0=[0.01; 0; 0];
x(:,1) = G_disc.x0;

% Inequality Constraints
x1_constraint = 0.3; % it is about 17 degrees
x2_constraint = 2; % 0.3/0.05 --> goes from 0.3 rad to 0 in one time step
x3_constraint = 400; % about 85 rad/sec (see datasheet)

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
 % Compute terminal set using MPT3
model = LTISystem(sys_DT);
model.x.min = -[x1_constraint; x2_constraint; x3_constraint];
model.x.max = [x1_constraint; x2_constraint; x3_constraint];
model.u.max = umax;
model.u.min = umin;
model.x.penalty = QuadFunction(weight.Q);
model.u.penalty = QuadFunction(weight.R);
Tset = model.LQRSet;
% PN = model.LQRPenalty;
% model.x.with('terminalSet');
% model.x.terminalSet = Tset;
% model.x.with('terminalPenalty');
% model.x.terminalPenalty = PN;
% % dim.N = 50;
% ctrl = MPCController(model, dim.N)
% x0 = [0.1; 0; 0];
% Nsim = 500;
% loop = ClosedLoop(ctrl, model);
% data = loop.simulate(x0, Nsim);
% figure(1)
% subplot(2,1,1)
% plot(1:Nsim, data.Y);
% hold on;
% title('outputs')
% subplot(2,1,2)
% plot(1:Nsim, data.U);
% title('inputs')

%% MPC
total_iterations = T_sim;
o = waitbar(0, 'Progress');
for k = 1 : T_sim
    Constraint = [];
    x_0 = x(:,k); s
    u_horizon = sdpvar(dim.nu*dim.N,1);  % [5,1] where 5 is horizon, 1, input dim          
    % MPC using quadprog
    
    A_in_QP = [A_in; Tset.A*Ae];
    b_in_QP = [(eps-eps_1*x_0); (Tset.b - Tset.A*be*x_0)];
    lb = umin*ones(dim.N,1);
    ub = umax*ones(dim.N,1);
    [u_horizon_QP,fval,exitflag,output]= quadprog(H, (h * x_0), A_in_QP, b_in_QP, [],[], lb, ub );
    exitflag
    % Contraints
%     Constraint = [Constraint, 
%                     umin <= u_horizon <= umax,
%                     A_in * u_horizon <= eps - eps_1*x_0,
%                     Tset.A * Ae * u_horizon <= Tset.b - Tset.A*be*x_0,
%                     ];                                        %define constraints
%     
    % A_in has size 126 x 20
    % Tset.A * Ae  has size 34 x 20 --> the augmented A has size 160, 20
    % size(eps - eps_1*x_0) = 126 x 1
    % size(Tset.b - Tset.A*be*x_0) = 34 x 1 -- > the augmented version = 160 ,1


    % Objective function

%     Objective = u_horizon' * H * u_horizon + (h * x_0)' * u_horizon;     
    
    % Optimize
%     ops = sdpsettings('verbose', 0, 'solver', 'quadprog');  % Specify 'quadprog' as the solver
%     ops =  sdpsettings('verbose',0);
% %     sol = optimize(Constraint,Objective, ops);                 %solve the problem
    
%     if sol.problem == 0
%         u_horizon = value(u_horizon);                          %assign the solution
        % Select the first input only
        u_rec(k)=u_horizon_QP(1);
        % Compute the state/output evolution
        x(:,k+1)=G_disc.A*x_0 + G_disc.B*u_rec(k);        

%         clear u_horizon
%     else
%         disp('Problem solving MPC optimization');
%         break
%     end
    
    waitbar(k /total_iterations , o, sprintf('Progress: %d%%', round(k / total_iterations * 100)));
end
close(o);

%% Plots

% MPC
time = (0:k) * 0.05;
figure(10)
stairs(time, x(1,:)),
hold on
stairs(time,x(2,:)),
stairs(time, x(3,:)),
xlabel('Time t [s]'), 
ylabel('States x MPC'), 
grid on;
legend('x_1','x_2', ' x_3', ' x_4');

time_lq = (0:k-1) * 0.05;
figure(20)
stairs(time_lq, u_rec),
xlabel('Time t [s]'), 
ylabel('Control input u MPC'), 
grid on;
legend('u');

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