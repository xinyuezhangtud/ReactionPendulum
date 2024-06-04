clear all
close all
clc
%%
% INSPIRATION TAKEN FROM MPC CODE 
a21 =  -10.4142;
a22 = -0.1793;
a24 = -0.0030;
a41 = 10.4142;
a42 = 0.1793;
a44 = -0.9433;
b2 = -3.0445;
b4 = 343.9576;
Ts = 0.05;

A = [0 1 0; 
     a21 a22 a24; 
     a41 a42 a44];
B = [0; b2; b4];
C = [1 0 0; 0 0 1];
D = [0;0];

% Discretize system using 'zoh'
sys_CT = ss(A,B,C,D);
sys_DT = c2d(sys_CT, Ts, 'zoh');

LTI.A = sys_DT.A; 
LTI.B = sys_DT.B;
LTI.C = sys_DT.C;
LTI.D = sys_DT.D;

v = [10^5 10^4 10^2];
weight.Q = diag(v);   % weight on output
weight.R = 10;
[P, K, eigvals] = dare(LTI.A, LTI.B, weight.Q, weight.R, [], []);
weight.P = P;

% R = 10, N = 40 works, otherise terminal ste is a bit not converged

% Simulation horizon and input constraints
T_sim = 500;
umin = -0.11; % 0.1 A as rated current
umax = -umin;  
dim.nx = size(LTI.A,1); % state dimension
dim.nu = size(LTI.B,2); % input dimension
dim.ny = size(LTI.C,1);
dim.N = 40;

% Define prediction model + cost
[S,T] = predmodgen(LTI,dim);

% Cost definition
Q_bar = blkdiag(kron(eye(dim.N), weight.Q), weight.P);
H = S' * Q_bar * S + kron(eye(dim.N),weight.R);   
h = S' * Q_bar * T;

% MPC formulation
x = zeros(dim.nx,T_sim+1);
u_rec = zeros(dim.nu,T_sim);
LTI.x0=[0.1; 0; 0];
x(:,1) = LTI.x0;

% Inequality Constraints
x1_constraint = 0.3; % it is about 17 degrees
x2_constraint = 6; % 0.3/0.05 --> goes from 0.3 rad to 0 in one time step
x3_constraint = 9; % about 85 rad/sec (see datasheet)

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
model.x.max = [x1_constraint; x2_constraint; x3_constraint];;
model.u.max = umax;
model.u.min = umin;
model.x.penalty = QuadFunction(weight.Q);
model.u.penalty = QuadFunction(weight.R);
Tset = model.LQRSet;
PN = model.LQRPenalty;
model.x.with('terminalSet');
model.x.terminalSet = Tset;
model.x.with('terminalPenalty');
model.x.terminalPenalty = PN;
% dim.N = 50;
ctrl = MPCController(model, dim.N)
x0 = [0.1; 0; 0];
Nsim = 500;
loop = ClosedLoop(ctrl, model);
data = loop.simulate(x0, Nsim);
figure(1)
subplot(2,1,1)
plot(1:Nsim, data.Y);
hold on;
title('outputs')
subplot(2,1,2)
plot(1:Nsim, data.U);
title('inputs')

%% MPC
total_iterations = T_sim;
o = waitbar(0, 'Progress');

for k = 1 : T_sim
    Constraint = [];
    x_0 = x(:,k); 
    u_horizon = sdpvar(dim.nu*dim.N,1);  % [5,1] where 5 is horizon, 1, input dim          
    
    % Contraints
    Constraint = [Constraint, 
                    umin <= u_horizon <= umax,
                    A_in * u_horizon <= eps - eps_1*x_0,
                    Tset.A*(Ae * u_horizon + be * x_0) <= Tset.b,
                    ];                                        %define constraints
    
    % Objective function
    Objective = u_horizon' * H * u_horizon + (h * x_0)' * u_horizon;     
    
    % Optimize
    ops =  sdpsettings('verbose',0);
    sol = optimize(Constraint,Objective, ops);                 %solve the problem
    
    if sol.problem == 0
        u_horizon = value(u_horizon);                          %assign the solution
        % Select the first input only
        u_rec(k)=u_horizon(1);
        % Compute the state/output evolution
        x(:,k+1)=LTI.A*x_0 + LTI.B*u_rec(k);        

        clear u_horizon
    else
        disp('Problem solving MPC optimization');
        break
    end
    
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

%%

function [l,Vf] = lyapunov(u,x0,Q,R,P,LTI)
x = zeros(4,length(u)+1);
x(:,1) = x0;
N = length(u);
for i = 1:N
    x(:,i+1) = LTI.A*x(:,i)+LTI.B*u(i);
end

Vf = x(:,end)'*P*x(:,end);
l = x(:,end)'*Q*x(:,end) + u(end)'*R*u(end);

end