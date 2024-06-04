clear all
close all
clc

%% INSPIRATION TAKEN FROM MPC CODE 
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

% Simulation horizon and input constraints
T_sim = 500;
umin = -0.11; % 0.1 A as rated current
umax = -umin;  
dim.nx = size(LTI.A,1); % state dimension
dim.nu = size(LTI.B,2); % input dimension
dim.ny = size(LTI.C,1);
dim.N = 50;

% Define prediction model + cost
[S,T] = predmodgen(LTI,dim);

% Cost definition
Q_bar = blkdiag(kron(eye(dim.N), weight.Q), weight.P);
H = S' * Q_bar * S + kron(eye(dim.N),weight.R);   
h = S' * Q_bar * T;

% OUTPUT FEEDBACK
% Definition of extended system dimension
LTI.Cd = [0; 0];  
dim.nd = size(LTI.Cd,2);      % disturbance dimension = 1
dime.nx= dim.nx + dim.nd;     % state dimension = 3 + 1 =4
dime.nu= dim.nu;              % input dimension = 1
dime.ny= dim.ny;              % output dimension = 2
dime.N= dim.N;                % Horizon, same as before

LTI.yref = [0; 0];
LTI.d = 0.05;
LTI.Bd =  [1; 0; 0];
LTIe.A=[LTI.A LTI.Bd; zeros(dim.nd,dim.nx) eye(dim.nd)];
LTIe.B=[LTI.B; zeros(dim.nd,dim.nu)];
LTIe.C=[LTI.C LTI.Cd];
LTIe.yref = LTI.yref;

% Cost definition extended 
weighte.Q = blkdiag(weight.Q, zeros(dim.nd));           % weight on output [x;d]*[Q 0; 0 0]*[x;d]
weighte.R = weight.R;                                   % weight on input = same as before
weighte.P = blkdiag(weight.P, zeros(dim.nd));           % weight on final state [x;d]*[P 0; 0 0]*[x;d]

[Se, Te] = predmodgen(LTIe, dime);
[He,he]=costgen(Se, Te, weighte, dime);  

% MPC define variables
x = zeros(dim.nx,T_sim+1);
u_rec = zeros(dim.nu,T_sim);
% Inequality Constraints
x1_constraint = 0.3;        % it is about 17 degrees
x2_constraint = 6;          % 0.3/0.05 --> goes from 0.3 rad to 0 in one time step
x3_constraint = 9;          % about 85 rad/sec (see datasheet)
up =[x1_constraint; x2_constraint; x3_constraint; inf];
lb = -up;
bounds = [up; -lb];
eps = repmat(bounds, dime.N+1,1);
% Constraints for extended system
fe = [eye(dime.nx); -eye(dime.nx)]; 
F = kron(eye(dime.N+1), fe);         % [(N+1)*2*dim.nx] x [(N+1)*dim.nx] 
eps_1 = [F*Te F*Se F];
A_in = F*Se;  
% Equality constraints (terminal)
Ae = Se(end-dime.nx+1:end,:);
be = Te(end-dime.nx+1:end,:);
% Initial condition
xehat(:,1)=[0.0 ;0 ;0 ;0.005];
y(:,1)=LTIe.C * xehat(:,1);

%% Observer
syse = ss(LTIe.A, LTIe.B, LTIe.C, [], Ts);
pol = pole(syse);
desired_observer_poles = [pol(1)*0.5, pol(2)*0.5, pol(3)*0.3, pol(4)*0.2];
L = place(LTIe.A', LTIe.C', desired_observer_poles)';

rank_r = rank([eye(3)-LTI.A -LTI.Bd;LTI.C LTI.Cd]); 
y_ref = repmat(LTI.yref(1), T_sim,1);

%% MPT3 used to compare plots and to compute LQR set
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
ctrl = MPCController(model, dim.N);
x0 = [0.0; 0; 0];
Nsim = T_sim;
% loop = ClosedLoop(ctrl, model);
% data = loop.simulate(x0, Nsim);

%% MPC
total_iterations = T_sim;
o = waitbar(0, 'Progress');

for k = 1 : T_sim
    Constraint = [];
    dhat = xehat(4,k);
    x_0 = xehat(:,k);
    u_horizon = sdpvar(dim.nu*dim.N,1);        
    eqconstraints = eqconstraintsgen(LTI, dim, dhat);

    [xr,ur]=optimalss(dim, eqconstraints);
    xer = [xr; dhat];    % Added thiiiiiiis

    xer_e = repmat(xer, (dim.N+1),1);
    uer_e = repmat(ur, (dim.N),1);
    termdisA = [Tset.A zeros(size(Tset.A,1),1)];
    Constraint = [Constraint, 
                    umin <= u_horizon <= umax,
                    A_in * u_horizon <= eps - F*Te*x_0 + F*Se*uer_e - F*xer_e,
                    termdisA*(Ae * u_horizon + be * x_0) <= Tset.b,
                    ];                                       
    
    % Objective function
    Objective = u_horizon' * He * u_horizon + (he*[x_0; xer; ur])' * u_horizon;     
    
    % Optimize
    ops =  sdpsettings('verbose',0);
    sol = optimize(Constraint,Objective, ops);                        % solve the problem
    
    if sol.problem == 0
        u_horizon = value(u_horizon);                                 % assign the solution
        u_rec(k) = u_horizon(1);
        % Update states estimation
        y(:,k) = LTIe.C * xehat(:,k);
        xehat(:,k+1) = LTIe.A * xehat(:,k) + LTIe.B * u_rec(k) + L * (y(:,k) - LTIe.C *  xehat(:,k));
        dhat_plot(k) = xehat(4,k);
        clear u_horizon
    else
        close(o);
        disp('Problem solving MPC optimization');
        break
    end
    
    waitbar(k /total_iterations , o, sprintf('Progress: %d%%', round(k / total_iterations * 100)));    
end
% close(o);
%% PLOTS
time = (0:T_sim-1) * 0.02;

e=y-kron(ones(1,T_sim),LTI.yref);
figure(10000)
plot(time,e(1,:))
hold on
plot(time,e(2,:))
hold on 
xlabel('k'), ylabel('Tracking error'), 
grid on;
legend('e_1','e_2', 'e_3');

figure(10)
stairs(time, y(1,:)),
hold on
stairs(time, y(2,:)),
xlabel('Time t [s]')
legend('y1', 'y2')

figure(20)
stairs(time, dhat_plot)
hold on
xlabel('Time t [s]'), 

legend('d hat1')

figure(21)
stairs(time, y(1,:))
hold on
stairs(time, kron(ones(1,T_sim),LTI.yref(1)))
xlabel('Time t [s]'), 

legend('y', ' yref')


%% Cost generation function

function [H,h]=costgen(S, T ,weight,dim)
    Qbar = blkdiag(kron(eye(dim.N),weight.Q),weight.P);
    Rbar =kron(eye(dim.N),weight.R);
    H = S'* Qbar* S+Rbar;   
    hx0 = S'*Qbar* T;
    hxref = - S'*Qbar*kron(ones(dim.N+1,1),eye(dim.nx));
    huref = -Rbar*kron(ones(dim.N,1),eye(dim.nu));
    h = [hx0 hxref huref];
end

%% Define prediction model 

% Generation of prediction model -> Prediction matrix from initial state
function [S,T]=predmodgen(LTI,dim)

    T = zeros(dim.nx*(dim.N+1),dim.nx); % nx = 4
    for k=0:dim.N
        T(k*dim.nx+1:(k+1)*dim.nx,:)=LTI.A^k;
    end

    % Prediction matrix from input
    S=zeros(dim.nx*(dim.N+1),dim.nu*(dim.N));
        for k=1:dim.N
            for i=0:k-1
                S(k*dim.nx+1:(k+1)*dim.nx,i*dim.nu+1:(i+1)*dim.nu)=LTI.A^(k-1-i)*LTI.B;
            end
        end
end

%% Equality constraints for OTS

function eqconstraints=eqconstraintsgen(LTI,dim,dtilde)

eqconstraints.A=[eye(dim.nx)-LTI.A -LTI.B; LTI.C(2,:)  zeros(dim.ny-1, dim.nu)];
% eqconstraints.b=[zeros(dim.nx,1); LTI.yref-dtilde];
eqconstraints.b=[LTI.Bd * dtilde; LTI.yref(2,:) - LTI.Cd(2,:) * dtilde];

end


%% OTS calculation

function[xr,ur] = optimalss(dim, eqconstraints)

H = blkdiag(zeros(dim.nx), eye(dim.nu));
h = zeros(dim.nx+dim.nu,1);

options1 = optimoptions(@quadprog); 
options1.OptimalityTolerance=1e-2;
options1.ConstraintTolerance=1e-2;
% options1.Display='off';

[xur,fval,exitflag,output] = quadprog(H,h,[],[],eqconstraints.A, eqconstraints.b,[],[],[],options1);

xr=xur(1:dim.nx);
ur=xur(dim.nx+1:end);

end
