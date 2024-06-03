clear all
close all
clc

% Parameter definiton THE STATES ARE: THETA, THETA_DOT, PHI_DOT
a21=  -10.4142;
a22= -0.1793;
a24= -0.0030;
a41= 10.4142;
a42= 0.1793;
a44= -0.9433;
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
LTI.Cd = [0;0]; % suppose No output disturbance
LTI.yref = [0; 0];

v = [10^4 10^2 10^2];
weight.Q = diag(v);   % weight on output
weight.R = 10;
[P, K, eigvals] = dare(LTI.A, LTI.B, weight.Q, weight.R, [], []);
weight.P = P;

% Compute terminal set using MPT3
model = LTISystem(sys_DT);
model.x.penalty = QuadFunction(weight.Q);
model.u.penalty = QuadFunction(weight.R);
Tset = model.LQRSet;

% Simulation horizon and input constraints
T_sim = 500;
umin = 0.11; % 0.1 A as rated current
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

% Output feedback DEFINITION
dim.nd = 1;
dime.nx= dim.nx + dim.nd; %state dimension
dime.nu= dim.nu;     %input dimension
dime.ny= dim.ny;     %output dimension
dime.N= dim.N;

% Extended system
LTI.Bd =  [1;0;0]; % -v(:,1) %[0.0035; 0.3490;0;0]  %LTI.B; %[1; 1; 1; 0];
LTIe.A=[LTI.A LTI.Bd; zeros(dim.nd,dim.nx) eye(dim.nd)];
LTIe.B=[LTI.B; zeros(dim.nd,dim.nu)];
LTIe.C=[LTI.C LTI.Cd];
LTIe.yref = LTI.yref;

% Cost definition extended 
weighte.Q = blkdiag(weight.Q, zeros(dim.nd));            % new weight Q ( no weight on disturbance)
weighte.R = weight.R;                                    
weighte.P = blkdiag(weight.P, zeros(dim.nd));            % new weight P (no weight on distrubace)

[Se, Te] = predmodgen(LTIe, dime);
% Cost definition extended MPC
[He, he] = costgen(Se, Te, weighte, dime);  

% MPC define variables
u_rec = zeros(dim.nu,T_sim);
total_iterations = T_sim;
o = waitbar(0, 'Progress');

% Inequality Constraints
x1_constraint = 0.3; % it is about 17 degrees
x2_constraint = 6; % 0.3/0.05 --> goes from 0.3 rad to 0 in one time step
x3_constraint = 9; % about 85 rad/sec (see datasheet)

up =[x1_constraint; x2_constraint; x3_constraint; inf]*100;
lb = -up;
bounds = [up; -lb];
eps = repmat(bounds, dime.N+1,1);
%extended
fe = [eye(dime.nx); -eye(dime.nx)];  % changed thiiiis
F = kron(eye(dime.N+1), fe);         % [(N+1)*2*dim.nx] x [(N+1)*dim.nx] 
eps_1 = [F*Te F*Se F];
A_in = F*Se;  

% Equality constraints (terminal)
Ae = Se(end-dime.nx+1:end,:);
be = Te(end-dime.nx+1:end,:);

% Initial condition
xehat(:,1)=[0.1; 0.0; 0; 0.0];
y(:,1)=LTIe.C * xehat(:,1);

% Observer
syse = ss(LTIe.A, LTIe.B, LTIe.C, [], Ts);
pol = pole(syse);
desired_observer_poles = [pol(1)*0.5,pol(2)*0.4,pol(3)*0.5,pol(4)*0.1];
L = place(LTIe.A', LTIe.C', desired_observer_poles)';
y_ref = repmat(LTI.yref(1), T_sim,1);


% MTP3
% mpc = MPCController(sys_DT);
model = LTISystem(sys_DT);
model.x.min = -[x1_constraint; x2_constraint; x3_constraint];
model.x.max = [x1_constraint; x2_constraint; x3_constraint];;
model.u.max = umax;
model.u.max = umin;
model.x.penalty = QuadFunction(weight.Q);
model.u.penalty = QuadFunction(weight.R);
Tset = model.LQRSet;
PN = model.LQRPenalty;
model.x.with('terminalSet');
model.x.terminalSet = Tset;
model.x.with('terminalPenalty');
model.x.terminalPenalty = PN;
N = 50;
ctrl = MPCController(model, N)
x0 = [0.1; 0; 0];
% u = mpc.evaluate(x0);
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

sys_new = ss(loop.system.A, loop.system.B, loop.system.C, loop.system.D, Ts);

t = 0:Ts:25;
u = zeros(size(t));
y = lsim(sys_new, u, t, [0.1 0.0 0.0]');

figure(2)
stairs(t, y(:,1))
hold on
stairs(t, y(:,2))
title('states')
legend('theta', 'theta dot', 'phi dot')

%% MPC
% Xn = importdata("Xn.mat");

for k = 1 : T_sim
    Constraint = [];
    dhat = xehat(4,k);
    x_0 = xehat(:,k);
    u_horizon = sdpvar(dime.nu * dime.N, 1);        
    eqconstraints = eqconstraintsgen(LTI, dim, dhat);
    
    [xr,ur]=optimalss(dim, eqconstraints);
    xer = [xr; dhat];    

    xer_e = repmat(xer, (dim.N+1),1);
    uer_e = repmat(ur, (dim.N),1);

    termdisA = [Tset.A zeros(size(Tset.A,1),1)];

    Constraint = [Constraint, 
                    umin <= u_horizon <= umax,
% %                     A_in * u_horizon <= eps - F*Te*x_0 + F*Se*uer_e - F*xer_e,
                    termdisA * (Ae * u_horizon + be * x_0) <= Tset.b,
                    ];                                       
%     eps_1 = [F*Te F*Se F];
    if k == 300
        break
    end
    % Objective function
    Objective = u_horizon' * He * u_horizon + (he*[x_0; xer; ur])' * u_horizon;     
%     Objective = [];
    % Optimize
    ops =  sdpsettings('verbose',0);
    sol = optimize(Constraint,Objective, ops);                      %solve the problem
    
    if sol.problem == 0
        u_horizon = value(u_horizon);                                  %assign the solution
        u_rec(k) = u_horizon(1);

%         Update states estimation
        y(:,k) = LTIe.C * xehat(:,k);
        xehat(:,k+1) = LTIe.A * xehat(:,k) + LTIe.B * u_rec(k);% + L * (y(:,k) - LTIe.C *  xehat(:,k));
%         y(:,k+1) = LTIe.C * xehat(:,k+1) + 0.1 * randn(dim.ny,1);
        dhat_plot(k) = xehat(4,k);
        clear u_horizon
    else
        close(o);
        disp('Problem solving MPC optimization');
        break
    end
   
    waitbar(k /total_iterations , o, sprintf('Progress: %d%%', round(k / total_iterations * 100)));
    
end

close(o);

%%
time = (0:k-1) * 0.05;

% e=y-kron(ones(1,T_sim),LTI.yref);
% figure(10000)
% plot(time,e(1,:))
% hold on
% plot(time,e(2,:))
% hold on 
% xlabel('k'), ylabel('Tracking error'), 
% grid on;
% legend('e_1','e_2', 'e_3');


figure(10)
stairs(y(1,:)),
hold on
stairs( y(2,:)),
xlabel('Time t [s]')
legend('y1', 'y2')

figure(20)
stairs( u_rec)
hold on
xlabel('Time t [s]'), 

% legend('d hat1')

figure(21)
stairs(xehat(1,:))
hold on
stairs( xehat(2,:))
stairs(xehat(3,:))


% stairs(time, kron(ones(1,T_sim),LTI.yref(1)))
xlabel('Time t [s]'), 


%% Define prediction model 

% Generation of prediction model -> Prediction matrix from initial state
function [S,T]=predmodgen(LTI,dim)

T = zeros(dim.nx*(dim.N+1),dim.nx); % nx = 4
    for k=0:dim.N
        T(k*dim.nx+1:(k+1)*dim.nx,:)=LTI.A^k;
    end

%Prediction matrix from input
S = zeros(dim.nx*(dim.N+1),dim.nu*(dim.N));
    for k=1:dim.N
        for i=0:k-1
            S(k*dim.nx+1:(k+1)*dim.nx,i*dim.nu+1:(i+1)*dim.nu)=LTI.A^(k-1-i)*LTI.B;
        end
    end

end

%% Cost generation extended system 
function [He,he, Qbar, Rbar]=costgen(S, T ,weight,dim)
    Qbar = blkdiag(kron(eye(dim.N),weight.Q),weight.P);
    Rbar =kron(eye(dim.N),weight.R);
    He = S'* Qbar* S+Rbar;   
    hx0 = S'*Qbar* T;
    hxref = - S'*Qbar*kron(ones(dim.N+1,1),eye(dim.nx));
    huref = -Rbar*kron(ones(dim.N,1),eye(dim.nu));
    he = [hx0 hxref huref];   
end

%% Equality constraints for OTS
function eqconstraints=eqconstraintsgen(LTI,dim,dtilde)

eqconstraints.A=[eye(dim.nx)-LTI.A -LTI.B; LTI.C(1,:)  zeros(dim.ny-1,dim.nu)];
%eqconstraints.b=[zeros(dim.nx,1); LTI.yref-dtilde]*;
eqconstraints.b=[LTI.Bd * dtilde; LTI.yref(1,:) - LTI.Cd(1) * dtilde];

end
%% OTS calculation

function[xr,ur] = optimalss(dim, eqconstraints)

H = blkdiag(zeros(dim.nx), eye(dim.nu));
h = zeros(dim.nx+dim.nu,1);

options1 = optimoptions(@quadprog); 
options1.OptimalityTolerance=1e-2;
options1.ConstraintTolerance=1e-2;

[xur,fval,exitflag,output] = quadprog(H,h,[],[],eqconstraints.A, eqconstraints.b,[],[],[],options1);

xr=xur(1:dim.nx);
ur=xur(dim.nx+1:end);
end