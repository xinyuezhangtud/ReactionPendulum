clear all
close all
clc
%% Compute sets
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
Ts = 0.05;
sys_DT = c2d(sys_CT,Ts,'zoh');

LTI.A=sys_DT.A; 
LTI.B=sys_DT.B;
LTI.C=sys_DT.C;
x0=[0;0;0;0];

weight.Q= eye(3);   % weight on output
weight.R=10;        % weight on input

Q =weight.Q;
R = weight.R;

Qx = eye(3);  % State weighting
Qn = Qx;                        % Terminal cost
R = 10;                          % Input weighting

A = LTI.A;
B = LTI.B;
C = LTI.C;
% Riccati Equation
[P,K,L] = idare(sys_DT.A, sys_DT.B, weight.Q, weight.R);
Qn = P;
K = -K;% Sign convention

nx=3; % state dimension
nu=1; % input dimension
ny = size(LTI.C,1);
N = 10;  % horizon

x1_constraint = 0.3; % it is about 17 degrees
x2_constraint = 6; % 0.3/0.05 --> goes from 0.3 rad to 0 in one time step
x3_constraint = 9; % about 85 rad/sec (see datasheet)

xlb = -[x1_constraint; x2_constraint; x3_constraint];
xub = [x1_constraint; x2_constraint; x3_constraint];

ulb = -0.11;
uub = 0.11;


% Calculate Xn and Xf (maximum LQR-invariant set) using normal penalty.
[Xn, V, Z] = findXn(A, B, K, N, xlb, xub, ulb, uub, 'lqr');
% K. V are vertices of these sets, and Z is the feasible input space
%% Projecting and Plotting in 2D

%plot 1 vs 2
dim = [1 2];
dimToPlot = [2 1];
Xn_projected = projectAllSets(Xn, dim);
plotProjectedSets(Xn_projected, dimToPlot);

%plot 1 vs 3
dim = [3 1];
dimToPlot = [2 1];
Xn_projected = projectAllSets(Xn, dim);
plotProjectedSets(Xn_projected, dimToPlot);

% %plot 1 vs 4
% dim = [3 2];
% dimToPlot = [4 1];
% Xn_projected = projectAllSets(Xn, dim);
% plotProjectedSets(Xn_projected, dimToPlot);

% %plot 2 vs 4
% dim = [3 1];
% dimToPlot = [4 2];
% Xn_projected = projectAllSets(Xn, dim);
% plotProjectedSets(Xn_projected, dimToPlot);

%plot 2 vs 3
dim = [2 3];
dimToPlot = [3 2];
Xn_projected = projectAllSets(Xn, dim);
plotProjectedSets(Xn_projected, dimToPlot);

% %plot 3 vs 4
% dim = [2 1];
% dimToPlot = [4 3];
% Xn_projected = projectAllSets(Xn, dim);
% plotProjectedSets(Xn_projected, dimToPlot);
%% Summary OF Functions 
    % findXn: Calculates the sequence of feasible sets (Xn) for the MPC controller, which are used to ensure the states remain within the constraints for each step in the control horizon.
    % 
    % hyperrectangle: Returns the halfspace representation of a hyperrectangle (i.e., a set defined by box constraints).
    % 
    % calcOinf: Calculates the maximum admissible set (also known as the positive invariant set) for the system under linear feedback control.
    % 
    % plotpoly: Given a polyhedron in a halfspace representation, this function calculates its vertices, sorts them, and plots the polygon.
    % 
    % solvelp_YALMIP: Solves a linear program using the YALMIP toolbox with the quadprog solver.
    % 
    % removeredundantcon: Removes redundant constraints from a set of linear inequalities.
    % 
    % halfspace2vertex: Converts a halfspace representation of a polyhedron to its vertex representation.
    % 
    % fullranksolve: Solves a system of linear equations if the system matrix has full rank.
    % 
    % polarsort: Sorts points in 2D based on their angle and radius from the origin.
    % 
    % computeX1: Computes the feasible set X1X1​ which is the set of all states that can be controlled into the terminal set XfXf​ in one step under the given constraints.
    % 
    % fouriermotzkin: A procedure used for projecting polyhedra (in halfspace representation) along a specified dimension. This is used in the computeX1 function for projection.
    % 
    % fmelim, removezerorows, and findinteriorpoint: These are helper functions for the Fourier-Motzkin elimination process and for finding a strictly feasible point inside a set of linear inequalities.
    %
    % projectAllSets
    %
    % projectSet
    %
    % plotProjectedSets







function [Xn, V, Z] = findXn(A, B, K, N, xlb, xub, ulb, uub, terminal)
% Calculates the various sets for the given system and parameters.
%
% Most inputs are self-explanatory. terminal should be either the string
% 'termeq' to use a terminal equality constraint, or 'lqr' to use the
% maximum LQR-invariant set.
%
% Output Xn is a cell array defining each of the \bbX_n sets. Xn{1} is the
% terminal set, Xn{2} is \bbX_1, etc. V is a cell array with each entry giving
% the extreme vertices of the corresponding Xn (as a 2 by N matrix). Z is a
% structure defining the feasible space.

% Define constraints for Z.
Nx = size(A, 1);
[Az, bz] = hyperrectangle([xlb; ulb], [xub; uub]);
Z = struct('G', Az(:,1:Nx), 'H', Az(:,(Nx + 1):end), 'psi', bz);

% Decide terminal constraint.
switch terminal
case 'termeq'
    % Equality constraint at the origin.
    Xf = [0; 0];
    
case 'lqr'
    % Build feasible region considering x \in X and Kx \in U.
    % Control input
    [A_U, b_U] = hyperrectangle(ulb, uub);
    A_lqr = A_U*K;
    b_lqr = b_U;
    
    % State input
    [A_X, b_X] = hyperrectangle(xlb, xub);
    Acon = [A_lqr; A_X];
    bcon = [b_lqr; b_X];

    % Use LQR-invariant set.
    Xf = struct();
    ApBK = A + B*K; % LQR evolution matrix.
    [Xf.A, Xf.b] = calcOinf(ApBK, Acon, bcon);
    [~, Xf.A, Xf.b] = removeredundantcon(Xf.A, Xf.b);

end

% Now do feasible sets computation.
figure();
hold('on');
colors = jet(N + 1);

Xn = cell(N + 1, 1);
Xn{1} = Xf;
names = cell(N + 1, 1);
names{1} = 'Xf';
V = cell(N + 1, 1);
for n = 1:(N + 1)
    % % Plot current set.
    if n > 1
        names{n} = sprintf('X%d', n - 1);
    end

    plotargs = {'-o', 'color', colors(n,:)};
    % 
    V{n} = plotpoly(Xn{n}, plotargs{:});
    if n == (N + 1)
        break
    end

%     Compute next set. Also need to prune constraints.
    nextXn = computeX1(Z, A, B, Xn{n});
    [~, nextXn.A, nextXn.b] = removeredundantcon(nextXn.A, nextXn.b);
    Xn{n + 1} = nextXn;
end

% legend(names{:});

end%function


function [A, b] = hyperrectangle(lb, ub)

% Returns halfspace representation of hyperrectangle with bounds lb and ub.

A = kron(eye(length(lb)), [1; -1]);
b = reshape([ub(:)'; -lb(:)'], 2*length(lb), 1);

% Any infinite or NaN bounds are ignored.
goodrows = ~isinf(b) & ~isnan(b);
A = A(goodrows,:);
b = b(goodrows);

end%function


function [AOinf, bOinf, tstar] = calcOinf(F, A, b)

% Calculates the maximum admissible set for x^+ = F*x subject to A*x <= b. Note
% that if F is unstable, this procedure will not work.
%
% tmax is the value of t to stop at, as an upper bound is not known a-priori.
% The default value is 100. If this bound is reached without termination, then
% tstar is set to inf.

% Arguments and sizes.

Nx = size(F, 1);

Nc = size(A, 1);

b = b(:);

% Start the algorithm.
Ft = eye(Nx);
AOinf = zeros(0, Nx);
bOinf = zeros(0, 1);
tstar = inf();

for t = 0:100
    % Add constraints for new time point.
    AOinf = [AOinf; A*Ft];
    bOinf = [bOinf; b];

    % Recalculate objective.
    Ft = F*Ft;
    fobj = A*Ft;

    % Maximize each component, stopping early if any bounds are violated.
    okay = true();
    for i = 1:Nc
        [obj, feas] = solvelp_YALMIP(fobj(i,:), AOinf, bOinf);
        if ~feas || obj > b(i)
            okay = false(); % N isn't high enough Need to keep going.
            continue
        end
    end

    % If everything was feasible, then we're done.
    if okay
        tstar = t;
        break
    end
end


end%function

function p = plotpoly(p, varargin)
% p = plotpoly(p, ...)
% p = plotpoly({A, b}, ...)
% p = plotpoly(struct('A', A, 'b', b), ...)
%
% Plots the polyhedron with vertices given in the 2 by N matrix p or given by
% the extreme points of A*x <= b. Any additional arguments are passed to plot.
%
% Returns the extreme points matrix p. If you only want this matrix, pass
% false as the only additional argument, and the plot will not be made.
if nargin() < 1
    error('Argument p is required!');
end

% First, find the vertices if {A, b} given.
if iscell(p)
    p = halfspace2vertex(p{1}, p{2})';
elseif isstruct(p)
    p = halfspace2vertex(p.A, p.b)';
end

% Next, sort the vertices.
ptilde = bsxfun(@rdivide, bsxfun(@plus, p, -mean(p, 2)), std(p, 0, 2));
x = ptilde(1,:);
y = ptilde(2,:);
[th, r] = cart2pol(x, y);
thneg = (th < 0);
th(thneg) = th(thneg) + 2*pi(); % Makes theta in [0, 2*pi] instaed of [-pi, pi].
[~, s] = sortrows([th', r']); % Sort on theta then r.
p = p(:,s);

% Duplicate first data point to give closed cycle.
p = p(:,[1:end,1]);

% Now plot.
if isempty(varargin) || ~isequal(varargin{1}, false())
    plot(p(1,:), p(2,:), varargin{:}, 'LineWidth', 2);
end

end%function


%%
function [obj, feas] = solvelp_YALMIP(f, A, b)

    options = sdpsettings('verbose',0,'solver','quadprog');
    
    x_ = sdpvar(length(f),1);                % define optimization variable

    Constraint = [A*x_ <= b];                  %define constraints

    Objective = f*x_;  %define cost function

    diagnostic = optimize(Constraint,Objective, options);  %solve the problem
    x_=value(x_);                  %assign the solution to uopt
    obj = -f*x_;
    feas = (diagnostic.problem == 0);

end


%%
function [nr, Anr, bnr, h, x0] = removeredundantcon(A, b, x0, tol, qhullopts)
    % [nr, Anr, bnr, h, x0] = removeredundantcon(A, b, [x0], [tol], qhullopts)
    %
    % Finds the non-redundant constraints for the polyhedron Ax <= b. nr is a
    % column vector with the non-redundant rows. If requested, Anr and bnr are the
    % non-redundant parts of A and b.
    %
    % If x0 is supplied, it must be in the strict interior of A. Otherwise an error
    % is thrown. Specifying a valid x0 will speed up the function.
    %
    % tol is used to decide how much on the interior we need to be. If not supplied,
    % the default value is 1e-8*max(1, norm(b)/N). Note that this may be too large
    % if b is poorly scaled.
    %
    % qhullopts is a string of options to pass to qhull. Defaults are described in
    % the documentation for convhulln (see `help convhulln`).
    %
    % h is the output of convhulln, which actually describes the facets, and x0
    % is a point on the interior of the polyhedron. Note that if the input
    % polyhedron is unbounded, h may have zeros in some entries corresponding to the
    % row of all zeros that needs to be added for the method to work.
    %
    % Note that this requires finding convex hull in N + 1 dimensions, where N
    % is the number of columns of A. Thus, this will get very slow if A has a lot
    % of columns.
    narginchk(2, 5);
    
    % Force b to column vector and check sizes.
    b = b(:);
    if isrow(b)
        b = b';
    end
    if size(A, 1) ~= length(b)
        error('A and b must have the same number of rows!');
    end
    if nargin() < 3
        x0 = [];
    end
    if nargin() < 4 || isempty(tol)
        tol = 1e-8*max(1, norm(b)/length(b));
    elseif tol <= 0
        error('tol must be strictly positive!')
    end
    if nargin() < 5
        if size(A, 2) <= 4
            qhullopts = {'Qt'};
        else
            qhullopts = {'Qt', 'Qx'};
        end
    end
    
    
    % Save copies before things get messed up.
    Anr = A;
    bnr = b;
    
    % First, get rid of any rows of A that are all zero.
    Anorms = max(abs(A), [], 2);
    badrows = (Anorms == 0);
    if any(b(badrows) < 0)
        error('A has infeasible trivial rows.')
    end
    A(badrows, :) = [];
    b(badrows) = [];
    goodrows = [0; find(~badrows)]; % Add zero for all zero row that gets added.
    
    % Need to find a point in the interior of the polyhedron.
    if isempty(x0)
        if all(b > 0)
            % If b is strictly positive, we know the origin works.
            x0 = zeros(size(A, 2), 1);
        else
            error('Must supply an interior point!');
        end
    else
        x0 = x0(:);
        if isrow(x0)
            x0 = x0';
        end
        if length(x0) ~= size(A,2)
            error('x0 must have as many entries as A has columns.')
        end
        if any(A*x0 >= b - tol)
            error('x0 is not in the strict interior of Ax <= b!')
        end
    end
    
    % Now, project the rows of P and find the convex hull.
    btilde = b - A*x0;
    if any(btilde <= 0)
        warning('Shifted b is not strictly positive. convhull will likely fail.')
    end
    Atilde = [zeros(1, size(A, 2)); bsxfun(@rdivide, A, btilde)];
    h = convhulln(Atilde, qhullopts);
    u = unique(h(:));
    
    nr = goodrows(u);
    if nr(1) == 0
        nr(1) = [];
    end
    h = goodrows(h);
    
    % Finally, grab the appropriate rows for Anr and bnr.
    Anr = Anr(nr, :);
    bnr = bnr(nr);

end%function

%%
function [V, nr] = halfspace2vertex(A, b, x0)
% [V, nr] = halfspace2vertex(A, b, [x0])
%
% Finds extreme points of polyhedron A*x <= b. Note that the polyhedron must
% have a point strictly on its interior.
%
% If provided, x0 must be a point on the interior of the polyhedron. If it is
% not given, one is found by solving a linear program.
%
% V is returned as an N by 2 matrix with each row giving an extreme point.
%
% Second output nr is a list of the non-redundant constraints of the polytope.

% Check inputs.
narginchk(2, 3);
Nc = size(A, 1);
Nx = size(A, 2);
if ~isvector(b) || length(b) ~= Nc
    error('b is the incorrect size!');
end
b = b(:); % Make sure b is a column vector.

% Sort out interior point.
if nargin() < 3
    if all(b > 0)
        % The origin is on the interior. Can rescale rows so that b = 1.
        x0 = zeros(Nx, 1);
        A = bsxfun(@rdivide, A, b);
        b = ones(size(b));
    else
        x0 = findinteriorpoint(A, b);
    end
elseif ~isvector(x0) || length(x0) ~= Nx
    error('Invalid size for x0!');
end
x0 = x0(:); % Make sure x0 is a column vector.

% Get non-redundant constraints from A and b.
[nr, ~, ~, k] = removeredundantcon(A, b, x0);

% Now loop through facets to find vertices.
V = zeros(size(k, 1), Nx);
keep = true(size(k, 1), 1);
for ix = 1:size(k, 1)
    F = A(k(ix,:),:);
    g = b(k(ix,:));
    [keep(ix), V(ix,:)] = fullranksolve(F, g);
end

V = V(keep,:);
[~, u] = unique(round(V*1e6), 'rows');
V = V(u,:);

% If in 2D, sort the vertices.
if Nx == 2
    V = polarsort(V);
end

end%function
%%
function [fullrank, x] = fullranksolve(A, b);
    % Checks whether the system is full rank and if so, solves it. If it is not
    % full rank, a vector of NaNs are returned.
    Nx = size(A, 1);
    [U, S, V] = svd(A);
    s = diag(S);
    tol = Nx*eps(s(1)); % Rank tolerance.
    fullrank = all(s >= tol);
    if fullrank
        x = V*diag(1./s)*U'*b;
    else
        x = NaN(Nx, 1);
    end
end%function
%%
function [p, s] = polarsort(p)
    % [p, s] = polarsort(p)
    %
    % Sorts the [n by 2] matrix p so that the points are counter-clockwise starting
    % at the theta = 0 axis. For ties in theta, sorts on radius.
    x = p(:,1);
    y = p(:,2);
    x = (x - mean(x))/std(x); % Normalize so that the origin is at the center.
    y = (y - mean(y))/std(y);
    [th, r] = cart2pol(x, y);
    [~, s] = sortrows([th, r]); % Sort on theta then r.
    p = p(s,:);
end%function
%%
function X1 = computeX1(Z, A, B, Xf)
% X1 = computeX1(Z, [A], [B], [Xf])
%
% Computes the feasible set X_1 for the system x^+ = Ax + Bu subject to
% constraints Gx + Hu <= psi and x^+ \in Xf.
%
% Z must be a struct with fields G, H, and psi.
%
% A and B are only necessary if the terminal constraint is given.
%
% Xf can be either a struct with fields A and b to define a polytope, or a
% single vector to define a point constraint. If not provided, it is assumed
% that Xf is the entire space.
%
% X1 is returned as a struct with fields A and b defining a set of inequality
% constraints.

% narginchk(minArgs,maxArgs) validates the number of input arguments in the call to the currently executing function. 
% narginchk throws an error if the number of inputs specified in the call is fewer than minArgs or greater than maxArgs. If the number of inputs is between minArgs and maxArgs (inclusive), then narginchk does nothing.

% Check arguments.
narginchk(1, 4);
if ~isstruct(Z) || ~all(isfield(Z, {'G', 'H', 'psi'}))
    error('Invalid input for Z!');
end
Nx = size(Z.G, 2);
Nu = size(Z.H, 2);
if nargin() >= 3
    sys = struct('A', A, 'B', B); % Save these.
end

% Preallocate constraint matrices.
A = [Z.G, Z.H];
b = Z.psi;
Aeq = zeros(0, Nx + Nu);
beq = zeros(0, 1);

% Do some more stuff if there is a terminal constraint.
if nargin() >= 4 && ~isempty(Xf)
    if ~isequal([size(sys.A), size(sys.B)], [Nx, Nx, Nx, Nu])
        error('Incorrect size for A or B!');
    end
    if isstruct(Xf)
        % Polyhedron.
        if ~all(isfield(Xf, {'A', 'b'}))
            error('Struct Xf must have fields A and b!');
        end
        A = [A; Xf.A*[sys.A, sys.B]];
        b = [b; Xf.b];
    elseif isvector(Xf) && length(Xf) == Nx
        % Terminal equality constraint.
        Aeq = [sys.A, sys.B];
        beq = Xf;
    else
        error('Invalid input for Xf!');
    end
end

% Now do the projection step.
for i = 1:Nu
    [A, b, Aeq, beq] = fouriermotzkin(A, b, Aeq, beq);
end
X1 = struct('A', A, 'b', b);

end%function

%%
function [A, b, Aeq, beq] = fouriermotzkin(A, b, Aeq, beq, ielim)
% [A, b, Aeq, beq] = fouriermotzkin(A, b, [Aeq], [beq], [ielim])
%
% Perform one step of Fourier-Motzkin elimination for the set defined by
%
%   {x \in R^n : A*x <= b, Aeq*x = beq}
%
% Note that the inequality constrants can be degenerate, although degeneracy
% will create significantly more constraints compared to explicit equality
% constraints.
%
% Optional argument ielim decides which column to eliminate. Default is the last
% column.

% Check arguments and get sizes.
narginchk(2, 5);
Nlt = size(A, 1);
if ~isvector(b) || length(b) ~= Nlt
    error('Invalid size for b!');
end
Nx = size(A, 2);
if nargin() < 3 || isempty(Aeq)
    Aeq = zeros(0, Nx);
    beq = zeros(0, 1);
    Neq = 0;
elseif nargin() < 4
    error('beq is required if Aeq is given!');
else
    if size(Aeq, 2) ~= Nx
        error('Aeq has the wrong number of columns!');
    end
    Neq = size(Aeq, 1);
    if ~isvector(beq) || length(beq) ~= Neq
        error('Invalid size for beq!');
    end
end
if nargin() < 5
   Nx
    ielim = Nx;
elseif ~isscalar(ielim) || ielim <= 0 || ielim > Nx || round(ielim) ~= ielim
    error('ielim must be a scalar positive integer less than Nx!');
end

% Now decide what to do. First, look for an equality constraint with the
% variable of interest.
[pivval, pivrow] = max(abs(Aeq(:,ielim)));
if isempty(pivrow) || pivval == 0
    % No suitable equality constraint found. Need to change inequality
    % constraints.
    [A, b] = fmelim(A, b, ielim);
    Aeq(:,ielim) = [];
else
    % An equality constraint is available. Handle that.
    a = Aeq(pivrow,:);
    c = a(ielim);

    % Make the pivots.
    Aeq = Aeq - Aeq(:,ielim)*a/c;
    beq = beq - Aeq(:,ielim)*beq(pivrow)/c;

    A = A - A(:,ielim)*a/c;
    b = b - A(:,ielim)*beq(pivrow)/c;

    % Get rid of the appropriate rows and columns.
    A(:,ielim) = [];
    Aeq(:,ielim) = [];
    Aeq(pivrow,:) = [];
    beq(pivrow) = [];
end

% Zap any rows that are all zeros.
[A, b] = removezerorows(A, b);
[Aeq, beq] = removezerorows(Aeq, beq);

end%function


% *****************************************************************************
% Helper Functions
% *****************************************************************************

function [Ae, be] = fmelim(A, b, ielim)
    % Performs one step of Fourier-Motzkin elimination for inequality
    % constraints.
    c = A(:,ielim);
    I0 = find(c == 0);
    Ip = find(c > 0);
    Im = find(c < 0);

    Nx = size(A, 2);
    Ne = length(I0) + length(Ip)*length(Im);

    E = [A(:,1:ielim - 1), A(:,ielim + 1:end), b];
    Ee = [E(I0,:); kron(c(Ip), E(Im,:)) - kron(E(Ip,:), c(Im))];

    Ae = Ee(:,1:end-1);
    be = Ee(:,end);
end%function
%%
function [A, b] = removezerorows(A, b)
    % Removes rows of A that are all zeros.
    keeprows = any(A, 2);
    A = A(keeprows,:);
    b = b(keeprows);
end%function
%%
function [x0, okay, feas, margin] = findinteriorpoint(A, b, Aeq, beq, tol, maxinside)
% [x0, okay, feas, margin] = findinteriorpoint(A, b, [Aeq], [beq], [tol])
%
% Find a strictly feasible point x0 so that A*x0 <= b - tol. If no such point
% can be found, okay is set to False.
%
% If there is at least a feasible point (but not necessarily on the interior),
% then feas is true, and x0 gives that point. If both okay and feas are false,
% then x0 is meaningless.
%
% margin gives the value e such that A*x0 <= b - e.
%
% The origin is always checked first. If it does not work, an LP is solved
% to find a valid point.
%
% tol is used to decide how much on the interior we need to be. If not
% supplied, the default value is 1e-8*max(1, norm(b)/N). Note that this may
% be too large if b is poorly scaled.
if nargin() < 2
    error('A and b are mandatory.')
elseif nargin() < 5
    tol = 1e-8*max(1, norm(b)/length(b));
end
if nargin() < 6
    maxinside = 1;
else
    maxinside = max(tol, maxinside);
end
[m, n] = size(A);
if nargin() < 4 || isempty(Aeq)
    Aeq = zeros(0, n);
    beq = [];
end
meq = size(Aeq, 1);
okay = false();

% Check whether the origin is on the inside.
if all(abs(beq) < tol) && all(b > tol)
    x0 = zeros(n, 1);
    okay = true();
    feas = true();
    margin = min(b);
end

% Try to use fminsearch if there are no equality constraints. Doesn't work
% well if the number of dimensions is very high, so we cap it at 10.
if ~okay && meq == 0 && m <= 10
    options = optimset('display', 'off');
    [x0, maxr] = fminsearch(@(x) max(max(A*x - b), -1e5*tol), A\b, options);
    okay = (maxr < -tol);
    feas = okay;
    margin = -maxr;
end

% Solve LP otherwise.
if ~okay
    c = [zeros(n, 1); -1];
    AA = [A, ones(m, 1)];
    AAeq = [Aeq, zeros(meq,1)];
    lb = [-inf(n, 1); 0];
    ub = [inf(n, 1); maxinside];
    if isOctave()
        ctype = [repmat('U', m, 1); repmat('S', meq, 1)];
        [xtilde, ~, err, extra] = glpk(c, [AA; AAeq], [b; beq], lb, ub,  ...
            ctype, repmat('C', n + 1, 1), 1, struct('msglev', 0));
        okay = (err == 0 && extra.status == 5);
    else
        options = optimoptions('linprog', 'display', 'off', 'algorithm', 'dual-simplex');
        [xtilde, ~, exitflag] = linprog(c, AA, b, AAeq, beq, lb, ub, [], options);
        okay = (exitflag == 1);
    end
    if isempty(xtilde)
        margin = -inf();
    else
        margin = xtilde(end); % Amount constraints could be tightened.
    end
    okay = okay && (margin >= tol);
    if isempty(xtilde)
        x0 = zeros(n, 1);
        okay = false();
    else
        x0 = xtilde(1:end-1);
    end

    % Check feasibility of x0.
    feas = all(A*x0 - b < tol);
    if feas && ~isempty(Aeq)
        feas = all(abs(Aeq*x0 - beq) < tol);
    end
    okay = okay && feas;
end
end%function

function Xn_projected = projectAllSets(Xn, dim)
    % Xn is a cell array where each cell contains a struct with fields 'A' and 'b'
    % This function projects each set in Xn onto the 1st and 2nd dimensions

    % intialize
    Xn_projected = cell(size(Xn));
    
    for i = 1:length(Xn)
        % Extract the current set's A and b
        A = Xn{i}.A;
        b = Xn{i}.b;
        
        % disp(['Size of A: ', mat2str(size(A))]);
        % disp(['Length of b: ', num2str(length(b))]);

        % Project the current set onto 2 dim
        [A_proj, b_proj] = projectSet(A, b, dim);

        % Store the projected set in the new cell array
        Xn_projected{i} = struct('A', A_proj, 'b', b_proj);
    end
end

function [A_proj, b_proj] = projectSet(A, b, dim)
    % empty equality constraints
    Aeq = [];
    beq = [];

    [A, b, Aeq, beq] = fouriermotzkin(A, b, Aeq, beq, dim(1));
    
    %Remove redundant constraints
    [~, A, b, ~] = removeredundantcon(A, b);

    [A, b, Aeq, beq] = fouriermotzkin(A, b, Aeq, beq, dim(2));

    %Remove redundant constraints
    [~, A, b, ~] = removeredundantcon(A, b);
    

    % The projected A matrix and b vector, now correctly reduced to the target dimensions
    A_proj = A;
    b_proj = b;
end

function plotProjectedSets(Xn_projected, dimToPlot)
    % This function plots all the projected sets contained in Xn_projected
    % with an emphasis on specific dimensions for each set.

    figure; % Create a new figure
    hold on; % Hold on to plot all the sets on the same figure
    colors = jet(length(Xn_projected)); % Generate a color map for the sets

    % Loop through each projected set and plot it with specific styling
    for i = 1:length(Xn_projected)
        p = Xn_projected{i}; % p is now a struct with fields 'A' and 'b'
        
            plotargs = {'-o', 'color', colors(i,:)};
            plotpoly(p, plotargs{:}); % Plot using the generic poly plotting function
        
    end

    hold off; 
    xlabel(sprintf('x%d', dimToPlot(2))); 
    ylabel(sprintf('x%d', dimToPlot(1)));
    title(sprintf('Projected Sets onto Dimension x%d vs x%d', dimToPlot(1), dimToPlot(2)));
    legends = ['Xf', arrayfun(@(n) sprintf('X%d', n), 1:(length(Xn_projected)-1), 'UniformOutput', false)];
    legend(legends); 
    grid on
end

