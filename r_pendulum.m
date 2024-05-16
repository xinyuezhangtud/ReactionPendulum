function [A,B,C,D, K, x0] = r_pendulum(par, T)

A = [0 1 0; -par(1)/par(2) 0 0; par(1)/par(2) 0 0];

B = [0; -par(4)/par(2); 
    par(4)*((1/par(2))+(1/par(3)))];

C = [1 0  0;
     0 0  1];
K = eye(3)*0; %[0;0;0];
D = [0; 0];
x0 = [par(5) 0 0]';
end


