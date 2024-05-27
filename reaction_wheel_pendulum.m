
function [dx, y] = reaction_wheel_pendulum(t,x,u, mgl, J_theta, b_theta, b_phi, k, J_phi, varargin)

    % Output equation
    y = [x(1);
         x(4)];

    % State equation
    dx = [x(2);
        -(mgl*cos(x(1))/J_theta) -(b_theta/J_theta)*x(2) + (b_phi/J_theta)*x(4) -(k/J_theta)*u;
        x(3);
         (mgl*cos(x(1))/J_theta) (b_theta/J_theta)*x(2) - (b_phi/J_theta)*x(4) + ((k/J_theta)+ (k/J_phi))*u;
         ];
end