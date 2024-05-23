load('square.mat')

phi_dot_with_outlier = [phi_dot.Data];

phi_dot_removal = filloutliers(phi_dot_with_outlier(1:100), 'linear');

phi_dot = vertcat(phi_dot_removal,phi_dot_with_outlier(101:1001));

figure(1)
plot(phi_dot_with_outlier)

figure(2)
plot(phi_dot)