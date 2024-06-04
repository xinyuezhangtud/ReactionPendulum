clear all
close all
clc

% load("deg_71.mat")
% load("deg_145.mat")
% load("pi_rad.mat")
% load("sixtythree_deg.mat")
% load("one_twentytwo_deg.mat")
% load("zero_rad.mat")

y = [ pi 0 -2.04204 2.1293 4.38078 2.53073]';
x = [3.14473 0 -2.18058 2.12434 4.37875 2.57171]';

rad = y;
deg = rad * 180 / pi 
c = polyfit(x,y,1);
slope = c(1);
offset = c(2);

yfit = slope*x + offset;
% yfit2 = slope*x - slope*offset;

figure(1)
scatter(x, y)
grid on
% plot(y)
figure(2)
scatter(x,y)
hold on
plot(x, yfit)
grid on
legend('Angle values', 'Linear fit', 'Interpreter', 'latex')
xlabel('Sensor values $$[rad]$$', 'Interpreter', 'latex')
ylabel('Manually veasured values $$[rad]$$', 'Interpreter', 'latex')

