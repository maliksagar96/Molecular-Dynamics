clear all, close all, clc


xAxis2 = importdata('xAxis.dat');
yAxis2 = importdata('yAxis.dat');
yAxis2 = yAxis2/sum(yAxis2);

v = xAxis2;
v = v.*v;
%y = ((6/pi)^(3/2))*(4*pi*v).*exp(-(6*v)) ;
y  = 1.25*(4*pi*v) * ((1/(2*pi))^3/2).*(exp(-v/2));
hold on
plot(xAxis2, y, 'Linewidth', 2)
plot(xAxis2, yAxis2, 'Linewidth', 2)
legend('Theoretical Curve', 'Simulation Data')
xlabel('Speed')
ylabel('Probablity Density')
title('Maxwell Boltzmann Distribution')