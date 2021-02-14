clear all, close all, clc

x1 = importdata('xAxis.dat');
y1 = importdata('yAxis.dat');
y1 = y1/sum(y1);
x2 = importdata('xAxis4.dat');
y2 = importdata('yAxis4.dat');
y2 = y2/sum(y2);

x3 = importdata('xAxis5.dat');
y3 = importdata('yAxis5.dat');
y3 = y3/sum(y3);

hold on
plot(x1, y1);
plot(x2, y2);
plot(x3, y3);

xlabel('Speed');
ylabel('Probability density');
title('Maxwell Boltzmann');
legend('KbT = 1' , 'KbT = 4', 'KbT = 9');