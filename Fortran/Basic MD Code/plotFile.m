clear all, close all, clc

KE = importdata('KE.dat');
PE = importdata('PE.dat');
TE = importdata('TE3.dat');
TE2 = importdata('TE2.dat');

hold on

% plot(KE(:,1), KE(:,2))
% plot(PE(:,1), PE(:,2))
plot(TE(:,1), TE(:,2))
plot(TE2(:,1), TE2(:,2))
hold off
xlabel('Time');
ylabel('Energy');
legend('dt = .005', 'dt = .0025');
