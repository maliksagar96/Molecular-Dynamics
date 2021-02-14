clear all, close all, clc


KE = importdata('KE.dat');
PE = importdata('PE.dat');
TE = importdata('TE.dat');

hold on

plot(KE(:,1), KE(:,2))
plot(PE(:,1), PE(:,2))
plot(TE(:,1), TE(:,2))
hold off
xlabel('Time');
ylabel('Energy');
legend('KE', 'PE', 'TE');
title('Thermostat')
