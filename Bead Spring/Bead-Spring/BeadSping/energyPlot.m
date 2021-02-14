clear all, close all, clc

energy = importdata('Energy.dat');


hold on

% plot(energy(:,1), energy(:,2), 'linewidth', 1.5);
% plot(energy(:,1), energy(:,3), 'linewidth', 1.5);
plot(energy(:,1), energy(:,4), 'linewidth', 1.5);
% 
% hold off
xlabel('Time');
ylabel('Energy');

% legend('PE', 'KE', 'TE');