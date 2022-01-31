
close all;
clear all;
clc;

% Define fontsize
fs = 28;

% Import digitized data
exp_lo = readmatrix('./plate_102012_data/exp_lower_q.csv');
exp_hi = readmatrix('./plate_102012_data/exp_higher_q.csv');
dowell = readmatrix('./plate_102012_data/102012_theory.csv');

% Construct ROM results
rom = [2, 1.46;   % 1.35 stable, trying 1.4
       3, 2.18;
       4, 2.91;
       5, 3.64];

flanto = [2, 2.12;
          3, 3.46;
          4, 4.74;
          5, 6.00];

% Plot results
figure;
lg1 = scatter(exp_lo(:,1), exp_lo(:,2), 'k'); hold on;
lg2 = scatter(exp_hi(:,1), exp_hi(:,2), 'k'); grid on;
lg3 = plot(dowell(:,1), dowell(:,2), 'b'); 
lg3.LineWidth = 2;
lg4 = plot(rom(:,1), rom(:,2), '-*r');
lg4.LineWidth = 2;
lg5 = plot(flanto(:,1), flanto(:,2), '-*g');
lg5.LineWidth = 2;
legend([lg1 lg3 lg4, lg5], 'Experimental', 'Dowell \& Voss', 'Euler-Lagrange', 'FLANTO2D', 'interpreter', 'latex', 'fontsize', 20, 'location', 'northwest');
xlim([1,6]); ylim([0,7]);
xlabel('$M$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$q_{\infty}$ [psi]', 'interpreter', 'latex', 'fontsize', fs);
title('Comparison of Experimental and Theoretical Results, 10-20-12 Plate', 'interpreter', 'latex', 'fontsize', fs);


% Plot results to be used in Turbo Expo presentation
figure;
lg1 = scatter(exp_lo(:,1), exp_lo(:,2), 'k'); hold on;
lg2 = scatter(exp_hi(:,1), exp_hi(:,2), 'k'); grid on;
lg3 = plot(dowell(:,1), dowell(:,2), 'b'); 
lg3.LineWidth = 2;
lg4 = plot(rom(:,1), rom(:,2), '-*r');
lg4.LineWidth = 2;
legend([lg1 lg3 lg4], 'Dowell \& Voss Experimental', 'Dowell \& Voss Theoretical', 'Euler-Lagrange', 'interpreter', 'latex', 'fontsize', 20, 'location', 'northwest');
xlim([1,6]); ylim([0,7]);
xlabel('$M$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$q_{\infty}$ [psi]', 'interpreter', 'latex', 'fontsize', fs);
title('Comparison of Experimental and Theoretical Panel Flutter Results', 'interpreter', 'latex', 'fontsize', fs);


% Plot results to be used in prelim proposal
figure;
lg1 = scatter(exp_lo(:,1), exp_lo(:,2), 45, 'k', 'filled'); hold on;
lg2 = scatter(exp_hi(:,1), exp_hi(:,2), 45, 'k', 'filled'); grid on;
lg3 = plot(dowell(:,1), dowell(:,2), 'b'); 
lg3.LineWidth = 2;
lg4 = plot(rom(:,1), rom(:,2), '-*r');
lg4.LineWidth = 2;
legend([lg1 lg3 lg4], 'Dowell \& Voss Experimental', 'Dowell \& Voss Theoretical', 'Euler-Lagrange ROM', 'interpreter', 'latex', 'fontsize', 26, 'location', 'northwest');
xlim([1,6]); ylim([0,7]);
a = get(gca, 'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',18);
xlabel('Mach Number', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Dynamic Pressure [psi]', 'interpreter', 'latex', 'fontsize', fs);
title('Comparison of Flutter Stability Boundary, 10-20-12 Plate', 'interpreter', 'latex', 'fontsize', fs);


% Plot results to be used in J. Turbomach. manuscript
figure;
lg1 = scatter(exp_lo(:,1), exp_lo(:,2) .* 6894.757, 45, 'k', 'filled'); hold on;  % convert from psi to Pa
lg2 = scatter(exp_hi(:,1), exp_hi(:,2) .* 6894.757, 45, 'k', 'filled'); grid on;  % convert from psi to Pa
lg3 = plot(dowell(:,1), dowell(:,2) .* 6894.757, 'b');  % convert from psi to Pa
lg3.LineWidth = 2;
lg4 = plot(rom(:,1), rom(:,2) .* 6894.757, '-*r');  % convert from psi to Pa
lg4.LineWidth = 2;
legend([lg1 lg3 lg4], 'Dowell \& Voss Experimental', 'Dowell \& Voss Theoretical', 'Euler-Lagrange ROM', 'interpreter', 'latex', 'fontsize', 26, 'location', 'northwest');
xlim([1,6]); %ylim([0,7]);
a = get(gca, 'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',18);
xlabel('Mach Number', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Dynamic Pressure [Pa]', 'interpreter', 'latex', 'fontsize', fs);
title('Comparison of Flutter Stability Boundary, 10-20-12 Plate', 'interpreter', 'latex', 'fontsize', fs);

