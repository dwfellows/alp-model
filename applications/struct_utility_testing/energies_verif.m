%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to attempt to verify energies_func
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all;
clear all;
clc;

addpath('../../struct_energies_utils');
fs = 26;

%% Numerical data

% Read in data for circular verification cases
% Structural steel material properties
E = 2e11; % Pa
nu = 0.3; % dimensionless
rho = 7850; %kg/m^3
mat_props = [E, nu, rho];

% Circular plate properties
R = 0.1;
H = 0.0001;
h = H/2;

vol = pi*(R^2)*H;

% Coarse circular verification
% Read in data
coarse_node_elem = readmatrix('./circular_coarse/coarse.csv');
coarse_ux_1 = readmatrix('./circular_coarse/coarse_ux1.txt');
coarse_uy_1 = readmatrix('./circular_coarse/coarse_uy1.txt');
coarse_uz_1 = readmatrix('./circular_coarse/coarse_uz1.txt');
coarse_ux_2 = readmatrix('./circular_coarse/coarse_ux2.txt');
coarse_uy_2 = readmatrix('./circular_coarse/coarse_uy2.txt');
coarse_uz_2 = readmatrix('./circular_coarse/coarse_uz2.txt');


% Compute avg tetrahedron edge length
coarse_elem_num = max(size(coarse_node_elem));
coarse_elem_avg_vol = vol / coarse_elem_num;
coarse_elem_h = (coarse_elem_avg_vol*6*sqrt(2))^(1/3);

% Simplify data
coarse_node_locs = coarse_ux_1(:,2:4);
coarse_ux_1 = coarse_ux_1(:,5);
coarse_uy_1 = coarse_uy_1(:,5);
coarse_uz_1 = coarse_uz_1(:,5);
coarse_ux_2 = coarse_ux_2(:,5);
coarse_uy_2 = coarse_uy_2(:,5);
coarse_uz_2 = coarse_uz_2(:,5);

% Normalize such that max displacement is equal to 1
coarse_u_1 = [coarse_ux_1, coarse_uy_1, coarse_uz_1];
coarse_u_1_norm = vecnorm(coarse_u_1, 2, 2);
coarse_ux_1 = coarse_ux_1 ./ max(coarse_u_1_norm);
coarse_uy_1 = coarse_uy_1 ./ max(coarse_u_1_norm);
coarse_uz_1 = coarse_uz_1 ./ max(coarse_u_1_norm);
coarse_psi_1 = [coarse_ux_1, coarse_uy_1, coarse_uz_1];

coarse_u_2 = [coarse_ux_2, coarse_uy_2, coarse_uz_2];
coarse_u_2_norm = vecnorm(coarse_u_2, 2, 2);
coarse_ux_2 = coarse_ux_2 ./ max(coarse_u_2_norm);
coarse_uy_2 = coarse_uy_2 ./ max(coarse_u_2_norm);
coarse_uz_2 = coarse_uz_2 ./ max(coarse_u_2_norm);
coarse_psi_2 = [coarse_ux_2, coarse_uy_2, coarse_uz_2];


% Mid circular verification
% Read in data
mid_node_elem = readmatrix('./circular_mid/mid.csv');
mid_ux_1 = readmatrix('./circular_mid/mid_ux1.txt');
mid_uy_1 = readmatrix('./circular_mid/mid_uy1.txt');
mid_uz_1 = readmatrix('./circular_mid/mid_uz1.txt');
mid_ux_2 = readmatrix('./circular_mid/mid_ux2.txt');
mid_uy_2 = readmatrix('./circular_mid/mid_uy2.txt');
mid_uz_2 = readmatrix('./circular_mid/mid_uz2.txt');

% Compute avg tetrahedron edge length
mid_elem_num = max(size(mid_node_elem));
mid_elem_avg_vol = vol / mid_elem_num;
mid_elem_h = (mid_elem_avg_vol*6*sqrt(2))^(1/3);

% Simplify data
mid_node_locs = mid_ux_1(:,2:4);
mid_ux_1 = mid_ux_1(:,5);
mid_uy_1 = mid_uy_1(:,5);
mid_uz_1 = mid_uz_1(:,5);
mid_ux_2 = mid_ux_2(:,5);
mid_uy_2 = mid_uy_2(:,5);
mid_uz_2 = mid_uz_2(:,5);

% Normalize such that max displacement is equal to 1
mid_u_1 = [mid_ux_1, mid_uy_1, mid_uz_1];
mid_u_1_norm = vecnorm(mid_u_1, 2, 2);
mid_ux_1 = mid_ux_1 ./ max(mid_u_1_norm);
mid_uy_1 = mid_uy_1 ./ max(mid_u_1_norm);
mid_uz_1 = mid_uz_1 ./ max(mid_u_1_norm);
mid_psi_1 = [mid_ux_1, mid_uy_1, mid_uz_1];

mid_u_2 = [mid_ux_2, mid_uy_2, mid_uz_2];
mid_u_2_norm = vecnorm(mid_u_2, 2, 2);
mid_ux_2 = mid_ux_2 ./ max(mid_u_2_norm);
mid_uy_2 = mid_uy_2 ./ max(mid_u_2_norm);
mid_uz_2 = mid_uz_2 ./ max(mid_u_2_norm);
mid_psi_2 = [mid_ux_2, mid_uy_2, mid_uz_2];


% Fine circular verification
% Read in data
fine_node_elem = readmatrix('./circular_fine/fine.csv');
fine_ux_1 = readmatrix('./circular_fine/fine_ux1.txt');
fine_uy_1 = readmatrix('./circular_fine/fine_uy1.txt');
fine_uz_1 = readmatrix('./circular_fine/fine_uz1.txt');
fine_ux_2 = readmatrix('./circular_fine/fine_ux2.txt');
fine_uy_2 = readmatrix('./circular_fine/fine_uy2.txt');
fine_uz_2 = readmatrix('./circular_fine/fine_uz2.txt');

% Compute avg tetrahedron edge length
fine_elem_num = max(size(fine_node_elem));
fine_elem_avg_vol = vol / fine_elem_num;
fine_elem_h = (fine_elem_avg_vol*6*sqrt(2))^(1/3);

% Simplify data
fine_node_locs = fine_ux_1(:,2:4);
fine_ux_1 = fine_ux_1(:,5);
fine_uy_1 = fine_uy_1(:,5);
fine_uz_1 = fine_uz_1(:,5);
fine_ux_2 = fine_ux_2(:,5);
fine_uy_2 = fine_uy_2(:,5);
fine_uz_2 = fine_uz_2(:,5);

% Normalize such that max displacement is equal to 1
fine_u_1 = [fine_ux_1, fine_uy_1, fine_uz_1];
fine_u_1_norm = vecnorm(fine_u_1, 2, 2);
fine_ux_1 = fine_ux_1 ./ max(fine_u_1_norm);
fine_uy_1 = fine_uy_1 ./ max(fine_u_1_norm);
fine_uz_1 = fine_uz_1 ./ max(fine_u_1_norm);
fine_psi_1 = [fine_ux_1, fine_uy_1, fine_uz_1];

fine_u_2 = [fine_ux_2, fine_uy_2, fine_uz_2];
fine_u_2_norm = vecnorm(fine_u_2, 2, 2);
fine_ux_2 = fine_ux_2 ./ max(fine_u_2_norm);
fine_uy_2 = fine_uy_2 ./ max(fine_u_2_norm);
fine_uz_2 = fine_uz_2 ./ max(fine_u_2_norm);
fine_psi_2 = [fine_ux_2, fine_uy_2, fine_uz_2];


% Compute coarse kinetic and potential energy
disp('Calculating coarse diagonal kinetic energy.');
ckId = kinetic_energy_diag(coarse_node_elem, coarse_node_locs, coarse_psi_1, mat_props);
disp('Calculating coarse cross kinetic energy.');
ckIc = kinetic_energy_cross(coarse_node_elem, coarse_node_locs, coarse_psi_1, coarse_psi_2, mat_props);

disp('Computing coarse diagonal potential energy.');
cpId = potential_energy_diag(coarse_node_elem, coarse_node_locs, coarse_psi_1, mat_props);
disp('Computing coarse cross potential energy.');
cpIc = potential_energy_cross(coarse_node_elem, coarse_node_locs, coarse_psi_1, coarse_psi_2, mat_props);


% Compute mid kinetic and potential energy
disp('Computing mid diagonal kinetic energy.');
mkId = kinetic_energy_diag(mid_node_elem, mid_node_locs, mid_psi_1, mat_props);
disp('Computing mid cross kinetic energy.');
mkIc = kinetic_energy_cross(mid_node_elem, mid_node_locs, mid_psi_1, mid_psi_2, mat_props);

disp('Computing mid diagonal potential energy.');
mpId = potential_energy_diag(mid_node_elem, mid_node_locs, mid_psi_1, mat_props);
disp('Computing mid cross potential energy.');
mpIc = potential_energy_cross(mid_node_elem, mid_node_locs, mid_psi_1, mid_psi_2, mat_props);


% Compute fine kinetic and potential energy
disp('Computing fine diagonal kinetic energy.');
fkId = kinetic_energy_diag(fine_node_elem, fine_node_locs, fine_psi_1, mat_props);
disp('Computing fine cross kinetic energy.');
fkIc = kinetic_energy_cross(fine_node_elem, fine_node_locs, fine_psi_1, fine_psi_2, mat_props);

disp('Computing fine diagonal potential energy.');
fpId = potential_energy_diag(fine_node_elem, fine_node_locs, fine_psi_1, mat_props);
disp('Computing fine cross potential energy.');
fpIc = potential_energy_cross(fine_node_elem, fine_node_locs, fine_psi_1, fine_psi_2, mat_props);


%% Kinetic energy analytical data
R = 0.1;
H = 0.0001;
h = H/2;

% Solve for frequencies given plate radius
fun = @(x) besselj(0,x)*besseli(1,x) + besseli(0,x)*besselj(1,x);
ka1 = fzero(fun,3.1);
ka2 = fzero(fun,6.3);
k1 = ka1 / R;
k2 = ka2 / R;

% Construct analytical expression for plate deformation
syms r

% Compute normalization factor
abar1 = 1 / (besselj(0,0) - (besselj(0,ka1)/besseli(0,ka1))*besseli(0,0));
abar2 = 1 / (besselj(0,0) - (besselj(0,ka2)/besseli(0,ka2))*besseli(0,0));

% Construct integrand
psi1 = abar1*(besselj(0,k1*r) - (besselj(0,ka1)/besseli(0,ka1))*besseli(0,k1*r));
psi2 = abar2*(besselj(0,k2*r) - (besselj(0,ka2)/besseli(0,ka2))*besseli(0,k2*r));
psi1_sq = psi1*psi1;

% Perform integration
akId = 2*pi*H*int(rho*psi1_sq*r, [0,R]);
akIc = 2*pi*H*int(rho*psi1*psi2*r, [0,R]);

% Numerical integration -- to check
%fun = @(r) (abar^2)*r.*(besselj(0,k*r) - (besselj(0,ka1)/besseli(0,ka1))*besseli(0,k*r)).^2;
%nint = integral(fun,0,R);

% Convert to numerical form
akIdvpa = vpa(akId);
akIcvpa = vpa(akIc);

% Construct error
ckIderr = abs(akIdvpa - ckId) / abs(akIdvpa);
mkIderr = abs(akIdvpa - mkId) / abs(akIdvpa);
fkIderr = abs(akIdvpa - fkId) / abs(akIdvpa);

ckIcerr = abs(akIcvpa - ckIc) / abs(akIcvpa);
mkIcerr = abs(akIcvpa - mkIc) / abs(akIcvpa);
fkIcerr = abs(akIcvpa - fkIc) / abs(akIcvpa);

c = exp(log(fkIderr) - 3*log(fine_elem_h));
yc = c.*[coarse_elem_h, mid_elem_h, fine_elem_h].^3;

% Plot diagonal kinetic energy results
figure;
lg1 = loglog([coarse_elem_h, mid_elem_h, fine_elem_h], [ckIderr, mkIderr, fkIderr], '-ob'); hold on;
lg1.LineWidth = 2;
lg2 = loglog([coarse_elem_h, mid_elem_h, fine_elem_h], yc, '--r'); hold on;
lg2.LineWidth = 2;
legend('Numerical', '$\mathcal{O} (h^3)$', 'interpreter', 'latex', 'fontsize', 24, 'location', 'northwest');
xlim([5.5e-4, 1e-3]); ylim([4e-3, 7e-2]);
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'fontsize',18);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',18);
grid on;
xlabel('Avg. Tetrahedron Side Length [m]', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Error: $|KE_{a} - KE_{n}| / |KE_{a}|$', 'interpreter', 'latex', 'fontsize', fs);
title('Diagonal Kinetic Energy Verification', 'interpreter', 'latex', 'fontsize', fs);

% Plot cross kinetic energy results

c = exp(log(fkIcerr) - 3*log(fine_elem_h));
yc = c.*[coarse_elem_h, mid_elem_h, fine_elem_h].^3;

figure;
lg1 = loglog([coarse_elem_h, mid_elem_h, fine_elem_h], [ckIcerr, mkIcerr, fkIcerr], '-ob'); hold on;
lg1.LineWidth = 2;
lg2 = loglog([coarse_elem_h, mid_elem_h, fine_elem_h], yc, '--r'); hold on;
lg2.LineWidth = 2;
legend('Numerical', '$\mathcal{O} (h^3)$', 'interpreter', 'latex', 'fontsize', 24, 'location', 'northwest');
%xlim([5.5e-4, 1e-3]); ylim([4e-3, 7e-2]);
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'fontsize',18);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',18);
grid on;
xlabel('Avg. Tetrahedron Side Length [m]', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Error: $|KE_{a} - KE_{n}| / |KE_{a}|$', 'interpreter', 'latex', 'fontsize', fs);
title('Cross Kinetic Energy Verification', 'interpreter', 'latex', 'fontsize', fs);



%% Potential energy analytical data
% Following P.A. Kelly, Ch. 6.9.3
syms w r

D = (E*H^3) / (12*(1-nu^2));

abar1 = 1 / (besselj(0,0) - (besselj(0,ka1)/besseli(0,ka1))*besseli(0,0));
abar2 = 1 / (besselj(0,0) - (besselj(0,ka2)/besseli(0,ka2))*besseli(0,0));
psi1 = abar1*(besselj(0,k1*r) - (besselj(0,ka1)/besseli(0,ka1))*besseli(0,k1*r));
psi2 = abar2*(besselj(0,k2*r) - (besselj(0,ka2)/besseli(0,ka2))*besseli(0,k2*r));
w1 = psi1; w2 = psi1 + psi2; w3 = psi2;

dw1dr = diff(w1,r);
d2w1dr2 = diff(dw1dr,r);
dw2dr = diff(w2,r);
d2w2dr2 = diff(dw2dr,r);
dw3dr = diff(w3,r);
d2w3dr2 = diff(dw3dr,r);

du1 = (d2w1dr2 + (1/r)*dw1dr)^2 - 2*(1-nu)*d2w1dr2*((1/r)*dw1dr);
du2 = (d2w2dr2 + (1/r)*dw2dr)^2 - 2*(1-nu)*d2w2dr2*((1/r)*dw2dr);
du3 = (d2w3dr2 + (1/r)*dw3dr)^2 - 2*(1-nu)*d2w3dr2*((1/r)*dw3dr);

U1 = pi*D*int(du1*r,[0,R]);
U2 = pi*D*int(du2*r,[0,R]);
U3 = pi*D*int(du3*r,[0,R]);
Udvpa = vpa(U1);
Ucvpa = vpa(U2) - vpa(U1) - vpa(U3);

% Construct error
cpIderr = abs(Udvpa - cpId) / abs(Udvpa);
mpIderr = abs(Udvpa - mpId) / abs(Udvpa);
fpIderr = abs(Udvpa - fpId) / abs(Udvpa);

cpIcerr = abs(Ucvpa - cpIc) / abs(Ucvpa);
mpIcerr = abs(Ucvpa - mpIc) / abs(Ucvpa);
fpIcerr = abs(Ucvpa - fpIc) / abs(Ucvpa);

c = exp(log(fpIderr) - 3*log(fine_elem_h));
yc = c.*[coarse_elem_h, mid_elem_h, fine_elem_h].^3;

% Plot diagonal potential energy results
figure;
lg1 = loglog([coarse_elem_h, mid_elem_h, fine_elem_h], [cpIderr, mpIderr, fpIderr], '-ob'); hold on;
lg1.LineWidth = 2;
lg2 = loglog([coarse_elem_h, mid_elem_h, fine_elem_h], yc, '--r'); hold on;
lg2.LineWidth = 2;
legend('Numerical', '$\mathcal{O} (h^3)$', 'interpreter', 'latex', 'fontsize', 24, 'location', 'northwest');
xlim([5.5e-4, 1e-3]); ylim([7e-2, 5]);
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'fontsize',18);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',18);
grid on;
xlabel('Avg. Tetrahedron Side Length [m]', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Error: $|PE_{a} - PE_{n}| / |PE_{a}|$', 'interpreter', 'latex', 'fontsize', fs);
title('Diagonal Strain Energy Verification', 'interpreter', 'latex', 'fontsize', fs);


c = exp(log(fpIcerr) - 3*log(fine_elem_h));
yc = c.*[coarse_elem_h, mid_elem_h, fine_elem_h].^3;

% Plot diagonal potential energy results
figure;
lg1 = loglog([coarse_elem_h, mid_elem_h, fine_elem_h], [cpIcerr, mpIcerr, fpIcerr], '-ob'); hold on;
lg1.LineWidth = 2;
lg2 = loglog([coarse_elem_h, mid_elem_h, fine_elem_h], yc, '--r'); hold on;
lg2.LineWidth = 2;
legend('Numerical', '$\mathcal{O} (h^3)$', 'interpreter', 'latex', 'fontsize', 24, 'location', 'northwest');
xlim([5.5e-4, 1e-3]); ylim([7e-2, 5]);
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'fontsize',18);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',18);
grid on;
xlabel('Avg. Tetrahedron Side Length [m]', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Error: $|PE_{a} - PE_{n}| / |PE_{a}|$', 'interpreter', 'latex', 'fontsize', fs);
title('Cross Strain Energy Verification', 'interpreter', 'latex', 'fontsize', fs);



