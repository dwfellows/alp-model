%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to attempt to verify energies_func
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all;
clear all;
clc;

addpath('../../struct_energies_utils');
fs = 28;

%% Numerical data

% Read in data for circular verification cases
% Structural steel material properties
E = 2e11; % Pa
nu = 0.3; % dimensionless
rho = 7850; %kg/m^3
mat_props = [E, nu, rho];

% Circular plate properties
R = 1;
H = 0.005;
h = H/2;

vol = pi*(R^2)*H;

% Calculate numerical kinetic and strain energies
elem_h_list = [];
k_list = [];
p_list = [];
for i = 1:6
  disp(strcat(['Case ', num2str(i)]));
  fileloc = strcat(['./hires_circular_data_2/c', num2str(i), '/']);
  node_elem = importdata(strcat([fileloc, 'elem_node.txt'])); node_elem = node_elem.data;
  ux = dlmread(strcat([fileloc, 'xdir.txt']));
  uy = dlmread(strcat([fileloc, 'ydir.txt']));
  uz = dlmread(strcat([fileloc, 'zdir.txt']));

  elem_num = max(size(node_elem));
  elem_avg_vol = vol / elem_num;
  elem_h = (elem_avg_vol*6*sqrt(2))^(1/3);
  elem_h_list = [elem_h_list, elem_h];

  node_locs = ux(:,2:4);
  u = [ux(:,5), uy(:,5), uz(:,5)];
  u_norm = vecnorm(u, 2,2);
  ux = ux(:,5) ./ max(u_norm);
  uy = uy(:,5) ./ max(u_norm);
  uz = uz(:,5) ./ max(u_norm);
  psi = [ux, uy, uz];
  disp('Computing kinetic energy');
  k_loc = kinetic_energy_diag(node_elem, node_locs, psi, mat_props);
  k_list = [k_list, k_loc];
  disp('Computing strain energy');
  p_loc = potential_energy_diag(node_elem, node_locs, psi, mat_props);
  p_list = [p_list, p_loc];

end 
elem_h_list = elem_h_list / R;

%% Kinetic energy analytical data
R = 1;
H = 0.005;
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
k_err = abs(akIdvpa - k_list) ./ abs(akIdvpa);

c = exp(log(k_err(end)) - 3*log(elem_h_list(end)));
yc = c.*elem_h_list.^3;

% Plot diagonal kinetic energy results
f1 = figure;
lg1 = loglog(elem_h_list, k_err, '-ob'); hold on;
lg1.LineWidth = 2;
lg2 = loglog(elem_h_list, yc, '--r'); hold on;
lg2.LineWidth = 2;
legend('Numerical', '$\mathcal{O} (h^3)$', 'interpreter', 'latex', 'fontsize', 30, 'location', 'northwest');
xlim([.009, .022]); ylim([8e-4, 1e-2]);
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'fontsize',18);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',18);
grid on;
xlabel('$h / R$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Error: $| T_{11} |_{a} - T_{11} |_{n} | / |T_{11} |_{a} |$', 'interpreter', 'latex', 'fontsize', fs);
title('$T_{11}$ Verification', 'interpreter', 'latex', 'fontsize', fs);
f1.Position = [100 100 900 700];

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
p_err = abs(Udvpa - (0.5).*p_list) ./ abs(Udvpa);
c = exp(log(p_err(end)) - 3*log(elem_h_list(end)));
yc = c.*elem_h_list.^3;

% Plot diagonal potential energy results
f3 = figure;
lg1 = loglog(elem_h_list, p_err, '-ob'); hold on;
lg1.LineWidth = 2;
lg2 = loglog(elem_h_list, yc, '--r'); hold on;
lg2.LineWidth = 2;
legend('Numerical', '$\mathcal{O} (h^3)$', 'interpreter', 'latex', 'fontsize', 30, 'location', 'northwest');
xlim([.009, .022]); ylim([8e-4, 3e-2]);
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'fontsize',18);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',18);
grid on;
xlabel('$h/R$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Error: $|U_{11} |_{a} - U_{11} |_{n}| / |U_{11} |_{a}|$', 'interpreter', 'latex', 'fontsize', fs);
title('$U_{11}$ Verification', 'interpreter', 'latex', 'fontsize', fs);
f3.Position = [100 100 900 700];


