
close all;
clear all;
clc;

format long

addpath('../../fluid_function_utils');
addpath('../../struct_energies_utils');

fs = 26;

disp('Defining parameters');
% Define 10-20-12 parameters
a = 8.5 * 0.0254; % m
b = 18.5 * 0.0254; % m
t = 0.012 * 0.0254; % m
area = a*b; % m^2
d = 2.5 * 0.0254; % m

% Construct deformation function
l1 = t; l2 = 1/l1;

% Define material properties
E = 7.31e10; % Pa
nu = 0.33; % [-]
rho = 2780; %[kg m^-3]
mat_props = [E, nu, rho];
D = (E*t^3)/(12*(1-nu^2));

% Define imposed fluid properties
M_inf = 5; % [-]
T_inf = 300; % K
q_inf = 370000 * 6894.76; % Pa
gam = 1.4;
R = 287;
lambda = (2*q_inf*a^3)/(sqrt(M_inf^2 - 1)*D);

% Compute dependent fluid properties
rho_inf = 2*q_inf/((M_inf^2)*gam*R*T_inf);
U_inf = M_inf*sqrt(gam*R*T_inf);
a_inf = sqrt(gam*R*T_inf);

% Define mean fluid values
rhol = @(x,y) rho_inf .* (x.^0).*(y.^0);
al = @(x,y) a_inf .* (x.^0).*(y.^0);
vl = @(x,y) [U_inf; 0; 0] .* (x.^0).*(y.^0);

% Define quadrature locations (integration points) -- utilize 3-pt rule as FEM data is second-order accurate
% Rule from Alex and is contained in FEM textbooks 
int_loc1 = [2/3, 1/6, 1/6];
int_loc2 = [1/6, 2/3, 1/6];
int_loc3 = [1/6, 1/6, 2/3];

disp('Importing structural data');
% Import first mode structural data
elem_node = importdata('./10_20_12_plate/mesh/elements.txt'); elem_node = elem_node.data;
node_elem = dlmread('./10_20_12_plate/mesh/node_elem.txt');
nodes = dlmread('./10_20_12_plate/mesh/nodes.txt');
structural_node_locs = nodes(:,2:4);
surface_node_data = dlmread('./10_20_12_plate/mesh/face_nodes.txt');
struct_faces = dlmread('./10_20_12_plate/mesh/surf_faces_trim.txt');

nmodes = 1;
disp('Importing second mode structural data and shape function');
mode_xdir_deform = dlmread('./10_20_12_plate/mode5/mode5_xdir_deform.txt');
mode_ydir_deform = dlmread('./10_20_12_plate/mode5/mode5_ydir_deform.txt');
mode_zdir_deform = dlmread('./10_20_12_plate/mode5/mode5_zdir_deform.txt');
def_shape_function = [mode_xdir_deform(:,5), mode_ydir_deform(:,5), mode_zdir_deform(:,5)];
def_shape_function_norms = vecnorm(def_shape_function'); max_def = l2*max(def_shape_function_norms);
def_shape_function = def_shape_function ./ max_def;

disp('Computing kinetic and strain energy of mode 1');
% Compute kinetic and strain energies
%T1 = kinetic_energy_diag(elem_node, structural_node_locs, def_shape_function_mode1, mat_props);
%U1 = potential_energy_diag(elem_node, structural_node_locs, def_shape_function_mode1, mat_props);
%T1 = 1.362405569524759e-09;
%U1 = 8.368843386290405e-05;

%disp('Computing kinetic and strain energy of mode 2');
%%T2 = kinetic_energy_diag(elem_node, structural_node_locs, def_shape_function_mode2, mat_props);
%%U2 = potential_energy_diag(elem_node, structural_node_locs, def_shape_function_mode2, mat_props);
T1 = 1.597669525012157e-09;
U1 = 6.799211788509672e-04;


% Compute fluid response coefficients
disp('Performing compliant surface integrations');
N = max(size(struct_faces));
disp(strcat(['Number of faces: ', num2str(N)]));
dF1dq1_arr = zeros(1,N);
tic;
%parfor (i = 1:N, 12)
for i = 1:N
  disp(strcat(['Face ', num2str(i)]));

  % Grab face element number and nodes on surface
  loc_elem = struct_faces(i,1);
  loc_elem_surf_nodes = struct_faces(i,2:11);

  % Get all element nodes and locations
  loc_elem_nodes = elem_node(loc_elem, :);
  loc_elem_nodelocs = structural_node_locs(loc_elem_nodes,:);
  loc_elem_surf_num = find(loc_elem_surf_nodes(1:4) == 0);
  loc_elem_surf = [1 2 3 4]; loc_elem_surf(loc_elem_surf_num) = [];

  psi_nodes = def_shape_function(loc_elem_nodes,:);

  %%% Go over each int loc and form integrand
  %%% Int loc 1
  % Get location of integration point in Cartesian (global) and barycentric (local) coords
  location1_bary = zeros(1,4); 
  location1_bary(loc_elem_surf) = int_loc1;
  location1 = tet_interp(loc_elem_nodelocs, location1_bary);

  % Interpolate deformation, normal and normal derivative at integration point
  psi1 = tet_interp(psi_nodes, location1_bary);
  n1 = determine_normal(loc_elem_nodelocs, location1, location1_bary, loc_elem_surf_num);
  dn1dq = determine_normal_deriv(loc_elem_nodelocs, psi_nodes, location1, location1_bary, loc_elem_surf_num);
  [dxdLi1, dxdLj1, ~] = def_derivs(loc_elem_nodelocs, psi_nodes, location1, location1_bary, loc_elem_surf_num);
  A1 = cross(dxdLi1, dxdLj1); A1 = norm(A1);

  xl = location1(1); yl = location1(2);
  rho1 = rhol(xl,yl);
  a1 = al(xl,yl);
  v1 = vl(xl,yl);


  %%% Int loc 2
  % Get location of integration point in Cartesian (global) and barycentric (local) coords
  location2_bary = zeros(1,4);
  location2_bary(loc_elem_surf) = int_loc2;
  location2 = tet_interp(loc_elem_nodelocs, location2_bary);

  % Interpolate deformation, normal, and normal derivative at integration point
  psi2 = tet_interp(psi_nodes, location2_bary);
  n2 = determine_normal(loc_elem_nodelocs, location2, location2_bary, loc_elem_surf_num);
  dn2dq = determine_normal_deriv(loc_elem_nodelocs, psi_nodes, location2, location2_bary, loc_elem_surf_num);
  [dxdLi2, dxdLj2, ~] = def_derivs(loc_elem_nodelocs, psi_nodes, location2, location2_bary, loc_elem_surf_num);
  A2 = cross(dxdLi2, dxdLj2); A2 = norm(A2);

  xl = location2(1); yl = location2(2);
  rho2 = rhol(xl,yl);
  a2 = al(xl,yl);
  v2 = vl(xl,yl);


  %%% Int loc 3
  % Get location of integration point in Cartesian (global) and barycentric (local) coords
  location3_bary = zeros(1,4);
  location3_bary(loc_elem_surf) = int_loc3;
  location3 = tet_interp(loc_elem_nodelocs, location3_bary);

  % Interpolate deformation, normal, and normal derivative at integration point
  psi3 = tet_interp(psi_nodes, location3_bary);
  n3 = determine_normal(loc_elem_nodelocs, location3, location3_bary, loc_elem_surf_num);
  dn3dq = determine_normal_deriv(loc_elem_nodelocs, psi_nodes, location3, location3_bary, loc_elem_surf_num);
  [dxdLi3, dxdLj3, ~] = def_derivs(loc_elem_nodelocs, psi_nodes, location3, location3_bary, loc_elem_surf_num);
  A3 = cross(dxdLi3, dxdLj3); A3 = norm(A3);

  xl = location3(1); yl = location3(2);
  rho3 = rhol(xl,yl);
  a3 = al(xl,yl);
  v3 = vl(xl,yl);

  if ((norm(n2 - [0 0 1])) > 1e-8)
    disp('hey'); pause;
  end

  % Compute integrands at integration points
  dF1dq1_int1 = -rho1*a1*(dot(v1,-1.*dn1dq')*dot(n1',psi1'))*A1;
  dF1dq1_int2 = -rho2*a2*(dot(v2,-1.*dn2dq')*dot(n2',psi2'))*A2;
  dF1dq1_int3 = -rho3*a3*(dot(v3,-1.*dn3dq')*dot(n3',psi3'))*A3;

  % Perform quadrature
  dF1dq1_loc = quadrature(dF1dq1_int1, dF1dq1_int2, dF1dq1_int3);

  dF1dq1_arr(i) = dF1dq1_loc;

end
toc

dFidqj = zeros(2,2);
dF1dq1 = sum(dF1dq1_arr);

disp('Results:');
disp(strcat(['T1: ', num2str(T1)]));
disp(strcat(['U1: ', num2str(U1)]));
disp(strcat(['dF1dq1: ', num2str(dF1dq1)]));

syms A
syms p
A(1) = T1*p^2 + (U1 - dF1dq1);
charpoly = det(A) == 0;
B = double(solve(charpoly));
disp(B);

%disp('Coefficients:')
%disp(strcat(['leading p^4:', num2str(T1*T2)]));
%disp(strcat(['leading p^2:', num2str(T1*(U2 - dFidqj(2,2)) + T2*(U1 - dFidqj(1,1)))]));
%disp(strcat(['leading p^0:', num2str((U1 - dFidqj(1,1))*(U2 - dFidqj(2,2)) - dFidqj(1,2)*dFidqj(2,1))]));
%
%disp('Coefficients w/ temp:')
%disp(strcat(['leading p^4:', num2str(T1*T2)]));
%disp(strcat(['leading p^2:', num2str(T1*(U2_temp - dFidqj(2,2)) + T2*(U1_temp - dFidqj(1,1)))]));
%disp(strcat(['leading p^0:', num2str((U1_temp - dFidqj(1,1))*(U2_temp - dFidqj(2,2)) - dFidqj(1,2)*dFidqj(2,1))]));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Triangular interpolation function for location
function [loc] = tri_loc(face_node_locs, interp_loc)
  % Linear triangle location from area coordinates
  % face_node_locs: Cartesian locations of each face node
  % interp_loc: location of point in triangle, area coordinates
  % loc: location of point in triangle, Cartesian coordinates

  loc = face_node_locs(1,:).*interp_loc(1) + face_node_locs(2,:).*interp_loc(2) + face_node_locs(3,:).*interp_loc(3);

end

% Triangular interpolation function for data
function [quant] = tri_interp(interp_loc, data);

  quant = data(1)*interp_loc(1) + data(2)*interp_loc(2) + data(3)*interp_loc(3);

end

% Third-order quadrature function
function [int] = quadrature(val1, val2, val3)

  % Weights obtained from Zienkiewics and modified to ensure that transformed surface integral provides A = 1/2.
  w1 = 1/6;
  w2 = 1/6;
  w3 = 1/6;

  int = w1*val1 + w2*val2 + w3*val3;

end


