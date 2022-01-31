
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
E = 68.9e9; % Pa
nu = 0.33; % [-]
rho = 2700; %[kg m^-3]
mat_props = [E, nu, rho];
D = (E*t^3)/(12*(1-nu^2));

% Define imposed fluid properties
M_inf = 5; % [-]
T_inf = 300; % K
q_inf = 50000 * 6894.76; % Pa
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

% Construct mode data using C-C beam modes
nmodes = 4;

g1 = 4.73; gp1 = 7.853; g2 = 10.996; gp2 = 14.137;
k1 = sin(g1/2)/sinh(g1/2); kp1 = -sin(gp1/2)/sinh(gp1/2); k2 = sin(g2/2)/sinh(g2/2); kp2 = -sin(gp2/2)/sinh(gp2/2);
exp1 = @(x,y) ( cos(g1.*(x./a - 1/2)) + k1.*cosh(g1.*(x./a - 1/2)) ).*( cos(g1.*(y./b - 1/2)) + k1.*cosh(g1.*(y./b - 1/2)) );
exp2 = @(x,y) ( sin(gp1.*(x./a - 1/2)) + kp1.*sinh(gp1.*(x./a - 1/2)) ).*( cos(g1.*(y./b - 1/2)) + k1.*cosh(g1.*(y./b - 1/2)) );
exp3 = @(x,y) ( cos(g2.*(x./a - 1/2)) + k2.*cosh(g2.*(x./a - 1/2)) ).*( cos(g1.*(y./b - 1/2)) + k1.*cosh(g1.*(y./b - 1/2)) );
exp4 = @(x,y) ( sin(gp2.*(x./a - 1/2)) + kp2.*sinh(gp2.*(x./a - 1/2)) ).*( cos(g1.*(y./b - 1/2)) + k1.*cosh(g1.*(y./b - 1/2)) );
N = max(size(structural_node_locs));

disp('Constructing first mode structural data and shape function');
mode1_zdir_deform = exp1(structural_node_locs(:,1), structural_node_locs(:,2));
def_shape_function = [zeros(N,2), mode1_zdir_deform];
def_shape_function_norms = vecnorm(def_shape_function'); max_def = l2*max(def_shape_function_norms);
def_shape_functions(:,:,1) = def_shape_function ./ max_def;

disp('Constructing second mode structural data and shape function');
mode2_zdir_deform = exp2(structural_node_locs(:,1), structural_node_locs(:,2));
def_shape_function = [zeros(N,2), mode2_zdir_deform];
def_shape_function_norms = vecnorm(def_shape_function'); max_def = l2*max(def_shape_function_norms);
def_shape_functions(:,:,2) = def_shape_function ./ max_def;

disp('Constructing third mode structural data and shape function');
mode3_zdir_deform = exp3(structural_node_locs(:,1), structural_node_locs(:,2));
def_shape_function = [zeros(N,2), mode3_zdir_deform];
def_shape_function_norms = vecnorm(def_shape_function'); max_def = l2*max(def_shape_function_norms);
def_shape_functions(:,:,3) = def_shape_function ./ max_def;

disp('Constructing fourth mode structural data and shape function');
mode4_zdir_deform = exp4(structural_node_locs(:,1), structural_node_locs(:,2));
def_shape_function = [zeros(N,2), mode4_zdir_deform];
def_shape_function_norms = vecnorm(def_shape_function'); max_def = l2*max(def_shape_function_norms);
def_shape_functions(:,:,4) = def_shape_function ./ max_def;

disp('Computing structural energies using C-C beam functions');
T(1:4) = (1e-8) .* [0.121932690447277   0.135019593908423   0.134417091026494   0.134451689894264];
U(1:4) = [3.738841031596209  13.517628191134563  28.086253120701834  48.198058132031640];
%omega = sqrt(U./T);
if ~exist('T')
  T = zeros(1, nmodes);
  U = zeros(1, nmodes);
  for i = 1:nmodes
    disp(strcat(['Computing structural coefficients of energy for mode ', num2str(i)]));
    T(i) = kinetic_energy_diag(elem_node, structural_node_locs, def_shape_functions(:,:,i), mat_props);
    U(i) = potential_energy_diag(elem_node, structural_node_locs, def_shape_functions(:,:,i), mat_props);
    omega(i) = sqrt(U(i)/T(i));
    disp(strcat(['T', num2str(i), ': ', num2str(T(i))]));
    disp(strcat(['U', num2str(i), ': ', num2str(U(i))]));
    disp(strcat(['omega', num2str(i), ': ', num2str(omega(i))]));
  end
end


% Compute fluid response coefficients
disp('Performing compliant surface integrations');
N = max(size(struct_faces));
disp(strcat(['Number of faces: ', num2str(N)]));
dFidqj = zeros(nmodes);
dFidqj_int1 = zeros(nmodes); dFidqj_int2 = zeros(nmodes); dFidqj_int3 = zeros(nmodes);
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

  psi_nodes = def_shape_functions(loc_elem_nodes,:,:);  % [10x3xnmodes]

  %%% Go over each int loc and form integrand
  %%% Int loc 1
  % Get location of integration point in Cartesian (global) and barycentric (local) coords
  location1_bary = zeros(1,4); 
  location1_bary(loc_elem_surf) = int_loc1;
  location1 = tet_interp(loc_elem_nodelocs, location1_bary);

  % Interpolate deformation, normal and normal derivative at integration point
  for j = 1:nmodes
    psi1_modes(j,:) = tet_interp(psi_nodes(:,:,j), location1_bary);
    dn1dq_modes(j,:) = determine_normal_deriv(loc_elem_nodelocs, psi_nodes(:,:,j), location1, location1_bary, loc_elem_surf_num);
  end
  n1 = determine_normal(loc_elem_nodelocs, location1, location1_bary, loc_elem_surf_num);
  [dxdLi1, dxdLj1, ~] = def_derivs(loc_elem_nodelocs, psi_nodes(:,:,1), location1, location1_bary, loc_elem_surf_num);
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
  for j = 1:nmodes
    psi2_modes(j,:) = tet_interp(psi_nodes(:,:,j), location2_bary);
    dn2dq_modes(j,:) = determine_normal_deriv(loc_elem_nodelocs, psi_nodes(:,:,j), location2, location2_bary, loc_elem_surf_num);
  end
  n2 = determine_normal(loc_elem_nodelocs, location2, location2_bary, loc_elem_surf_num);
  [dxdLi2, dxdLj2, ~] = def_derivs(loc_elem_nodelocs, psi_nodes(:,:,1), location2, location2_bary, loc_elem_surf_num);
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
  for j = 1:nmodes
    psi3_modes(j,:) = tet_interp(psi_nodes(:,:,j), location3_bary);
    dn3dq_modes(j,:) = determine_normal_deriv(loc_elem_nodelocs, psi_nodes(:,:,j), location3, location3_bary, loc_elem_surf_num);
  end
  n3 = determine_normal(loc_elem_nodelocs, location3, location3_bary, loc_elem_surf_num);
  [dxdLi3, dxdLj3, ~] = def_derivs(loc_elem_nodelocs, psi_nodes(:,:,1), location3, location3_bary, loc_elem_surf_num);
  A3 = cross(dxdLi3, dxdLj3); A3 = norm(A3);

  xl = location3(1); yl = location3(2);
  rho3 = rhol(xl,yl);
  a3 = al(xl,yl);
  v3 = vl(xl,yl);

  if ((norm(n2 - [0 0 1])) > 1e-8)
    disp('hey'); pause;
  end

  % Compute integrands at integration points
%  dF1dq1_int1 = -rho1*a1*(dot(v1,-1.*dn1dq_mode1')*dot(n1',psi1_mode1'))*A1;
%  dF1dq2_int1 = -rho1*a1*(dot(v1,-1.*dn1dq_mode2')*dot(n1',psi1_mode1'))*A1;
%  dF2dq1_int1 = -rho1*a1*(dot(v1,-1.*dn1dq_mode1')*dot(n1',psi1_mode2'))*A1;
%  dF2dq2_int1 = -rho1*a1*(dot(v1,-1.*dn1dq_mode2')*dot(n1',psi1_mode2'))*A1;
%
%  dF1dq1_int2 = -rho2*a2*(dot(v2,-1.*dn2dq_mode1')*dot(n2',psi2_mode1'))*A2;
%  dF1dq2_int2 = -rho2*a2*(dot(v2,-1.*dn2dq_mode2')*dot(n2',psi2_mode1'))*A2;
%  dF2dq1_int2 = -rho2*a2*(dot(v2,-1.*dn2dq_mode1')*dot(n2',psi2_mode2'))*A2;
%  dF2dq2_int2 = -rho2*a2*(dot(v2,-1.*dn2dq_mode2')*dot(n2',psi2_mode2'))*A2;
%
%  dF1dq1_int3 = -rho3*a3*(dot(v3,-1.*dn3dq_mode1')*dot(n3',psi3_mode1'))*A3;
%  dF1dq2_int3 = -rho3*a3*(dot(v3,-1.*dn3dq_mode2')*dot(n3',psi3_mode1'))*A3;
%  dF2dq1_int3 = -rho3*a3*(dot(v3,-1.*dn3dq_mode1')*dot(n3',psi3_mode2'))*A3;
%  dF2dq2_int3 = -rho3*a3*(dot(v3,-1.*dn3dq_mode2')*dot(n3',psi3_mode2'))*A3;

  % Compute value at interpolation points
  for j = 1:nmodes
    for k = 1:nmodes
      dFidqj_int1(j,k) = -(2*q_inf/(U_inf*sqrt(M_inf^2 - 1)))*(dot(v1,-1.*dn1dq_modes(k,:)')*dot(n1',psi1_modes(j,:)))*A1;
      dFidqj_int2(j,k) = -(2*q_inf/(U_inf*sqrt(M_inf^2 - 1)))*(dot(v2,-1.*dn2dq_modes(k,:)')*dot(n2',psi2_modes(j,:)))*A2;
      dFidqj_int3(j,k) = -(2*q_inf/(U_inf*sqrt(M_inf^2 - 1)))*(dot(v3,-1.*dn3dq_modes(k,:)')*dot(n3',psi3_modes(j,:)))*A3;
    end
  end
  
  % Perform quadrature
  for j = 1:nmodes
    for k = 1:nmodes
      dFidqj(j,k) = dFidqj(j,k) + quadrature(dFidqj_int1(j,k), dFidqj_int2(j,k), dFidqj_int3(j,k));
    end
  end


end
toc

% Construct matrix
syms A
syms p
for i = 1:nmodes
  for j = 1:nmodes
    if (i==j)
      A(i,i) = T(i)*p^2 + (U(i) - dFidqj(i,j));
    else
      A(i,j) = -dFidqj(i,j);
    end
  end
end
charpoly = (det(A) == 0);
B = double(solve(charpoly));
disp(B);



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


