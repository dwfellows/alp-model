
close all;
clear all;
clc;

addpath('../../../fluid_function_utils');

%% Import fluid data -- all quantities, from surface
% 1 node num | 2 X [m] | 3 Y [m] | 4 Z [m] | 5 p_a [Pa] | 6 A [m^2] | 7 \rho [kg m^-3] | 8 F_x [N] | 9 F_y [N] 
% 10 F_z [N] | 11 a [m s^-1] | 12 Mach # | 13 Mach # Stn | 14 normal_x | 15 normal_y | 16 normal_z | 17 v_u stn [m s^1]
% 18 v_v stn [m s^-1] | 19 v_w stn [m s^-1] | 20 v_u [m s^-1] | 21 v_v [m s^-1] | 22 v_w [m s^-1] | 23 X [m]
% 24 Y [m] | 25 Z [m]
%fluid_data = dlmread('../turbine_blade_data/free_slip/rotating_inner_wall_fluid_data.csv');
%face_conn = dlmread('../turbine_blade_data/free_slip/rotating_inner_wall_face_conn.csv');
%fluid_rho = fluid_data(:,7);
%fluid_p = fluid_data(:,5);
%fluid_a = fluid_data(:,11);
%fluid_vu = fluid_data(:,20);
%fluid_vv = fluid_data(:,21);
%fluid_vw = fluid_data(:,22);
%fluid_node_locs = fluid_data(:,2:4);

% Construct fluid data using circular nodes
fluid_node_surf_circular = dlmread('../fem_data/circular_fluid/fluid_data.csv');
face_conn = dlmread('../fem_data/circular_fluid/face_conn.csv');
fluid_node_locs = fluid_node_surf_circular(:,2:4);

% Fix current fluid node locations and enforce method of manufactured solutions fluid variables
R = 0.1; h = 0.0001;
dev = 0;
N = max(size(fluid_node_locs));

fluid_p = zeros(N,1); fluid_rho = zeros(N,1); fluid_a = zeros(N,1);
fluid_vu = zeros(N,1); fluid_vv = zeros(N,1); fluid_vw = zeros(N,1);

for j = 1:N
  x = fluid_node_locs(j,1);
  y = fluid_node_locs(j,2);
  z = fluid_node_locs(j,3);
  r = sqrt(x^2 + y^2);
  fluid_node_locs(j,3) = dev*(R^2 - r^2) + z;

  fluid_p(j) = (r/R)^2 + 1000;
  fluid_rho(j) = (r/R)^2 + 1;
  fluid_a(j) = sqrt(1.4*287*300);
  fluid_vu(j) = 1; fluid_vv(j) = 1;
  fluid_vw(j) = 0;
end

%% Import structural data -- all quantities
l = 3;
if (l==1)
  elem_node = dlmread('../fem_data/circular_coarse/coarse.csv');
  node_elem = dlmread('../fem_data/circular_coarse/coarse_node_elem.csv');
  xdir_deform = dlmread('../fem_data/circular_coarse/coarse_xdir_deformation.txt');
  ydir_deform = dlmread('../fem_data/circular_coarse/coarse_ydir_deformation.txt');
  zdir_deform = dlmread('../fem_data/circular_coarse/coarse_zdir_deformation.txt');
  def_shape_function = [xdir_deform(:,5), ydir_deform(:,5), zdir_deform(:,5)];
  def_shape_function_norms = vecnorm(def_shape_function'); max_def = max(def_shape_function_norms);
  def_shape_function = def_shape_function ./ max_def;
  structural_node_locs = xdir_deform(:,2:4);
  surface_node_data = dlmread('../fem_data/circular_coarse/coarse_surface.txt');
elseif (l==2)
  elem_node = dlmread('../fem_data/circular_mid/mid.csv');
  node_elem = dlmread('../fem_data/circular_mid/mid_node_elem.csv');
  xdir_deform = dlmread('../fem_data/circular_mid/mid_xdir_deformation.txt');
  ydir_deform = dlmread('../fem_data/circular_mid/mid_ydir_deformation.txt');
  zdir_deform = dlmread('../fem_data/circular_mid/mid_zdir_deformation.txt');
  def_shape_function = [xdir_deform(:,5), ydir_deform(:,5), zdir_deform(:,5)];
  def_shape_function_norms = vecnorm(def_shape_function'); max_def = max(def_shape_function_norms);
  def_shape_function = def_shape_function ./ max_def;
  structural_node_locs = xdir_deform(:,2:4);
  surface_node_data = dlmread('../fem_data/circular_mid/mid_surface.txt');
else
  elem_node = dlmread('../fem_data/circular_fine/fine.csv');
  node_elem = dlmread('../fem_data/circular_fine/fine_node_elem.csv');
  xdir_deform = dlmread('../fem_data/circular_fine/fine_xdir_deformation.txt');
  ydir_deform = dlmread('../fem_data/circular_fine/fine_ydir_deformation.txt');
  zdir_deform = dlmread('../fem_data/circular_fine/fine_zdir_deformation.txt');
  def_shape_function = [xdir_deform(:,5), ydir_deform(:,5), zdir_deform(:,5)];
  def_shape_function_norms = vecnorm(def_shape_function'); max_def = max(def_shape_function_norms);
  def_shape_function = def_shape_function ./ max_def;
  structural_node_locs = xdir_deform(:,2:4);
  surface_node_data = dlmread('../fem_data/circular_fine/fine_surface.txt');
end

k = 3.1962/R;
abar = 1 ./ (besselj(0,0) - (besselj(0,k.*R)/besseli(0,k.*R)).*besseli(0,0));
psi3_anal = @(r) abar.*(besselj(0,k.*r) - (besselj(0,k.*R)/besseli(0,k.*R))*besseli(0,k.*r));
for k = 1:max(size(def_shape_function))
  x = structural_node_locs(k,1);
  y = structural_node_locs(k,2);
  r = sqrt(x^2 + y^2);
  def_shape_function(k,:) = [0, 0, psi3_anal(r)];
end

% Fix current structural node locations
%M = max(size(structural_node_locs));
%for j = 1:M
%  x = structural_node_locs(j,1);
%  y = structural_node_locs(j,2);
%  z = structural_node_locs(j,3);
%  r = sqrt(x^2 + y^2);
%  structural_node_locs(j,3) = dev*(R^2 - r^2) + z;
%end
%
%M = max(size(surface_node_data));
%for j = 1:M
%  x = surface_node_data(j,2);
%  y = surface_node_data(j,3);
%  z = surface_node_data(j,4);
%  r = sqrt(x^2 + y^2);
%  surface_node_data(j,4) = dev*(R^2 - r^2) + z;
%end

% Compute face areas
N = max(size(face_conn));
face_areas = zeros(N,1);

% Define quadrature locations (integration points)
int_loc1 = [1/3, 1/3, 1/3];
int_loc2 = [0.6, 0.2, 0.2];
int_loc3 = [0.2, 0.6, 0.2];
int_loc4 = [0.2, 0.2, 0.6];

%% Compute dF1dq
disp(strcat(['Number of faces: ', num2str(N)]));
dF1dq = 0;
for i = 1:1
  disp(strcat(['Face ', num2str(i)]));

  % Compute face area
  face_nodes = face_conn(i,:);
  face_nodes = face_nodes + 1; % convert from zero-index to one-index
  face_node_locs = fluid_node_locs(face_nodes', :);
  local_area = norm(cross(face_node_locs(2,:) - face_node_locs(1,:), face_node_locs(3,:) - face_node_locs(1,:))) / 2;
  face_areas(i) = local_area;

  % Interpolate fluid quantities at integration points
  rho1 = tri_interp(int_loc1, fluid_rho(face_nodes'));
  rho2 = tri_interp(int_loc2, fluid_rho(face_nodes'));
  rho3 = tri_interp(int_loc3, fluid_rho(face_nodes'));
  rho4 = tri_interp(int_loc4, fluid_rho(face_nodes'));

%  p1 = tri_interp(int_loc1, fluid_p(face_nodes'));
%  p2 = tri_interp(int_loc2, fluid_p(face_nodes'));
%  p3 = tri_interp(int_loc3, fluid_p(face_nodes'));
%  p4 = tri_interp(int_loc4, fluid_p(face_nodes'));

  a1 = tri_interp(int_loc1, fluid_a(face_nodes'));
  a2 = tri_interp(int_loc2, fluid_a(face_nodes'));
  a3 = tri_interp(int_loc3, fluid_a(face_nodes'));
  a4 = tri_interp(int_loc4, fluid_a(face_nodes'));

  vu1 = tri_interp(int_loc1, fluid_vu(face_nodes'));
  vu2 = tri_interp(int_loc2, fluid_vu(face_nodes'));
  vu3 = tri_interp(int_loc3, fluid_vu(face_nodes'));
  vu4 = tri_interp(int_loc4, fluid_vu(face_nodes'));

  vv1 = tri_interp(int_loc1, fluid_vv(face_nodes'));
  vv2 = tri_interp(int_loc2, fluid_vv(face_nodes'));
  vv3 = tri_interp(int_loc3, fluid_vv(face_nodes'));
  vv4 = tri_interp(int_loc4, fluid_vv(face_nodes'));

  vw1 = tri_interp(int_loc1, fluid_vw(face_nodes'));
  vw2 = tri_interp(int_loc2, fluid_vw(face_nodes'));
  vw3 = tri_interp(int_loc3, fluid_vw(face_nodes'));
  vw4 = tri_interp(int_loc4, fluid_vw(face_nodes'));

  v1 = [vu1; vv1; vw1];
  v2 = [vu2; vv2; vw2];
  v3 = [vu3; vv3; vw3];
  v4 = [vu4; vv4; vw4];

  % Interpolate structural quantities at integration points
  loc1 = tri_loc(face_node_locs, int_loc1);
  loc1_elem = determine_elem(elem_node, node_elem, surface_node_data, structural_node_locs, loc1);
  loc1_elem_nodes = elem_node(loc1_elem, :);
  loc1_elem_nodelocs = structural_node_locs(loc1_elem_nodes,:);
  struc_loc1_bary = cart2bary(loc1_elem_nodelocs, loc1);
  psi1_nodes = def_shape_function(loc1_elem_nodes,:);
  psi1 = tet_interp(psi1_nodes, struc_loc1_bary);
  n1 = determine_normal(loc1_elem_nodelocs, loc1, struc_loc1_bary);
  dn1dq = determine_normal_deriv(loc1_elem_nodelocs, psi1_nodes, loc1, struc_loc1_bary);

  loc2 = tri_loc(face_node_locs, int_loc2);
  loc2_elem = determine_elem(elem_node, node_elem, surface_node_data, structural_node_locs, loc2);
  loc2_elem_nodes = elem_node(loc2_elem, :);
  loc2_elem_nodelocs = structural_node_locs(loc2_elem_nodes,:);
  struc_loc2_bary = cart2bary(loc2_elem_nodelocs, loc2);
  psi2_nodes = def_shape_function(loc2_elem_nodes,:);
  psi2 = tet_interp(psi2_nodes, struc_loc2_bary);
  n2 = determine_normal(loc2_elem_nodelocs, loc2, struc_loc2_bary);
  dn2dq = determine_normal_deriv(loc2_elem_nodelocs, psi2_nodes, loc2, struc_loc2_bary);

  loc3 = tri_loc(face_node_locs, int_loc3);
  loc3_elem = determine_elem(elem_node, node_elem, surface_node_data, structural_node_locs, loc3);
  loc3_elem_nodes = elem_node(loc3_elem, :);
  loc3_elem_nodelocs = structural_node_locs(loc3_elem_nodes,:);
  struc_loc3_bary = cart2bary(loc3_elem_nodelocs, loc3);
  psi3_nodes = def_shape_function(loc3_elem_nodes,:);
  psi3 = tet_interp(psi3_nodes, struc_loc3_bary);
  n3 = determine_normal(loc3_elem_nodelocs, loc3, struc_loc3_bary);
  dn3dq = determine_normal_deriv(loc3_elem_nodelocs, psi3_nodes, loc3, struc_loc3_bary);

  loc4 = tri_loc(face_node_locs, int_loc4);
  loc4_elem = determine_elem(elem_node, node_elem, surface_node_data, structural_node_locs, loc4);
  loc4_elem_nodes = elem_node(loc4_elem, :);
  loc4_elem_nodelocs = structural_node_locs(loc4_elem_nodes,:);
  struc_loc4_bary = cart2bary(loc4_elem_nodelocs, loc4);
  psi4_nodes = def_shape_function(loc4_elem_nodes,:);
  psi4 = tet_interp(psi4_nodes, struc_loc4_bary);
  n4 = determine_normal(loc4_elem_nodelocs, loc4, struc_loc4_bary);
  dn4dq = determine_normal_deriv(loc4_elem_nodelocs, psi4_nodes, loc4, struc_loc4_bary);

  % Compute integrands at integration points
  dF1dq_int1 = -rho1*a1*(dot(v1,dn1dq')*dot(n1',psi1') + dot(v1,n1')*dot(dn1dq',psi1'));
  dF1dq_int2 = -rho2*a2*(dot(v2,dn2dq')*dot(n2',psi2') + dot(v2,n2')*dot(dn2dq',psi2'));
  dF1dq_int3 = -rho3*a3*(dot(v3,dn3dq')*dot(n3',psi3') + dot(v3,n3')*dot(dn3dq',psi3'));
  dF1dq_int4 = -rho4*a4*(dot(v4,dn4dq')*dot(n4',psi4') + dot(v4,n4')*dot(dn4dq',psi4'));

  % Perform quadrature
  dF1dq_loc = quadrature(dF1dq_int1, dF1dq_int2, dF1dq_int3, dF1dq_int4, local_area);
  dF1dq = dF1dq + dF1dq_loc;

end

locs = [loc1; loc2; loc3; loc4];

rp = sqrt(locs(:,1).^2 + locs(:,2).^2);
tp = atan2(locs(:,2), locs(:,1));

% Integrate analytically at q = 0
clear r
k = 3.1962/R;
rhol = @(r) (r./R).^2 +1;
al = sqrt(1.4*287*300);
vl = [1;1;0];
dmdr = @(r) -2.*dev.*r;
d2mdrdq = @(r) -k.*besselj(1,k.*r) - k.*(besselj(0,k.*R)/besseli(0,k.*R))*besseli(1,k.*r);
dn1dq_anal = @(r,t) (-d2mdrdq(r).*cos(t).*(sqrt((dmdr(r).^2)+1)-(dmdr(r).^2)./sqrt((dmdr(r).^2)+1)))./((dmdr(r).^2)+1);
dn2dq_anal = @(r,t) (-d2mdrdq(r).*sin(t).*(sqrt((dmdr(r).^2)+1)-(dmdr(r).^2)./sqrt((dmdr(r).^2)+1)))./((dmdr(r).^2)+1);
dn3dq_anal = @(r,t) (-(dmdr(r).*d2mdrdq(r))./sqrt((dmdr(r).^2)+1))./((dmdr(r).^2)+1);
n1_anal = @(r,t) (-dmdr(r).*cos(t))./(sqrt((dmdr(r).^2)+1)); 
n2_anal = @(r,t) (-dmdr(r).*sin(t))./(sqrt((dmdr(r).^2)+1)); 
n3_anal = @(r,t) 1./sqrt((dmdr(r).^2)+1);
abar = 1 / (besselj(0,0) - (besselj(0,k.*R)/besseli(0,k.*R))*besseli(0,0));
psi3_anal = @(r) abar.*(besselj(0,k.*r) - (besselj(0,k.*R)/besseli(0,k.*R))*besseli(0,k.*r));

fun = @(r,t) rhol(r) .* al .* ((dn1dq_anal(r,t) + dn2dq_anal(r,t)).*(n3_anal(r,t).*psi3_anal(h)) + (n1_anal(r,t) + n2_anal(r,t)).*(dn3dq_anal(r,t).*psi3_anal(r))) .* r;

dF1dq_anal = -1 * integral2(fun, 0, R, 0, 2*pi);

% Investigate local error within structural element
%L2_temp = [0, 1, 0, linspace(0.1,0.3,10), linspace(0,0.5,10), linspace(0.1,0.3,10), linspace(0,1,10), zeros(1,10)];
%L3_temp = [1, 0, 0, linspace(0,0.5,10), linspace(0.1,0.3,10), linspace(0.1,0.3,10), zeros(1,10), linspace(0,1,10)];
%L4_temp = 1 - L2_temp - L3_temp;
%psi_err = zeros(length(L2_temp),1);
%loc = zeros(length(L2_temp),3);
%for i = 1:length(L2_temp)
%  L1 = 0;
%  L2 = L2_temp(i);
%  L3 = L3_temp(i);
%  L4 = L4_temp(i);
%  loc_bary = [L1, L2, L3, L4];
%  loc(i,:) = tet_interp(loc1_elem_nodelocs, loc_bary);
%  x = loc(i,1);
%  y = loc(i,2);
%  r_loc = sqrt(x^2 + y^2);
%  psi_loc = tet_interp(psi1_nodes, loc_bary);
%  psi_loc_anal = psi3_anal(r_loc);
%  psi_err(i) = psi_loc_anal - psi_loc(3);
%end
%figure; hold on;
%scatter3(loc1_elem_nodelocs(2:4,1), loc1_elem_nodelocs(2:4,2), zeros(3,1), 'k');
%scatter3(loc(:,1), loc(:,2), psi_err(:), 24, psi_err(:));
%colormap('jet'); colorbar;

function [loc] = tri_loc(face_node_locs, interp_loc)
  % Linear triangle location from area coordinates
  % face_node_locs: Cartesian locations of each face node
  % interp_loc: location of point in triangle, area coordinates
  % loc: location of point in triangle, Cartesian coordinates

  loc = face_node_locs(1,:).*interp_loc(1) + face_node_locs(2,:).*interp_loc(2) + face_node_locs(3,:).*interp_loc(3);

end

function [quant] = tri_interp(interp_loc, data);

  quant = data(1)*interp_loc(1) + data(2)*interp_loc(2) + data(3)*interp_loc(3);

end

function [int] = quadrature(val1, val2, val3, val4, area)

  % Weights obtained from Zienkiewics and modified to ensure that transformed surface integral provides A = 1/2.
  w1 = -27/96;
  w2 = 25/96;
  w3 = 25/96;
  w4 = 25/96;

  int = w1*val1 + w2*val2 + w3*val3 + w4*val4;
  int = 2*area*int;

end

