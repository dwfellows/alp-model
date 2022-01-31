
% Perform quadrature about structural element faces

close all;
clear all;
clc;

addpath('../../../fluid_function_utils');

%% Import fluid data -- all quantities, from surface
% 1 node num | 2 X [m] | 3 Y [m] | 4 Z [m] | 5 p_a [Pa] | 6 A [m^2] | 7 \rho [kg m^-3] | 8 F_x [N] | 9 F_y [N] 
% 10 F_z [N] | 11 a [m s^-1] | 12 Mach # | 13 Mach # Stn | 14 normal_x | 15 normal_y | 16 normal_z | 17 v_u stn [m s^1]
% 18 v_v stn [m s^-1] | 19 v_w stn [m s^-1] | 20 v_u [m s^-1] | 21 v_v [m s^-1] | 22 v_w [m s^-1]
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
fluid_node_faces = dlmread('../fem_data/circular_fluid/fluid_node_faces.txt');

% Define plate parameters
R = 0.1; h = 0.0001;
area = pi*R^2;

% Construct deformation function
l1 = 0.0005; l2 = 1/l1;
%l1 = 1; l2 = 1/l1;
besselfuns = @(l) besselj(0,l)*besseli(1,l) + besseli(0,l)*besselj(1,l);
lambda = fsolve(besselfuns, 3.1962);
k = lambda/R;
abar = l1 ./ (besselj(0,0) - (besselj(0,k.*R)/besseli(0,k.*R)).*besseli(0,0));
psi3_anal = @(f) abar.*(besselj(0,k.*f) - (besselj(0,k.*R)/besseli(0,k.*R))*besseli(0,k.*f));

% Define analytical values
%rhol = @(r,t) ((r./R).^3).*(sin(wrapTo2Pi(t)./2)) +1;
rhol = @(r,t) ((r./R).^3).*(sin(wrapTo2Pi(t))) + 1;
%rhol = @(r,t) (t.^0).*((r/R).^2) + 1
%rhol = @(r,t) 1.*(r.^0).*(t.^0) + 1;
al = @(r,t) sqrt(1.4*287*300) .* ( t.^0 ) .* (r.^0);
vl = @(r,t) [1;1;0] .* (r.^0) .* (t.^0);
dmdr = @(r) 0*(r.^0);
d2mdrdq = @(r) abar.*(-k.*besselj(1,k.*r) - k.*(besselj(0,k.*R)/besseli(0,k.*R))*besseli(1,k.*r));
dn1dq_anal = @(r,t) (-d2mdrdq(r).*cos(t).*(sqrt((dmdr(r).^2)+1)-(dmdr(r).^2)./sqrt((dmdr(r).^2)+1)))./((dmdr(r).^2)+1);
dn2dq_anal = @(r,t) (-d2mdrdq(r).*sin(t).*(sqrt((dmdr(r).^2)+1)-(dmdr(r).^2)./sqrt((dmdr(r).^2)+1)))./((dmdr(r).^2)+1);
dn3dq_anal = @(r,t) (-(dmdr(r).*d2mdrdq(r))./sqrt((dmdr(r).^2)+1))./((dmdr(r).^2)+1);
n1_anal = @(r,t) (-dmdr(r).*cos(t))./(sqrt((dmdr(r).^2)+1)); 
n2_anal = @(r,t) (-dmdr(r).*sin(t))./(sqrt((dmdr(r).^2)+1)); 
n3_anal = @(r,t) 1./sqrt((dmdr(r).^2)+1);

% Define quadrature locations (integration points) -- utilize 3-pt rule as FEM data is second-order accurate
% Rule from Alex and is contained in FEM textbooks
int_loc1 = [2/3, 1/6, 1/6];
int_loc2 = [1/6, 2/3, 1/6];
int_loc3 = [1/6, 1/6, 2/3];

dF1dq_approx = [];
elem_h_list = [];

%% Import structural data -- all quantities
for l = [6 7 1 4 2 5 3]
  if (l==1)
    elem_node = dlmread('../fem_data/circular_coarse/coarse.csv');
    node_elem = dlmread('../fem_data/circular_coarse/coarse_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_coarse/coarse_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_coarse/coarse_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_coarse/coarse_zdir_deformation.txt');
    structural_node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_coarse/coarse_surface.txt');
    struct_faces = dlmread('../fem_data/circular_coarse/struc_data/coarse_faces.txt');
  elseif (l==2)
    elem_node = dlmread('../fem_data/circular_mid/mid.csv');
    node_elem = dlmread('../fem_data/circular_mid/mid_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_mid/mid_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_mid/mid_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_mid/mid_zdir_deformation.txt');
    structural_node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_mid/mid_surface.txt');
    struct_faces = dlmread('../fem_data/circular_mid/struc_data/mid_faces.txt');
  elseif (l==3)
    elem_node = dlmread('../fem_data/circular_fine/fine.csv');
    node_elem = dlmread('../fem_data/circular_fine/fine_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_fine/fine_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_fine/fine_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_fine/fine_zdir_deformation.txt');
    structural_node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_fine/fine_surface.txt');
    struct_faces = dlmread('../fem_data/circular_fine/struc_data/fine_faces.txt');
  elseif (l==4)
    elem_node = importdata('../fem_data/circular_cm/cm_elements.txt');
    elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/circular_cm/cm_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_cm/cm_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_cm/cm_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_cm/cm_zdir_deformation.txt');
    structural_node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_cm/cm_surf_nodes.txt');
    struct_faces = dlmread('../fem_data/circular_cm/struc_data/cm_faces.txt');
  elseif (l==5)
    elem_node = importdata('../fem_data/circular_fm/fm_elements.txt');
    elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/circular_fm/fm_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_fm/fm_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_fm/fm_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_fm/fm_zdir_deformation.txt');
    structural_node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_fm/fm_surface.txt');
    struct_faces = dlmread('../fem_data/circular_fm/struc_data/fm_faces.txt');
  elseif (l==6)
    elem_node = importdata('../fem_data/new_circ_data/circ_c1/circ_c1_elems.txt');
    elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/new_circ_data/circ_c1/circ_c1_node_elem.txt');
    xdir_deform = dlmread('../fem_data/new_circ_data/circ_c1/xdir_deform.txt');
    ydir_deform = dlmread('../fem_data/new_circ_data/circ_c1/ydir_deform.txt');
    zdir_deform = dlmread('../fem_data/new_circ_data/circ_c1/zdir_deform.txt');
    structural_node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/new_circ_data/circ_c1/circ_c1_face_nodes.txt');
    struct_faces = dlmread('../fem_data/new_circ_data/circ_c1/circ_c1_faces.txt');
  elseif (l==7)
    elem_node = importdata('../fem_data/new_circ_data/circ_c2/circ_c2_elems.txt');
    elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/new_circ_data/circ_c2/circ_c2_node_elem.txt');
    xdir_deform = dlmread('../fem_data/new_circ_data/circ_c2/xdir_deform.txt');
    ydir_deform = dlmread('../fem_data/new_circ_data/circ_c2/ydir_deform.txt');
    zdir_deform = dlmread('../fem_data/new_circ_data/circ_c2/zdir_deform.txt');
    structural_node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/new_circ_data/circ_c2/circ_c2_face_nodes.txt');
    struct_faces = dlmread('../fem_data/new_circ_data/circ_c2/circ_c2_faces.txt');
  end

  % Construct FEM data
  def_shape_function = [xdir_deform(:,5), ydir_deform(:,5), zdir_deform(:,5)];
  def_shape_function_norms = vecnorm(def_shape_function'); max_def = l2*max(def_shape_function_norms);
  def_shape_function = def_shape_function ./ max_def;

%  % Impose analytical data directly
%  for j = 1:max(size(xdir_deform))
%    x = structural_node_locs(j,1); y = structural_node_locs(j,2);
%    r = sqrt(x^2 + y^2); t = atan2(y,x); t = wrapTo2Pi(t);
%    def_shape_function(j,:) = [0 0 psi3_anal(r)];
%  end

  %% Construct flow functions
  N = max(size(fluid_node_locs));
  fluid_rho = zeros(N,1); fluid_a = zeros(N,1);
  fluid_vu = zeros(N,1); fluid_vv = zeros(N,1); fluid_vw = zeros(N,1);
  for j = 1:max(size(fluid_node_locs))
    x = fluid_node_locs(j,1); y = fluid_node_locs(j,2); 
    rl = sqrt(x^2 + y^2); tl = atan2(y,x);
    fluid_rho(j) = rhol(rl, wrapTo2Pi(tl));
    fluid_a(j) = al(rl, wrapTo2Pi(tl));
    temp_v = vl(rl,tl);
    fluid_vu(j) = temp_v(1); fluid_vv(j) = temp_v(2);
    fluid_vw(j) = temp_v(3);
  end

  % Compute characteristic size of element
  N = max(size(struct_faces));
  elem_a = area / N;
  elem_h_list = [elem_h_list, sqrt(elem_a*4/sqrt(3))];

  %% Compute dF1dq on local mesh
  disp(strcat(['Mesh Case Number: ', num2str(l)]));
  disp(strcat(['Number of faces: ', num2str(N)]));
  dF1dq = 0;
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
  
    %%% Go over each int loc and form integrand
    %%% Int loc 1
    % Get location of integration point in Cartesian (global) and barycentric (local) coords
    location1_bary = zeros(1,4); 
    location1_bary(loc_elem_surf) = int_loc1;
    location1 = tet_interp(loc_elem_nodelocs, location1_bary);
  
    rl = sqrt(location1(1)^2 + location1(2)^2);
    tl = atan2(location1(2), location1(1));
  
    % Interpolate deformation, normal and normal derivative at integration point
    psi_nodes = def_shape_function(loc_elem_nodes,:);
    psi1 = tet_interp(psi_nodes, location1_bary);
    n1 = determine_normal(loc_elem_nodelocs, location1, location1_bary, loc_elem_surf_num);
    dn1dq = determine_normal_deriv(loc_elem_nodelocs, psi_nodes, location1, location1_bary, loc_elem_surf_num);
    [dxdLi1, dxdLj1, ~] = def_derivs(loc_elem_nodelocs, psi_nodes, location1, location1_bary, loc_elem_surf_num);
    A1 = cross(dxdLi1, dxdLj1); A1 = norm(A1);
  
    % Determine fluid face containing integration point and face coordinates of structural point in the face
    [loc1_fluid_face, loc1_struct_bary] = determine_fluid_face(face_conn, fluid_node_surf_circular, fluid_node_faces, location1);
    loc1_fluid_face_nodes = face_conn(loc1_fluid_face, :) + 1;
  
    % Construct local fluid face nodes for interpolation of rho, a, velocity
    loc1_rho_nodes = fluid_rho(loc1_fluid_face_nodes);
    loc1_a_nodes = fluid_a(loc1_fluid_face_nodes);
    loc1_vu_nodes = fluid_vu(loc1_fluid_face_nodes);
    loc1_vv_nodes = fluid_vv(loc1_fluid_face_nodes);
    loc1_vw_nodes = fluid_vw(loc1_fluid_face_nodes);
  
    % Perform triangular interpolation
    rho1 = tri_interp(loc1_struct_bary, loc1_rho_nodes);
    a1 = tri_interp(loc1_struct_bary, loc1_a_nodes);
    vu1 = tri_interp(loc1_struct_bary, loc1_vu_nodes);
    vv1 = tri_interp(loc1_struct_bary, loc1_vv_nodes);
    vw1 = tri_interp(loc1_struct_bary, loc1_vw_nodes);
    v1 = [vu1; vv1; vw1];
  
    %%% Int loc 2
    % Get location of integration point in Cartesian (global) and barycentric (local) coords
    location2_bary = zeros(1,4);
    location2_bary(loc_elem_surf) = int_loc2;
    location2 = tet_interp(loc_elem_nodelocs, location2_bary);
  
    rl = sqrt(location2(1)^2 + location2(2)^2);
    tl = atan2(location2(2), location2(1));
  
    % Interpolate deformation, normal, and normal derivative at integration point
    psi2 = tet_interp(psi_nodes, location2_bary);
    n2 = determine_normal(loc_elem_nodelocs, location2, location2_bary, loc_elem_surf_num);
    dn2dq = determine_normal_deriv(loc_elem_nodelocs, psi_nodes, location2, location2_bary, loc_elem_surf_num);
    [dxdLi2, dxdLj2, ~] = def_derivs(loc_elem_nodelocs, psi_nodes, location2, location2_bary, loc_elem_surf_num);
    A2 = cross(dxdLi2, dxdLj2); A2 = norm(A2);
  
    % Determine fluid face containing integration point and face coordinates of structural point in the face
    [loc2_fluid_face, loc2_struct_bary] = determine_fluid_face(face_conn, fluid_node_surf_circular, fluid_node_faces, location2);
    loc2_fluid_face_nodes = face_conn(loc2_fluid_face, :) + 1;
  
    % Construct local fluid face nodes for interpolation of rho, a, velocity
    loc2_rho_nodes = fluid_rho(loc2_fluid_face_nodes);
    loc2_a_nodes = fluid_a(loc2_fluid_face_nodes);
    loc2_vu_nodes = fluid_vu(loc2_fluid_face_nodes);
    loc2_vv_nodes = fluid_vv(loc2_fluid_face_nodes);
    loc2_vw_nodes = fluid_vw(loc2_fluid_face_nodes);
  
    % Perform triangular interpolation
    rho2 = tri_interp(loc2_struct_bary, loc2_rho_nodes);
    a2 = tri_interp(loc2_struct_bary, loc2_a_nodes);
    vu2 = tri_interp(loc2_struct_bary, loc2_vu_nodes);
    vv2 = tri_interp(loc2_struct_bary, loc2_vv_nodes);
    vw2 = tri_interp(loc2_struct_bary, loc2_vw_nodes);
    v2 = [vu2; vv2; vw2];
  
    %%% Int loc 3
    % Get location of integration point in Cartesian (global) and barycentric (local) coords
    location3_bary = zeros(1,4);
    location3_bary(loc_elem_surf) = int_loc3;
    location3 = tet_interp(loc_elem_nodelocs, location3_bary);
  
    rl = sqrt(location3(1)^2 + location3(2)^2);
    tl = atan2(location3(2), location3(1));
  
    % Interpolate deformation, normal, and normal derivative at integration point
    psi3 = tet_interp(psi_nodes, location3_bary);
    n3 = determine_normal(loc_elem_nodelocs, location3, location3_bary, loc_elem_surf_num);
    dn3dq = determine_normal_deriv(loc_elem_nodelocs, psi_nodes, location3, location3_bary, loc_elem_surf_num);
    [dxdLi3, dxdLj3, ~] = def_derivs(loc_elem_nodelocs, psi_nodes, location3, location3_bary, loc_elem_surf_num);
    A3 = cross(dxdLi3, dxdLj3); A3 = norm(A3);
  
    % Determine fluid face containing integration point and face coordinates of structural point in the face
    [loc3_fluid_face, loc3_struct_bary] = determine_fluid_face(face_conn, fluid_node_surf_circular, fluid_node_faces, location3);
    loc3_fluid_face_nodes = face_conn(loc3_fluid_face, :) + 1;
  
    % Construct local fluid face nodes for interpolation of rho, a, velocity
    loc3_rho_nodes = fluid_rho(loc3_fluid_face_nodes);
    loc3_a_nodes = fluid_a(loc3_fluid_face_nodes);
    loc3_vu_nodes = fluid_vu(loc3_fluid_face_nodes);
    loc3_vv_nodes = fluid_vv(loc3_fluid_face_nodes);
    loc3_vw_nodes = fluid_vw(loc3_fluid_face_nodes);
  
    % Perform triangular interpolation
    rho3 = tri_interp(loc3_struct_bary, loc3_rho_nodes);
    a3 = tri_interp(loc3_struct_bary, loc3_a_nodes);
    vu3 = tri_interp(loc3_struct_bary, loc3_vu_nodes);
    vv3 = tri_interp(loc3_struct_bary, loc3_vv_nodes);
    vw3 = tri_interp(loc3_struct_bary, loc3_vw_nodes);
    v3 = [vu3; vv3; vw3];
  
  
    % Compute integrands at integration points
    dF1dq_int1 = -rho1*a1*(dot(v1,dn1dq')*dot(n1',psi1') + dot(v1,n1')*dot(dn1dq',psi1'))*A1;
    dF1dq_int2 = -rho2*a2*(dot(v2,dn2dq')*dot(n2',psi2') + dot(v2,n2')*dot(dn2dq',psi2'))*A2;
    dF1dq_int3 = -rho3*a3*(dot(v3,dn3dq')*dot(n3',psi3') + dot(v3,n3')*dot(dn3dq',psi3'))*A3;
  
    % Perform quadrature
    dF1dq_loc = quadrature(dF1dq_int1, dF1dq_int2, dF1dq_int3);
    dF1dq = dF1dq + dF1dq_loc;
  
  end

  dF1dq_approx = [dF1dq_approx, dF1dq];

end

% Compute solution analytically
fun = @(r,t) rhol(r,t) .* al(r,t) .* ((dn1dq_anal(r,t) + dn2dq_anal(r,t)).*(n3_anal(r,t).*psi3_anal(h)) + (n1_anal(r,t) + n2_anal(r,t)).*(dn3dq_anal(r,t).*psi3_anal(r))) .* r;
dF1dq_anal = -1 * integral2(fun, 0, R, 0, 2*pi);

% Compute difference metric
if (abs(dF1dq_anal) > 1e-10)
  diffs = abs(dF1dq_approx - dF1dq_anal) ./ abs(dF1dq_anal);
else
  diffs = abs(dF1dq_approx - dF1dq_anal);
end

% Plot errors on log-log plot
c = diffs(end) / (elem_h_list(end)^2);
yc = c.*elem_h_list.^2;

figure;
loglog(elem_h_list, diffs, '-ob'); hold on;
loglog(elem_h_list, yc, '--r'); grid on;
legend('Numerical', '$\mathcal{O}(h)$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 18);
xlabel('$h$', 'interpreter', 'latex', 'fontsize', 18);
ylabel('$|| I_n - I_a ||$', 'interpreter', 'latex', 'fontsize', 18);
title('Numerical Approximation of Imposed $\frac{dF_1}{dq_1} (q_1^e = 0)$', 'interpreter', 'latex', 'fontsize', 18);

% Compare self convergence
if (abs(dF1dq_approx(end)) > 1e-10)
  diffs = abs(dF1dq_approx(1:end-1) - dF1dq_approx(end)) / abs(dF1dq_approx(end));
else
  diffs = abs(dF1dq_approx(1:end-1) - dF1dq_approx(end));
end

c = diffs(end) / elem_h_list(end-1);
yc = c.*elem_h_list(1:end-1);

figure;
loglog(elem_h_list(1:end-1), diffs, '-ob'); hold on;
loglog(elem_h_list(1:end-1), yc, '--r'); grid on;
legend('Numerical', '$\mathcal{O}(h)$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 18);
xlabel('$h$', 'interpreter', 'latex', 'fontsize', 18);
ylabel('$|| I_n - I_{n,f} || / || I_{n,f} ||$', 'interpreter', 'latex', 'fontsize', 18);
title('FEM Self Convergence of Approximation of Imposed $\frac{dF_1}{dq_1} (q_1^e = 0)$', 'interpreter', 'latex', 'fontsize', 18);


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

