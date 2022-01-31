

% Driver function to compute stability of given high-pressure turbine configurations

close all;
clear all;
clc;

addpath('../../fluid_function_utils');
addpath('../../struct_energies_utils');

%% Import fluid data -- all quantities, hp_turbine default surface
disp('Reading fluid data...');
fluid_data = dlmread('./computation_data/fluid_solution_data/lp_wall_mesh_data.csv');
fluid_node_locs = fluid_data(:,2:4);

% Clear fluid data to regain some memory
clear fluid_data;


%% Read in structural data
disp('Reading structural data...');
elem_node = importdata('./computation_data/structural_solution_data/lp_elems.txt'); elem_node = elem_node.data;
node_elem = dlmread('./computation_data/structural_solution_data/lp_node_elem.txt');
nodes = dlmread('./computation_data/structural_solution_data/lp_nodes.txt');
structural_node_locs = nodes(:,2:4);  
surface_node_data = dlmread('./computation_data/structural_solution_data/lp_face_nodes.txt');
surface_node_locs = surface_node_data(:,2:4);
struct_faces = dlmread('./computation_data/structural_solution_data/lp_faces.txt');


%% Align coordinate systems
disp('Aligning coordinate systems')

% Define structural reference points
structural_mesh_info = dlmread('./computation_data/structural_solution_data/mesh_ref_pts_lpt.csv');
structural_origin = structural_mesh_info(1,:);
structural_vec1_pt = structural_mesh_info(2,:);
structural_vec2_pt = structural_mesh_info(3,:);

% Define fluid reference points -- in order origin, vec1_pt, vec2_pt
fluid_mesh_info = dlmread('./computation_data/fluid_solution_data/mesh_ref_pts_lpt.csv');
fluid_origin = fluid_mesh_info(1,:);
fluid_vec1_pt = fluid_mesh_info(2,:);
fluid_vec2_pt = fluid_mesh_info(3,:);


% Utilize common point as origin
fluid_node_locs = fluid_node_locs - fluid_origin;
structural_node_locs = structural_node_locs - structural_origin;
surface_node_locs = surface_node_locs - structural_origin;


% Perform fluid rotation -- align top of turbine with z-axis
fluid_vec1 = fluid_vec1_pt - fluid_origin;
fluid_vec2 = fluid_vec2_pt - fluid_origin;
fluid_surf_norm = cross(fluid_vec1, fluid_vec2);
fluid_surf_norm = fluid_surf_norm ./ norm(fluid_surf_norm);

a = fluid_surf_norm;
b = [0,0,1];
R = rotation_matrix_align(a,b);   Rf = R;
points_temp = R*(fluid_node_locs');
fluid_node_locs = points_temp';

fluid_vec1 = R*(fluid_vec1');
fluid_vec1 = fluid_vec1';
fluid_vec2 = R*(fluid_vec2');
fluid_vec2 = fluid_vec2';
fluid_vec1_pt = fluid_vec1;
fluid_vec2_pt = fluid_vec2;

% Perform structural rotation -- align top of turbine with z-axis
structural_vec1 = structural_vec1_pt - structural_origin;
structural_vec2 = structural_vec2_pt - structural_origin;
structural_surf_norm = cross(structural_vec2, structural_vec1);
structural_surf_norm = structural_surf_norm ./ norm(structural_surf_norm);

a = structural_surf_norm;
b = [0,1,0];
R = rotation_matrix_align(a,b);
points_temp = R*(structural_node_locs');
structural_node_locs = points_temp';
points_temp = R*(surface_node_locs');
surface_node_locs = points_temp';

structural_vec1 = R*(structural_vec1');
structural_vec1 = structural_vec1';
structural_vec2 = R*(structural_vec2');
structural_vec2 = structural_vec2';
structural_vec1_pt = structural_vec1;
structural_vec2_pt = structural_vec2;
structural_surf_norm = R*(structural_surf_norm');
structural_surf_norm = structural_surf_norm';


a = structural_surf_norm;
b = [0,0,-1];
R = rotation_matrix_align(a,b);
points_temp = R*(structural_node_locs');
structural_node_locs = points_temp';
points_temp = R*(surface_node_locs');
surface_node_locs = points_temp';

structural_vec1 = R*(structural_vec1');
structural_vec1 = structural_vec1';
structural_vec2 = R*(structural_vec2');
structural_vec2 = structural_vec2';
structural_vec1_pt = structural_vec1;
structural_vec2_pt = structural_vec2;

% Scale structural points by a factor to make them lie on fluid points
struct_norm = norm(structural_vec1_pt);
fluid_norm = norm(fluid_vec1_pt);
scale_factor = fluid_norm / struct_norm;
structural_node_locs = scale_factor .* structural_node_locs;
surface_node_locs = scale_factor .* surface_node_locs;
structural_vec1_pt = scale_factor .* structural_vec1_pt;
structural_vec2_pt = scale_factor .* structural_vec2_pt;

% Align structural orientation with fluid orientation
a = structural_vec1 ./ norm(structural_vec1);
b = fluid_vec1 ./ norm(fluid_vec1);
R = rotation_matrix_align(a,b);
points_temp = R*(structural_node_locs');
structural_node_locs = points_temp';
points_temp = R*(surface_node_locs');
surface_node_locs = points_temp';

surface_node_data = [surface_node_data(:,1), surface_node_locs];



%% Choose deformation function based on choice of run priority to analyze

% Define material properties
import_mat_props = dlmread('./computation_data/structural_solution_data/mat_props_lpt.csv');
E = import_mat_props(1); % [Pa]
nu = import_mat_props(2); % [-]
rho = import_mat_props(3); % [kg m^-3]
mat_props = [E, nu, rho];

full_mode_data = readtable('./computation_data/run_data_lpt.dat');
run_priority = input('Which run priority (10/12/15/16) is to be analyzed?');
mode_type = ' ';
while (mode_type == ' ')
  if ((run_priority == 10) | (run_priority == 12) | (run_priority == 15) | (run_priority == 16))
    subset = full_mode_data(run_priority,:);
    mode_type = table2array(subset(1,1)); mode_type = mode_type{1};
    fluid_data_fname = table2array(subset(1,2)); fluid_data_fname = fluid_data_fname{1};
    fluid_data = dlmread(fluid_data_fname);
    T1 = table2array(subset(1,3));
    U1 = table2array(subset(1,4));
    deform_data_fname = table2array(subset(1,5)); deform_data_fname = deform_data_fname{1};

    disp('Computing structural energy coefficients for relevant mode shape.');
    deform_data = importdata(deform_data_fname);
    xdir_deform = deform_data(:,1:5);
    ydir_deform = [deform_data(:,1:4), deform_data(:,6)];
    zdir_deform = [deform_data(:,1:4), deform_data(:,7)];
  else
    disp('Cannot be analyzed. Enter a different number.');
    run_priority = input('Which run priority (10/12/15/16) is to be analyzed?');
  end
  if (mode_type ~= ' ')
    disp('Reading in pre-interpolated fluid data at structural integration points.');
  end
end


% Construct deformation function to have defined amplitude
l1 = 0.001; l2 = 1/l1;
def_shape_function = [xdir_deform(:,5), ydir_deform(:,5), zdir_deform(:,5)];
def_shape_function_norms = vecnorm(def_shape_function'); max_def = l2*max(def_shape_function_norms);
def_shape_function = def_shape_function ./ max_def;


%% Define quadrature and perform surface integrations
int_loc1 = [2/3, 1/6, 1/6];
int_loc2 = [1/6, 2/3, 1/6];
int_loc3 = [1/6, 1/6, 2/3];

disp('Beginning integration.');
N = max(size(struct_faces));
dF1dq = 0;
G11 = 0;
disp(strcat(['Number of faces: ', num2str(N)]));
int_point_counter = 1;
disp('Press any button to begin.'); pause;
for i = 1:N
  disp(strcat(['Face ', num2str(i)]));

  % Grab face element number and nodes on surface
  loc_elem = struct_faces(i,1);
  loc_elem_surf_nodes = struct_faces(i,2:11);

  % Get all element nodes and locations
  loc_elem_nodes = elem_node(loc_elem, :);
  loc_elem_nodelocs = structural_node_locs(loc_elem_nodes,:);
  if (nnz(loc_elem_surf_nodes(1:4)) == 3)
    loc_elem_surf_num = find(loc_elem_surf_nodes(1:4) == 0);
  else
    if (nnz(loc_elem_surf_nodes(5:7)) == 3)
      loc_elem_surf_num = 4;
    else
      if (loc_elem_surf_nodes(5) > 0)
        loc_elem_surf_num = 3;
      elseif (loc_elem_surf_nodes(6) > 0)
        loc_elem_surf_num = 1;
      else
        loc_elem_surf_num = 2;
      end
    end
  end
  loc_elem_surf = [1 2 3 4]; loc_elem_surf(loc_elem_surf_num) = [];

  %%% Go over each int loc and form integrand
  %%% Int loc 1
  % Get location of integration point in Cartesian (global) and barycentric (local) coords
  location1_bary = zeros(1,4); 
  location1_bary(loc_elem_surf) = int_loc1;
  location1 = tet_interp(loc_elem_nodelocs, location1_bary);

  % Interpolate deformation, normal and normal derivative at integration point
  psi_nodes = def_shape_function(loc_elem_nodes,:);
  psi1 = tet_interp(psi_nodes, location1_bary);
  n1 = determine_normal(loc_elem_nodelocs, location1, location1_bary, loc_elem_surf_num);
  dn1dq = determine_normal_deriv(loc_elem_nodelocs, psi_nodes, location1, location1_bary, loc_elem_surf_num);
  [dxdLi1, dxdLj1, ~] = def_derivs(loc_elem_nodelocs, psi_nodes, location1, location1_bary, loc_elem_surf_num);
  A1 = cross(dxdLi1, dxdLj1); A1 = norm(A1);

  % Read in pre-interpolated fluid data
  rho1 = fluid_data(int_point_counter, 2);
  a1 = fluid_data(int_point_counter, 4);
  vu1 = fluid_data(int_point_counter, 5);
  vv1 = fluid_data(int_point_counter, 6);
  vw1 = fluid_data(int_point_counter, 7);
  v1 = [vu1; vv1; vw1]; 
  v1 = Rf*v1; % Rotate fluid vector data (to match rotation of mesh) -- IS THIS STILL NECESSARY?
  int_point_counter = int_point_counter + 1;


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

  % Read in pre-interpolated fluid data
  rho2 = fluid_data(int_point_counter, 2);
  a2 = fluid_data(int_point_counter, 4);
  vu2 = fluid_data(int_point_counter, 5);
  vv2 = fluid_data(int_point_counter, 6);
  vw2 = fluid_data(int_point_counter, 7);
  v2 = [vu2; vv2; vw2]; 
  v2 = Rf*v2; % Rotate fluid vector data (to match rotation of mesh) -- IS THIS STILL NECESSARY?
  int_point_counter = int_point_counter + 1;


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

  % Read in pre-interpolated fluid data
  rho3 = fluid_data(int_point_counter, 2);
  a3 = fluid_data(int_point_counter, 4);
  vu3 = fluid_data(int_point_counter, 5);
  vv3 = fluid_data(int_point_counter, 6);
  vw3 = fluid_data(int_point_counter, 7);
  v3 = [vu3; vv3; vw3]; 
  v3 = Rf*v3; % Rotate fluid vector data (to match rotation of mesh) -- IS THIS STILL NECESSARY
  int_point_counter = int_point_counter + 1;


  % Compute integrands at integration points
  dF1dq_int1 = -rho1*a1*(dot(v1,-1.*dn1dq')*dot(n1',psi1'))*A1;
  dF1dq_int2 = -rho2*a2*(dot(v2,-1.*dn2dq')*dot(n2',psi2'))*A2;
  dF1dq_int3 = -rho3*a3*(dot(v2,-1.*dn3dq')*dot(n3',psi3'))*A3;

  G11_int1 = -rho1*a1*(dot(n1',psi1')*dot(n1',psi1'))*A1;
  G11_int2 = -rho2*a2*(dot(n2',psi2')*dot(n2',psi2'))*A2;
  G11_int3 = -rho3*a3*(dot(n3',psi3')*dot(n3',psi3'))*A3;

  % Perform quadrature
  dF1dq_loc = quadrature(dF1dq_int1, dF1dq_int2, dF1dq_int3);
  G11_loc = quadrature(G11_int1, G11_int2, G11_int3);
  dF1dq = dF1dq + dF1dq_loc;
  G11 = G11 + G11_loc;

end

% Compute stability
syms A
syms p
A(1,1) = T1*p^2 + (-G11)*p + (U1 - dF1dq);
charpoly = det(A) == 0;
B = double(solve(charpoly));
disp(B);

omega_n = sqrt(U1 / T1);
omega_n_mod = sqrt((U1 - dF1dq) / T1);
f_n_mod = (1/(2*pi))*omega_n_mod;
zeta = (-G11) / (2*T1*omega_n);
omega_d = omega_n*sqrt(1-zeta^2);
f_d = (1/(2*pi))*omega_d;
disp(strcat(['Modified natural frequency, f_n_mod: ', num2str(f_n_mod)]));
disp(strcat(['Damped frequency, f_d: ', num2str(f_d)]));
disp(strcat(['Damping ratio, zeta: ', num2str(zeta)]));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration and Projection Routines

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

function [R] = rotation_matrix_align(a,b)

  % Rotation matrix which aligns a onto b

  v = cross(a,b); c = dot(a,b);
  vx = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
  R = eye(3) + vx + vx*vx.*(1/(1+c));

end
