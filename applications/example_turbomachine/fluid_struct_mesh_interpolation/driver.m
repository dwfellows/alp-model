
close all;
clear all;
clc;

% Import mesh interpolation utilities
addpath('../../../mesh_interpolation_utils');
addpath('../../../fluid_function_utils');

%% Construct information arrays
disp('Constructing information arrays.');

% Construct structural info array
% Note -- All files, except "hp_integration_pt_locs", must exist before executing this driver function
% The integration pt locations will be derived from this data and placed in the directory described below
% by this driver function if (perform_construct_int_locs == 1).
struct_info = ["../computation_data/structural_solution_data/hp_elems.txt", ...
               "../computation_data/structural_solution_data/hp_node_elem.txt", ...
               "../computation_data/structural_solution_data/hp_nodes.txt", ...
               "../computation_data/structural_solution_data/hp_face_nodes.txt", ...
               "../computation_data/structural_solution_data/hp_faces.txt", ...
               "../computation_data/structural_solution_data/hp_integration_pt_locs.txt", ...
               "../computation_data/structural_solution_data/mesh_ref_pts_hpt.csv"];

% Construct fluid info array
% Note -- All files, except "hp_point_face_data", must exist before executing this driver function
% The fluid point-face connectivity matrix will be derived from this data and placed in the directory described
% below by this driver function if (perform_face_to_pt == 1).
fluid_info = ["../computation_data/fluid_solution_data/hp_face_data.csv", ...
              "../computation_data/fluid_solution_data/hp_pt_data.csv", ...
              "../computation_data/fluid_solution_data/hp_point_face_data.csv", ...
              "../computation_data/fluid_solution_data/mesh_ref_pts_hpt.csv"];


%% Perform routines necessary for fluid-structure mesh interpolation

perform_face_to_pt = 1;
if (perform_face_to_pt == 1)
  disp('Deriving fluid face-point connectivity array.');
  message = face_to_pt(fluid_info);
end

perform_construct_int_locs = 1;
if (perform_construct_int_locs == 1)
  disp('Deriving integration point locations on structural mesh.');
  message = construct_integration_locs(struct_info, fluid_info);
end

perform_mesh_interpolation = 1;
if (perform_mesh_interpolation == 1)
  disp('Performing mesh-to-mesh interpolation of fluid mesh onto structural mesh.');
  message = fluid_mesh_interpolation(struct_info, fluid_info);
end

