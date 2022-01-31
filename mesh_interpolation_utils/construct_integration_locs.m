
% Perform quadrature about structural element faces

%addpath('../../tools');

function [message] = construct_integration_locs(struct_info, fluid_info)

  message = 0;

  %%%% Read in structural data
  disp('Reading structural data...');
  elem_node = importdata(struct_info(1)); elem_node = elem_node.data;
  node_elem = dlmread(struct_info(2));
  nodes = dlmread(struct_info(3));
  structural_node_locs = nodes(:,2:4); structural_node_locs_orig = structural_node_locs;
  surface_node_data = dlmread(struct_info(4));
  struct_faces = dlmread(struct_info(5));
  
  
  %%%% Read in fluid data
  fluid_solution_data = dlmread(fluid_info(2));
  surf_points = fluid_solution_data(:,2:4);
  
  %%%% Align fluid and structural meshes
  % Define structural reference points
  structural_mesh_info = dlmread(struct_info(7));
  structural_origin = structural_mesh_info(1,:);
  structural_vec1_pt = structural_mesh_info(2,:);
  structural_vec2_pt = structural_mesh_info(3,:);
  
  % Transform points in structural mesh to points in fluid mesh
  fluid_mesh_info = dlmread(fluid_info(4));
  fluid_origin = fluid_mesh_info(1,:);
  fluid_vec1_pt = fluid_mesh_info(2,:);
  fluid_vec2_pt = fluid_mesh_info(3,:);
  
  
  % Utilize common point as origin
  surf_points = surf_points - fluid_origin;
  structural_node_locs = structural_node_locs - structural_origin;
  
  % Perform fluid rotation -- align top of turbine with z-axis
  fluid_vec1 = fluid_vec1_pt - fluid_origin;
  fluid_vec2 = fluid_vec2_pt - fluid_origin;
  fluid_surf_norm = cross(fluid_vec1, fluid_vec2);
  fluid_surf_norm = fluid_surf_norm ./ norm(fluid_surf_norm);
  
  a = fluid_surf_norm;
  b = [0,0,1];
  R = rotation_matrix_align(a,b);
  points_temp = R*(surf_points');
  surf_points = points_temp';
  
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
  b = [0,0,1];
  R = rotation_matrix_align(a,b);
  points_temp = R*(structural_node_locs');
  structural_node_locs = points_temp';
  
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
  struct_points = scale_factor .* structural_node_locs;
  structural_vec1_pt = scale_factor .* structural_vec1_pt;
  structural_vec2_pt = scale_factor .* structural_vec2_pt;
  
  % Align structural orientation with fluid orientation
  a = structural_vec1 ./ norm(structural_vec1);
  b = fluid_vec1 ./ norm(fluid_vec1);
  R = rotation_matrix_align(a,b);
  points_temp = R*(structural_node_locs');
  structural_node_locs = points_temp';
  
  
  %%%%%% Define quadrature
  % Define quadrature locations (integration points) -- utilize 3-pt rule as FEM data is second-order accurate
  % Rule from Alex and is contained in FEM textbooks
  int_loc1 = [2/3, 1/6, 1/6];
  int_loc2 = [1/6, 2/3, 1/6];
  int_loc3 = [1/6, 1/6, 2/3];
  
  %%%%%%%% Store integration point locations and output -- for struct_faces
  N = max(size(struct_faces));
  integration_locs = zeros(N,18);
  for i = 1:N
    disp(strcat(['Face ', num2str(i)]));
   
    % Grab face element number and nodes on surface
    loc_elem = struct_faces(i,1);
    loc_elem_surf_nodes = struct_faces(i,2:11);
  
    % Get all element nodes and locations
    loc_elem_nodes = elem_node(loc_elem, :);
    loc_elem_nodelocs = structural_node_locs(loc_elem_nodes,:);  % Transformed coordinates
    loc_elem_nodelocs_orig = structural_node_locs_orig(loc_elem_nodes,:);  % Untransformed coordinates
    if (nnz(loc_elem_surf_nodes(1:4)) == 3)
      loc_elem_surf_num = find(loc_elem_surf_nodes(1:4) == 0);
    else
      if (nnz(loc_elem_surf_nodes(5:7)) == 3)
        loc_elem_surf_num = 4;
      else
        if (loc_elem_surf_nodes(5) > 0)
          loc_elem_surf_num = 3;
        elseif (loc_elem_surf_nodes(6) == 0)
          loc_elem_surf_num = 1;
        else
          loc_elem_surf_num = 2;
        end
      end
    end
    loc_elem_surf = [1 2 3 4]; loc_elem_surf(loc_elem_surf_num) = [];
   
    %%% Int loc 1
    % Get location of integration point in Cartesian (global) and barycentric (local) coords
    location1_bary = zeros(1,4); 
    location1_bary(loc_elem_surf) = int_loc1;
    location1 = tet_interp(loc_elem_nodelocs, location1_bary);  % Transformed coordinates
    location1_orig = tet_interp(loc_elem_nodelocs_orig, location1_bary);  % Untransformed coordinates
    n1 = determine_normal(loc_elem_nodelocs, location1, location1_bary, loc_elem_surf_num);
  
    %%% Int loc 2
    % Get location of integration point in Cartesian (global) and barycentric (local) coords
    location2_bary = zeros(1,4); 
    location2_bary(loc_elem_surf) = int_loc2;
    location2 = tet_interp(loc_elem_nodelocs, location2_bary);  % Transformed coordinates
    location2_orig = tet_interp(loc_elem_nodelocs_orig, location2_bary);  % Untransformed coordinates
    n2 = determine_normal(loc_elem_nodelocs, location2, location2_bary, loc_elem_surf_num);
  
    %%% Int loc 3
    % Get location of integration point in Cartesian (global) and barycentric (local) coords
    location3_bary = zeros(1,4); 
    location3_bary(loc_elem_surf) = int_loc3;
    location3 = tet_interp(loc_elem_nodelocs, location3_bary);  % Transformed coordinates
    location3_orig = tet_interp(loc_elem_nodelocs_orig, location3_bary);  % Untransformed coordinates
    n3 = determine_normal(loc_elem_nodelocs, location2, location3_bary, loc_elem_surf_num);
  
    integration_locs(i,:) = [location1, n1, location2, n2, location3, n3];  % Transformed coordinates
    integration_locs(i,:) = [location1_orig, n1, location2_orig, n2, location3_orig, n3];
  end
  dlmwrite('../computation_data/structural_solution_data/hp_integration_pt_locs.txt', integration_locs);
  disp('Integration points determined and written.');
  message = 1;

end

