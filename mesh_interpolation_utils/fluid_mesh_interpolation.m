% Perform mesh-to-mesh interpolation

function [message] = fluid_mesh_interpolation(struct_info, fluid_info)

  message = 0;

  % Read in structural data
  disp('Reading structural data...');
  elem_node = importdata(struct_info(1)); elem_node = elem_node.data;
  node_elem = dlmread(struct_info(2));
  nodes = dlmread(struct_info(3));
  structural_node_locs = nodes(:,2:4); structural_node_locs_orig = structural_node_locs;
  surface_node_data = dlmread(struct_info(4));
  struct_faces = dlmread(struct_info(5));
  
  % Read in fluid data
  fluid_face_conn = dlmread(fluid_info(1));
  fluid_solution_data = dlmread(fluid_info(2));
  fluid_point_face_data = dlmread(fluid_info(3));

  
  surf_points = fluid_solution_data(:,2:4);
  surf_normals = fluid_solution_data(:,9:11);
  rho_surf = fluid_solution_data(:,6);
  p_surf = fluid_solution_data(:,12);
  vel_rot_surf = fluid_solution_data(:,14:16);
  M_rot_surf = fluid_solution_data(:,8);
  U_rot_surf = fluid_solution_data(:,13);
  a_rot_surf = fluid_solution_data(:,7);
  
  % Fluid solution data is in following configuration:
  % 1 node # | 2 x | 3 y | 4 z | 5 A | 6 rho | 7 a | 8 M | 9 n_x | 10 n_y | 11 n_z | 12 p | 13 v_mag 
  % | 14 u | 15 v | 16 w
  
  %%%%%%% Align Fluid and Structural Grids
  
  % Read in structural data -- contains both integration point locations and the normal vectors at those locations
  % The way the scripts are currently written, the integration points must be transformed but the normals are computed post-transformation
  
  struct_points_data = dlmread(struct_info(6));
  struct_points = [struct_points_data(:,1), struct_points_data(:,2), struct_points_data(:,3);
                   struct_points_data(:,7), struct_points_data(:,8), struct_points_data(:,9);
                   struct_points_data(:,13), struct_points_data(:,14), struct_points_data(:,15)];
  
  struct_point_normals = [struct_points_data(:,4), struct_points_data(:,5), struct_points_data(:,6);
                          struct_points_data(:,10), struct_points_data(:,11), struct_points_data(:,12);
                          struct_points_data(:,16), struct_points_data(:,17), struct_points_data(:,18)];
  
  struct_points_fluidref = [struct_points_data(:,1), struct_points_data(:,2), struct_points_data(:,3);
                            struct_points_data(:,7), struct_points_data(:,8), struct_points_data(:,9);
                            struct_points_data(:,13), struct_points_data(:,14), struct_points_data(:,15)];
  
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
  struct_points = struct_points - structural_origin;
  
  
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
  points_temp = R*(surf_normals');
  surf_normals = points_temp';
  
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
  points_temp = R*(struct_points');
  struct_points = points_temp';
  
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
  struct_points = scale_factor .* struct_points;
  structural_vec1_pt = scale_factor .* structural_vec1_pt;
  structural_vec2_pt = scale_factor .* structural_vec2_pt;
  
  % Align structural orientation with fluid orientation
  a = structural_vec1 ./ norm(structural_vec1);
  b = fluid_vec1 ./ norm(fluid_vec1);
  R = rotation_matrix_align(a,b);
  points_temp = R*(struct_points');
  struct_points = points_temp';
  
  
  %%%%% Perform point projection and interpolation
  N = max(size(struct_points));
  interp_pressure = zeros(N,4);
  interp_density = zeros(N,4);
  for i = 1:N
  
    % Pull off integration point
    disp(i);
    struct_point = struct_points(i,:);
    struct_point_normal = struct_point_normals(i,:);
   
    % Find closest point on fluid mesh and containing cells
    loc_fluid_point = find_closest_fluid_point(struct_point, surf_points);
    loc_fluid_cells = fluid_point_face_data(loc_fluid_point, :);
  
    % Pick closest containing cell based on which normals are most closely aligned
    [closest_fluid_cell, computed_normal] = find_aligned_fluid_cell(struct_point_normal, loc_fluid_cells, fluid_face_conn, surf_points, surf_normals);
  
    % Project integration point into the cell
    closest_fluid_cell_nodes = fluid_face_conn(closest_fluid_cell, :) + 1;  % Nodes are zero-indexed
    closest_fluid_cell_nodelocs = surf_points(closest_fluid_cell_nodes, :);
    struct_point_proj = point_projection(struct_point, closest_fluid_cell_nodelocs);
  
    % Perform interpolation
    p_locs = p_surf(closest_fluid_cell_nodes);
    p_int = point_interpolation(struct_point_proj, closest_fluid_cell_nodelocs, p_locs);
  
    rho_locs = rho_surf(closest_fluid_cell_nodes);
    rho_int = point_interpolation(struct_point_proj, closest_fluid_cell_nodelocs, rho_locs);
  
    a_locs = a_rot_surf(closest_fluid_cell_nodes);
    a_int = point_interpolation(struct_point_proj, closest_fluid_cell_nodelocs, a_locs);
  
    velx_locs = vel_rot_surf(closest_fluid_cell_nodes, 1);
    vely_locs = vel_rot_surf(closest_fluid_cell_nodes, 2);
    velz_locs = vel_rot_surf(closest_fluid_cell_nodes, 3);
  
    velx_int = point_interpolation(struct_point_proj, closest_fluid_cell_nodelocs, velx_locs);
    vely_int = point_interpolation(struct_point_proj, closest_fluid_cell_nodelocs, vely_locs);
    velz_int = point_interpolation(struct_point_proj, closest_fluid_cell_nodelocs, velz_locs);
  
    % Save quantities for plotting
    interp_pressure(i,:) = [struct_point, p_int];
    interp_density(i,:) = [struct_point, rho_int];
    interp_quants(i,1) = 1;
    interp_quants(i,2) = rho_int;
    interp_quants(i,3) = p_int;
    interp_quants(i,4) = a_int;
    interp_quants(i,5) = velx_int;
    interp_quants(i,6) = vely_int;
    interp_quants(i,7) = velz_int;
  
  end 
  
  % Write interpolated quantities
  dlmwrite('./normal_vect_interp_quants.csv', interp_quants,'precision',8);
  message = 1;

end

