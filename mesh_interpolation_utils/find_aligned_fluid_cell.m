function [closest_fluid_cell, computed_normal] = find_aligned_fluid_cell(struct_point_normal, loc_fluid_cells, fluid_face_conn, surf_points, surf_normals)

  % Return fluid cell with normal most aligned to that of structural integration point

  % Derive fluid mesh nodes belonging to each fluid cell
  end_ind = find(loc_fluid_cells==0, 1, 'first');
  loc_fluid_cells = loc_fluid_cells(1:end_ind-1);
  loc_fluid_cells_nodes = fluid_face_conn(loc_fluid_cells, :) + 1;  % Nodes are zero-indexed

  % Construct error in alignment at each cell
  Nt = max(size(loc_fluid_cells));
  alignment_errors = zeros(Nt,1);
  cell_normals = zeros(Nt,3);
  for j = 1:Nt
    % Construct normal on cell
    fluid_mesh_nodes = loc_fluid_cells_nodes(j,:);
    fluid_mesh_nodelocs = surf_points(fluid_mesh_nodes, :);
    n = plane_normal(fluid_mesh_nodelocs);
    loc_fluid_mesh_normal = surf_normals(fluid_mesh_nodes(1),:);
    align_metric = dot(n, loc_fluid_mesh_normal);

    % Ensure outward pointing normal
    if (align_metric < 0)
      n = -1.*n;
    end

    % Save normal
    cell_normals(j,:) = n;

    % Construct error metric
    diff = n - struct_point_normal;
    alignment_errors(j) = norm(diff);

  end

  % Find cell with smallest alignment error
  [min_val, min_ind] = min(alignment_errors);

  % Return cell number with closest alignment and local normal
  closest_fluid_cell = loc_fluid_cells(min_ind);
  computed_normal = cell_normals(min_ind, :);

end

