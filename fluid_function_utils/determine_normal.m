function [normal] = determine_normal(struc_elem_locs, fluid_node_loc, fluid_node_loc_barycentric, bary_coord_num)

  % Find which barycentric coordinate is zero, if not specified
  if (nargin == 3)
    [~,bary_coord_num] = min(abs(fluid_node_loc_barycentric));
  end
  normal = bary_deriv_normal(struc_elem_locs, fluid_node_loc, fluid_node_loc_barycentric, bary_coord_num);  % [x,y,z]

end
