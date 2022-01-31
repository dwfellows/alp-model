function [dx_dLi, dx_dLj, order] = def_derivs(struc_elem_locs, def_shape_function, fluid_node_loc, fluid_node_loc_barycentric, bary_coord_num)

  % Find which surface we are on and provide normal
  if (bary_coord_num == 1)
    % Compute local derivatives with respect to transformed coordinates
    dx_dL3 = struc_elem_locs(2,:).*(-3 + 4*fluid_node_loc_barycentric(3) + 4*fluid_node_loc_barycentric(4)) ...
           + struc_elem_locs(3,:).*(4*fluid_node_loc_barycentric(3) - 1) ...
           + 4.*struc_elem_locs(6,:).*(1 - 2*fluid_node_loc_barycentric(3) - fluid_node_loc_barycentric(4)) ...
           + 4.*struc_elem_locs(9,:).*(-fluid_node_loc_barycentric(4)) + 4.*struc_elem_locs(10,:).*fluid_node_loc_barycentric(4);
    dx_dL4 = struc_elem_locs(2,:).*(-3 + 4*fluid_node_loc_barycentric(3) + 4*fluid_node_loc_barycentric(4)) ...
           + struc_elem_locs(4,:).*(4*fluid_node_loc_barycentric(4) - 1) ...
           + 4.*struc_elem_locs(6,:).*(-fluid_node_loc_barycentric(3)) ...
           + 4.*struc_elem_locs(9,:).*(1 - fluid_node_loc_barycentric(3) - 2*fluid_node_loc_barycentric(4)) ...
           + 4.*struc_elem_locs(10,:).*fluid_node_loc_barycentric(3);

    normal = cross(dx_dL3, dx_dL4);
    normal = normal ./ norm(normal);

    % Move in normal direction and see if it is inward or outward pointing
    vec_test = struc_elem_locs(1,:) - fluid_node_loc;
    dp = dot(normal, vec_test);
    if (dp > 0)
      dx_dLi = dx_dL4;
      dx_dLj = dx_dL3;
      order = [4,3];
    else
      dx_dLi = dx_dL3;
      dx_dLj = dx_dL4;
      order = [3,4];
    end

  elseif (bary_coord_num == 2)
    dx_dL3 = struc_elem_locs(1,:).*(-3 + 4*fluid_node_loc_barycentric(3) + 4*fluid_node_loc_barycentric(4)) ...
           + struc_elem_locs(3,:).*(4*fluid_node_loc_barycentric(3) - 1) ...
           + 4.*struc_elem_locs(7,:).*(1-2*fluid_node_loc_barycentric(3) - fluid_node_loc_barycentric(4)) + 4.*struc_elem_locs(8,:).*(-fluid_node_loc_barycentric(4)) ...
           + 4.*struc_elem_locs(10,:).*fluid_node_loc_barycentric(4);
    dx_dL4 = struc_elem_locs(1,:).*(-3 + 4*fluid_node_loc_barycentric(3) + 4*fluid_node_loc_barycentric(4)) ...
           + struc_elem_locs(4,:).*(4*fluid_node_loc_barycentric(4) - 1) ...
           + 4.*struc_elem_locs(7,:).*(-fluid_node_loc_barycentric(3)) ...
           + 4.*struc_elem_locs(8,:).*(1 - fluid_node_loc_barycentric(3) - 2*fluid_node_loc_barycentric(4)) ...
           + 4.*struc_elem_locs(10,:).*fluid_node_loc_barycentric(3);

    normal = cross(dx_dL3, dx_dL4);
    normal = normal ./ norm(normal);

    % Move in normal direction and see if it is inward or outward pointing
    vec_test = struc_elem_locs(2,:) - fluid_node_loc;
    dp = dot(normal, vec_test);
    if (dp > 0)
      dx_dLi = dx_dL4;
      dx_dLj = dx_dL3;
      order = [4,3];
    else
      dx_dLi = dx_dL3;
      dx_dLj = dx_dL4;
      order = [3,4];
    end

  elseif (bary_coord_num == 3)
    dx_dL2 = struc_elem_locs(1,:).*(-3 + 4*fluid_node_loc_barycentric(2) + 4*fluid_node_loc_barycentric(4)) ...
           + struc_elem_locs(2,:).*(4*fluid_node_loc_barycentric(2) - 1) ...
           + 4.*struc_elem_locs(5,:).*(1 - 2*fluid_node_loc_barycentric(2) - fluid_node_loc_barycentric(4)) ...
           + 4.*struc_elem_locs(8,:).*(-fluid_node_loc_barycentric(4)) + 4.*struc_elem_locs(9,:).*fluid_node_loc_barycentric(4);
    dx_dL4 = struc_elem_locs(1,:).*(-3 + 4*fluid_node_loc_barycentric(2) + 4*fluid_node_loc_barycentric(4)) ...
           + struc_elem_locs(4,:).*(4*fluid_node_loc_barycentric(4) - 1) ...
           + 4.*struc_elem_locs(5,:).*(-fluid_node_loc_barycentric(2)) ...
           + 4.*struc_elem_locs(8,:).*(1 - fluid_node_loc_barycentric(2) - 2*fluid_node_loc_barycentric(4)) ...
           + 4.*struc_elem_locs(9,:).*fluid_node_loc_barycentric(2);

    normal = cross(dx_dL2, dx_dL4);
    normal = normal ./ norm(normal);

    % Move in normal direction and see if it is inward or outward pointing
    vec_test = struc_elem_locs(3,:) - fluid_node_loc;
    dp = dot(normal, vec_test);
    if (dp > 0)
      dx_dLi = dx_dL4;
      dx_dLj = dx_dL2;
      order = [4,2];
    else
      dx_dLi = dx_dL2;
      dx_dLj = dx_dL4;
      order = [2,4];
    end

  else
    dx_dL2 = struc_elem_locs(1,:).*(-3 + 4*fluid_node_loc_barycentric(2) + 4*fluid_node_loc_barycentric(3)) ...
           + struc_elem_locs(2,:).*(4*fluid_node_loc_barycentric(2) - 1) ...
           + 4.*struc_elem_locs(5,:).*(1 - 2*fluid_node_loc_barycentric(2) - fluid_node_loc_barycentric(3)) ...
           + 4.*struc_elem_locs(6,:).*fluid_node_loc_barycentric(3) + 4.*struc_elem_locs(7,:).*(-fluid_node_loc_barycentric(3));
    dx_dL3 = struc_elem_locs(1,:).*(-3 + 4*fluid_node_loc_barycentric(2) + 4*fluid_node_loc_barycentric(3)) ...
           + struc_elem_locs(3,:).*(4*fluid_node_loc_barycentric(3) - 1) + 4.*struc_elem_locs(5,:).*(-fluid_node_loc_barycentric(2)) ...
           + 4.*struc_elem_locs(6,:).*fluid_node_loc_barycentric(2) + 4.*struc_elem_locs(7,:).*(1 - fluid_node_loc_barycentric(2) - 2*fluid_node_loc_barycentric(3));

    normal = cross(dx_dL2, dx_dL3);
    normal = normal ./ norm(normal);

    % Move in normal direction and see if it is inward or outward pointing
    vec_test = struc_elem_locs(4,:) - fluid_node_loc;
    dp = dot(normal, vec_test);
    if (dp > 0)
      dx_dLi = dx_dL3;
      dx_dLj = dx_dL2;
      order = [3,2];
    else
      dx_dLi = dx_dL2;
      dx_dLj = dx_dL3;
      order = [2,3];
    end

  end

end
