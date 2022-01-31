function [normal_deriv] = determine_normal_deriv(struc_elem_locs, def_shape_function, fluid_node_loc, fluid_node_loc_barycentric, bary_coord_num)

  if (nargin==4)
    % Find which barycentric coordinate is zero, if needed
    [~,bary_coord_num] = min(abs(fluid_node_loc_barycentric));
  end
  [dx0_dLi, dx0_dLj, order] = def_derivs(struc_elem_locs, def_shape_function, fluid_node_loc, fluid_node_loc_barycentric, bary_coord_num);
  [dpsi_dLi, dpsi_dLj] = psi_derivs(def_shape_function, fluid_node_loc, fluid_node_loc_barycentric, bary_coord_num, order);

  def_derivs_cross = cross(dx0_dLi, dx0_dLj);
  def_derivs_cross_dq = cross(dpsi_dLi, dx0_dLj) + cross(dx0_dLi, dpsi_dLj);

  dn1_dq = (dpsi_dLi(2)*dx0_dLj(3) + dx0_dLi(2)*dpsi_dLj(3)) - (dpsi_dLi(3)*dx0_dLj(2) + dx0_dLi(3)*dpsi_dLj(2));
  dn2_dq = (dpsi_dLi(3)*dx0_dLj(1) + dx0_dLi(3)*dpsi_dLj(1)) - (dpsi_dLi(1)*dx0_dLj(3) + dx0_dLi(1)*dpsi_dLj(3));
  dn3_dq = (dpsi_dLi(1)*dx0_dLj(2) + dx0_dLi(1)*dpsi_dLj(2)) - (dpsi_dLi(2)*dx0_dLj(1) + dx0_dLi(2)*dpsi_dLj(1));
  def_derivs_cross_norm_dq = (def_derivs_cross(1)*dn1_dq + def_derivs_cross(2)*dn2_dq + def_derivs_cross(3)*dn3_dq) / sqrt(def_derivs_cross(1)^2 + def_derivs_cross(2)^2 + def_derivs_cross(3)^2);

  normal_deriv = (norm(def_derivs_cross).*def_derivs_cross_dq - def_derivs_cross.*def_derivs_cross_norm_dq)./(norm(def_derivs_cross)^2);

end
