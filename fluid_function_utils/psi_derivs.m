function [dpsi_dLi, dpsi_dLj] = psi_derivs(dsf, fluid_node_loc, fluid_node_loc_barycentric, bary_coord_num, order)

  % dsf == "Deformation shape function" i.e. psi_1, component-wise

  % Find which surface we are on and acquire correct derivatives
  if (bary_coord_num == 1)
    % Compute local derivatives with respect to transformed coordinates
    dpsi_dL3 = dsf(2,:).*(-3 + 4*fluid_node_loc_barycentric(3) + 4*fluid_node_loc_barycentric(4)) ...
               + dsf(3,:).*(4*fluid_node_loc_barycentric(3) - 1) ...
               + 4*dsf(6,:).*(1 - 2*fluid_node_loc_barycentric(3) - fluid_node_loc_barycentric(4)) ...
               + 4*dsf(9,:).*(-fluid_node_loc_barycentric(4)) + 4*dsf(10,:).*fluid_node_loc_barycentric(4);

    dpsi_dL4 = dsf(2,:).*(-3 + 4*fluid_node_loc_barycentric(3) + 4*fluid_node_loc_barycentric(4)) ...
               + dsf(4,:).*(4*fluid_node_loc_barycentric(4) - 1) ...
               + 4*dsf(6,:).*(-fluid_node_loc_barycentric(3)) ...
               + 4*dsf(9,:).*(1 - fluid_node_loc_barycentric(3) - 2*fluid_node_loc_barycentric(4)) ...
               + 4*dsf(10,:).*fluid_node_loc_barycentric(3);

    % Match precedent from surface normal computations
    if (order(1) == 3)
      dpsi_dLi = dpsi_dL3;
      dpsi_dLj = dpsi_dL4;
    else
      dpsi_dLi = dpsi_dL4;
      dpsi_dLj = dpsi_dL3;
    end

  elseif (bary_coord_num == 2);
    % Compute local derivativeds with respect to transformed coordinates
   dpsi_dL3 = dsf(1,:).*(-3 + 4*fluid_node_loc_barycentric(3) + 4*fluid_node_loc_barycentric(4)) ...
           + dsf(3,:).*(4*fluid_node_loc_barycentric(3) - 1) ...
           + 4.*dsf(7,:).*(1-2*fluid_node_loc_barycentric(3) - fluid_node_loc_barycentric(4)) ...
           + 4.*dsf(8,:).*(-fluid_node_loc_barycentric(4)) ...
           + 4.*dsf(10,:).*fluid_node_loc_barycentric(4);

    dpsi_dL4 = dsf(1,:).*(-3 + 4*fluid_node_loc_barycentric(3) + 4*fluid_node_loc_barycentric(4)) ...
           + dsf(4,:).*(4*fluid_node_loc_barycentric(4) - 1) ...
           + 4.*dsf(7,:).*(-fluid_node_loc_barycentric(3)) ...
           + 4.*dsf(8,:).*(1 - fluid_node_loc_barycentric(3) - 2*fluid_node_loc_barycentric(4)) ...
           + 4.*dsf(10,:).*fluid_node_loc_barycentric(3);

    % Match precedent from surface normal computations
    if (order(1) == 3)
      dpsi_dLi = dpsi_dL3;
      dpsi_dLj = dpsi_dL4;
    else
      dpsi_dLi = dpsi_dL4;
      dpsi_dLj = dpsi_dL3;
    end

  elseif (bary_coord_num == 3)
    % Compute local derivatives with respect to transformed coordinates
    dpsi_dL2 = dsf(1,:).*(-3 + 4*fluid_node_loc_barycentric(2) + 4*fluid_node_loc_barycentric(4)) ...
           + dsf(2,:).*(4*fluid_node_loc_barycentric(2) - 1) ...
           + 4.*dsf(5,:).*(1 - 2*fluid_node_loc_barycentric(2) - fluid_node_loc_barycentric(4)) ...
           + 4.*dsf(8,:).*(-fluid_node_loc_barycentric(4)) + 4.*dsf(9,:).*fluid_node_loc_barycentric(4);

    dpsi_dL4 = dsf(1,:).*(-3 + 4*fluid_node_loc_barycentric(2) + 4*fluid_node_loc_barycentric(4)) ...
           + dsf(4,:).*(4*fluid_node_loc_barycentric(4) - 1) ...
           + 4.*dsf(5,:).*(-fluid_node_loc_barycentric(2)) ...
           + 4.*dsf(8,:).*(1 - fluid_node_loc_barycentric(2) - 2*fluid_node_loc_barycentric(4)) ...
           + 4.*dsf(9,:).*fluid_node_loc_barycentric(2);

    % Match precedent from surface normal computations
    if (order(1) == 2)
      dpsi_dLi = dpsi_dL2;
      dpsi_dLj = dpsi_dL4;
    else
      dpsi_dLi = dpsi_dL4;
      dpsi_dLj = dpsi_dL2;
    end

  else
    % Compute local derivatives with respect to transformed coordinates
    dpsi_dL2 = dsf(1,:).*(-3 + 4*fluid_node_loc_barycentric(2) + 4*fluid_node_loc_barycentric(3)) ...
           + dsf(2,:).*(4*fluid_node_loc_barycentric(2) - 1) ...
           + 4.*dsf(5,:).*(1 - 2*fluid_node_loc_barycentric(2) - fluid_node_loc_barycentric(3)) ...
           + 4.*dsf(6,:).*fluid_node_loc_barycentric(3) + 4.*dsf(7,:).*(-fluid_node_loc_barycentric(3));

    dpsi_dL3 = dsf(1,:).*(-3 + 4*fluid_node_loc_barycentric(2) + 4*fluid_node_loc_barycentric(3)) ...
           + dsf(3,:).*(4*fluid_node_loc_barycentric(3) - 1) + 4.*dsf(5,:).*(-fluid_node_loc_barycentric(2)) ...
           + 4.*dsf(6,:).*fluid_node_loc_barycentric(2) ...
           + 4.*dsf(7,:).*(1 - fluid_node_loc_barycentric(2) - 2*fluid_node_loc_barycentric(3));

    % Match precedent from surface normal computations
    if (order(1) == 2)
      dpsi_dLi = dpsi_dL2;
      dpsi_dLj = dpsi_dL3;
    else
      dpsi_dLi = dpsi_dL3;
      dpsi_dLj = dpsi_dL2;
    end

  end

end
