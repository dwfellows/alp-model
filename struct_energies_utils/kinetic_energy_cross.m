function [kI] = kinetic_energy_cross(elem_nodes, node_loc, psi_i, psi_j, mat_props)

  % Numerically computes integrals arising from 
  % "cross" component of 
  % kinetic energy calculation -- \int_V \rho_m \psi_i \cdot \psi_j dV

  % Inputs:
  % elem_nodes: nelm x 10 array, giving all the nodes belonging to the element in the ith row
  % node_loc: nnode x 3 array, giving the physical location of each node
  % psi_i: [x,y,z]-directional deformation components from mode i
  % psi_j: [x,y,z]-directional deformation components from mode j
  % mat_props: [E, \nu, \rho_m]

  % Pull off number of elements
  nelm = max(size(elem_nodes));
  
  % Pull out material density
  rhom = mat_props(3);

  % Compute volume integral for kinetic energy
  kI = 0;
  
  % 4-pt quadrature rule: O(h^3) and so preserves third-order quadratic tetrahedron accuracy
  % Values taken from Zienkiewicz Ch.5. Alpha/beta are locations, weight is same for all locs
  alpha = 0.585410196624968;
  beta = 0.138196601125010;
  weight = 1/24; % to preserve volume of tetrahedron in local coordinate system
  
  for i = 1:nelm
  
    % Pull off node locations
    global_coords = node_loc(elem_nodes(i,1:10), 1:3);

    % Pull out nodal displacement values
    psi_i_nodes = psi_i(elem_nodes(i,1:10), :);  % [uI; ...; uR], 10x3 matrix
    psi_j_nodes = psi_j(elem_nodes(i,1:10), :); 
    psi_ij_dot_nodes = dot(psi_i_nodes, psi_j_nodes, 2);
  
    % Interpolate displacement values at quadrature points
    psi1_ij = shape_interpolation([alpha,beta,beta,beta], psi_ij_dot_nodes);
    psi2_ij = shape_interpolation([beta,alpha,beta,beta], psi_ij_dot_nodes);
    psi3_ij = shape_interpolation([beta,beta,alpha,beta], psi_ij_dot_nodes);
    psi4_ij = shape_interpolation([beta,beta,beta,alpha], psi_ij_dot_nodes);
  
    % Compute Jacobians at each location
    [detJ1, J1, J1inv, L1] = jacobian(beta, beta, beta, global_coords);
    [detJ2, J2, J2inv, L2] = jacobian(alpha, beta, beta, global_coords);
    [detJ3, J3, J3inv, L3] = jacobian(beta, alpha, beta, global_coords);
    [detJ4, J4, J4inv, L4] = jacobian(beta, beta, alpha, global_coords);

    p1 = rhom*psi1_ij;
    p2 = rhom*psi2_ij;
    p3 = rhom*psi3_ij;
    p4 = rhom*psi4_ij;

    % Integrate
    vol_elm = weight*(detJ1*p1 + detJ2*p2 + detJ3*p3 + detJ4*p4);
  
    % Combine with existing integral tracking
    kI = kI + vol_elm;
  
  end
  
end

