function [kI] = kinetic_energy_diag(elem_nodes, node_loc, psi_i, mat_props)

  % Numerically computes integrals arising from 
  % "diagonal" component of 
  % kinetic energy calculation -- \int_V \rho_m \psi_i \cdot \psi_i dV

  % Inputs:
  % elem_nodes: nelm x 10 array, giving all the nodes belonging to the element in the ith row
  % node_loc: nnode x 3 array, giving the physical location of each node
  % psi_i: [x,y,z]-directional deformation components
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
  
    % Interpolate displacement values at quadrature points
    psi1 = shape_interpolation([alpha,beta,beta,beta], psi_i_nodes);
    psi2 = shape_interpolation([beta,alpha,beta,beta], psi_i_nodes);
    psi3 = shape_interpolation([beta,beta,alpha,beta], psi_i_nodes);
    psi4 = shape_interpolation([beta,beta,beta,alpha], psi_i_nodes);

    psi1_sq = dot(psi1,psi1);
    psi2_sq = dot(psi2,psi2);
    psi3_sq = dot(psi3,psi3);
    psi4_sq = dot(psi4,psi4);
  
    % Compute Jacobians at each location
    [detJ1, J1, J1inv, L1] = jacobian(beta, beta, beta, global_coords);
    [detJ2, J2, J2inv, L2] = jacobian(alpha, beta, beta, global_coords);
    [detJ3, J3, J3inv, L3] = jacobian(beta, alpha, beta, global_coords);
    [detJ4, J4, J4inv, L4] = jacobian(beta, beta, alpha, global_coords);

    p1 = rhom*psi1_sq;
    p2 = rhom*psi2_sq;
    p3 = rhom*psi3_sq;
    p4 = rhom*psi4_sq;

    % Integrate
    vol_elm = weight*(detJ1*p1 + detJ2*p2 + detJ3*p3 + detJ4*p4);
  
    % Combine with existing integral tracking
    kI = kI + vol_elm;
  
  end
  
end

