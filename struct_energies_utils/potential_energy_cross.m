function [pI] = potential_energy_cross(elem_nodes, node_loc, psi_i, psi_j, mat_props)

  % Numerically computes strain energy arising from
  % potential energy calculation -- \int_V \sigma : \epsilon dV
  % = \int_V \lambda (\epsilon_{11} + \epsilon_{22} + \epsilon_{33})^2 dV ...
  %    + \int_V \mu (\epsilon_{11}^2 + \epsilon_{22}^2 + \epsilon_{33}^2) dV ...
  %    + \int_V \mu (2*\epsilon_{12}^2 + 2*\epsilon_{23}^2 + 2*\epsilon_{13}^2) dV
  % with "cross" mode terms

  % Inputs:
  % elem_nodes: nelm x 10 array, giving all the nodes belonging to the ith element (these are ordered as node I -- node R)
  % node_loc: nnodes x 3 array, giving the physical location of each node
  % psi_i: [x,y,z]-directional deformation components of mode i
  % psi_j: [x,y,z]-directional deformation components of mode j
  % mat_props: [Young's modulus, Poisson's ratio]

  % Pull off number of elements
  nelm = max(size(elem_nodes));

  % Store potential energy
  pI = 0;
  
  % Define Young's modulus and Poisson's ratio
  E = mat_props(1);   % Pascals
  nu = mat_props(2);  % Dimensionless
  
  % Define Lame's constant and shear modulus
  lambda = (E*nu)/((1+nu)*(1-2*nu));
  mu = E/(2*(1+nu));

  % 4-pt quadrature rule: O(h^3) and so preserves third-order quadratic tetrahedron accuracy
  % Values taken from Zienkiewicz Ch.5. Alpha/beta are locations, weight is same for all locs
  alpha = 0.58541020;
  beta = 0.13819660;
  weight = 1/24;  % to preserve volume of tetrahedron in local coordinate system

  % Compute volume integral for potential energy
  for i = 1:nelm
  
    % Pull off corner nodes
    global_coords = node_loc(elem_nodes(i,1:10), 1:3);

    % Compute jacobians at each integration point
    [detJ1, J1, J1inv, L1] = jacobian(beta, beta, beta, global_coords);
    [detJ2, J2, J2inv, L2] = jacobian(alpha, beta, beta, global_coords);
    [detJ3, J3, J3inv, L3] = jacobian(beta, alpha, beta, global_coords);
    [detJ4, J4, J4inv, L4] = jacobian(beta, beta, alpha, global_coords);

    % Pull out component displacement values
    psi_i_nodes = psi_i(elem_nodes(i,1:10), :); % [uiI; ...; uiR] 10x3
    uxNi = psi_i_nodes(:,1)';
    uyNi = psi_i_nodes(:,2)';
    uzNi = psi_i_nodes(:,3)';
    qi = reshape([uxNi' uyNi' uzNi']', [], 1);

    psi_j_nodes = psi_j(elem_nodes(i,1:10), :); % [ujI; ...; ujR] 10x3
    uxNj = psi_j_nodes(:,1)';
    uyNj = psi_j_nodes(:,2)';
    uzNj = psi_j_nodes(:,3)';
    qj = reshape([uxNj' uyNj' uzNj']', [], 1);

    % Acquire strains at integration points computed separately
    [eps_i_xx1, eps_i_yy1, eps_i_zz1, eps_i_xy1, eps_i_yz1, eps_i_xz1] = strain(L1, J1inv, qi);
    [eps_i_xx2, eps_i_yy2, eps_i_zz2, eps_i_xy2, eps_i_yz2, eps_i_xz2] = strain(L2, J2inv, qi);
    [eps_i_xx3, eps_i_yy3, eps_i_zz3, eps_i_xy3, eps_i_yz3, eps_i_xz3] = strain(L3, J3inv, qi);
    [eps_i_xx4, eps_i_yy4, eps_i_zz4, eps_i_xy4, eps_i_yz4, eps_i_xz4] = strain(L4, J4inv, qi);

    [eps_j_xx1, eps_j_yy1, eps_j_zz1, eps_j_xy1, eps_j_yz1, eps_j_xz1] = strain(L1, J1inv, qj);
    [eps_j_xx2, eps_j_yy2, eps_j_zz2, eps_j_xy2, eps_j_yz2, eps_j_xz2] = strain(L2, J2inv, qj);
    [eps_j_xx3, eps_j_yy3, eps_j_zz3, eps_j_xy3, eps_j_yz3, eps_j_xz3] = strain(L3, J3inv, qj);
    [eps_j_xx4, eps_j_yy4, eps_j_zz4, eps_j_xy4, eps_j_yz4, eps_j_xz4] = strain(L4, J4inv, qj);

    p1_normal = 2*mu*(eps_i_xx1*eps_j_xx1 + eps_i_yy1*eps_j_yy1 + eps_i_zz1*eps_j_zz1);
    p2_normal = 2*mu*(eps_i_xx2*eps_j_xx2 + eps_i_yy2*eps_j_yy2 + eps_i_zz2*eps_j_zz2);
    p3_normal = 2*mu*(eps_i_xx3*eps_j_xx3 + eps_i_yy3*eps_j_yy3 + eps_i_zz3*eps_j_zz3);
    p4_normal = 2*mu*(eps_i_xx4*eps_j_xx4 + eps_i_yy4*eps_j_yy4 + eps_i_zz4*eps_j_zz4);

    p1_shear = 8*mu*(eps_i_xy1*eps_j_xy1 + eps_i_yz1*eps_j_yz1 + eps_i_xz1*eps_j_xz1);
    p2_shear = 8*mu*(eps_i_xy2*eps_j_xy2 + eps_i_yz2*eps_j_yz2 + eps_i_xz2*eps_j_xz2);
    p3_shear = 8*mu*(eps_i_xy3*eps_j_xy3 + eps_i_yz3*eps_j_yz3 + eps_i_xz3*eps_j_xz3);
    p4_shear = 8*mu*(eps_i_xy4*eps_j_xy4 + eps_i_yz4*eps_j_yz4 + eps_i_xz4*eps_j_xz4);

    p1_dil = 2*lambda*((eps_i_xx1 + eps_i_yy1 + eps_i_zz1)*(eps_j_xx1 + eps_j_yy1 + eps_j_zz1));
    p2_dil = 2*lambda*((eps_i_xx2 + eps_i_yy2 + eps_i_zz2)*(eps_j_xx2 + eps_j_yy2 + eps_j_zz2));
    p3_dil = 2*lambda*((eps_i_xx3 + eps_i_yy3 + eps_i_zz3)*(eps_j_xx3 + eps_j_yy3 + eps_j_zz3));
    p4_dil = 2*lambda*((eps_i_xx4 + eps_i_yy4 + eps_i_zz4)*(eps_j_xx4 + eps_j_yy4 + eps_j_zz4));


    p1 = p1_normal + p1_shear + p1_dil;
    p2 = p2_normal + p2_shear + p2_dil;
    p3 = p3_normal + p3_shear + p3_dil;
    p4 = p4_normal + p4_shear + p4_dil;


    % Integrate
    vol_elm = weight*(detJ1*p1+detJ2*p2+detJ3*p3+detJ4*p4);

    % Combine with existing integral tracking
    pI = pI + vol_elm;
  
  end

end

