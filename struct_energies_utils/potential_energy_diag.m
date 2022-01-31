function [pI] = potential_energy_diag(elem_nodes, node_loc, psi_i, mat_props)

  % Numerically computes strain energy arising from
  % potential energy calculation -- \int_V \sigma : \epsilon dV
  % = \int_V \lambda (\epsilon_{11} + \epsilon_{22} + \epsilon_{33})^2 dV ...
  %    + \int_V \mu (\epsilon_{11}^2 + \epsilon_{22}^2 + \epsilon_{33}^2) dV ...
  %    + \int_V \mu (2*\epsilon_{12}^2 + 2*\epsilon_{23}^2 + 2*\epsilon_{13}^2) dV
  % with "diagonal" mode terms

  % Inputs:
  % elem_nodes: nelm x 10 array, giving all the nodes belonging to the ith element (these are ordered as node I -- node R)
  % node_loc: nnodes x 3 array, giving the physical location of each node
  % psi_i: [x,y,z]-directional deformation components
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
  alpha = 0.585410196624968;
  beta = 0.138196601125010;
  weight = 1/24;  % to preserve volume of tetrahedron in local coordinate system

  % Compute volume integral for potential energy
  vol_elms = zeros(nelm,1);
  detJ_ratio = zeros(nelm,1);
  for i = 1:nelm
  
    % Pull off corner nodes
    global_coords = node_loc(elem_nodes(i,1:10), 1:3);

    % Compute jacobians at each integration point
    [detJ1, J1, J1inv, L1] = jacobian(beta, beta, beta, global_coords);
    [detJ2, J2, J2inv, L2] = jacobian(alpha, beta, beta, global_coords);
    [detJ3, J3, J3inv, L3] = jacobian(beta, alpha, beta, global_coords);
    [detJ4, J4, J4inv, L4] = jacobian(beta, beta, alpha, global_coords);

    % TEMP (4/11/21): check Jacobian ratios
    max_detJ_loc = max([detJ1, detJ2, detJ3, detJ4]);
    min_detJ_loc = min([detJ1, detJ2, detJ3, detJ4]);
    detJ_ratio(i) = max_detJ_loc / min_detJ_loc;

    % Pull out component displacement values
    psi_i_nodes = psi_i(elem_nodes(i,1:10), :); % [uI; ...; uR] 10x3
    uxN = psi_i_nodes(:,1)';
    uyN = psi_i_nodes(:,2)';
    uzN = psi_i_nodes(:,3)';
    q = reshape([uxN' uyN' uzN']', [], 1);

    % Acquire strains at integration points by computing derivatives of nodal
    % displacement data. For diagonal terms, the strains in the element are directly
    % defined by the derivatives of the (single mode) displacement data
    [eps_xx1, eps_yy1, eps_zz1, eps_xy1, eps_yz1, eps_xz1] = strain(L1, J1inv, q);
    [eps_xx2, eps_yy2, eps_zz2, eps_xy2, eps_yz2, eps_xz2] = strain(L2, J2inv, q);
    [eps_xx3, eps_yy3, eps_zz3, eps_xy3, eps_yz3, eps_xz3] = strain(L3, J3inv, q);
    [eps_xx4, eps_yy4, eps_zz4, eps_xy4, eps_yz4, eps_xz4] = strain(L4, J4inv, q);

    p1 = (lambda/2)*(eps_xx1 + eps_yy1 + eps_zz1)^2 ...
       + mu*(eps_xx1^2 + eps_yy1^2 + eps_zz1^2 ...
       + 2*eps_xy1^2 + 2*eps_yz1^2 + 2*eps_xz1^2);

    p2 = (lambda/2)*(eps_xx2 + eps_yy2 + eps_zz2)^2 ...
       + mu*(eps_xx2^2 + eps_yy2^2 + eps_zz2^2 ...
       + 2*eps_xy2^2 + 2*eps_yz2^2 + 2*eps_xz2^2);

    p3 = (lambda/2)*(eps_xx3 + eps_yy3 + eps_zz3)^2 ...
       + mu*(eps_xx3^2 + eps_yy3^2 + eps_zz3^2 ...
       + 2*eps_xy3^2 + 2*eps_yz3^2 + 2*eps_xz3^2);

    p4 = (lambda/2)*(eps_xx4 + eps_yy4 + eps_zz4)^2 ...
       + mu*(eps_xx4^2 + eps_yy4^2 + eps_zz4^2 ...
       + 2*eps_xy4^2 + 2*eps_yz4^2 + 2*eps_xz4^2);


    % Integrate
%    vol_elm = weight*(detJ1*p1+detJ2*p2+detJ3*p3+detJ4*p4);   % integral of leading term for U(q_j)
    vol_elm = 2*weight*(detJ1*p1+detJ2*p2+detJ3*p3+detJ4*p4);  % leading term after taking dU/dq_j
    vol_elms(i) = vol_elm;

    % Combine with existing integral tracking
    pI = pI + vol_elm;
  
  end


end

