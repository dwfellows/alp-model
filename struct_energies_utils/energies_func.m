%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55%%
% An initial script to compute kinetic and potential energies for method of assumed modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Main routines

function [kI] = kinetic_energy(elem_nodes, ux, uy, uz)

  % Numerically computes integrals arising from 
  % kinetic energy calculation -- \int_V \psi^2 dV

  % Pull off number of elements
  nelm = max(size(elem_nodes));
  
  % Construct total displacement components
  u = sqrt(ux.^2 + uy.^2 + uz.^2);

  % Normalize total displacement shape function
  u = u ./ max(u);

  % Compute volume integral for kinetic energy
  kI = 0;
  
  % 4-pt quadrature rule: O(h^3) and so preserves third-order quadratic tetrahedron accuracy
  % Values taken from Zienkiewicz Ch.5. Alpha/beta are locations, weight is same for all locs
  alpha = 0.58541020;
  beta = 0.13819660;
  weight = 0.25;
  
  for i = 1:nelm
  
    % Pull off corner nodes
    x1 = node_loc(elem_nodes(i,1), 1:3);
    x2 = node_loc(elem_nodes(i,2), 1:3);
    x3 = node_loc(elem_nodes(i,3), 1:3);
    x4 = node_loc(elem_nodes(i,4), 1:3);
  
    % Compute jacobian -- follow implementation by Orio
    J = [x2(1) - x1(1), x3(1) - x1(1), x4(1) - x1(1);
         x2(2) - x1(2), x3(2) - x1(2), x4(2) - x1(2);
         x2(3) - x1(3), x3(3) - x1(3), x4(3) - x1(3)];
  
    detJ = det(J);
  
    % Pull out nodal displacement values
    uI = u(elem_nodes(i,1));
    uJ = u(elem_nodes(i,2));
    uK = u(elem_nodes(i,3));
    uL = u(elem_nodes(i,4));
    uM = u(elem_nodes(i,5));
    uN = u(elem_nodes(i,6));
    uO = u(elem_nodes(i,7));
    uP = u(elem_nodes(i,8));
    uQ = u(elem_nodes(i,9));
    uR = u(elem_nodes(i,10));
    nodal_disp = [uI, uJ, uK, uL, uM, uN, uO, uP, uQ, uR];
  
    % Interpolate displacement values at quadrature points
    psi1 = shape_interpolation([alpha,beta,beta,beta], nodal_disp);
    psi2 = shape_interpolation([beta,alpha,beta,beta], nodal_disp);
    psi3 = shape_interpolation([beta,beta,alpha,beta], nodal_disp);
    psi4 = shape_interpolation([beta,beta,beta,alpha], nodal_disp);
  
    % Integrate
    vol_elm = weight*detJ*(psi1^2 + psi2^2 + psi3^2 + psi4^2);
  
    % Combine with existing integral tracking
    kI = kI + vol_elm;
  
  end
  
end

function [pI] = potential_energy(elem_nodes, ux, uy, uz, mat_props)

  % Numerically computes strain energy arising from
  % potential energy calculation -- \int_V \sigma : \epsilon dV

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
  
  % Compute volume integral for potential energy
  for i = 1:nelm
  
    % Pull off corner nodes
    x1  = node_loc(elem_nodes(i,1), 1:3);
    x2  = node_loc(elem_nodes(i,2), 1:3);
    x3  = node_loc(elem_nodes(i,3), 1:3);
    x4  = node_loc(elem_nodes(i,4), 1:3);
  
    % Compute jacobian -- follow implementation by Orio
    J = [x2(1) - x1(1), x3(1) - x1(1), x4(1) - x1(1);
         x2(2) - x1(2), x3(2) - x1(2), x4(2) - x1(2);
         x2(3) - x1(3), x3(3) - x1(3), x4(3) - x1(3)];
  
    detJ = det(J);
  
    % Compute volume coordinate components
    % Volume
    V = [1, x1(1), x1(2), x1(3);
         1, x2(1), x2(2), x2(3);
         1, x3(1), x3(2), x3(3);
         1, x4(1), x4(2), x4(3)];
    V = (1/6)*det(V);
  
    % L1
    a1 = [x2(1), x2(2), x2(3);
          x3(1), x3(2), x3(3);
          x4(1), x4(2), x4(3)];
    a1 = det(a1);
  
    b1 = [1, x2(2), x2(3);
          1, x3(2), x3(3);
          1, x4(2), x4(3)];
    b1 = -1*det(b1);
  
    c1 = [x2(1), 1, x2(3);
          x3(1), 1, x3(3);
          x4(1), 1, x4(3)];
    c1 = -1*det(c1);
  
    d1 = [x2(1), x2(2), 1;
          x3(1), x3(2), 1;
          x4(1), x4(2), 1];
    d1 = -1*det(d1);
  
    % L2
    a2 = [x3(1), x3(2), x3(3);
          x4(1), x4(2), x4(3);
          x1(1), x1(2), x1(3)];
    a2 = det(a2);
  
    b2 = [1, x3(2), x3(3);
          1, x4(2), x4(3);
          1, x1(2), x1(3)];
    b2 = -1*det(b2);
  
    c2 = [x3(1), 1, x3(3);
          x4(1), 1, x4(3);
          x1(1), 1, x1(3)];
    c2 = -1*det(c2);
  
    d2 = [x3(1), x3(2), 1;
          x4(1), x4(2), 1;
          x1(1), x1(2), 1];
    d2 = -1*det(d2);
  
    % L3
    a3 = [x4(1), x4(2), x4(3);
          x1(1), x1(2), x1(3);
          x2(1), x2(2), x2(3)];
    a3 = det(a3);
  
    b3 = [1, x4(2), x4(3);
          1, x1(2), x1(3);
          1, x2(2), x2(3)];
    b3 = -1*det(b3);
  
    c3 = [x4(1), 1, x4(3);
          x1(1), 1, x1(3);
          x2(1), 1, x2(3)];
    c3 = -1*det(c3);
  
    d3 = [x4(1), x4(2), 1;
          x1(1), x1(2), 1;
          x2(1), x2(2), 1];
    d3 = -1*det(d3);
  
    % L4
    a4 = [x1(1), x1(2), x1(3);
          x2(1), x2(2), x2(3);
          x3(1), x3(2), x3(3)];
    a4 = det(a4);
  
    b4 = [1, x1(2), x1(3);
          1, x2(2), x2(3);
          1, x3(2), x3(3)];
    b4 = -1*det(b4);
  
    c4 = [x1(1), 1, x1(3);
          x2(1), 1, x2(3);
          x3(1), 1, x3(3)];
    c4 = -1*det(c4);
  
    d4 = [x1(1), x1(2), 1;
          x2(1), x2(2), 1;
          x3(1), x3(2), 1];
    d4 = -1*det(d4);
  
    % Assemble into array
    coeffs = [V, a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3, a4, b4, c4, d4];
  
    % Pull out component displacement values
    uI = ux(elem_nodes(i,1));
    uJ = ux(elem_nodes(i,2));
    uK = ux(elem_nodes(i,3));
    uL = ux(elem_nodes(i,4));
    uM = ux(elem_nodes(i,5));
    uN = ux(elem_nodes(i,6));
    uO = ux(elem_nodes(i,7));
    uP = ux(elem_nodes(i,8));
    uQ = ux(elem_nodes(i,9));
    uR = ux(elem_nodes(i,10));
  
    uxN = [uI, uJ, uK, uL, uM, uN, uO, uP, uQ, uR];
  
    vI = uy(elem_nodes(i,1));
    vJ = uy(elem_nodes(i,2));
    vK = uy(elem_nodes(i,3));
    vL = uy(elem_nodes(i,4));
    vM = uy(elem_nodes(i,5));
    vN = uy(elem_nodes(i,6));
    vO = uy(elem_nodes(i,7));
    vP = uy(elem_nodes(i,8));
    vQ = uy(elem_nodes(i,9));
    vR = uy(elem_nodes(i,10));
  
    uyN = [vI, vJ, vK, vL, vM, vN, vO, vP, vQ, vR];
  
    wI = uz(elem_nodes(i,1));
    wJ = uz(elem_nodes(i,2));
    wK = uz(elem_nodes(i,3));
    wL = uz(elem_nodes(i,4));
    wM = uz(elem_nodes(i,5));
    wN = uz(elem_nodes(i,6));
    wO = uz(elem_nodes(i,7));
    wP = uz(elem_nodes(i,8));
    wQ = uz(elem_nodes(i,9));
    wR = uz(elem_nodes(i,10));
  
    uzN = [wI, wJ, wK, wL, wM, wN, wO, wP, wQ, wR];
  
    % Acquire derivatives of displacements at integration points
    % du/dx
    ux_dx1 = def_deriv(uxN, 1, [alpha, beta, beta, beta], coeffs);
    ux_dx2 = def_deriv(uxN, 1, [beta, alpha, beta, beta], coeffs);
    ux_dx3 = def_deriv(uxN, 1, [beta, beta, alpha, beta], coeffs);
    ux_dx4 = def_deriv(uxN, 1, [beta, beta, beta, alpha], coeffs);
    
    % du/dy
    ux_dy1 = def_deriv(uxN, 2, [alpha, beta, beta, beta], coeffs);
    ux_dy2 = def_deriv(uxN, 2, [beta, alpha, beta, beta], coeffs);
    ux_dy3 = def_deriv(uxN, 2, [beta, beta, alpha, beta], coeffs);
    ux_dy4 = def_deriv(uxN, 2, [beta, beta, beta, alpha], coeffs);
   
    % du/dz
    ux_dz1 = def_deriv(uxN, 3, [alpha, beta, beta, beta], coeffs);
    ux_dz2 = def_deriv(uxN, 3, [beta, alpha, beta, beta], coeffs);
    ux_dz3 = def_deriv(uxN, 3, [beta, beta, alpha, beta], coeffs);
    ux_dz4 = def_deriv(uxN, 3, [beta, beta, beta, alpha], coeffs);
  
    % dv/dx
    uy_dx1 = def_deriv(uyN, 1, [alpha, beta, beta, beta], coeffs);
    uy_dx2 = def_deriv(uyN, 1, [beta, alpha, beta, beta], coeffs);
    uy_dx3 = def_deriv(uyN, 1, [beta, beta, alpha, beta], coeffs);
    uy_dx4 = def_deriv(uyN, 1, [beta, beta, beta, alpha], coeffs);
  
    % dv/dy
    uy_dy1 = def_deriv(uyN, 2, [alpha, beta, beta, beta], coeffs);
    uy_dy2 = def_deriv(uyN, 2, [beta, alpha, beta, beta], coeffs);
    uy_dy3 = def_deriv(uyN, 2, [beta, beta, alpha, beta], coeffs);
    uy_dy4 = def_deriv(uyN, 2, [beta, beta, beta, alpha], coeffs);
  
    % dv/dz
    uy_dz1 = def_deriv(uyN, 3, [alpha, beta, beta, beta], coeffs);
    uy_dz2 = def_deriv(uyN, 3, [beta, alpha, beta, beta], coeffs);
    uy_dz3 = def_deriv(uyN, 3, [beta, beta, alpha, beta], coeffs);
    uy_dz4 = def_deriv(uyN, 3, [beta, beta, beta, alpha], coeffs);
  
    % dw/dx
    uz_dx1 = def_deriv(uzN, 1, [alpha, beta, beta, beta], coeffs);
    uz_dx2 = def_deriv(uzN, 1, [beta, alpha, beta, beta], coeffs);
    uz_dx3 = def_deriv(uzN, 1, [beta, beta, alpha, beta], coeffs);
    uz_dx4 = def_deriv(uzN, 1, [beta, beta, beta, alpha], coeffs);
  
    % dw/dy
    uz_dy1 = def_deriv(uzN, 2, [alpha, beta, beta, beta], coeffs);
    uz_dy2 = def_deriv(uzN, 2, [beta, alpha, beta, beta], coeffs);
    uz_dy3 = def_deriv(uzN, 2, [beta, beta, alpha, beta], coeffs);
    uz_dy4 = def_deriv(uzN, 2, [beta, beta, beta, alpha], coeffs);
  
    % dw/dz
    uz_dz1 = def_deriv(uzN, 3, [alpha, beta, beta, beta], coeffs);
    uz_dz2 = def_deriv(uzN, 3, [beta, alpha, beta, beta], coeffs);
    uz_dz3 = def_deriv(uzN, 3, [beta, beta, alpha, beta], coeffs);
    uz_dz4 = def_deriv(uzN, 3, [beta, beta, beta, alpha], coeffs);
  
    % Compute integrands at integration points
    p1 = (lambda/2)*(ux_dx1 + uy_dy1 + uz_dz1)^2 ...
       + mu*(ux_dx1^2 + uy_dy1^2 + uz_dz1^2 ...
       + 0.5*(ux_dy1^2 + 2*ux_dy1*uy_dx1 + uy_dx1^2) ...
       + 0.5*(uy_dz1^2 + 2*uy_dz1*uz_dy1 + uz_dy1^2) ...
       + 0.5*(ux_dz1^2 + 2*ux_dz1*uz_dx1 + uz_dx1^2));
  
    p2 = (lambda/2)*(ux_dx2 + uy_dy2 + uz_dz2)^2 ...
       + mu*(ux_dx2^2 + uy_dy2^2 + uz_dz2^2 ...
       + 0.5*(ux_dy2^2 + 2*ux_dy2*uy_dx2 + uy_dx2^2) ...
       + 0.5*(uy_dz2^2 + 2*uy_dz2*uz_dy2 + uz_dy2^2) ...
       + 0.5*(ux_dz2^2 + 2*ux_dz2*uz_dx2 + uz_dx2^2));
  
    p3 = (lambda/2)*(ux_dx3 + uy_dy3 + uz_dz3)^2 ...
       + mu*(ux_dx3^2 + uy_dy3^2 + uz_dz3^2 ...
       + 0.5*(ux_dy3^2 + 2*ux_dy3*uy_dx3 + uy_dx3^2) ...
       + 0.5*(uy_dz3^2 + 2*uy_dz3*uz_dy3 + uz_dy3^2) ...
       + 0.5*(ux_dz3^2 + 2*ux_dz3*uz_dx3 + uz_dx3^2));
  
    p4 = (lambda/2)*(ux_dx4 + uy_dy4 + uz_dz4)^2 ...
       + mu*(ux_dx4^2 + uy_dy4^2 + uz_dz4^2 ...
       + 0.5*(ux_dy4^2 + 2*ux_dy4*uy_dx4 + uy_dx4^2) ...
       + 0.5*(uy_dz4^2 + 2*uy_dz4*uz_dy4 + uz_dy4^2) ...
       + 0.5*(ux_dz4^2 + 2*ux_dz4*uz_dx4 + uz_dx4^2));
  
    % Integrate
    vol_elm = weight*detJ*(p1+p2+p3+p4);
  
    % Combine with existing integral tracking
    pI = pI + vol_elm;
  
  end

end


%%%%%%% Subroutines utilized in above routines

function [physical_loc] = tet2phys(tet_loc, x1, x2, x3, x4)

  L1 = tet_loc(1);
  L2 = tet_loc(2);
  L3 = tet_loc(3);
  L4 = tet_loc(4);

  x = L1*x1(1) + L2*x2(1) + L3*x3(1) + L4*x4(1);
  y = L1*x1(2) + L2*x2(2) + L3*x3(2) + L4*x4(2);
  z = L1*x1(3) + L2*x2(3) + L3*x3(3) + L4*x4(3);

  physical_loc = [x,y,z];

end

function [interpolation_val] = shape_interpolation(tet_loc, nodal_displacement)

  % Pull off tet_locs
  L1 = tet_loc(1);
  L2 = tet_loc(2);
  L3 = tet_loc(3);
  L4 = tet_loc(4);

  % Pull off nodal displacements
  uI = nodal_displacement(1);
  uJ = nodal_displacement(2);
  uK = nodal_displacement(3);
  uL = nodal_displacement(4);
  uM = nodal_displacement(5);
  uN = nodal_displacement(6);
  uO = nodal_displacement(7);
  uP = nodal_displacement(8);
  uQ = nodal_displacement(9);
  uR = nodal_displacement(10);

  % Construct interpolation using shape function from ANSYS Mechanical Theory Ref (and Zienkiewics Ch.4)
  interpolation_val = uI*(2*L1-1)*L1 + uJ*(2*L2-1)*L2 + uK*(2*L3-1)*L3 + uL*(2*L4-1)*L4 + ...
                      4*uM*L1*L2 + 4*uN*L2*L3 + 4*uO*L1*L3 + 4*uP*L1*L4 + 4*uQ*L2*L4 + ...
                      4*uR*L3*L4;

end


function [dLidxj] = vcd(i,j,coeffs)
  % vcd -- Volume coordinate derivative
  % L_i = (1/6V)(a_i + b_i x + c_i y + d_i z)
  % dL_i / dx_j = (1/6V) e_j; e_1 = b_i, e_2 = c_i, e_3 = d_i

  V = coeffs(1);
  a1 = coeffs(2);
  b1 = coeffs(3);
  c1 = coeffs(4);
  d1 = coeffs(5);
  a2 = coeffs(6);
  b2 = coeffs(7);
  c2 = coeffs(8);
  d2 = coeffs(9);
  a3 = coeffs(10);
  b3 = coeffs(11);
  c3 = coeffs(12);
  d3 = coeffs(13);
  a4 = coeffs(14);
  b4 = coeffs(15);
  c4 = coeffs(16);
  d4 = coeffs(17);

  if (i==1)  %Return dL1 / dx_j
    if (j==1)
      dLidxj = b1;
    elseif (j==2)
      dLidxj = c1;
    elseif (j==3)
      dLidxj = d1;
    end
  elseif (i==2) % Return dL2 / dx_j
    if (j==1)
      dLidxj = b2;
    elseif (j==2)
      dLidxj = c2;
    elseif (j==3)
      dLidxj = d2;
    end
  elseif (i==3) % Return dL3 / dx_j
    if (j==1)
      dLidxj = b3;
    elseif (j==2)
      dLidxj = c3;
    elseif (j==3)
      dLidxj = d3;
    end
  else % Return dL4 / dx_j
    if (j==1)
      dLidxj = b4;
    elseif (j==2)
      dLidxj = c4;
    elseif (j==3)
      dLidxj = d4;
    end
  end

  % Augment by (1/(6V))
  dLidxj = dLidxj / (6*V);

end

function [dui_dxj] = def_deriv(ui, j, vol_coords, coeffs)

  % Return dui_dxj.
  % Inputs:
  % ui: the nodal data for deformation component in question
  % j: the derivative to take respect to. 1: x, 2: y, 3: z
  % vol_coords: volume coordinates to evaluate derivative at
  % coeffs: volume coordinate coeffs

  dui_dxj = ui(1)*vcd(1,j,coeffs)*(4*vol_coords(1)-1) + ui(2)*vcd(2,j,coeffs)*(4*vol_coords(2)-1) ...
          + ui(3)*vcd(3,j,coeffs)*(4*vol_coords(3)-1) + ui(4)*vcd(4,j,coeffs)*(4*vol_coords(4)-1) ...
          + 4*ui(5)*(vol_coords(1)*vcd(2,j,coeffs) + vol_coords(2)*vcd(1,j,coeffs)) ...
          + 4*ui(6)*(vol_coords(2)*vcd(3,j,coeffs) + vol_coords(3)*vcd(2,j,coeffs)) ...
          + 4*ui(7)*(vol_coords(1)*vcd(3,j,coeffs) + vol_coords(3)*vcd(1,j,coeffs)) ...
          + 4*ui(8)*(vol_coords(4)*vcd(1,j,coeffs) + vol_coords(1)*vcd(4,j,coeffs)) ...
          + 4*ui(9)*(vol_coords(4)*vcd(2,j,coeffs) + vol_coords(2)*vcd(4,j,coeffs)) ...
          + 4*ui(10)*(vol_coords(4)*vcd(3,j,coeffs) + vol_coords(3)*vcd(4,j,coeffs));

end


