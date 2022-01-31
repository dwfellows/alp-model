function [detJ, J, Jinv, L] = jacobian(xi, eta, zeta, global_coords)

  % Constructs the Jacobian utilizing 10-node tet shape functions
  % and returns its determinant

  % N1 = (2*L1 - 1)*L1
  dN1dxi = 1-4*(1-xi-eta-zeta);
  dN1deta = 1-4*(1-xi-eta-zeta);
  dN1dzeta = 1-4*(1-xi-eta-zeta);

  % N2 = (2*L2 - 1)*L2
  dN2dxi = 4*xi-1;
  dN2deta = 0;
  dN2dzeta = 0;

  % N3 = (2*L3 - 1)*L3
  dN3dxi = 0;
  dN3deta = 4*eta-1;
  dN3dzeta = 0;

  % N4 = (2*L4 - 1)*L4
  dN4dxi = 0;
  dN4deta = 0;
  dN4dzeta = 4*zeta-1;

  % N5 = 4*L1*L2
  dN5dxi = 4*(1-2*xi-eta-zeta);
  dN5deta = -4*xi;
  dN5dzeta = -4*xi;

  % N6 = 4*L2*L3
  dN6dxi = 4*eta;
  dN6deta = 4*xi;
  dN6dzeta = 0;

  % N7 = 4*L1*L3
  dN7dxi = -4*eta;
  dN7deta = 4*(1-xi-2*eta-zeta);
  dN7dzeta = -4*eta;

  % N8 = 4*L1*L4
  dN8dxi = -4*zeta;
  dN8deta = -4*zeta;
  dN8dzeta = 4*(1-xi-eta-2*zeta);

  % N9 = 4*L2*L4
  dN9dxi = 4*zeta;
  dN9deta = 0;
  dN9dzeta = 4*xi;

  % N10 = 4*L3*L4
  dN10dxi = 0;
  dN10deta = 4*zeta;
  dN10dzeta = 4*eta;

  L = [dN1dxi, dN2dxi, dN3dxi, dN4dxi, dN5dxi, dN6dxi, dN7dxi, dN8dxi, dN9dxi, dN10dxi;
       dN1deta, dN2deta, dN3deta, dN4deta, dN5deta, dN6deta, dN7deta, dN8deta, dN9deta, dN10deta;
       dN1dzeta, dN2dzeta, dN3dzeta, dN4dzeta, dN5dzeta, dN6dzeta, dN7dzeta, dN8dzeta, dN9dzeta, dN10dzeta];
 
  % Construct Jacobian: [
  J = L*global_coords;

  detJ = det(J);

  Jinv = inv(J);

end
