function [dui_dxj] = def_deriv(ui, j, vol_coords, Jinv, coeffs)

  % Return dui_dxj.
  % Inputs:
  % ui: the nodal data for deformation component in question
  % j: the derivative to take respect to. 1: x, 2: y, 3: z
  % vol_coords: volume coordinates to evaluate derivative at
  % coeffs: volume coordinate coeffs

  xi = vol_coords(2);
  eta = vol_coords(3);
  zeta = vol_coords(4);

  dui_dxi = ui(1)*(4*(1-xi-eta-zeta)*(-1) + 1) + ui(2)*(4*xi-1) + 4*ui(5)*(1-xi-eta-zeta-xi) ...
          + 4*ui(6)*eta + 4*ui(7)*(-eta) + 4*ui(8)*(-zeta) + 4*ui(9)*(zeta);
  dui_deta = ui(1)*(4*(1-xi-eta-zeta)*(-1) + 1) + ui(3)*(4*eta-1) + 4*ui(5)*(-xi) ...
           + 4*ui(6)*(xi) + 4*ui(7)*(-eta + (1-xi-eta-zeta)) + 4*ui(8)*(-zeta) + 4*ui(10)*(zeta);
  dui_dzeta = ui(1)*(4*(1-xi-eta-zeta)+1) + ui(4)*(4*zeta-1) + 4*ui(5)*(-xi) + 4*ui(7)*(-eta) ...
            + 4*ui(8)*(-zeta + (1-xi-eta-zeta)) + 4*ui(9)*(xi) + 4*ui(10)*(eta);

  dui = Jinv*[dui_dxi; dui_deta; dui_dzeta];

  if (j==1)
    dui_dxj = dui(1);
  else if (j==2)
    dui_dxj = dui(2);
  else
    dui_dxj = dui(3);
  end


%  dui_dxj = ui(1)*vcd(1,j,coeffs)*(4*vol_coords(1)-1) + ui(2)*vcd(2,j,coeffs)*(4*vol_coords(2)-1) ...
%          + ui(3)*vcd(3,j,coeffs)*(4*vol_coords(3)-1) + ui(4)*vcd(4,j,coeffs)*(4*vol_coords(4)-1) ...
%          + 4*ui(5)*(vol_coords(1)*vcd(2,j,coeffs) + vol_coords(2)*vcd(1,j,coeffs)) ...
%          + 4*ui(6)*(vol_coords(2)*vcd(3,j,coeffs) + vol_coords(3)*vcd(2,j,coeffs)) ...
%          + 4*ui(7)*(vol_coords(1)*vcd(3,j,coeffs) + vol_coords(3)*vcd(1,j,coeffs)) ...
%          + 4*ui(8)*(vol_coords(4)*vcd(1,j,coeffs) + vol_coords(1)*vcd(4,j,coeffs)) ...
%          + 4*ui(9)*(vol_coords(4)*vcd(2,j,coeffs) + vol_coords(2)*vcd(4,j,coeffs)) ...
%          + 4*ui(10)*(vol_coords(4)*vcd(3,j,coeffs) + vol_coords(3)*vcd(4,j,coeffs));

end

