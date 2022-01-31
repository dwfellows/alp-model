function [quantity_interp] = point_interpolation(point_loc, cell_point_locs, quantity)

  [L1, L2, L3] = barycentric(point_loc, cell_point_locs);
  quantity_interp = L1*quantity(1) + L2*quantity(2) + L3*quantity(3);

end

